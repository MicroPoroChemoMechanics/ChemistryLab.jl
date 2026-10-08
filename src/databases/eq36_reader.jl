# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# ── Reading a data0 file of EQ3/6 ────────────────────────────────────────────
#
# A data0 file is a database of reactions, as a PHREEQC database is, laid out in
# blocks separated by lines of `+---`: a block of parameters, the hard core
# diameters, the elements, then the species, each with its composition, its
# dissociation reaction and its log K on a grid of eight temperatures
# [DavelerWolery1992; chapter 4](@cite). It is read into the same representation
# as a PHREEQC database, its energies counted from its strict
# basis species, which are at zero.

# A log K on the grid 0-25-60-100 / 150-200-250-300 °C, "500" marking a lack of
# data. EQPT, which prepares a data0 file for EQ3/6, fits an interpolating
# polynomial through the valid points of each part of the grid, 0-100 °C and
# 100-300 °C, of degree three at most below and four above, a constant where a
# part has a single point [DavelerWolery1992; § 3.1](@cite).
struct GridLogK <: AbstractLogK
    low::Vector{Float64}                  # 0-100 °C, coefficients in °C
    high::Vector{Float64}                 # 100-300 °C
    split::Float64                        # °C
end

function _interpolating(t, v)
    isempty(t) && return Float64[]
    V = [ti^k for ti in t, k in 0:(length(t) - 1)]
    return V \ v
end

function GridLogK(temperatures::AbstractVector, values::AbstractVector)
    length(temperatures) == 8 || throw(ArgumentError("a log K grid has eight temperatures; got $(length(temperatures))"))
    valid(r) = [i for i in r if values[i] != 500.0]
    lo, hi = valid(1:4), valid(4:8)
    isempty(lo) && isempty(hi) && throw(ArgumentError("a log K grid without any value"))
    low = _interpolating(Float64.(temperatures[lo]), Float64.(values[lo]))
    high = _interpolating(Float64.(temperatures[hi]), Float64.(values[hi]))
    # A part without data takes the other's polynomial, the species being
    # valid only where it has data (EQPT suppresses it there).
    isempty(low) && (low = high)
    isempty(high) && (high = low)
    return GridLogK(low, high, Float64(temperatures[4]))
end

function _log10K(k::GridLogK, T)
    tc = _celsius(T)
    p = _plain(tc) <= k.split ? k.low : k.high
    L, dL, d2L = zero(tc), zero(tc), zero(tc)
    for (j, a) in enumerate(p)
        n = j - 1
        L += a * tc^n
        n >= 1 && (dL += n * a * tc^(n - 1))
        n >= 2 && (d2L += n * (n - 1) * a * tc^(n - 2))
    end
    return L, dL, d2L
end

# A reaction line: one or two terms, a coefficient then a name that may hold
# single spaces (`acetic acid(aq)`), as the fixed format 2(1x,f10.4,2x,a24) sets
# them out.
function _eq36_terms(line)
    out = Tuple{Float64, String}[]
    for m in eachmatch(r"(-?\d+\.\d*)\s{1,3}(\S(?:\S| (?! ))*)", line)
        push!(out, (parse(Float64, m.captures[1]), String(m.captures[2])))
    end
    return out
end

"""
    read_eq36_database(path) -> (df_elements, df_substances, df_reactions)

Read a data0 file of EQ3/6, of the `com` archetype (an extended Debye-Hückel
activity model): its parameters, its elements, and its species, auxiliary basis
species, minerals, liquids and gases, each with its composition, its reaction and
its log K on the temperature grid [DavelerWolery1992; chapter 4](@cite).

The species are resolved against the strict basis species as those of a PHREEQC
database are against its master species, and laid out in the same three tables
as [`read_phreeqc_database`](@ref) lays them out, from which
[`build_species`](@ref) builds: their energies are counted from the strict basis
species, which are at zero, a gauge recorded under `:gauge`. Their
composition is the one each block states, checked against the balance of its
reaction. A log K is evaluated between the temperatures of the grid as EQPT
evaluates it, by an interpolating polynomial through the valid points of each
part of the grid. The solid solutions and the Pitzer parameters of the `hmw`
archetype are noted and not read.

No data0 file is distributed with the package or downloaded by it; a copy is
installed by hand ([`install_database`](@ref)).
"""
function read_eq36_database(path::AbstractString)
    file = basename(path)
    digest = bytes2hex(open(sha256, path))
    raw = [(n, rstrip(l)) for (n, l) in enumerate(eachline(path))]
    notes = String[]
    first_line = isempty(raw) ? "" : strip(raw[1][2])
    startswith(lowercase(first_line), "data0") ||
        throw(ArgumentError("$file: the first line is `$first_line`, not the header of a data0 file"))
    archetype = length(first_line) >= 9 ? lowercase(first_line[7:9]) : "com"
    archetype in ("hmw", "pit") && push!(notes, "$file: archetype `$archetype`, whose Pitzer parameters are not read")
    # Comment lines begin with `*` in the first column.
    lines = [(n, l) for (n, l) in raw if !startswith(l, "*")]
    text(i) = strip(lines[i][2])
    i = 1
    findnext_line(pred, from) = findnext(j -> pred(lowercase(text(j))), eachindex(lines), from)

    # The parameters block: the grid, then A, B, B-dot, the coefficients of CO2
    # and the log K of the reaction of the electron, each introduced by its name.
    numbers_after(k, n) = begin
        out = Float64[]
        j = k + 1
        while length(out) < n && j <= length(lines)
            v = tryparse.(Float64, split(text(j)))
            all(!isnothing, v) && !isempty(v) && append!(out, Float64.(v))
            j += 1
        end
        out
    end
    k = findnext_line(l -> startswith(l, "temperatures"), 1)
    k === nothing && throw(ArgumentError("$file: no grid of temperatures"))
    temperatures = numbers_after(k, 8)
    parameters = Dict{String, Vector{Float64}}()
    for (key, n) in (("pressures", 8), ("debye huckel a", 8), ("debye huckel b", 8), ("bdot", 8), ("cco2", 5), ("log k for eh reaction", 8))
        j = findnext_line(l -> startswith(l, key), k)
        j === nothing || (parameters[key] = numbers_after(j, n))
    end

    # The blocks: a header line names a superblock, a block a species.
    superblocks = ("basis species", "auxiliary basis species", "aqueous species", "solids", "liquids", "gases", "solid solutions", "references", "stop.")
    blocks = Vector{Vector{Tuple{Int, String}}}()
    current = Tuple{Int, String}[]
    section = ""
    sections = String[]
    started = false
    for (n, l) in lines
        s = strip(l)
        if startswith(s, "+--")
            isempty(current) || (push!(blocks, current); push!(sections, section))
            current = Tuple{Int, String}[]
            continue
        end
        low = lowercase(s)
        if low in superblocks || low == "elements" || startswith(low, "bdot parameters")
            isempty(current) || (push!(blocks, current); push!(sections, section))
            current = Tuple{Int, String}[]
            section = low
            started = true
            continue
        end
        started && !isempty(s) && push!(current, (n, String(l)))
    end
    isempty(current) || (push!(blocks, current); push!(sections, section))

    sizes = Dict{String, Tuple{Float64, Int}}()
    records = NamedTuple[]
    basis = String[]
    for (blk, sec) in zip(blocks, sections)
        if startswith(sec, "bdot parameters")
            for (_, l) in blk
                m = match(r"^(.*?)\s+(-?\d+\.\d*)\s+(-?\d+)\s*$", l)
                m === nothing || (sizes[_phreeqc_name(strip(m.captures[1]))] = (parse(Float64, m.captures[2]), parse(Int, m.captures[3])))
            end
            continue
        end
        sec in ("basis species", "auxiliary basis species", "aqueous species", "solids", "liquids", "gases") || continue
        line0, head = blk[1]
        # The name fills the first 24 columns; what follows is a formula.
        name = _phreeqc_name(strip(first(head, min(24, length(head)))))
        charge, volume = 0.0, nothing
        comp = Dict{Symbol, Float64}()
        terms = Tuple{Float64, String}[]
        grid = Float64[]
        j = 2
        while j <= length(blk)
            l = blk[j][2]
            low = lowercase(l)
            if (m = match(r"charge\s*=\s*(-?[\d.]+)", low)) !== nothing
                charge = parse(Float64, m.captures[1])
            elseif (m = match(r"v0prtr\s*=\s*(-?[\d.eE+]+)", low)) !== nothing
                volume = parse(Float64, m.captures[1])
            elseif (m = match(r"^\s*(\d+)\s+chemical elements", low)) !== nothing
                nel = parse(Int, m.captures[1])
                while length(comp) < nel && j < length(blk)
                    j += 1
                    for t in eachmatch(r"(-?\d+\.\d*)\s+([a-z][a-z]?)\b", lowercase(blk[j][2]))
                        e = Symbol(uppercasefirst(t.captures[2]))
                        comp[e] = parse(Float64, t.captures[1])
                    end
                end
            elseif (m = match(r"^\s*(\d+)\s+species in (?:data0 )?reaction", low)) !== nothing
                nt = parse(Int, m.captures[1])
                while length(terms) < nt && j < length(blk)
                    j += 1
                    append!(terms, [(c, _phreeqc_name(s)) for (c, s) in _eq36_terms(blk[j][2])])
                end
            elseif occursin(r"log k grid", low)
                while length(grid) < 8 && j < length(blk)
                    j += 1
                    v = tryparse.(Float64, split(blk[j][2]))
                    all(!isnothing, v) && append!(grid, Float64.(v))
                end
            end
            j += 1
        end
        sec == "basis species" && push!(basis, name)
        push!(records, (; name, section = sec, line = line0, charge, volume, comp, terms, grid))
    end
    return _eq36_database(file, String(path), digest, archetype, temperatures, parameters, sizes, basis, records, notes)
end

function _eq36_database(file, path, digest, archetype, temperatures, parameters, sizes, basis, records, notes)
    logks = AbstractLogK[]
    by_name = Dict{String, Any}()
    for rec in records
        if haskey(by_name, rec.name)
            push!(notes, "$file:$(rec.line): $(rec.name) is defined again; the definition at line $(by_name[rec.name][1].line) is replaced")
        end
        k = 0
        if !(rec.name in basis)
            if length(rec.grid) != 8
                push!(notes, "$file:$(rec.line): $(rec.name) has no log K grid of eight values; not read")
                continue
            end
            push!(logks, GridLogK(temperatures, rec.grid))
            k = length(logks)
        end
        by_name[rec.name] = (rec, k)
    end
    resolved = Dict{String, Any}()
    visiting = Set{String}()
    function resolve(name)
        haskey(resolved, name) && return resolved[name]
        haskey(by_name, name) || return (resolved[name] = nothing)
        name in visiting && return nothing
        push!(visiting, name)
        rec, k = by_name[name]
        out = if name in basis
            (Dict{Symbol, Float64}(rec.comp), rec.charge, Dict{Int, Float64}())
        else
            # The reaction lists the species itself, with a negative coefficient.
            _resolve_defined(name, rec.terms, k, logks, resolve)
        end
        delete!(visiting, name)
        return (resolved[name] = out)
    end
    species, phases = _ReactionEntry[], _ReactionEntry[]
    for (rec, _) in sort!(collect(values(by_name)); by = x -> x[1].line)
        r = resolve(rec.name)
        if r === nothing
            push!(notes, "$file:$(rec.line): $(rec.name) rests on a species the file does not define; not read")
            continue
        end
        atoms, z, formation = r
        stated = Dict{Symbol, Float64}(e => n for (e, n) in rec.comp if !iszero(n))
        if !isempty(stated) && (keys(stated) != keys(atoms) || any(abs(stated[e] - atoms[e]) > 1.0e-6 for e in keys(atoms)) || abs(z - rec.charge) > 1.0e-6)
            push!(notes, "$file:$(rec.line): the reaction of $(rec.name) does not balance against its stated composition; not read")
            continue
        end
        unknown = [e for e in keys(atoms) if !haskey(elements.bysymbol, e)]
        if !isempty(unknown)
            push!(notes, "$file:$(rec.line): $(rec.name) is made of $(join(string.(unknown), ", ")), not of chemical elements; not read")
            continue
        end
        if rec.section in ("solids", "liquids", "gases") || endswith(rec.name, "(g)")
            push!(phases, _ReactionEntry(rec.name, :phase, atoms, z, formation, nothing, nothing, false, rec.volume, nothing, rec.line, rec.name))
        else
            å = get(sizes, rec.name, nothing)
            push!(species, _ReactionEntry(rec.name, :solution, atoms, z, formation, nothing, å === nothing ? nothing : å[1], false, nothing, nothing, rec.line, rec.name))
        end
    end
    masters = Dict{String, String}(b => b for b in basis)
    gauge = "strict basis species of $file ($(digest[1:12]))"
    params = (; archetype, temperatures, parameters, sizes)
    reactions = NamedTuple[
        (; symbol = rec.name, kind = Symbol(replace(rec.section, " " => "_")), equation = join(("$(c) $(n)" for (c, n) in rec.terms), " + "), logk = logks[k], line = rec.line)
            for (rec, k) in values(by_name) if k > 0
    ]
    sort!(reactions; by = r -> r.line)
    data = _ReactionData(
        file, path, "$file sha256 $(digest[1:12])", gauge, :eq36, masters, species, phases, logks,
        reactions, params, Set{String}(), notes,
    )
    return _reaction_tables(data)
end
