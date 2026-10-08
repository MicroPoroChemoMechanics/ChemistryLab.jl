# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# ── Reading a thermo dataset of The Geochemist's Workbench ───────────────────
#
# A thermo dataset (`.tdat`, formerly `thermo.dat`) is a database of reactions,
# as a PHREEQC database is: its basis species carry their elemental composition,
# and every other species, mineral and gas a dissociation reaction into species
# of the dataset with the log K of that reaction, as a function of temperature
# [BethkeFarrell2026; chapter 3](@cite). It is read into the same representation
# as a PHREEQC database, its energies counted from its basis
# species, which are at zero.
#
# What the reference manual describes is read: the sections in order (elements,
# basis species, redox couples, aqueous species, the free electron, minerals,
# solid solutions, gases, oxides), the log K as a table at the principal
# temperatures or as the polynomial of format "jan19", and the header variables.
# The solid solutions, the oxide components and the virial coefficients are noted
# and not read.

# A log K given at the principal temperatures, "500" marking a lack of data: the
# dataset's applications fit the values present to a polynomial of degree four in
# the temperature in °C [BethkeFarrell2026; § 3.2.1](@cite), which is evaluated
# here at every temperature (they take the tabulated value at a principal
# temperature itself, which with `span = on` they do not, and the polynomial is
# the continuous one). The fit is by least squares, the number of values present
# exceeding five as a rule; with five or fewer, the polynomial of lower degree
# through them.
struct TabulatedLogK <: AbstractLogK
    poly::Vector{Float64}                 # coefficients in °C, constant first
    range::Tuple{Float64, Float64}        # °C, the principal temperatures with data
end

function TabulatedLogK(temperatures::AbstractVector, values::AbstractVector)
    keep = [i for i in eachindex(values) if values[i] != 500.0]
    isempty(keep) && throw(ArgumentError("a log K table without any value"))
    t, v = Float64.(temperatures[keep]), Float64.(values[keep])
    degree = min(4, length(t) - 1)
    V = [ti^k for ti in t, k in 0:degree]
    return TabulatedLogK(V \ v, (minimum(t), maximum(t)))
end

function _log10K(k::TabulatedLogK, T)
    tc = _celsius(T)
    L, dL, d2L = zero(tc), zero(tc), zero(tc)
    for (j, a) in enumerate(k.poly)
        n = j - 1
        L += a * tc^n
        n >= 1 && (dL += n * a * tc^(n - 1))
        n >= 2 && (d2L += n * (n - 1) * a * tc^(n - 2))
    end
    return L, dL, d2L
end

# The polynomial of format "jan19" [BethkeFarrell2026; § 3.2.2](@cite),
#     a + b (T − Tr) + c (T² − Tr²) + d (1/T − 1/Tr) + e (1/T² − 1/Tr²) + f ln(T/Tr),
# T in kelvin and Tr = 298.15 K.
struct GWBPolynomialLogK <: AbstractLogK
    coeffs::NTuple{6, Float64}
    range::Tuple{Float64, Float64}        # K
end

function _log10K(k::GWBPolynomialLogK, T)
    a, b, c, d, e, f = k.coeffs
    Tr = _T25
    L = a + b * (T - Tr) + c * (T^2 - Tr^2) + d * (1 / T - 1 / Tr) + e * (1 / T^2 - 1 / Tr^2) + f * log(T / Tr)
    dL = b + 2c * T - d / T^2 - 2e / T^3 + f / T
    d2L = 2c + 2d / T^3 + 6e / T^4 - f / T^2
    return L, dL, d2L
end

# The lines of a dataset, blank lines removed, with their numbers; a comment
# line (`*` in the first column) is kept as `*` alone, for the comments that
# separate the tables of the header, and skipped everywhere else.
function _gwb_lines(path)
    out = Tuple{Int, String}[]
    for (n, raw) in enumerate(eachline(path))
        startswith(raw, "*") && (push!(out, (n, "*")); continue)
        s = strip(raw)
        isempty(s) || push!(out, (n, String(s)))
    end
    return out
end

# A cursor over the lines, read word by word as the programs read them.
mutable struct _GWBCursor
    lines::Vector{Tuple{Int, String}}
    i::Int
end
function _skip_comments!(c::_GWBCursor)
    while c.i <= length(c.lines) && c.lines[c.i][2] == "*"
        c.i += 1
    end
    return c
end
_peek(c::_GWBCursor) = (_skip_comments!(c); c.i <= length(c.lines) ? c.lines[c.i][2] : nothing)
_line(c::_GWBCursor) = (_skip_comments!(c); c.i <= length(c.lines) ? c.lines[c.i][1] : 0)
function _next!(c::_GWBCursor)
    _skip_comments!(c)
    c.i > length(c.lines) && throw(ArgumentError("the dataset ends in the middle of an entry"))
    c.i += 1
    return c.lines[c.i - 1][2]
end
_numbers_on(line) = (v = tryparse.(Float64, split(line)); any(isnothing, v) ? nothing : Float64.(v))

# The numbers on the following lines, `n` of them.
function _numbers!(c::_GWBCursor, n)
    out = Float64[]
    while length(out) < n
        line = _next!(c)
        v = _numbers_on(line)
        v === nothing && throw(ArgumentError("line $(c.lines[c.i - 1][1]): `$line` where $n numbers were expected"))
        append!(out, v)
    end
    length(out) == n || throw(ArgumentError("line $(c.lines[c.i - 1][1]): $(length(out)) numbers where $n were expected"))
    return out
end

# A temperature expansion: a table at the principal temperatures, or the
# polynomial "a= … f=" with an optional range line.
function _gwb_expansion!(c::_GWBCursor, temperatures)
    line = _peek(c)
    if line !== nothing && startswith(line, "a=")
        words = String[]
        while (l = _peek(c)) !== nothing && occursin(r"^[a-f]=", l)
            append!(words, split(replace(_next!(c), "=" => "= ")))
        end
        coeffs = zeros(6)
        for (j, key) in enumerate(("a=", "b=", "c=", "d=", "e=", "f="))
            k = findfirst(==(key), words)
            k === nothing || (coeffs[j] = parse(Float64, words[k + 1]))
        end
        range = (_celsius_to_kelvin(first(temperatures)), _celsius_to_kelvin(last(temperatures)))
        if (l = _peek(c)) !== nothing && occursin(r"^T(min|max)[KC]=", l)
            w = split(replace(_next!(c), "=" => "= "))
            vals = [parse(Float64, w[k + 1]) for k in eachindex(w) if endswith(w[k], "=")]
            range = occursin("TminC", l) ? (_celsius_to_kelvin(vals[1]), _celsius_to_kelvin(vals[2])) : (vals[1], vals[2])
        end
        return GWBPolynomialLogK(Tuple(coeffs), range)
    end
    return TabulatedLogK(temperatures, _numbers!(c, length(temperatures)))
end

# The header of an entry, from its name to the line that counts its reaction or
# its elements: the fields the reference manual names.
function _gwb_entry_header!(c::_GWBCursor)
    text = String[]
    while (l = _peek(c)) !== nothing && !occursin(r"^\d+\s+(species in reaction|elements in species)", l)
        push!(text, _next!(c))
    end
    h = join(text, " ")
    field(re) = (m = match(re, h); m === nothing ? nothing : parse(Float64, m.captures[1]))
    return (;
        charge = something(field(r"charge=\s*(-?[\d.]+)"), 0.0),
        ion_size = field(r"ion size=\s*(-?[\d.]+)"),
        volume = field(r"mole vol\.=\s*(-?[\d.eE+]+)"),
        Tc = field(r"Tcrit=\s*([\d.eE+]+)"),
        Pc = field(r"Pcrit=\s*([\d.eE+]+)"),
        omega = field(r"omega=\s*(-?[\d.eE+]+)"),
    )
end

# The `N species in reaction` line and the reaction, three terms a line; names
# with spaces are quoted.
function _gwb_reaction!(c::_GWBCursor)
    head = _next!(c)
    m = match(r"^(\d+)\s+species in reaction", head)
    m === nothing && throw(ArgumentError("line $(c.lines[c.i - 1][1]): `$head` where the count of a reaction was expected"))
    n = parse(Int, m.captures[1])
    terms = Tuple{Float64, String}[]
    while length(terms) < n
        for t in eachmatch(r"(-?[\d.]+(?:[eE][-+]?\d+)?)\s+(\"[^\"]+\"|\S+)", _next!(c))
            push!(terms, (parse(Float64, t.captures[1]), _phreeqc_name(strip(t.captures[2], '"'))))
        end
    end
    return terms
end

const _GWB_SPECIES_SECTIONS = ("basis species", "redox couples", "aqueous species", "free electron", "minerals", "gases")

"""
    read_gwb_database(path) -> (df_elements, df_substances, df_reactions)

Read a thermo dataset of The Geochemist's Workbench (`.tdat`, format "jan19" or
later, and the earlier formats whose sections are the same): its elements, basis
species, redox couples, aqueous species, free electron, minerals and gases, each
with its reaction and its log K, as a table at the principal temperatures or as
a polynomial [BethkeFarrell2026; chapter 3](@cite).

The species are resolved against the basis species as those of a PHREEQC
database are against its master species, and laid out in the same three tables
as [`read_phreeqc_database`](@ref) lays them out, from which
[`build_species`](@ref) builds: their energies are counted from the basis
species, which are at zero, a gauge recorded under `:gauge`. A log K table is
evaluated as the dataset's applications evaluate it between principal
temperatures, by the polynomial of degree four fitted to it. What is read and
not used is listed in `notes`: the solid solutions, the oxide components, the
virial coefficients.

No thermo dataset is distributed with the package or downloaded by it: the
datasets of The Geochemist's Workbench come under the terms of their publishers,
and a copy is installed by hand ([`install_database`](@ref)).
"""
function read_gwb_database(path::AbstractString)
    file = basename(path)
    digest = bytes2hex(open(sha256, path))
    c = _GWBCursor(_gwb_lines(path), 1)
    notes = String[]
    header = Dict{String, String}()
    # The initial lines, up to the first table of numbers.
    while (l = _peek(c)) !== nothing && _numbers_on(l) === nothing
        m = match(r"^([a-z ]+):\s*(.*)$", _next!(c))
        m === nothing || (header[strip(m.captures[1])] = strip(m.captures[2]))
    end
    # The principal temperatures, four a line, up to the comment or the first
    # value that does not continue the increasing list.
    temperatures = Float64[]
    _skip_comments!(c)
    while c.i <= length(c.lines) && c.lines[c.i][2] != "*" && (v = _numbers_on(c.lines[c.i][2])) !== nothing
        !isempty(temperatures) && first(v) <= last(temperatures) && break
        append!(temperatures, v)
        c.i += 1
    end
    isempty(temperatures) && throw(ArgumentError("$file: no principal temperatures; not a thermo dataset of The Geochemist's Workbench"))
    (l = _peek(c)) !== nothing && occursin(r"^T(min|max)C=", l) && _next!(c)
    # The header variables, up to the elements: pressure, Debye-Hückel A and B,
    # then, for the B-dot model, B-dot and the coefficients for CO2 and water.
    variables = AbstractLogK[]
    while (l = _peek(c)) !== nothing && match(r"^\d+\s+elements$", l) === nothing
        occursin(r"^P(min|max)=", l) ? _next!(c) : push!(variables, _gwb_expansion!(c, temperatures))
    end
    _peek(c) === nothing && throw(ArgumentError("$file: no section of elements; not a thermo dataset of The Geochemist's Workbench"))

    element_names = Dict{String, Symbol}()
    basis = Dict{String, Tuple{Dict{Symbol, Float64}, Float64}}()
    records = NamedTuple[]
    while (l = _peek(c)) !== nothing
        m = match(r"^(\d+)\s+(.+)$", l)
        if m === nothing
            push!(notes, "$file:$(_line(c)): `$l` outside any section; not read")
            _next!(c)
            continue
        end
        section = strip(m.captures[2])
        start = _line(c)
        _next!(c)
        while (l = _peek(c)) !== nothing && !startswith(l, "-end-")
            line_no = _line(c)
            if section == "elements"
                e = match(r"^(.+?)\s*\(\s*(\S+)\s*\)", _next!(c))
                e === nothing || (element_names[strip(e.captures[1])] = Symbol(e.captures[2]))
            elseif section in _GWB_SPECIES_SECTIONS
                name = _phreeqc_name(strip(first(split(_next!(c), "formula="))))
                h = _gwb_entry_header!(c)
                if section == "basis species"
                    k = match(r"^(\d+)\s+elements in species", _next!(c))
                    k === nothing && throw(ArgumentError("$file:$line_no: basis species $name gives no composition"))
                    comp = Dict{Symbol, Float64}()
                    while length(comp) < parse(Int, k.captures[1])
                        for t in eachmatch(r"(-?[\d.]+)\s+(\S+)", _next!(c))
                            comp[Symbol(t.captures[2])] = parse(Float64, t.captures[1])
                        end
                    end
                    basis[name] = (comp, h.charge)
                    push!(records, (; name, section, line = line_no, header = h, terms = Tuple{Float64, String}[], logk = nothing))
                else
                    terms = _gwb_reaction!(c)
                    logk = _gwb_expansion!(c, temperatures)
                    push!(records, (; name, section, line = line_no, header = h, terms, logk))
                end
                (l2 = _peek(c)) !== nothing && startswith(l2, "Vm(appelo)") && _next!(c)
            else
                _next!(c)
            end
        end
        section in ("elements", _GWB_SPECIES_SECTIONS...) ||
            push!(notes, "$file:$start: section `$section` is not read")
        _peek(c) === nothing || _next!(c)        # -end-
    end
    return _gwb_database(file, String(path), digest, header, temperatures, variables, element_names, basis, records, notes)
end

# The species of a dataset resolved against its basis species: a basis species
# by its composition, at zero; any other by the balance of its dissociation
# reaction, `s = Σ cᵢ Pᵢ`, and the log K of its formation from the basis.
function _gwb_database(file, path, digest, header, temperatures, variables, element_names, basis, records, notes)
    logks = AbstractLogK[]
    by_name = Dict{String, Any}()
    for rec in records
        haskey(by_name, rec.name) &&
            push!(notes, "$file:$(rec.line): $(rec.name) is defined again; the definition at line $(by_name[rec.name][1].line) is replaced")
        k = rec.logk === nothing ? 0 : (push!(logks, rec.logk); length(logks))
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
        out = if haskey(basis, name)
            comp, z = basis[name]
            (Dict{Symbol, Float64}(comp), z, Dict{Int, Float64}())
        else
            terms = vcat([(-1.0, "\$self")], rec.terms)
            r = _resolve_defined("\$self", terms, k, logks, resolve)
            r === nothing ? nothing : (r[1], r[2], r[3])
        end
        delete!(visiting, name)
        return (resolved[name] = out)
    end
    species, phases = _ReactionEntry[], _ReactionEntry[]
    for rec in values(by_name) |> collect |> v -> sort(v; by = x -> x[1].line)
        r = rec[1]
        res = resolve(r.name)
        if res === nothing
            push!(notes, "$file:$(r.line): $(r.name) rests on a species the dataset does not define; not read")
            continue
        end
        atoms, z, formation = res
        unknown = [e for e in keys(atoms) if !haskey(elements.bysymbol, e)]
        if !isempty(unknown)
            push!(notes, "$file:$(r.line): $(r.name) is made of $(join(string.(unknown), ", ")), not of chemical elements; not read")
            continue
        end
        h = r.header
        if r.section in ("minerals", "gases")
            # The critical pressure in atm, as the entries of the other formats hold it.
            critical = (h.Tc === nothing || h.Pc === nothing) ? nothing :
                (h.Tc, ustrip(us"Pa", h.Pc * u"bar") / ustrip(us"Pa", _DQConstants.atm), something(h.omega, 0.0))
            push!(phases, _ReactionEntry(r.name, :phase, atoms, z, formation, nothing, nothing, false, h.volume, critical, r.line, r.name))
        else
            push!(species, _ReactionEntry(r.name, :solution, atoms, z, formation, nothing, h.ion_size, false, nothing, nothing, r.line, r.name))
        end
    end
    masters = Dict{String, String}(name => name for name in keys(basis))
    gauge = "basis species of $file ($(digest[1:12]))"
    parameters = (; format = get(header, "dataset format", "unstated"), activity_model = get(header, "activity model", "unstated"), temperatures, variables)
    reactions = NamedTuple[
        (; symbol = rec.name, kind = Symbol(replace(rec.section, " " => "_")), equation = join(("$(c) $(n)" for (c, n) in rec.terms), " + "), logk = rec.logk, line = rec.line)
            for rec in records if rec.logk !== nothing
    ]
    data = _ReactionData(
        file, path, "$file sha256 $(digest[1:12])", gauge, :gwb, masters, species, phases, logks,
        reactions, parameters, Set{String}(), notes,
    )
    return _reaction_tables(data)
end
