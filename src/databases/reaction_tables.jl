# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# ── Databases of reactions ────────────────────────────────────────────────────
#
# A PHREEQC database, a thermo dataset of The Geochemist's Workbench and a data0
# file of EQ3/6 define their species the same way: a few basis (master) species
# carry their composition, and every other species, mineral and gas is defined by
# a reaction among species of the database and the log K of that reaction as a
# function of temperature. They are read into one representation (`_ReactionData`)
# and laid out in the same three tables (`_reaction_tables`), from which
# `build_species` builds: what the formats share is written here once, the
# formats themselves in `phreeqc_reader.jl`, `gwb_reader.jl` and `eq36_reader.jl`.
# The readers of a database of formation properties (ThermoFun, Reaktoro) share
# the naming rule and the builder of a species from a row.

# ── Provenance ───────────────────────────────────────────────────────────────

# What a value read from `path` says of where it comes from: the file and the
# start of its SHA-256 (`_digest`, computed once per file and session).
_short_digest(path) = first(_digest(path), 12)
_source_tag(path) = "$(basename(path)) sha256 $(_short_digest(path))"

# ── Species names ────────────────────────────────────────────────────────────

# A species name as PHREEQC means it: a run of signs is a charge (`SO4--` is
# `SO4-2`), and a charge of one is written with its sign alone. The names of
# GWB, EQ3/6 and of the PITZER and SIT blocks are read through it too.
function _phreeqc_name(name::AbstractString)
    s = String(strip(name))
    m = match(r"^(.*?)(\++|-+)$", s)
    if m !== nothing && length(m.captures[2]) > 1 && !isempty(m.captures[1])
        return m.captures[1] * m.captures[2][1:1] * string(length(m.captures[2]))
    end
    m = match(r"^(.*[^+-])([+-])1$", s)
    m === nothing || return m.captures[1] * m.captures[2]
    return s
end

# The ChemistryLab symbol of a dissolved species from the name a database gives
# it: the `(aq)` with which GWB, EQ3/6 and Reaktoro name a neutral solute is
# dropped, a neutral solute gets `@`, water is the solvent `H2O@`, and an ion
# keeps its name. Every reader applies this one rule.
function _solute_symbol(name::AbstractString, charge = _name_charge(name))
    base = replace(name, r"\(aq\)$"i => "")
    lowercase(base) == "h2o" && return "H2O@"
    return iszero(charge) ? base * "@" : String(base)
end

# The class of a species of a table, from its state and its symbol.
_class_of(state::AggregateState, symbol::AbstractString) =
    state == AS_AQUEOUS ? (symbol == "H2O@" ? SC_AQSOLVENT : SC_AQSOLUTE) :
    state == AS_GAS ? SC_GASFLUID : SC_COMPONENT

# Whether a row is among the species named (by the database's name or by the
# ChemistryLab symbol), all of them when none is.
_listed(row, names) = names === nothing || row.name in names || row.symbol in names

# The species of a row of a table, without its thermodynamics: its composition
# in integers when every count is one.
function _row_species(row)
    atoms = all(isinteger, values(row.atoms)) ? Dict{Symbol, Int}(e => Int(n) for (e, n) in row.atoms) :
        Dict{Symbol, Float64}(e => n for (e, n) in row.atoms)
    z = row.charge
    return Species(
        atoms, isinteger(z) ? Int(z) : z;
        name = row.name, symbol = row.symbol, aggregate_state = row.aggregate_state, class = row.class,
    )
end

# ── Equilibrium constants ────────────────────────────────────────────────────

"""
    AbstractLogK

The log K of one reaction of a database of reactions, as a function of
temperature, in the form its format gives: `_log10K(k, T)` returns it with its
first two derivatives.
"""
abstract type AbstractLogK end

# The constants a log K adds to its own (PHREEQC's `-add_logk`); none for the
# other formats.
_added(::AbstractLogK) = ()

# The standard properties of a reaction from its log K and the first two
# temperature derivatives of it, `r = (L, L′, L″)` (van 't Hoff):
#     ΔG = −RT ln 10 L,   ΔH = RT² ln 10 L′,   ΔS = R ln 10 (L + T L′),
#     ΔCp = R ln 10 (2T L′ + T² L″).
_logK_gibbs(r, T) = -_R_LN10 * T * r[1]
_logK_enthalpy(r, T) = _R_LN10 * T^2 * r[2]
_logK_entropy(r, T) = _R_LN10 * (r[1] + T * r[2])
_logK_heat_capacity(r, T) = _R_LN10 * (2T * r[2] + T^2 * r[3])

# A polynomial of the temperature in °C, coefficients constant first, with its
# first two derivatives: the log K of GWB tables and of EQ3/6 grids.
function _celsius_polynomial(coefficients, tc)
    L, dL, d2L = zero(tc), zero(tc), zero(tc)
    for (j, a) in enumerate(coefficients)
        n = j - 1
        L += a * tc^n
        n >= 1 && (dL += n * a * tc^(n - 1))
        n >= 2 && (d2L += n * (n - 1) * a * tc^(n - 2))
    end
    return L, dL, d2L
end

# The coefficients, constant first, of the polynomial of degree `degree` that
# fits the values `v` at the temperatures `t` (°C): by least squares, through the
# points when there are `degree + 1` of them.
_polynomial_fit(t, v, degree = length(t) - 1) = [ti^k for ti in t, k in 0:degree] \ v

# The numbers of a line, or `nothing` when a word of it is not a number: the
# tables of GWB and EQ3/6.
_numbers_on(line) = (v = tryparse.(Float64, split(line)); any(isnothing, v) ? nothing : Float64.(v))

"""
    FormationLogK

The log K of formation of a species of a database of reactions from its master
species, as the column `formation` of its table of substances holds it: a
combination `Σ cₖ log Kₖ(T)` of the constants of the database, each in the form
its format gives ([`AbstractLogK`](@ref)).
"""
struct FormationLogK
    terms::Vector{Tuple{Float64, AbstractLogK}}
end

function _log10K(f::FormationLogK, T)
    L = dL = d2L = zero(T) * 0.0
    for (c, k) in f.terms
        l, d, d2 = _log10K(k, T)
        L += c * l
        dL += c * d
        d2L += c * d2
    end
    return L, dL, d2L
end

# A combination of constants, by index, as one `FormationLogK`: the constants a
# PHREEQC constant adds (`-add_logk`) carried with their coefficients, depth
# first.
function _formation(logks, f::Dict{Int, Float64})
    terms = Tuple{Float64, AbstractLogK}[]
    for (j, c) in sort!(collect(f); by = first)
        _add_term!(terms, logks, j, c)
    end
    return FormationLogK(terms)
end
function _add_term!(terms, logks, j, c)
    k = logks[j]
    push!(terms, (c, k))
    for (i, ci) in _added(k)
        _add_term!(terms, logks, i, c * ci)
    end
    return terms
end

# ── The database ─────────────────────────────────────────────────────────────

# A species or a phase of a database of reactions, resolved: its composition and
# charge by the balance of its reaction, and its log K of formation from the
# master species as a combination `Σ c_k log K_k(T)` of the constants of the
# database (`formation`, index => coefficient). `gamma` holds the `(å, b)` of
# PHREEQC's `-gamma`, `ion_size` the size å of PHREEQC's `-llnl_gamma` or of the
# ion size of the other formats, `co2_gamma` whether a neutral species takes the
# activity coefficient of CO2 (`-co2_llnl_gamma`), `volume` the molar volume of
# a phase (cm³/mol), `critical` the `(T_c, P_c, ω)` of a gas (K, Pa).
struct _ReactionEntry
    name::String
    kind::Symbol                          # :solution or :phase
    atoms::Dict{Symbol, Float64}
    charge::Float64
    formation::Dict{Int, Float64}
    gamma::Union{Nothing, Tuple{Float64, Float64}}
    ion_size::Union{Nothing, Float64}
    co2_gamma::Bool
    volume::Union{Nothing, Float64}
    critical::Union{Nothing, NTuple{3, Float64}}
    line::Int
    formula::String                       # the name of a solute, the formula of a phase
end

# A database of reactions as read, before it is laid out in tables
# (`_reaction_tables`): PHREEQC, The Geochemist's Workbench or EQ3/6.
struct _ReactionData
    name::String
    path::String
    source::String
    gauge::String
    format::Symbol
    masters::Dict{String, String}         # element, or element(valence) => master species
    species::Vector{_ReactionEntry}
    phases::Vector{_ReactionEntry}
    logks::Vector{AbstractLogK}
    reactions::Vector{NamedTuple}         # (symbol, kind, equation, logk, line), as written
    parameters::Any                       # the activity model's, as the format gives them
    keywords::Set{String}                 # the blocks of a PHREEQC file
    notes::Vector{String}
end

# ── Resolution ───────────────────────────────────────────────────────────────

# The composition, the charge and the formation of a species, as resolved.
const _Resolved = Tuple{Dict{Symbol, Float64}, Float64, Dict{Int, Float64}}

# The species of a database resolved on demand against those they rest on:
# `define(name, resolver)` resolves one species, calling the resolver for each
# species of its reaction, and returns `nothing` when one of them is not defined.
# Each species is resolved once, and a species met again while it is being
# resolved (a cycle of reactions) resolves to `nothing`.
struct _Resolver{F}
    define::F
    resolved::Dict{String, Union{Nothing, _Resolved}}
    visiting::Set{String}
end
_Resolver(define) = _Resolver(define, Dict{String, Union{Nothing, _Resolved}}(), Set{String}())

function (r::_Resolver)(name::AbstractString)
    haskey(r.resolved, name) && return r.resolved[name]
    name in r.visiting && return nothing
    push!(r.visiting, name)
    out = r.define(name, r)
    delete!(r.visiting, name)
    return (r.resolved[name] = out)
end

# The defined species of a reaction `Σ νᵢ Aᵢ = 0` (products positive): its
# composition and charge by the balances, and its log K of formation from the
# masters, `L_d = (log K − Σ_{i≠d} νᵢ Lᵢ) / ν_d`, from `Σ νᵢ Gᵢ = −RT ln 10 log K`
# and `Gᵢ = −RT ln 10 Lᵢ`.
function _resolve_defined(defined, terms, k, logks, resolve)
    ν_d = 0.0
    atoms = Dict{Symbol, Float64}()
    z = 0.0
    formation = Dict{Int, Float64}(k => 1.0)
    for (ν, s) in terms
        if s == defined
            ν_d += ν
            continue
        end
        r = resolve(s)
        r === nothing && return nothing
        a, zi, f = r
        for (e, n) in a
            atoms[e] = get(atoms, e, 0.0) - ν * n
        end
        z -= ν * zi
        for (j, c) in f
            formation[j] = get(formation, j, 0.0) - ν * c
        end
    end
    iszero(ν_d) && return nothing
    # The balances are sums of decimal coefficients: a whole or simple number
    # comes back within rounding of itself, and is put back on it.
    clean(x) = abs(x - round(x)) < 1.0e-9 ? round(x) : round(x; digits = 12)
    for e in keys(atoms)
        atoms[e] = clean(atoms[e] / ν_d)
    end
    filter!(p -> !iszero(p.second), atoms)
    for j in keys(formation)
        formation[j] /= ν_d
    end
    filter!(p -> !iszero(p.second), formation)
    return atoms, clean(z / ν_d), formation
end

# Whether a resolved species is made of something other than chemical elements
# (a pseudo-element of PHREEQC, an element a GWB or EQ3/6 file does not tie to
# one), and the note saying it is not read. ChemistryLab conserves chemical
# elements and solves the redox equilibrium, so these are not read.
function _not_elements(atoms, file, line, name, notes)
    unknown = [e for e in keys(atoms) if !haskey(elements.bysymbol, e)]
    isempty(unknown) && return false
    push!(notes, "$file:$line: $name is made of $(join(string.(unknown), ", ")), not of chemical elements; not read")
    return true
end

# ── Tables ───────────────────────────────────────────────────────────────────

# A database of reactions laid out in tables, as `read_thermofun_database` lays
# out a database of formation properties: the master species, the substances
# (one row per species or phase, its log K of formation in the column
# `formation`) and the reactions as written. What concerns the database as a
# whole is metadata of the table of substances.
function _reaction_tables(d::_ReactionData)
    rows = NamedTuple[]
    cm3 = ustrip(us"m^3", 1.0u"cm^3")
    for sp in vcat(d.species, d.phases)
        lowercase(sp.name) == "e-" && continue
        sym = sp.kind === :phase ? sp.name : _solute_symbol(sp.name, sp.charge)
        is_gas = sp.kind === :phase && (sp.critical !== nothing || endswith(lowercase(sp.name), "(g)"))
        state = sp.kind === :solution ? AS_AQUEOUS : (is_gas ? AS_GAS : AS_CRYSTAL)
        push!(
            rows, (;
                symbol = sym, name = sp.name, formula = sp.formula, aggregate_state = state, class = _class_of(state, sym),
                charge = sp.charge, atoms = sp.atoms, formation = _formation(d.logks, sp.formation),
                gamma = something(sp.gamma, missing), ion_size = something(sp.ion_size, missing),
                co2_gamma = sp.co2_gamma,
                molar_volume = sp.volume === nothing ? missing : sp.volume * cm3,
                T_c = sp.critical === nothing ? missing : sp.critical[1],
                P_c = sp.critical === nothing ? missing : sp.critical[2],
                omega = sp.critical === nothing ? missing : sp.critical[3],
                gauge = d.gauge, line = sp.line,
            ),
        )
    end
    df_substances = DataFrame(Tables.dictrowtable(rows))
    df_elements = DataFrame(element = collect(keys(d.masters)), master = collect(values(d.masters)))
    sort!(df_elements, :element)
    df_reactions = isempty(d.reactions) ?
        DataFrame(symbol = String[], kind = Symbol[], equation = String[], logk = AbstractLogK[], line = Int[]) :
        DataFrame(Tables.dictrowtable(d.reactions))
    for (key, value) in (
            "format" => String(d.format), "path" => d.path, "source" => d.source, "gauge" => d.gauge,
            "notes" => d.notes, "parameters" => d.parameters, "keywords" => sort!(collect(d.keywords)),
        )
        metadata!(df_substances, key, value; style = :note)
    end
    return df_elements, df_substances, df_reactions
end

# One species of a table of a database of reactions. Its standard Gibbs energy is
# `−RT ln 10 · L(T)`, `L` its log K of formation from the master species, which
# are at zero; its enthalpy, entropy and heat capacity follow from the
# temperature dependence of `L`.
function _table_species(row)
    s = _row_species(row)
    f = row.formation
    refs = _NF_DEFAULT_REFS
    s[:ΔₐG⁰] = NumericFunc((T, P) -> _logK_gibbs(_log10K(f, T), T), (:T, :P), refs, u"J/mol")
    s[:ΔₐH⁰] = NumericFunc((T, P) -> _logK_enthalpy(_log10K(f, T), T), (:T, :P), refs, u"J/mol")
    s[:S⁰] = NumericFunc((T, P) -> _logK_entropy(_log10K(f, T), T), (:T, :P), refs, u"J/(mol*K)")
    s[:Cp⁰] = NumericFunc((T, P) -> _logK_heat_capacity(_log10K(f, T), T), (:T, :P), refs, u"J/(mol*K)")
    if row.class == SC_AQSOLVENT
        M = ustrip(us"kg/mol", s[:M])
        s[:V⁰] = NumericFunc((T, P) -> M / _hgk_density(T, P), (:T, :P), refs, u"m^3/mol")
    elseif !ismissing(row.molar_volume) && row.aggregate_state == AS_CRYSTAL
        v = row.molar_volume
        s[:V⁰] = NumericFunc((T, P) -> v + zero(T), (:T, :P), refs, u"m^3/mol")
    end
    s[:gauge] = row.gauge
    ismissing(row.T_c) || (s[:critical_constants] = [row.T_c, row.P_c, row.omega])
    return s
end

# The species of a table of a database of reactions, those named (by the
# database's name or ChemistryLab's symbol) or all.
_build_reaction_species(df::AbstractDataFrame, names = nothing) =
    Species[_table_species(r) for r in eachrow(df) if _listed(r, names)]
