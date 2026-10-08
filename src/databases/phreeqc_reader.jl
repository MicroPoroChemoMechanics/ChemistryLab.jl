# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using SHA: sha256

# ── Reading a PHREEQC database ───────────────────────────────────────────────
#
# A PHREEQC database is a list of keyword data blocks (Parkhurst and Appelo 2013,
# "Description of Data Input"). The species it defines carry no energy of their
# own: each is defined by a reaction and the log K of that reaction, as a
# function of temperature, relative to the master species of
# SOLUTION_MASTER_SPECIES, which are defined by identity reactions of log K 0.
#
# The energies ChemistryLab minimizes follow by setting the master species to
# zero at every temperature: the standard Gibbs energy of a species is then
# -RT ln K of its formation from the masters. That is a gauge. Adding to every
# energy a combination of its element amounts and its charge leaves the
# equilibrium unchanged under the balances, and the masters, one per element plus
# H+, H2O and e-, span the element and charge space; but the energies mean
# something only beside the other species of the same database, and a system
# that mixes them with the energies of another database is refused
# (`theory/formation_quantities.md`).
#
# The composition of a species is taken from its reaction, by the balance of
# each element and of charge, and not from its name: a name is a label to
# PHREEQC (`CaSO4:2H2O`, `Fe(OH)3(a)`), and the reaction is what its energy
# belongs to.

# The keyword data blocks of PHREEQC 3, as the manual lists them, and the later
# MEAN_GAMMAS. A block ends at the next keyword; a database ends at its first
# END, after which PHREEQC reads nothing.
const _PHREEQC_KEYWORDS = Set(
    [
        "ADVECTION", "CALCULATE_VALUES", "COPY", "DATABASE", "DELETE", "DUMP", "END",
        "EQUILIBRIUM_PHASES", "EXCHANGE", "EXCHANGE_MASTER_SPECIES", "EXCHANGE_SPECIES",
        "GAS_PHASE", "INCLUDE\$", "INCREMENTAL_REACTIONS", "INVERSE_MODELING", "ISOTOPES",
        "ISOTOPE_ALPHAS", "ISOTOPE_RATIOS", "KINETICS", "KNOBS", "LLNL_AQUEOUS_MODEL_PARAMETERS",
        "MEAN_GAMMAS", "MIX", "NAMED_EXPRESSIONS", "PHASES", "PITZER", "PRINT", "RATES",
        "REACTION", "REACTION_PRESSURE", "REACTION_TEMPERATURE", "RUN_CELLS", "SAVE",
        "SELECTED_OUTPUT", "SIT", "SOLID_SOLUTIONS", "SOLUTION", "SOLUTION_MASTER_SPECIES",
        "SOLUTION_SPECIES", "SOLUTION_SPREAD", "SURFACE", "SURFACE_MASTER_SPECIES",
        "SURFACE_SPECIES", "TITLE", "TRANSPORT", "USE", "USER_GRAPH", "USER_PRINT", "USER_PUNCH",
    ]
)

function _is_phreeqc_keyword(token::AbstractString)
    u = uppercase(token)
    return u in _PHREEQC_KEYWORDS || endswith(u, "_RAW") || endswith(u, "_MODIFY")
end

"""
    PhreeqcBlock

One keyword data block of a PHREEQC file: its keyword, the line it starts on, and
its logical lines as `(line number, text)`, comments removed and the lines a
semicolon joins split apart.
"""
struct PhreeqcBlock
    keyword::String
    line::Int
    lines::Vector{Tuple{Int, String}}
end

"""
    phreeqc_blocks(path) -> Vector{PhreeqcBlock}

The keyword data blocks of a PHREEQC file, in order, up to its first `END`.
"""
function phreeqc_blocks(path::AbstractString)
    blocks = PhreeqcBlock[]
    current = nothing
    for (n, raw) in enumerate(eachline(path))
        text = strip(first(split(raw, '#'; limit = 2)))
        isempty(text) && continue
        token = first(split(text))
        if _is_phreeqc_keyword(token)
            uppercase(token) == "END" && break
            current = PhreeqcBlock(uppercase(token), n, Tuple{Int, String}[])
            push!(blocks, current)
            continue
        end
        current === nothing && continue
        # A semicolon puts several logical lines on one: `-T_c 154.6; -P_c 49.8`.
        # BASIC code (RATES) keeps its own semicolons.
        pieces = current.keyword == "RATES" ? [text] : split(text, ';')
        for piece in pieces
            p = strip(piece)
            isempty(p) || push!(current.lines, (n, String(p)))
        end
    end
    return blocks
end

# ── Equations ────────────────────────────────────────────────────────────────

# A species name as PHREEQC means it: a run of signs is a charge (`SO4--` is
# `SO4-2`), and a charge of one is written with its sign alone.
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

# The charge a PHREEQC name carries.
function _phreeqc_name_charge(name::AbstractString)
    s = _phreeqc_name(name)
    m = match(r"([+-])(\d+(?:\.\d+)?)$", s)
    m === nothing || return (m.captures[1] == "+" ? 1.0 : -1.0) * parse(Float64, m.captures[2])
    endswith(s, "+") && return 1.0
    endswith(s, "-") && return -1.0
    return 0.0
end

const _COEFFICIENT = r"^((?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)(.*)$"

# One side of an equation, `[(coefficient, species)]`. Terms are separated by
# `+`, or by `-` for a term subtracted (`1.000Am+3 - 0.500H2O`, in sit.dat); a
# sign may stand alone or carry the coefficient that follows it (`HS- +2.0000
# O2`, in llnl.dat); a coefficient stands alone or is written against its
# species (`2H2O`). A species name never starts with a sign, and one ending with
# a sign (`OH-`) is a single token.
function _phreeqc_side(text::AbstractString)
    out = Tuple{Float64, String}[]
    coefficient = nothing
    sign = 1.0
    for token in split(text)
        t = String(token)
        t == "+" && (sign = 1.0; continue)
        t == "-" && (sign = -1.0; continue)
        if (startswith(t, "+") || startswith(t, "-")) && length(t) > 1 && (isdigit(t[2]) || t[2] == '.')
            sign = t[1] == '-' ? -1.0 : 1.0
            t = t[2:end]
        end
        v = tryparse(Float64, t)
        if v !== nothing
            coefficient = v
            continue
        end
        m = match(_COEFFICIENT, t)
        c, name = (m !== nothing && !isempty(m.captures[2])) ?
            (something(coefficient, 1.0) * parse(Float64, m.captures[1]), m.captures[2]) :
            (something(coefficient, 1.0), t)
        push!(out, (sign * c, _phreeqc_name(name)))
        coefficient, sign = nothing, 1.0
    end
    return out
end

"""
    phreeqc_equation(text) -> (lhs, rhs)

The two sides of a PHREEQC equation, each a vector of `(coefficient, species)`.
"""
function phreeqc_equation(text::AbstractString)
    depth, cut = 0, 0
    for (i, c) in pairs(text)
        c == '[' && (depth += 1)
        c == ']' && (depth -= 1)
        if c == '=' && depth == 0
            cut == 0 || throw(ArgumentError("more than one equal sign: $text"))
            cut = i
        end
    end
    cut == 0 && throw(ArgumentError("no equal sign: $text"))
    return _phreeqc_side(text[1:prevind(text, cut)]), _phreeqc_side(text[nextind(text, cut):end])
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

"""
    PhreeqcLogK

The log K of one reaction as a PHREEQC database states it: at 25 °C, with an
enthalpy of reaction for the van 't Hoff equation, or an analytical expression
of the temperature, which takes precedence over both (Parkhurst and Appelo
2013, `SOLUTION_SPECIES` and `PHASES`). `add` holds the named expressions the
constant adds (`-add_logk`), as `(index, coefficient)`.
"""
struct PhreeqcLogK <: AbstractLogK
    log_k::Float64
    delta_h::Float64                      # kJ/mol
    analytic::Union{Nothing, NTuple{6, Float64}}
    add::Vector{Tuple{Int, Float64}}
end

const _T25 = T_STANDARD
_added(k::PhreeqcLogK) = k.add

# log10 K at T (kelvin) and its first two temperature derivatives, from the
# reaction's own data, without the named expressions it adds.
function _log10K(k::PhreeqcLogK, T)
    if k.analytic !== nothing
        A1, A2, A3, A4, A5, A6 = k.analytic
        L = A1 + A2 * T + A3 / T + A4 * log10(T) + A5 / T^2 + A6 * T^2
        dL = A2 - A3 / T^2 + A4 / (T * log(10)) - 2A5 / T^3 + 2A6 * T
        d2L = 2A3 / T^3 - A4 / (T^2 * log(10)) + 6A5 / T^4 + 2A6
        return L, dL, d2L
    end
    c = 1000 * k.delta_h / (R_GAS * log(10))
    return k.log_k - c * (1 / T - 1 / _T25), c / T^2, -2c / T^3
end

# ── The database ─────────────────────────────────────────────────────────────

# A species or a phase of a database of reactions, resolved: its composition and
# charge by the balance of its reaction, and its log K of formation from the
# master species as a combination `Σ c_k log K_k(T)` of the constants of the
# database (`formation`, index => coefficient). `gamma` holds the `(å, b)` of
# PHREEQC's `-gamma`, `ion_size` the size å of PHREEQC's `-llnl_gamma` or of the
# ion size of the other formats, `co2_gamma` whether a neutral species takes the
# activity coefficient of CO2 (`-co2_llnl_gamma`), `volume` the molar volume of
# a phase (cm³/mol), `critical` the `(T_c, P_c, ω)` of a gas (K, atm).
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

# The options of a species or a phase, in every spelling the manual lists:
# case does not matter, the leading dash is optional, and an abbreviation is the
# one the manual gives in brackets.
const _PHREEQC_SPECIES_OPTIONS = Dict(
    :log_k => ("log_k", "logk", "-log_k", "-logk", "-l"),
    :delta_h => ("delta_h", "deltah", "-delta_h", "-deltah", "-d"),
    :analytic => (
        "-analytic", "analytic", "analytical_expression", "-analytical_expression",
        "-analytical", "analytical", "a_e", "ae", "-a_e", "-ae", "-a",
    ),
    :gamma => ("-gamma", "-g"),
    :llnl_gamma => ("-llnl_gamma", "llnl_gamma", "-ll"),
    :co2_llnl_gamma => ("-co2_llnl_gamma", "co2_llnl_gamma", "-co"),
    :vm => ("-vm", "vm"),
    :dw => ("-dw", "dw"),
    :millero => ("-millero", "millero", "-mi"),
    :mole_balance => ("-mole_balance", "mole_balance", "-mass_balance", "mass_balance", "-mo", "-mass"),
    :no_check => ("-no_check", "no_check", "-n"),
    :add_logk => ("-add_logk", "add_logk"),
    :add_constant => ("-add_constant", "add_constant"),
    :t_c => ("-t_c", "t_c"),
    :p_c => ("-p_c", "p_c"),
    :omega => ("-omega", "omega"),
    :cd_music => ("-cd_music", "cd_music"),
    :erm_ddl => ("-erm_ddl", "erm_ddl"),
    :activity_water => ("-activity_water", "activity_water"),
    :tracer_diffusion => ("-tracer_diffusion",),
    :dalfg => ("-dalfg",),
)
function _species_option(token::AbstractString)
    t = lowercase(token)
    for (opt, spellings) in _PHREEQC_SPECIES_OPTIONS
        t in spellings && return opt
    end
    return startswith(t, "-") ? :unknown : nothing
end

# kJ/mol from a `delta_h` and the unit that may follow it.
function _delta_h_kJ(value::Float64, unit)
    unit === nothing && return value
    u = lowercase(unit)
    startswith(u, "kcal") && return CALORIE * value
    startswith(u, "cal") && return CALORIE * value / 1000
    startswith(u, "kj") && return value
    startswith(u, "j") && return 1.0e-3 * value
    return value
end

# A record under construction: its equation, constant and options.
mutable struct _Record
    name::String
    equation::String
    line::Int
    log_k::Float64
    delta_h::Float64
    analytic::Union{Nothing, NTuple{6, Float64}}
    add::Vector{Tuple{String, Float64}}
    gamma::Union{Nothing, Tuple{Float64, Float64}}
    llnl_gamma::Union{Nothing, Float64}
    vm::Union{Nothing, Float64}
    critical::Vector{Union{Nothing, Float64}}
    mole_balance::Union{Nothing, String}
    no_check::Bool
    co2_gamma::Bool
    refused::Vector{String}
end
_Record(name, equation, line) = _Record(
    name, equation, line, 0.0, 0.0, nothing, Tuple{String, Float64}[], nothing, nothing,
    nothing, Union{Nothing, Float64}[nothing, nothing, nothing], nothing, false, false, String[],
)

function _numbers(parts)
    out = Float64[]
    for p in parts
        v = tryparse(Float64, p)
        v === nothing && break
        push!(out, v)
    end
    return out
end

# Applies one option line to `rec`; returns a note when the line is not used.
function _apply_option!(rec::_Record, opt::Symbol, parts, n, file)
    args = parts[2:end]
    v = _numbers(args)
    if opt === :log_k
        isempty(v) || (rec.log_k = v[1])
    elseif opt === :delta_h
        isempty(v) || (rec.delta_h = _delta_h_kJ(v[1], length(args) > 1 ? args[2] : nothing))
    elseif opt === :analytic
        isempty(v) || (rec.analytic = ntuple(i -> i <= length(v) ? v[i] : 0.0, 6))
    elseif opt === :gamma
        length(v) >= 2 || return "$file:$n: the `-gamma` of $(rec.name) does not parse as two numbers; not read"
        rec.gamma = (v[1], v[2])
    elseif opt === :llnl_gamma
        isempty(v) || (rec.llnl_gamma = v[1])
    elseif opt === :vm
        isempty(v) || (rec.vm = v[1])
    elseif opt === :t_c
        isempty(v) || (rec.critical[1] = v[1])
    elseif opt === :p_c
        isempty(v) || (rec.critical[2] = v[1])
    elseif opt === :omega
        isempty(v) || (rec.critical[3] = v[1])
    elseif opt === :add_logk
        # The name of the expression, then its coefficient (1 by default).
        isempty(args) && return "$file:$n: `-add_logk` names no expression; not read"
        c = length(args) >= 2 ? something(tryparse(Float64, args[2]), 1.0) : 1.0
        push!(rec.add, (String(args[1]), c))
    elseif opt === :add_constant
        isempty(v) || push!(rec.add, ("\$constant", v[1]))
    elseif opt === :mole_balance
        isempty(args) || (rec.mole_balance = String(args[1]))
    elseif opt === :no_check
        rec.no_check = true
    elseif opt === :co2_llnl_gamma
        rec.co2_gamma = true
    elseif opt in (:dw, :millero, :activity_water, :tracer_diffusion, :dalfg)
        # Transport and density data, not thermodynamics.
        return nothing
    elseif opt in (:cd_music, :erm_ddl)
        push!(rec.refused, "$file:$n: `$(parts[1])` (a surface electrostatic option) is not read")
    else
        return "$file:$n: option `$(parts[1])` of $(rec.name) is not one the PHREEQC manual documents; not read"
    end
    return nothing
end

# The records of SOLUTION_SPECIES (or another *_SPECIES block): each one starts
# with an equation, and the species it defines is the first on the right.
function _species_records(block::PhreeqcBlock, file, notes)
    recs = _Record[]
    for (n, text) in block.lines
        parts = split(text)
        opt = _species_option(parts[1])
        if opt in (nothing, :unknown) && occursin('=', text)
            _, rhs = phreeqc_equation(text)
            isempty(rhs) && (push!(notes, "$file:$n: an equation with nothing to its right; not read"); continue)
            push!(recs, _Record(rhs[1][2], text, n))
        elseif !isempty(recs) && opt !== nothing
            note = _apply_option!(recs[end], opt, parts, n, file)
            note === nothing || push!(notes, note)
        else
            push!(notes, "$file:$n: `$text` belongs to no species; not read")
        end
    end
    return recs
end

# The records of PHASES: a name alone on its line, then its dissolution
# equation, whose first term on the left is the phase.
function _phase_records(block::PhreeqcBlock, file, notes)
    recs = _Record[]
    pending = nothing
    for (n, text) in block.lines
        parts = split(text)
        opt = _species_option(parts[1])
        if opt in (nothing, :unknown) && occursin('=', text)
            if pending === nothing
                push!(notes, "$file:$n: an equation with no phase name before it; not read")
                continue
            end
            push!(recs, _Record(pending[2], text, pending[1]))
            pending = nothing
        elseif opt === nothing
            pending = (n, String(parts[1]))
        elseif !isempty(recs)
            note = _apply_option!(recs[end], opt, parts, n, file)
            note === nothing || push!(notes, note)
        end
    end
    return recs
end

# NAMED_EXPRESSIONS: a name, then the options of a log K.
function _named_records(block::PhreeqcBlock, file, notes)
    recs = _Record[]
    for (n, text) in block.lines
        parts = split(text)
        opt = _species_option(parts[1])
        if opt === nothing
            push!(recs, _Record(String(parts[1]), "", n))
        elseif !isempty(recs)
            note = _apply_option!(recs[end], opt, parts, n, file)
            note === nothing || push!(notes, note)
        end
    end
    return recs
end

# LLNL_AQUEOUS_MODEL_PARAMETERS: the temperatures and, at each, the A and B of
# the Debye-Hückel equation and the B-dot, and the coefficients of the activity
# coefficient of CO2.
function _llnl_parameters(block::PhreeqcBlock)
    current, out = nothing, Dict{Symbol, Vector{Float64}}()
    for (_, text) in block.lines
        parts = split(text)
        # An option, and not a negative number (`-1.0312` of `-co2_coefs`).
        if startswith(parts[1], "-") && tryparse(Float64, parts[1]) === nothing
            current = Symbol(lowercase(lstrip(parts[1], '-')))
            out[current] = _numbers(parts[2:end])
        elseif current !== nothing
            append!(out[current], _numbers(parts))
        end
    end
    return (; (k => v for (k, v) in out)...)
end

"""
    read_phreeqc_database(path) -> (df_elements, df_substances, df_reactions)

Read a PHREEQC database into three tables, as [`read_thermofun_database`](@ref)
reads a ThermoFun file: its master species (`df_elements`: `element`,
`master`), its species and phases (`df_substances`), and its reactions as
written, each with its constant (`df_reactions`: `symbol`, `kind`, `equation`,
`logk`, `line`).

Every species is resolved against the master species, its composition by the
balance of its reaction and its energy by the log K of its formation from them.
`df_substances` has one row per species or phase: `symbol` (ChemistryLab's),
`name` (PHREEQC's), `formula`, `aggregate_state`, `class`, `charge`, `atoms`,
`formation` (a [`FormationLogK`](@ref)), the activity parameters `gamma` (the
`(å, b)` of `-gamma`), `ion_size` (`-llnl_gamma`) and `co2_gamma`, the
`molar_volume` of a phase (m³/mol), the critical constants `T_c` (K), `P_c` (Pa)
and `omega` of a gas, `gauge` and `line`; [`build_species`](@ref) builds from it,
and [`database_activity_model`](@ref) finds in it the activity model of the
database. Its metadata (`metadata(df_substances, "notes")`) hold `format`,
`path`, `source`, `gauge`, the `parameters` of the LLNL model, the `keywords` of
the blocks, and the `notes`: what is read and not used, with its line —
transport data (`-dw`, `-Millero`), an option the manual does not document, a
species whose equation does not balance (`-mole_balance`, `-no_check`), and the
blocks that are not thermodynamic data (`RATES`).

Any PHREEQC database is read the same way. The ones PHREEQC distributes are
obtained by `datapath`: `datapath("phreeqc.dat")`, `"llnl.dat"`,
`"minteq.v4.dat"`, `"wateq4f.dat"`, `"sit.dat"`, `"pitzer.dat"`.
"""
function read_phreeqc_database(path::AbstractString)
    file = basename(path)
    digest = bytes2hex(open(sha256, path))
    blocks = phreeqc_blocks(path)
    notes = String[]

    masters = Dict{String, String}()
    pseudo = Dict{String, String}()       # master species => pseudo-element it carries
    srecs, precs, nrecs = _Record[], _Record[], _Record[]
    llnl = nothing
    for b in blocks
        if b.keyword == "SOLUTION_MASTER_SPECIES"
            for (_, text) in b.lines
                parts = split(text)
                length(parts) >= 2 || continue
                element, master = String(parts[1]), _phreeqc_name(parts[2])
                masters[element] = master
                _is_pseudo_element(element, parts) && (pseudo[master] = element)
            end
        elseif b.keyword == "SOLUTION_SPECIES"
            append!(srecs, _species_records(b, file, notes))
        elseif b.keyword == "PHASES"
            append!(precs, _phase_records(b, file, notes))
        elseif b.keyword == "NAMED_EXPRESSIONS"
            append!(nrecs, _named_records(b, file, notes))
        elseif b.keyword == "LLNL_AQUEOUS_MODEL_PARAMETERS"
            llnl = _llnl_parameters(b)
        elseif b.keyword in (
                "PITZER", "SIT", "EXCHANGE_MASTER_SPECIES", "EXCHANGE_SPECIES",
                "SURFACE_MASTER_SPECIES", "SURFACE_SPECIES",
            )
            # Read by the readers of the activity and sorption models.
            continue
        elseif b.keyword in ("RATES", "ISOTOPES", "ISOTOPE_RATIOS", "ISOTOPE_ALPHAS", "CALCULATE_VALUES", "MEAN_GAMMAS")
            push!(notes, "$file:$(b.line): $(b.keyword) is not thermodynamic data; not read")
        else
            push!(notes, "$file:$(b.line): $(b.keyword) is an input block, not database data; not read")
        end
    end

    # A species or a phase defined again replaces the earlier definition.
    srecs = _last_definitions(srecs, file, notes)
    precs = _last_definitions(precs, file, notes)

    # The constants: the named expressions first, then one per record.
    logks = AbstractLogK[]
    named = Dict{String, Int}()
    function logk_of(rec)
        add = Tuple{Int, Float64}[]
        for (nm, c) in rec.add
            if nm == "\$constant"
                push!(logks, PhreeqcLogK(c, 0.0, nothing, Tuple{Int, Float64}[]))
                push!(add, (length(logks), 1.0))
            elseif haskey(named, nm)
                push!(add, (named[nm], c))
            else
                push!(rec.refused, "$file:$(rec.line): `-add_logk $nm` names no expression of NAMED_EXPRESSIONS")
            end
        end
        push!(logks, PhreeqcLogK(rec.log_k, rec.delta_h, rec.analytic, add))
        return length(logks)
    end
    for rec in nrecs
        named[rec.name] = logk_of(rec)
    end

    # The species, by name, each with its record.
    by_name = Dict{String, Tuple{_Record, Int}}()
    for rec in srecs
        by_name[rec.name] = (rec, logk_of(rec))
    end

    # Composition and charge of a master species, from its formula.
    function master_composition(name)
        name == "e-" && return Dict{Symbol, Float64}(), -1.0
        # A pseudo-element: kept as itself, which `_pseudo_element` refuses.
        haskey(pseudo, name) && return Dict{Symbol, Float64}(Symbol("pseudo-element ", pseudo[name]) => 1.0), _phreeqc_name_charge(name)
        atoms = Dict{Symbol, Float64}(Symbol(k) => Float64(v) for (k, v) in parse_formula(name))
        return atoms, _phreeqc_name_charge(name)
    end

    resolved = Dict{String, Union{Nothing, Tuple{Dict{Symbol, Float64}, Float64, Dict{Int, Float64}}}}()
    visiting = Set{String}()
    function resolve(name)
        haskey(resolved, name) && return resolved[name]
        if !haskey(by_name, name)
            name == "e-" && return (resolved[name] = (Dict{Symbol, Float64}(), -1.0, Dict{Int, Float64}()))
            return (resolved[name] = nothing)
        end
        name in visiting && return nothing
        push!(visiting, name)
        rec, k = by_name[name]
        lhs, rhs = phreeqc_equation(rec.equation)
        terms = vcat([(-c, s) for (c, s) in lhs], [(c, s) for (c, s) in rhs])
        out = if length(lhs) == 1 && length(rhs) == 1 && lhs[1][2] == rhs[1][2]
            # An identity reaction: a master species, at zero by the gauge.
            atoms, z = master_composition(name)
            (atoms, z, Dict{Int, Float64}())
        else
            _resolve_defined(name, terms, k, logks, resolve)
        end
        delete!(visiting, name)
        return (resolved[name] = out)
    end

    species = _ReactionEntry[]
    for rec in srecs
        r = resolve(rec.name)
        if r === nothing
            push!(notes, "$file:$(rec.line): $(rec.name) rests on a species the database does not define; not read")
            continue
        end
        _refusal(rec, r, file, notes) && continue
        _pseudo_element(rec, r, file, notes) && continue
        atoms, z, formation = r
        push!(species, _ReactionEntry(rec.name, :solution, atoms, z, formation, rec.gamma, rec.llnl_gamma, rec.co2_gamma, nothing, nothing, rec.line, rec.name))
    end

    phases = _ReactionEntry[]
    phase_logk = Dict{String, Int}()
    for rec in precs
        lhs, rhs = phreeqc_equation(rec.equation)
        isempty(lhs) && (push!(notes, "$file:$(rec.line): phase $(rec.name) has no formula on the left; not read"); continue)
        k = logk_of(rec)
        phase_logk[rec.name] = k
        # The phase is the first term on the left; everything else is aqueous.
        terms = vcat([(-lhs[1][1], "\$phase")], [(-c, s) for (c, s) in lhs[2:end]], [(c, s) for (c, s) in rhs])
        r = _resolve_defined("\$phase", terms, k, logks, resolve)
        if r === nothing
            push!(notes, "$file:$(rec.line): phase $(rec.name) rests on a species the database does not define; not read")
            continue
        end
        _refusal(rec, r, file, notes; formula = lhs[1][2]) && continue
        _pseudo_element(rec, r, file, notes) && continue
        atoms, z, formation = r
        critical = all(x -> x !== nothing, rec.critical) ? (rec.critical[1], rec.critical[2], rec.critical[3]) : nothing
        push!(phases, _ReactionEntry(rec.name, :phase, atoms, z, formation, nothing, nothing, false, rec.vm, critical, rec.line, lhs[1][2]))
    end

    gauge = "PHREEQC master species of $file ($(digest[1:12]))"
    # The reactions as written, each with its constant.
    reactions = NamedTuple[]
    for rec in nrecs
        push!(reactions, (; symbol = rec.name, kind = :named, equation = "", logk = logks[named[rec.name]], line = rec.line))
    end
    for rec in srecs
        haskey(by_name, rec.name) &&
            push!(reactions, (; symbol = rec.name, kind = :solution, equation = rec.equation, logk = logks[by_name[rec.name][2]], line = rec.line))
    end
    for rec in precs
        haskey(phase_logk, rec.name) &&
            push!(reactions, (; symbol = rec.name, kind = :phase, equation = rec.equation, logk = logks[phase_logk[rec.name]], line = rec.line))
    end
    data = _ReactionData(
        file, String(path), "$file sha256 $(digest[1:12])", gauge, :phreeqc, masters, species, phases, logks,
        reactions, llnl, Set(b.keyword for b in blocks), notes,
    )
    return _reaction_tables(data)
end

function _last_definitions(recs, file, notes)
    last = Dict{String, Int}()
    for (k, rec) in enumerate(recs)
        if haskey(last, rec.name)
            push!(notes, "$file:$(rec.line): $(rec.name) is defined again; the definition at line $(recs[last[rec.name]].line) is replaced")
        end
        last[rec.name] = k
    end
    return [rec for (k, rec) in enumerate(recs) if last[rec.name] == k]
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

# Whether a record is refused, and the note saying why: a species whose
# equation does not balance by design (`-mole_balance`), one the database
# exempts from the check (`-no_check`) whose name says another composition, and
# one carrying an option that changes its thermodynamics and is not read.
function _refusal(rec, r, file, notes; formula = rec.name)
    if !isempty(rec.refused)
        append!(notes, rec.refused)
        push!(notes, "$file:$(rec.line): $(rec.name) not read")
        return true
    end
    if rec.mole_balance !== nothing
        push!(notes, "$file:$(rec.line): $(rec.name) is given another composition than its equation's (`-mole_balance $(rec.mole_balance)`); not read")
        return true
    end
    if rec.no_check
        named = try
            Dict{Symbol, Float64}(Symbol(k) => Float64(v) for (k, v) in parse_formula(formula))
        catch
            nothing
        end
        atoms = r[1]
        if named === nothing || keys(named) != keys(atoms) || any(abs(named[e] - atoms[e]) > 1.0e-9 for e in keys(atoms))
            push!(notes, "$file:$(rec.line): the equation of $(rec.name) does not balance (`-no_check`); not read")
            return true
        end
    end
    return false
end

# Whether a line of SOLUTION_MASTER_SPECIES defines a pseudo-element: an element
# a database defines to keep a dissolved gas out of redox equilibrium (`Mtg`,
# `Ntg`, `Hdg`, `Oxg`, `Sg` in phreeqc.dat). Its name is not that of a chemical
# element, or it is one (`Sg` is seaborgium) whose gram formula weight the line
# states 3 % or more away from the element's (34.08 for `Sg`, the H2S it stands
# for; a radioactive element's own mass depends on the isotope taken, 147 for Pm
# in llnl.dat against 145). Alkalinity and the electron are not elements either, and are not
# species.
function _is_pseudo_element(element::AbstractString, parts)
    occursin('(', element) && return false            # a valence state of an element
    element in ("E", "Alkalinity") && return false
    startswith(element, "[") && return true           # an isotope, `[18O]`
    sym = Symbol(element)
    haskey(elements.bysymbol, sym) || return true
    gfw = length(parts) >= 5 ? tryparse(Float64, parts[5]) : nothing
    gfw === nothing && return false
    mass = Float64(elements.bysymbol[sym].atomic_mass.val)
    return abs(gfw - mass) > 0.03 * mass
end

# A species made of a pseudo-element. ChemistryLab conserves chemical elements
# and solves the redox equilibrium, so these are not read.
function _pseudo_element(rec, r, file, notes)
    unknown = [e for e in keys(r[1]) if !haskey(elements.bysymbol, e)]
    isempty(unknown) && return false
    push!(notes, "$file:$(rec.line): $(rec.name) is made of $(join(string.(unknown), ", ")), not of chemical elements; not read")
    return true
end

# ── Tables ───────────────────────────────────────────────────────────────────

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
# PHREEQC constant adds (`-add_logk`) carried with their coefficients.
function _formation(logks, f::Dict{Int, Float64})
    terms = Tuple{Float64, AbstractLogK}[]
    function add!(j, c)
        k = logks[j]
        push!(terms, (c, k))
        for (i, ci) in _added(k)
            add!(i, c * ci)
        end
        return nothing
    end
    for (j, c) in sort!(collect(f); by = first)
        add!(j, c)
    end
    return FormationLogK(terms)
end

# The ChemistryLab symbol of a species of a database of reactions: `@` marks a
# neutral solute, as the other readers of this package write it, and replaces the
# `(aq)` with which GWB and EQ3/6 name one; water is the solvent.
function _phreeqc_solute_symbol(name, charge = _phreeqc_name_charge(name))
    base = replace(name, r"\(aq\)$"i => "")
    lowercase(base) == "h2o" && return "H2O@"
    return iszero(charge) ? base * "@" : name
end

# A database of reactions laid out in tables, as `read_thermofun_database` lays
# out a database of formation properties: the master species, the substances
# (one row per species or phase, its log K of formation in the column
# `formation`) and the reactions as written. What concerns the database as a
# whole is metadata of the table of substances.
function _reaction_tables(d::_ReactionData)
    rows = NamedTuple[]
    cm3 = ustrip(us"m^3", 1.0u"cm^3")
    atm = ustrip(us"Pa", _DQConstants.atm)
    for sp in vcat(d.species, d.phases)
        lowercase(sp.name) == "e-" && continue
        sym = sp.kind === :phase ? sp.name : _phreeqc_solute_symbol(sp.name, sp.charge)
        is_gas = sp.kind === :phase && (sp.critical !== nothing || endswith(lowercase(sp.name), "(g)"))
        state = sp.kind === :solution ? AS_AQUEOUS : (is_gas ? AS_GAS : AS_CRYSTAL)
        cls = sp.kind === :solution ? (sym == "H2O@" ? SC_AQSOLVENT : SC_AQSOLUTE) : (is_gas ? SC_GASFLUID : SC_COMPONENT)
        push!(
            rows, (;
                symbol = sym, name = sp.name, formula = sp.formula, aggregate_state = state, class = cls,
                charge = sp.charge, atoms = sp.atoms, formation = _formation(d.logks, sp.formation),
                gamma = something(sp.gamma, missing), ion_size = something(sp.ion_size, missing),
                co2_gamma = sp.co2_gamma,
                molar_volume = sp.volume === nothing ? missing : sp.volume * cm3,
                T_c = sp.critical === nothing ? missing : sp.critical[1],
                P_c = sp.critical === nothing ? missing : sp.critical[2] * atm,
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
function _reaction_species(row)
    atoms = all(isinteger, values(row.atoms)) ? Dict{Symbol, Int}(e => Int(n) for (e, n) in row.atoms) :
        Dict{Symbol, Float64}(e => n for (e, n) in row.atoms)
    z = row.charge
    s = Species(
        atoms, isinteger(z) ? Int(z) : z;
        name = row.name, symbol = row.symbol, aggregate_state = row.aggregate_state, class = row.class,
    )
    f = row.formation
    k = R_GAS * log(10)
    G(T, P) = -k * T * first(_log10K(f, T))
    H(T, P) = (r = _log10K(f, T); k * T^2 * r[2])
    S(T, P) = (H(T, P) - G(T, P)) / T
    Cp(T, P) = (r = _log10K(f, T); k * (2T * r[2] + T^2 * r[3]))
    refs = _NF_DEFAULT_REFS
    s[:ΔₐG⁰] = NumericFunc(G, (:T, :P), refs, u"J/mol")
    s[:ΔₐH⁰] = NumericFunc(H, (:T, :P), refs, u"J/mol")
    s[:S⁰] = NumericFunc(S, (:T, :P), refs, u"J/(mol*K)")
    s[:Cp⁰] = NumericFunc(Cp, (:T, :P), refs, u"J/(mol*K)")
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
function _build_reaction_species(df::AbstractDataFrame, names = nothing)
    keep(r) = names === nothing || r.name in names || r.symbol in names
    return Species[_reaction_species(r) for r in eachrow(df) if keep(r)]
end

# ── Activity model ───────────────────────────────────────────────────────────

"""
    database_activity_model(df_substances) -> AbstractActivityModel

The activity model PHREEQC applies with the database whose table of substances
is `df_substances` ([`read_phreeqc_database`](@ref)), which its blocks decide:
Pitzer's with a `PITZER` block ([`PitzerActivityModel`](@ref), the parameters of
the block, a pair it leaves out taken as zero), the specific ion interaction
theory with a `SIT` block ([`SITActivityModel`](@ref)), the aqueous model of
Lawrence Livermore with `LLNL_AQUEOUS_MODEL_PARAMETERS`
([`LLNLActivityModel`](@ref)), and the WATEQ and Davies equations otherwise
([`TruesdellJonesActivityModel`](@ref)).

Each is built with the parameters the database gives its species, its
Debye–Hückel parameters following temperature, and the activity of water PHREEQC
computes: the osmotic coefficient of the model for Pitzer and SIT, and
`a_w = 1 − 0.017 Σᵢ mᵢ` for the other two [ParkhurstAppelo1999; Eq. 25](@cite).
The species are keyed by their ChemistryLab symbols, those of
[`build_species`](@ref).

# Example

```julia
_, substances, _ = read_phreeqc_database(datapath("llnl.dat"))
model = database_activity_model(substances)          # an LLNLActivityModel
```
"""
function database_activity_model(df::AbstractDataFrame)
    format = metadata(df, "format", "unstated")
    format == "phreeqc" || throw(
        ArgumentError(
            "the activity model of a database of format `$format` is not built from it: " *
                (format in ("gwb", "eq36") ? "its conventions for neutral species and water rest on formulas its format documentation does not give; " : "") *
                "choose an activity model (HKFActivityModel, LLNLActivityModel, …)",
        )
    )
    keywords = Set(metadata(df, "keywords", String[]))
    path = metadata(df, "path")
    solutes = [r for r in eachrow(df) if r.aggregate_state == AS_AQUEOUS]
    if "PITZER" in keywords
        return PitzerActivityModel(;
            parameters = build_pitzer_parameters(path; format = :phreeqc), temperature_dependent = true,
            missing_pairs = :zero,
        )
    elseif "SIT" in keywords
        return SITActivityModel(; parameters = build_sit_parameters(path), temperature_dependent = true, water = :osmotic)
    end
    gamma = Dict{String, Tuple{Float64, Float64}}(r.symbol => r.gamma for r in solutes if !ismissing(r.gamma))
    l = metadata(df, "parameters", nothing)
    if l !== nothing
        for key in (:temperatures, :dh_a, :dh_b, :bdot, :co2_coefs)
            haskey(l, key) || throw(ArgumentError("LLNL_AQUEOUS_MODEL_PARAMETERS gives no `-$key`"))
        end
        sizes = Dict{String, Float64}(r.symbol => r.ion_size for r in solutes if !ismissing(r.ion_size))
        co2 = [r.symbol for r in solutes if r.co2_gamma]
        return LLNLActivityModel(;
            temperatures = l.temperatures, A = l.dh_a, B = l.dh_b, Bdot = l.bdot, co2 = l.co2_coefs,
            sizes, co2_species = co2, gamma, water = :phreeqc,
        )
    end
    return TruesdellJonesActivityModel(; parameters = gamma, temperature_dependent = true, water = :phreeqc)
end

# ── Parameters and master species, for the readers that need only these ─────

"""
    phreeqc_gamma_parameters(path) -> Dict{String, Tuple{Float64, Float64}}

The WATEQ activity-coefficient parameters `(å, b)` that the `-gamma` option of a
PHREEQC database gives its aqueous species, by species symbol, for
[`TruesdellJonesActivityModel`](@ref).

Every `SOLUTION_SPECIES` block is read, the master species as well as the
others, without resolving the species against the master species: a fragment
of a database is read as well as a whole one. The species a reaction defines is
the first one to the right of its equal sign, as the PHREEQC manual requires,
and its name becomes a ChemistryLab symbol by the rule of
[`build_species`](@ref): a neutral species gets `@` appended (`CO2` is `CO2@`),
and `H2O` is the solvent.

A second `-gamma` for a species replaces the first, as PHREEQC reads the options
of a species in order. A `-gamma` that does not give two numbers is reported
with its line and skipped.
"""
function phreeqc_gamma_parameters(path::AbstractString)
    file = basename(path)
    notes = String[]
    params = Dict{String, Tuple{Float64, Float64}}()
    for b in phreeqc_blocks(resolve_data_path(path))
        b.keyword == "SOLUTION_SPECIES" || continue
        for rec in _species_records(b, file, notes)
            rec.gamma === nothing || (params[_phreeqc_solute_symbol(rec.name)] = rec.gamma)
        end
    end
    for note in notes
        occursin("`-gamma`", note) && @warn "phreeqc_gamma_parameters: $note"
    end
    return params
end

"""
    extract_primary_species(path) -> DataFrame

The primary species of a PHREEQC database: the master species its
`SOLUTION_SPECIES` define by an identity reaction (`Ca+2 = Ca+2`), every other
species being a combination of them, and the electron, written `Zz`, last.

One row per species, with the columns `species` (the PHREEQC name), `symbol`,
`formula`, `aggregate_state`, `atoms`, `charge` and `gamma` (the `(å, b)` of
`-gamma`, empty without one).
"""
function extract_primary_species(path::AbstractString)
    _, substances, _ = read_phreeqc_database(resolve_data_path(path))
    rows = NamedTuple{(:species, :symbol, :gamma), Tuple{String, String, Vector{Float64}}}[]
    for r in eachrow(substances)
        r.aggregate_state == AS_AQUEOUS && isempty(r.formation.terms) || continue
        push!(rows, (species = r.name, symbol = r.symbol, gamma = ismissing(r.gamma) ? Float64[] : [r.gamma...]))
    end
    push!(rows, (species = "Zz", symbol = "Zz", gamma = Float64[]))
    df = DataFrame(rows)
    df.formula = copy(df.symbol)
    df.aggregate_state .= "AS_AQUEOUS"
    df.atoms = parse_formula.(df.symbol)
    df.charge = [s == "Zz" ? 1 : extract_charge(s) for s in df.symbol]
    return df
end
