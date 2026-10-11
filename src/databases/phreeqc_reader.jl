# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

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
phreeqc_blocks(path::AbstractString) = _phreeqc_blocks(eachline(path))

# The scanner every reader of a PHREEQC file goes through: the databases, their
# PITZER and SIT blocks, and the sorption models. A block runs from its keyword
# to the next, whatever the indentation of its lines. With `comments = true` a
# logical line keeps the comment of its physical line, where a sorption model
# writes the source and the uncertainty of a constant (`ref:`, `error:`).
function _phreeqc_blocks(lines; comments::Bool = false)
    blocks = PhreeqcBlock[]
    current = nothing
    for (n, raw) in enumerate(lines)
        body, comment = _split_comment(raw)
        text = strip(body)
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
            isempty(p) && continue
            push!(current.lines, (n, comments && !isempty(comment) ? "$p # $comment" : String(p)))
        end
    end
    return blocks
end

# A line as its text and its comment, which follows the first `#`.
function _split_comment(line::AbstractString)
    i = findfirst('#', line)
    return i === nothing ? (line, "") : (line[1:prevind(line, i)], strip(line[nextind(line, i):end]))
end

# ── Equations ────────────────────────────────────────────────────────────────

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
    PhreeqcLogK

The log K of one reaction as a PHREEQC database states it: at 25 °C, with an
enthalpy of reaction for the van 't Hoff equation, or an analytical expression
of the temperature, which takes precedence over both (Parkhurst and Appelo
2013, `SOLUTION_SPECIES` and `PHASES`). `add` holds the named expressions the
constant adds (`-add_logk`), as `(index, coefficient)`.
"""
struct PhreeqcLogK{T <: Real} <: AbstractLogK
    log_k::T
    delta_h::T                            # kJ/mol
    analytic::Union{Nothing, NTuple{6, T}}
    add::Vector{Tuple{Int, Float64}}
end
function PhreeqcLogK(log_k::Real, delta_h::Real, analytic, add)
    T = promote_type(typeof(log_k), typeof(delta_h), (analytic === nothing ? () : map(typeof, analytic))...)
    return PhreeqcLogK{T}(log_k, delta_h, analytic === nothing ? nothing : map(T, analytic), add)
end
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
    c = 1000 * k.delta_h / _R_LN10
    return k.log_k - c * (1 / T - 1 / T_STANDARD), c / T^2, -2c / T^3
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
# Each option is a method of `_option!`, which receives the words after the
# option (`args`) and the numbers that lead them (`v`).
function _apply_option!(rec::_Record, opt::Symbol, parts, n, file)
    args = parts[2:end]
    return _option!(Val(opt), rec, args, _numbers(args), (; name = parts[1], n, file))
end

_option!(::Val{:log_k}, rec, args, v, at) = (isempty(v) || (rec.log_k = v[1]); nothing)
_option!(::Val{:delta_h}, rec, args, v, at) =
    (isempty(v) || (rec.delta_h = _delta_h_kJ(v[1], length(args) > 1 ? args[2] : nothing)); nothing)
_option!(::Val{:analytic}, rec, args, v, at) =
    (isempty(v) || (rec.analytic = ntuple(i -> i <= length(v) ? v[i] : 0.0, 6)); nothing)
function _option!(::Val{:gamma}, rec, args, v, at)
    length(v) >= 2 || return "$(at.file):$(at.n): the `-gamma` of $(rec.name) does not parse as two numbers; not read"
    rec.gamma = (v[1], v[2])
    return nothing
end
_option!(::Val{:llnl_gamma}, rec, args, v, at) = (isempty(v) || (rec.llnl_gamma = v[1]); nothing)
_option!(::Val{:vm}, rec, args, v, at) = (isempty(v) || (rec.vm = v[1]); nothing)
_option!(::Val{:t_c}, rec, args, v, at) = (isempty(v) || (rec.critical[1] = v[1]); nothing)
_option!(::Val{:p_c}, rec, args, v, at) = (isempty(v) || (rec.critical[2] = v[1]); nothing)
_option!(::Val{:omega}, rec, args, v, at) = (isempty(v) || (rec.critical[3] = v[1]); nothing)
function _option!(::Val{:add_logk}, rec, args, v, at)
    # The name of the expression, then its coefficient (1 by default).
    isempty(args) && return "$(at.file):$(at.n): `-add_logk` names no expression; not read"
    c = length(args) >= 2 ? something(tryparse(Float64, args[2]), 1.0) : 1.0
    push!(rec.add, (String(args[1]), c))
    return nothing
end
_option!(::Val{:add_constant}, rec, args, v, at) = (isempty(v) || push!(rec.add, ("\$constant", v[1])); nothing)
_option!(::Val{:mole_balance}, rec, args, v, at) = (isempty(args) || (rec.mole_balance = String(args[1])); nothing)
_option!(::Val{:no_check}, rec, args, v, at) = (rec.no_check = true; nothing)
_option!(::Val{:co2_llnl_gamma}, rec, args, v, at) = (rec.co2_gamma = true; nothing)
# Transport and density data, not thermodynamics.
_option!(::Union{Val{:dw}, Val{:millero}, Val{:activity_water}, Val{:tracer_diffusion}, Val{:dalfg}}, rec, args, v, at) = nothing
# A surface electrostatic option changes the thermodynamics of the species: the
# species is not read.
_option!(::Union{Val{:cd_music}, Val{:erm_ddl}}, rec, args, v, at) =
    (push!(rec.refused, "$(at.file):$(at.n): `$(at.name)` (a surface electrostatic option) is not read"); nothing)
_option!(::Val, rec, args, v, at) =
    "$(at.file):$(at.n): option `$(at.name)` of $(rec.name) is not one the PHREEQC manual documents; not read"

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
    blocks = phreeqc_blocks(path)
    notes = String[]
    db = _phreeqc_records(blocks, file, notes)

    # The constants: the named expressions first, then one per species, in the
    # order of the file.
    c = _PhreeqcConstants(AbstractLogK[], Dict{String, Int}())
    for rec in db.named
        c.named[rec.name] = _constant!(c, rec, file)
    end
    by_name = Dict{String, Tuple{_Record, Int}}()
    for rec in db.species
        by_name[rec.name] = (rec, _constant!(c, rec, file))
    end
    resolve = _Resolver((name, r) -> _phreeqc_define(name, by_name, db.pseudo, c.logks, r))

    species = _ReactionEntry[]
    for rec in db.species
        r = resolve(rec.name)
        if r === nothing
            push!(notes, "$file:$(rec.line): $(rec.name) rests on a species the database does not define; not read")
            continue
        end
        _refusal(rec, r, file, notes) && continue
        _not_elements(r[1], file, rec.line, rec.name, notes) && continue
        atoms, z, formation = r
        push!(species, _ReactionEntry(rec.name, :solution, atoms, z, formation, rec.gamma, rec.llnl_gamma, rec.co2_gamma, nothing, nothing, rec.line, rec.name))
    end

    phases = _ReactionEntry[]
    phase_logk = Dict{String, Int}()
    for rec in db.phases
        lhs, rhs = phreeqc_equation(rec.equation)
        isempty(lhs) && (push!(notes, "$file:$(rec.line): phase $(rec.name) has no formula on the left; not read"); continue)
        k = _constant!(c, rec, file)
        phase_logk[rec.name] = k
        # The phase is the first term on the left; everything else is aqueous.
        terms = vcat([(-lhs[1][1], "\$phase")], [(-x, s) for (x, s) in lhs[2:end]], [(x, s) for (x, s) in rhs])
        r = _resolve_defined("\$phase", terms, k, c.logks, resolve)
        if r === nothing
            push!(notes, "$file:$(rec.line): phase $(rec.name) rests on a species the database does not define; not read")
            continue
        end
        _refusal(rec, r, file, notes; formula = lhs[1][2]) && continue
        _not_elements(r[1], file, rec.line, rec.name, notes) && continue
        atoms, z, formation = r
        # The critical pressure in pascals; PHREEQC gives it in atm.
        critical = all(!isnothing, rec.critical) ? (rec.critical[1], rec.critical[2] * _ONE_ATM, rec.critical[3]) : nothing
        push!(phases, _ReactionEntry(rec.name, :phase, atoms, z, formation, nothing, nothing, false, rec.vm, critical, rec.line, lhs[1][2]))
    end

    # The reactions as written, each with its constant.
    reactions = NamedTuple[]
    for rec in db.named
        push!(reactions, (; symbol = rec.name, kind = :named, equation = "", logk = c.logks[c.named[rec.name]], line = rec.line))
    end
    for rec in db.species
        haskey(by_name, rec.name) &&
            push!(reactions, (; symbol = rec.name, kind = :solution, equation = rec.equation, logk = c.logks[by_name[rec.name][2]], line = rec.line))
    end
    for rec in db.phases
        haskey(phase_logk, rec.name) &&
            push!(reactions, (; symbol = rec.name, kind = :phase, equation = rec.equation, logk = c.logks[phase_logk[rec.name]], line = rec.line))
    end
    gauge = "PHREEQC master species of $file ($(_short_digest(path)))"
    data = _ReactionData(
        file, String(path), _source_tag(path), gauge, :phreeqc, db.masters, species, phases, c.logks,
        reactions, db.llnl, Set(b.keyword for b in blocks), notes,
    )
    return _reaction_tables(data)
end

# The records of the blocks of a database that hold thermodynamic data, each kind
# in the order of the file, a species or a phase defined again replacing the
# earlier definition; and a note for each block that is not read.
function _phreeqc_records(blocks, file, notes)
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
    species = _last_definitions(srecs, file, notes)
    phases = _last_definitions(precs, file, notes)
    return (; masters, pseudo, species, phases, named = nrecs, llnl)
end

# The constants of a database, as they are defined: `logks` by index, and the
# index of each named expression.
struct _PhreeqcConstants
    logks::Vector{AbstractLogK}
    named::Dict{String, Int}
end

# The constant of a record, after those it adds (`-add_constant`; a named
# expression it adds, `-add_logk`, is one defined before it); its index.
function _constant!(c::_PhreeqcConstants, rec, file)
    add = Tuple{Int, Float64}[]
    for (nm, x) in rec.add
        if nm == "\$constant"
            push!(c.logks, PhreeqcLogK(x, 0.0, nothing, Tuple{Int, Float64}[]))
            push!(add, (length(c.logks), 1.0))
        elseif haskey(c.named, nm)
            push!(add, (c.named[nm], x))
        else
            push!(rec.refused, "$file:$(rec.line): `-add_logk $nm` names no expression of NAMED_EXPRESSIONS")
        end
    end
    push!(c.logks, PhreeqcLogK(rec.log_k, rec.delta_h, rec.analytic, add))
    return length(c.logks)
end

# A species of SOLUTION_SPECIES resolved against the master species: a master by
# its formula, its identity reaction putting it at zero (the gauge), any other by
# its reaction. The electron is a master of its own.
function _phreeqc_define(name, by_name, pseudo, logks, resolve)
    if !haskey(by_name, name)
        name == "e-" && return (Dict{Symbol, Float64}(), -1.0, Dict{Int, Float64}())
        return nothing
    end
    rec, k = by_name[name]
    lhs, rhs = phreeqc_equation(rec.equation)
    if length(lhs) == 1 && length(rhs) == 1 && lhs[1][2] == rhs[1][2]
        atoms, z = _master_composition(name, pseudo)
        return (atoms, z, Dict{Int, Float64}())
    end
    terms = vcat([(-x, s) for (x, s) in lhs], [(x, s) for (x, s) in rhs])
    return _resolve_defined(name, terms, k, logks, resolve)
end

# The composition and charge of a master species, from its formula. A
# pseudo-element is kept as itself, which `_not_elements` refuses.
function _master_composition(name, pseudo)
    name == "e-" && return Dict{Symbol, Float64}(), -1.0
    haskey(pseudo, name) && return Dict{Symbol, Float64}(Symbol("pseudo-element ", pseudo[name]) => 1.0), _name_charge(name)
    atoms = Dict{Symbol, Float64}(Symbol(k) => Float64(v) for (k, v) in parse_formula(name))
    return atoms, _name_charge(name)
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
            rec.gamma === nothing || (params[_solute_symbol(rec.name)] = rec.gamma)
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
