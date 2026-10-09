# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using DataFrames
using DynamicQuantities
using JSON
using OrderedCollections
using ProgressMeter
using Tables
using TOML

"""
    HKF_SI_CONVERSIONS

Hardcoded conversion factors from SUPCRT [Johnson1992](@cite) (cal, bar) to SI
(J, Pa) for `eos_hkf_coeffs`. The JSON unit metadata for a3 and a4 omit the
`/bar` the SUPCRT units carry, so the conversion factors are set here
explicitly.

| Symbol | SUPCRT unit             | SI unit             | Factor     |
|--------|-------------------------|---------------------|------------|
| a1     | cal/(mol·bar)           | J/(mol·Pa)          | 4.184e-5   |
| a2     | cal/mol                 | J/mol               | 4.184      |
| a3     | (cal·K)/(mol·bar)       | (J·K)/(mol·Pa)      | 4.184e-5   |
| a4     | cal·K/mol               | J·K/mol             | 4.184      |
| c1     | cal/(mol·K)             | J/(mol·K)           | 4.184      |
| c2     | cal·K/mol               | J·K/mol             | 4.184      |
| wref   | cal/mol                 | J/mol               | 4.184      |
"""
const HKF_SI_CONVERSIONS = OrderedDict{Symbol, Float64}(
    :a1 => CALORIE / ustrip(us"Pa", 1.0u"bar"),
    :a2 => CALORIE,
    :a3 => CALORIE / ustrip(us"Pa", 1.0u"bar"),
    :a4 => CALORIE,
    :c1 => CALORIE,
    :c2 => CALORIE,
    :wref => CALORIE,
)

# Banners and progress bars are for someone watching a terminal. Written to a file,
# or captured into the documentation, they are three boxes and two bars per
# database loaded, between the reader and the result; they are drawn only when the
# stream is a terminal.
_banner(title, color) = stdout isa Base.TTY &&
    print_title(title; crayon = Crayon(; foreground = color), style = :box, indent = "")
_progress(n) = stderr isa Base.TTY ? Progress(n) : nothing
_tick!(::Nothing) = nothing
_tick!(p) = next!(p)

"""
    read_thermofun_database(filename::AbstractString) -> (DataFrame, DataFrame, DataFrame)

Read a ThermoFun database from a JSON file.

# Arguments

  - `filename`: path to the JSON database file.

# Returns

  - `df_elements`: DataFrame of chemical elements.
  - `df_substances`: DataFrame of chemical substances (species). A substance the
    database defines by a reaction carries that reaction's record in the column
    `defining_reaction`, from which [`build_species`](@ref) computes it.
  - `df_reactions`: DataFrame of chemical reactions.
"""
function read_thermofun_database(filename)
    path = resolve_data_path(filename)
    # `display_data_path`, not `path`: a banner must not carry a build machine's
    # directories.
    _banner("Loading database: $(display_data_path(path))", :green)
    data = JSON.parsefile(path)
    df_substances = DataFrame(Tables.dictrowtable(data["substances"]))
    df_reactions = DataFrame(Tables.dictrowtable(data["reactions"]))
    # A substance defined by a reaction names it by its symbol. A symbol given to
    # two records keeps the first, as the two copies Cemdata18 carries of some of
    # its reactions differ only in listing each reactant twice.
    if hasproperty(df_substances, :reaction)
        by_symbol = Dict{String, Any}()
        for r in data["reactions"]
            get!(by_symbol, String(r["symbol"]), r)
        end
        df_substances.defining_reaction = Any[
            (ismissing(r) || r === nothing) ? missing : get(by_symbol, String(r), missing)
                for r in df_substances.reaction
        ]
    end
    df_elements = DataFrame(Tables.dictrowtable(data["elements"]))
    return df_elements, df_substances, df_reactions
end

# Validate the entire syntax tree before calling uparse: its parser can evaluate
# arbitrary call heads. Only unit arithmetic and named constants are data here.
function is_unit_expression(ex)
    if ex isa Real
        return isfinite(ex)
    elseif ex isa Symbol
        return true # uparse itself resolves symbols only from its unit registry.
    elseif ex isa Expr && ex.head == :call && length(ex.args) >= 2
        op = first(ex.args)
        arity = length(ex.args) - 1
        allowed = if op in (:+, :*)
            arity >= 1
        elseif op == :-
            arity in (1, 2)
        elseif op in (:/, ://, :^)
            arity == 2
        elseif op in (:sqrt, :√, :cbrt, :∛)
            arity == 1
        else
            false
        end
        return allowed && all(is_unit_expression, ex.args[2:end])
    elseif ex isa Expr && ex.head == :. && length(ex.args) == 2
        return ex.args[1] == :Constants && ex.args[2] isa QuoteNode &&
            ex.args[2].value isa Symbol
    end
    return false
end

"""
    extract_unit(v, default_unit=u"1") -> AbstractQuantity

Parse unit arithmetic (including powers, roots, and `Constants` names).
Reject executable syntax before calling `uparse`. Returns `default_unit` for
unsupported expressions, unknown units, or malformed input.
"""
function extract_unit(v, default_unit = u"1")
    return try
        is_unit_expression(Meta.parse(v)) || return default_unit
        uparse(v)
    catch
        default_unit
    end
end

# Classification labels are enum names, never Julia expressions. Preserve the
# existing undefined fallback for missing, malformed, and unsupported labels.
function extract_classification(value, fallback::T) where {T <: Enum}
    return try
        label = only(values(value))
        index = findfirst(x -> string(x) == label, instances(T))
        isnothing(index) ? fallback : instances(T)[index]
    catch
        fallback
    end
end

"""
    extract_value(row, field; verbose=false, default_unit=u"1", with_units=true) -> Union{AbstractQuantity, Number, Missing}

Extract a scalar value (optionally with units) from a nested ThermoFun DataFrame `row`.
Returns `missing` when the field is absent, missing, or cannot be parsed.
"""
function extract_value(
        row, field::Symbol; verbose = false, default_unit = u"1", with_units = true
    )
    if haskey(row, field) && !ismissing(row[field]) && haskey(row[field], :values)
        try
            val = only(row[field].values)
            if with_units
                if iszero(val)
                    val *= default_unit
                elseif haskey(row[field], :units)
                    vunit = only(get(row[field], :units, [""]))
                    val *= extract_unit(vunit, default_unit)
                else
                    val *= default_unit
                end
            end
            if verbose
                println("$(row.symbol) => $field=$val")
            end
            return val
        catch
            return missing
        end
    else
        return missing
    end
end

correct_volume_unit(v::AbstractQuantity) = uamount(v) != -1 ? v / 1u"mol" : v

correct_volume_unit(v) = v

"""
    _reference_cp_interval(methods, Tref) -> method or nothing

The heat-capacity method of a ThermoFun substance that applies at the reference
temperature `Tref` (K).

A substance whose heat capacity changes form at a phase transition, such as
quartz or hematite, lists one `cp_ft_equation` per temperature interval. The
thermodynamic functions are anchored at `Tref`, so the interval containing it is
the one retained; a method without temperature limits applies everywhere, and
when no interval contains `Tref` the first one listed is kept. The functions are
anchored there and followed into the other intervals by
`complete_thermo_functions!` (see `_cp_intervals`).
"""
function _reference_cp_interval(methods, Tref::Real)
    cps = [
        m for m in methods if
            only(values(m.method)) == "cp_ft_equation" && haskey(m, :m_heat_capacity_ft_coeffs)
    ]
    isempty(cps) && return nothing
    i = findfirst(cps) do m
        lim = get(m, :limitsTP, nothing)
        lim === nothing && return true
        get(lim, :lowerT, -Inf) <= Tref <= get(lim, :upperT, Inf)
    end
    return cps[something(i, 1)]
end

"""
    _cp_intervals(methods) -> Vector{NamedTuple}

The heat-capacity intervals of a ThermoFun substance, sorted by temperature:
`lower` and `upper` (K), `coeffs` (the `aᵢ` with their units, as
`cp_ft_equation` takes them) and `transition`, `nothing` or the `(T, dS, dH)` of
the phase transition the record places at the top of the interval
(`m_phase_trans_props`). A method without temperature limits makes the list
empty: it applies everywhere.
"""
function _cp_intervals(methods)
    out = NamedTuple[]
    for m in methods
        only(values(m.method)) == "cp_ft_equation" || continue
        haskey(m, :m_heat_capacity_ft_coeffs) || continue
        lim = get(m, :limitsTP, nothing)
        lim === nothing && return NamedTuple[]
        coeffs = m.m_heat_capacity_ft_coeffs
        units = extract_unit.(coeffs.units)
        params = [
            Symbol("a", subscriptnumber(i - 1)) => float(coeffs.values[i] * units[i])
                for i in 1:min(length(coeffs.values), length(units))
        ]
        tr = nothing
        if haskey(m, :m_phase_trans_props)
            names = String.(m.m_phase_trans_props.names)
            v = Float64.(m.m_phase_trans_props.values)
            get_(n) = (k = findfirst(==(n), names); k === nothing ? 0.0 : v[k])
            tr = (T = get_("Temperature"), dS = get_("dS"), dH = get_("dH"))
        end
        push!(out, (lower = Float64(get(lim, :lowerT, -Inf)), upper = Float64(get(lim, :upperT, Inf)), coeffs = params, transition = tr))
    end
    sort!(out; by = x -> x.lower)
    return out
end

"""
    temperature_range(s::AbstractSpecies) -> (lower, upper)

The temperatures, in K, over which the record of `s` declares its heat capacity:
the `limitsTP` of its `cp_ft_equation` methods in a ThermoFun database, from the
lowest bound of its first interval to the highest of its last. `(-Inf, Inf)`
for a species whose record declares none, which includes every solute described
by the HKF equations [Helgeson1981](@cite), whose domain is a density of water
rather than a temperature (see [Where the HKF equations
hold](@ref sec-theory-hkf-domain)), and the solvent.

Nothing refuses a temperature outside the range: the heat capacity is carried
past it by the nearest interval's function. The ranges are narrow in places.
Cemdata18 declares its AFm phases to 50 °C, its ettringites to 60 °C and its
clinker phases from 25 °C only, so that a calculation at 80 °C extrapolates
ettringite and one at 5 °C extrapolates alite; a warning at every such
evaluation would fire in most calculations below 25 °C, and the function is
there for the calculation that needs to know.

```jldoctest
julia> s = Species("CaCO3"; aggregate_state = AS_CRYSTAL);

julia> temperature_range(s)
(-Inf, Inf)
```
"""
function temperature_range(s::AbstractSpecies)
    haskey(properties(s), :T_range) || return (-Inf, Inf)
    r = s[:T_range]
    return (Float64(r[1]), Float64(r[2]))
end

"""
    complete_species_with_thermo_model!(species, row; verbose=false)

Populate thermodynamic reference values and build thermodynamic functions on `species`
from a ThermoFun substance DataFrame `row`. Mutates `species.properties` in place.
"""
# ThermoFun's standard reference state is 298.15 K and 1 bar, and a record that
# omits it is referred to that state: one substance of the slop98 organic
# database, `Eth@`, carries no `Tst`. Read as missing, it made the whole database
# unreadable.
const _THERMOFUN_TST = T_STANDARD
const _THERMOFUN_PST = P_STANDARD
_reference_value(row, key, default) = (
    v = hasproperty(row, key) ? getproperty(row, key) : missing;
    (ismissing(v) || v === nothing) ? default : v
)

# The methods of a ThermoFun record computed as their publications define them:
# the volume and the Landau transition of a mineral of Holland and Powell,
# `mv_eos_murnaghan_hp98` and `landau_holland_powell98` (`_holland_powell98!`), and
# a dissolved gas of Akinfiev and Diamond, `solute_aknifiev_diamond03`
# (`_akinfiev_diamond!`).
#
# The methods of a ThermoFun record that need nothing of their own here:
#   - `standard_entropy_cp_integration` extrapolates the reference values with
#     the heat capacity held constant, which is what a record without a
#     heat-capacity function gets (`complete_thermo_functions!`);
#   - `water_diel_*` is the permittivity of the solvent, which the package takes
#     from its own model of water (`water_electro_props_jn`);
#   - `fluid_comp_redlich_kwong_hp91`, the compensated Redlich-Kwong fluid of
#     Holland and Powell, leaves the standard Gibbs energy that of the ideal gas,
#     as ThermoFun applies it, and moves the departure from the ideal gas into
#     the fugacity, which here is the gas phase's equation of state
#     ([`peng_robinson`](@ref)) [HollandPowell1991](@cite).
# Any other method is not computed, and the species carrying it is refused by
# `build_species`.
const _THERMOFUN_EQUIVALENT_METHODS = (
    "standard_entropy_cp_integration", "water_diel_jnort91_reaktoro", "fluid_comp_redlich_kwong_hp91",
)

function complete_species_with_thermo_model!(species, row; verbose = false)
    Tst = _reference_value(row, :Tst, _THERMOFUN_TST)
    Tref = Tst * u"K"
    Pref = _reference_value(row, :Pst, _THERMOFUN_PST) * u"Pa"
    species.Tref = Tref
    species.Pref = Pref
    values0 = [
        :Cp⁰ => extract_value(
            row, :sm_heat_capacity_p; verbose = verbose, default_unit = u"J/K/mol"
        ),
        :ΔₐH⁰ => extract_value(row, :sm_enthalpy; verbose = verbose, default_unit = u"J/mol"),
        :S⁰ =>
            extract_value(row, :sm_entropy_abs; verbose = verbose, default_unit = u"J/K/mol"),
        :ΔₐG⁰ =>
            extract_value(row, :sm_gibbs_energy; verbose = verbose, default_unit = u"J/mol"),
        :V⁰ => correct_volume_unit(extract_value(row, :sm_volume; verbose = verbose, default_unit = u"J/bar")),
    ]
    # A solid whose molar volume is recorded as exactly zero has none defined:
    # Cemdata18 prints "not defined" for the amorphous and microcrystalline
    # Fe(OH)3, and its ThermoFun file writes 0. Kept, the zero would make a
    # precipitated hydroxide occupy no volume; left out, the species is reported
    # by `missing_molar_volumes` instead.
    v⁰ = last(values0[end])
    if aggregate_state(species) == AS_CRYSTAL && v⁰ isa Union{Number, DynamicQuantities.AbstractQuantity} && iszero(ustrip(v⁰))
        values0 = values0[1:(end - 1)]
    end
    species[:thermo_params] = [values0; :T => Tref; :P => Pref]
    # The energies are of formation from the elements (`_refuse_mixed_gauges`).
    species[:gauge] = "formation from the elements"
    TPMethods = row.TPMethods
    if !ismissing(TPMethods)
        cp_interval = _reference_cp_interval(TPMethods, Tst)
        # Every heat-capacity interval, with the transition at its top when the
        # record gives one, so that the functions follow T past the interval
        # that holds Tref (`complete_thermo_functions!`).
        intervals = _cp_intervals(TPMethods)
        # The temperatures the record declares its heat capacity for, kept so
        # that a calculation can be checked against them (`temperature_range`).
        isempty(intervals) || (species[:T_range] = [first(intervals).lower, last(intervals).upper])
        # Held as a function returning the list, a property being a number, a
        # function, a string or a vector of numbers or of pairs.
        length(intervals) > 1 && (species[:cp_intervals] = () -> intervals)
        murnaghan, landau = false, nothing
        for method in TPMethods
            method_type = only(values(method.method))
            if method_type == "cp_ft_equation" && method === cp_interval
                species[:thermo_method] = "cp_ft_equation"
                coeffs = method.m_heat_capacity_ft_coeffs
                vals = coeffs.values
                units = extract_unit.(coeffs.units)
                params = [
                    Symbol("a", subscriptnumber(i - 1)) => float(vals[i] * units[i]) for
                        i in 1:min(length(vals), length(units))
                ]
                species[:thermo_params] = [params; species[:thermo_params]]

            elseif method_type == "solute_hkf88_reaktoro" && haskey(method, :eos_hkf_coeffs)
                species[:thermo_method] = "solute_hkf88_reaktoro"
                coeffs = method.eos_hkf_coeffs
                vals = float.(coeffs.values)
                names = [:a1, :a2, :a3, :a4, :c1, :c2, :wref]
                hkf_params = [
                    names[i] => vals[i] * HKF_SI_CONVERSIONS[names[i]] for
                        i in 1:min(length(vals), length(names))
                ]
                z = float(get(row, :formula_charge, 0))
                push!(hkf_params, :z => z)
                species[:thermo_params] = [hkf_params; species[:thermo_params]]

            elseif method_type in ("mv_constant", "mv_pvnrt")
                species[:V_method] = method_type
            elseif startswith(method_type, "water_eos")
                # The solvent's standard state is that of liquid water by the
                # equation of state of Haar, Gallagher and Kell, at any T and P,
                # anchored on the record at its reference
                # (`_solvent_from_water_eos!`).
                species[:V_method] = "water_eos"
            elseif method_type in _THERMOFUN_EQUIVALENT_METHODS
                # Computed the way the package computes a record without a
                # method of its own; see `_THERMOFUN_EQUIVALENT_METHODS`.
                nothing
            elseif method_type == "mv_eos_murnaghan_hp98"
                murnaghan = true
            elseif method_type == "landau_holland_powell98" && haskey(method, :m_landau_phase_trans_props)
                landau = Float64.(method.m_landau_phase_trans_props.values)
            elseif method_type == "solute_aknifiev_diamond03" && haskey(method, :eos_akinfiev_diamond_coeffs)
                # It needs the water of the database: refused until
                # `_species_from_row` applies it (`_akinfiev_diamond!`).
                species[:ad03] = Float64.(method.eos_akinfiev_diamond_coeffs.values[1:3])
                _refuse_method!(species, method_type)
            elseif !(method_type == "cp_ft_equation" || method_type == "solute_hkf88_reaktoro")
                _refuse_method!(species, method_type)
            end
        end
        if murnaghan || landau !== nothing
            p = _hp98_parameters(row, landau, Tst, ustrip(us"Pa", Pref))
            if p !== nothing
                _holland_powell98!(species, p)
            elseif landau === nothing
                # Without a bulk modulus, the volume of the record is held
                # constant, as ThermoFun holds it (Gibbsite and Boehmite in aq17).
                species[:V_method] = "mv_constant"
            else
                _refuse_method!(species, "landau_holland_powell98")
            end
        end
    end
    return species
end

_refuse_method!(species, method) = (
    species[:refused_method] = haskey(properties(species), :refused_method) ?
        species[:refused_method] * "`, `" * method : method
)

"""
    build_species(df_substances::AbstractDataFrame, list_symbols=nothing; verbose=false) -> Vector{Species}

Build Species objects from a substance DataFrame: that of a ThermoFun file
([`read_thermofun_database`](@ref)), of a database of reactions
([`read_phreeqc_database`](@ref), [`read_gwb_database`](@ref),
[`read_eq36_database`](@ref)) or of a database of Reaktoro
([`read_reaktoro_database`](@ref)). The substances named in `list_symbols` are
looked up by their symbol, or by their name in the database.

# Arguments

  - `df_substances`: DataFrame containing substance data.
  - `list_symbols`: optional list of symbols to filter (default: nothing, process all).
  - `verbose`: if true, print details during processing (default: false).

# Returns

  - Vector of `Species`.

# Substances defined by a reaction

A substance the database defines by a reaction (`defining_reaction`, see
[`read_thermofun_database`](@ref)) gets its standard properties at any
temperature from that reaction, as ThermoFun computes them: with the reaction
``\\sum_i \\nu_i \\mathrm{A}_i = 0`` and its properties ``\\Delta_r X`` from its
``\\log K(T)``,
``X_s = (\\Delta_r X - \\sum_{i \\ne s} \\nu_i X_i) / \\nu_s`` for ``G``, ``H``,
``S`` and ``C_p``, and for the molar volume ``V`` of a solute; a solid keeps the
molar volume its record gives. The other species of the reaction are built for
the purpose when they are not in `list_symbols`.

# Methods of aq17

The volume of a mineral of Holland and Powell (`mv_eos_murnaghan_hp98`, the
Murnaghan equation with a thermal expansion) and its order-disorder transition
(`landau_holland_powell98`, tricritical Landau theory) are computed as
[HollandPowell1998](@cite) defines them, the volume integrated from the
reference pressure; a mineral whose record gives no bulk modulus keeps the
volume of its record, its Gibbs energy moving by ``V^\\circ (P - P^\\circ)``, as
ThermoFun computes Gibbsite and Boehmite. A dissolved gas of Akinfiev and
Diamond (`solute_aknifiev_diamond03`) is computed as [AkinfievDiamond2003](@cite)
defines it, with the water and the water vapor of the same table. ThermoFun
integrates the volume of a mineral from zero pressure and leaves it out at the
reference state: away from that state its energies exceed these by the molar
volume times one bar, a few J/mol.

# Methods that are not computed

A substance whose record computes it by a method the package does not implement
(a method ThermoFun does not list, or a dissolved gas without the water of its
database) is built for its reference state only: its standard properties are
those of its record at the reference temperature, and within one bar of the
reference pressure; anywhere else they raise an error naming the method. A
warning lists such substances.
"""
function build_species(
        df_substances::AbstractDataFrame, list_symbols = nothing; verbose = false
    )
    # The tables of a database of reactions and of a database of Reaktoro.
    hasproperty(df_substances, :formation) && return _build_reaction_species(df_substances, list_symbols)
    hasproperty(df_substances, :standard_model) && return _build_reaktoro_species(df_substances, list_symbols)
    local_df_substances = if isnothing(list_symbols)
        df_substances
    else
        @view df_substances[df_substances.symbol .∈ Ref(list_symbols), :]
    end
    keylist = String[]
    species_list = Species[]
    refused = Tuple{String, String}[]
    _banner("Building species", :blue)
    progress = _progress(nrow(local_df_substances))
    for row in eachrow(local_df_substances)
        if verbose
            println(row[:symbol])
        end
        species = _species_from_row(row, df_substances; verbose = verbose)
        _tick!(progress)
        if haskey(properties(species), :refused_method)
            push!(refused, (row.symbol, species[:refused_method]))
            _reference_state_only!(species)
        end
        key = row.symbol
        if key in keylist
            @warn("Symbol $key is used for multiple species")
        end
        push!(keylist, key)
        push!(species_list, species)
    end
    isempty(refused) || @warn "build_species: $(length(refused)) substance(s) computed by a method ChemistryLab does not implement, usable at their reference state only: " *
        join(("$s (`$m`)" for (s, m) in refused), ", ")
    _define_by_reactions!(species_list, df_substances; verbose = verbose)
    return species_list
end

# The species of a row of `df`; a dissolved gas of Akinfiev and Diamond takes the
# water of `df` (`_akinfiev_diamond!`).
function _species_from_row(row, df = nothing; verbose = false)
    species = Species(
        row.formula;
        name = row.name,
        symbol = row.symbol,
        aggregate_state = extract_classification(get(row, :aggregate_state, missing), AS_UNDEF),
        class = extract_classification(get(row, :class_, missing), SC_UNDEF),
    )
    complete_species_with_thermo_model!(species, row; verbose = verbose)
    haskey(properties(species), :ad03) && _akinfiev_diamond!(species, df)
    return species
end

const _ONE_BAR = ustrip(us"Pa", 1.0u"bar")

# The functions of a substance computed by a method the package does not
# implement, restricted to its reference state, where the method leaves the
# record's values: the reference temperature (unless `temperature = false`, for a
# method that moves the properties with pressure only), and one bar about the
# reference pressure. Anywhere else they raise an error naming the method.
function _reference_state_only!(s; temperature::Bool = true)
    complete_thermo_functions!(s)
    sym, method = symbol(s), s[:refused_method]
    Tr, Pr = ustrip(us"K", s.Tref), ustrip(us"Pa", s.Pref)
    function check(T, P)
        ((!temperature || abs(_plain(T) - Tr) <= 1.0e-6) && abs(_plain(P) - Pr) <= _ONE_BAR) || throw(
            ArgumentError(
                "$sym: its record computes it by the method `$method`, which ChemistryLab does not implement; " *
                    "its properties are known " * (temperature ? "at its reference state only, $(Tr) K and " : "at ") *
                    "$(Pr) Pa within 1 bar, not at $(_plain(T)) K and $(_plain(P)) Pa. " *
                    "Leave it out of the system, or take it from another database.",
            )
        )
        return nothing
    end
    for key in (:ΔₐG⁰, :ΔₐH⁰, :S⁰, :Cp⁰, :V⁰)
        haskey(properties(s), key) || continue
        f = s[key]
        f isa AbstractFunc || continue
        s[key] = NumericFunc((T, P) -> (check(T, P); f(T = T, P = P)), (:T, :P), (T = s.Tref, P = s.Pref), f.unit)
    end
    return s
end


# Computes the substances of `species_list` that the database defines by a
# reaction, building the other species of each reaction from `df` as needed,
# themselves defined by a reaction or not.
function _define_by_reactions!(species_list, df; verbose = false)
    hasproperty(df, :defining_reaction) || return species_list
    rows = Dict(String(r.symbol) => r for r in eachrow(df))
    cache = Dict{String, Species}(symbol(s) => s for s in species_list)
    done, visiting = Set{String}(), Set{String}()
    function define!(s)
        sym = symbol(s)
        sym in done && return s
        row = get(rows, sym, nothing)
        rec = row === nothing ? missing : row.defining_reaction
        if !ismissing(rec)
            sym in visiting && throw(ArgumentError("$sym is defined by a reaction that rests on itself"))
            push!(visiting, sym)
            _define_by_reaction!(s, rec, species_of)
            delete!(visiting, sym)
        end
        push!(done, sym)
        return s
    end
    function species_of(sym)
        s = get(cache, sym, nothing)
        if s === nothing
            row = get(rows, sym, nothing)
            row === nothing && throw(ArgumentError("a reaction of the database names $sym, which it does not define"))
            s = _species_from_row(row, df; verbose = verbose)
            haskey(properties(s), :refused_method) && _reference_state_only!(s)
            cache[sym] = s
        end
        return define!(s)
    end
    foreach(define!, species_list)
    return species_list
end

# ── Substances defined by a reaction ─────────────────────────────────────────
#
# A ThermoFun database may define a substance by a reaction rather than by a
# heat-capacity function of its own: its record gives the standard properties at
# the reference temperature, and `reaction` names the reaction whose log K, with
# the properties of the other species of that reaction, gives them at any
# temperature. PSI/Nagra 12/07 defines 440 substances that way, Cemdata18 17
# (the zeolites, Fe(OH)3, MgSiO3@, S-2, ...). Extrapolated from the reference
# with its own heat capacity, such a substance follows neither its reaction nor
# ThermoFun, which computes it from the reaction.

# The log K of a ThermoFun reaction at T (kelvin) and its first two derivatives:
# the seven terms of `logk_fpt_function`,
#     log K = A₀ + A₁T + A₂/T + A₃ ln T + A₄/T² + A₅T² + A₆/√T,
# or, for a record without them, the value at the reference with the enthalpy
# and the heat capacity of reaction there, held constant (van 't Hoff's
# equation, integrated with a constant heat capacity of reaction).
function _thermofun_log10K(rec)
    coeffs = nothing
    for m in something(get(rec, "TPMethods", nothing), [])
        c = get(m, "logk_ft_coeffs", nothing)
        c === nothing || (coeffs = Float64.(c["values"]))
    end
    if coeffs !== nothing && any(!iszero, coeffs)
        length(coeffs) > 7 && any(!iszero, coeffs[8:end]) && throw(
            ArgumentError("reaction $(rec["symbol"]): logk_fpt_function has nonzero coefficients past the seventh, which no published form gives")
        )
        A = ntuple(i -> i <= length(coeffs) ? coeffs[i] : 0.0, 7)
        return function (T)
            L = A[1] + A[2] * T + A[3] / T + A[4] * log(T) + A[5] / T^2 + A[6] * T^2 + A[7] / sqrt(T)
            dL = A[2] - A[3] / T^2 + A[4] / T - 2A[5] / T^3 + 2A[6] * T - A[7] / (2 * T^1.5)
            d2L = 2A[3] / T^3 - A[4] / T^2 + 6A[5] / T^4 + 2A[6] + 3A[7] / (4 * T^2.5)
            return L, dL, d2L
        end
    end
    value(key) = (v = get(rec, key, nothing); v === nothing ? 0.0 : Float64(only(v["values"])))
    L0, H, Cp = value("logKr"), value("drsm_enthalpy"), value("drsm_heat_capacity_p")
    Tr = Float64(get(rec, "Tst", _THERMOFUN_TST))
    k = R_GAS * log(10)
    return function (T)
        L = L0 - H / k * (1 / T - 1 / Tr) + Cp / k * (Tr / T - 1 + log(T / Tr))
        dL = H / (k * T^2) + Cp / k * (1 / T - Tr / T^2)
        d2L = -2H / (k * T^3) + Cp / k * (2Tr / T^3 - 1 / T^2)
        return L, dL, d2L
    end
end

# The reactants of a reaction record, by symbol: a record that lists a species
# twice (the second copy of some reactions of Cemdata18 does) counts it once, and
# two different coefficients for one species are refused.
function _thermofun_reactants(rec)
    out = OrderedDict{String, Float64}()
    for r in rec["reactants"]
        s, c = String(r["symbol"]), Float64(r["coefficient"])
        if haskey(out, s) && out[s] != c
            throw(ArgumentError("reaction $(rec["symbol"]) gives $s two coefficients, $(out[s]) and $c"))
        end
        out[s] = c
    end
    return out
end

# Defines the standard properties of `s` by the reaction `rec`, from those of the
# other species of the reaction, which `species_of(symbol)` returns. With the
# reaction written Σ νᵢ Aᵢ = 0 (products positive) and its properties ΔᵣX:
#
#     X_s = (ΔᵣX − Σ_{i≠s} νᵢ Xᵢ) / ν_s ,
#
# for G, H, S and Cp, and V for a solute, where
#     ΔᵣG = −RT ln 10 log K(T) + ΔᵣV (P − P°),  ΔᵣH = RT² ln 10 d log K/dT + ΔᵣV (P − P°),
#     ΔᵣS = R ln 10 (log K + T d log K/dT),     ΔᵣCp = R ln 10 (2T d log K/dT + T² d² log K/dT²),
# and ΔᵣV the constant volume of reaction (`dr_volume_constant`). The molar
# volume of a solid or a gas is the one its record gives, or none when the
# record gives zero (Cemdata18 for Fe(OH)3): it is what the volume of the phase
# is computed from, which a volume of reaction does not give.
function _define_by_reaction!(s, rec, species_of)
    reactants = _thermofun_reactants(rec)
    me = symbol(s)
    haskey(reactants, me) || throw(ArgumentError("reaction $(rec["symbol"]) does not contain $me, which it is said to define"))
    ν_s = reactants[me]
    others = [(species_of(k), ν) for (k, ν) in reactants if k != me]
    logK = _thermofun_log10K(rec)
    ΔV = let v = get(rec, "drsm_volume", nothing)
        v === nothing ? 0.0 : ustrip(us"m^3/mol", Float64(only(v["values"])) * u"J/(bar*mol)")
    end
    Pr = Float64(get(rec, "Pst", _THERMOFUN_PST))
    k = R_GAS * log(10)
    ΔG(T, P) = -k * T * first(logK(T)) + ΔV * (P - Pr)
    ΔH(T, P) = (r = logK(T); k * T^2 * r[2] + ΔV * (P - Pr))
    ΔS(T, P) = (r = logK(T); k * (r[1] + T * r[2]))
    ΔCp(T, P) = (r = logK(T); k * (2T * r[2] + T^2 * r[3]))
    ΔVf(T, P) = ΔV + zero(T)
    refs = (T = s.Tref, P = s.Pref)
    for (key, Δ, unit) in (
            (:ΔₐG⁰, ΔG, u"J/mol"), (:ΔₐH⁰, ΔH, u"J/mol"), (:S⁰, ΔS, u"J/(mol*K)"),
            (:Cp⁰, ΔCp, u"J/(mol*K)"), (:V⁰, ΔVf, u"m^3/mol"),
        )
        key === :V⁰ && aggregate_state(s) != AS_AQUEOUS && continue
        all(o -> haskey(first(o), key), others) || continue
        fs = [(first(o)[key], last(o)) for o in others]
        f = function (T, P)
            acc = Δ(T, P)
            for (g, ν) in fs
                acc -= ν * g(T = T, P = P)
            end
            return acc / ν_s
        end
        s[key] = NumericFunc(f, (:T, :P), refs, unit)
    end
    s[:defining_reaction] = String(rec["symbol"])
    return s
end

function build_species(filename, list_symbols = nothing; verbose = false)
    _, df_substances, _ = read_thermofun_database(filename)
    return build_species(df_substances, list_symbols; verbose = verbose)
end

"""
    complete_reaction_with_thermo_model!(reaction, row; verbose=false)

Populate thermodynamic reference values and build thermodynamic functions on `reaction`
from a ThermoFun reaction DataFrame `row`. Mutates `reaction.properties` in place.
"""
function complete_reaction_with_thermo_model!(reaction, row; verbose = false)
    Tref = _reference_value(row, :Tst, _THERMOFUN_TST) * u"K"
    Pref = _reference_value(row, :Pst, _THERMOFUN_PST) * u"Pa"
    reaction.Tref = Tref
    reaction.Pref = Pref
    values0 = [
        :ΔᵣCp⁰ => extract_value(
            row, :drsm_heat_capacity_p; verbose = verbose, default_unit = u"J/K/mol"
        ),
        :ΔᵣH⁰ => extract_value(row, :drsm_enthalpy; verbose = verbose, default_unit = u"J/mol"),
        # ThermoFun writes the entropy of reaction as `drsm_entropy`; an older
        # spelling, `drsm_entropy_abs`, is read when it is the one given.
        :ΔᵣS⁰ => coalesce(
            extract_value(row, :drsm_entropy; verbose = verbose, default_unit = u"J/K/mol"),
            extract_value(row, :drsm_entropy_abs; verbose = verbose, default_unit = u"J/K/mol"),
        ),
        :ΔᵣG⁰ =>
            extract_value(row, :drsm_gibbs_energy; verbose = verbose, default_unit = u"J/mol"),
        :ΔᵣV⁰ => correct_volume_unit(extract_value(row, :drsm_volume; verbose = verbose, default_unit = u"J/bar")),
        :logKr => extract_value(row, :logKr; verbose = verbose, default_unit = u"1"),
    ]
    reaction[:thermo_params] = [values0; :T => Tref; :P => Pref]
    TPMethods = row.TPMethods
    if !ismissing(TPMethods)
        for method in TPMethods
            method_type = only(values(method.method))
            if method_type == "logk_fpt_function" && haskey(method, :logk_ft_coeffs)
                reaction[:logk_method] = "logk_fpt_function"
                coeffs = method.logk_ft_coeffs
                vals = coeffs.values
                # The seventh term is A₆/√T, as ThermoFun evaluates it.
                units = dimension.([1, u"1/K", u"K", 1, u"K^2", u"1/K^2", u"√K"])
                params = [
                    Symbol("A", subscriptnumber(i - 1)) =>
                        float(Quantity(vals[i], units[i])) for
                        i in 1:min(length(vals), length(units))
                ]
                reaction[:thermo_params] = [params; reaction[:thermo_params]]
            elseif method_type == "dr_volume_constant"
                reaction[:V_method] = "dr_volume_constant"
            end
        end
    end
    return reaction
end

"""
    build_reactions(df_reactions::AbstractDataFrame, dict_species=Dict(), list_symbols=nothing; verbose=false) -> Vector{Reaction}

Build Reaction objects from a reaction DataFrame.

# Arguments

  - `df_reactions`: DataFrame containing reaction data.
  - `species_list`: vector of existing `Species` objects to use in reactions.
  - `list_symbols`: optional list of reaction symbols to filter (default: nothing, process all).
  - `verbose`: if true, print details during processing (default: false).

# Returns

  - Vector of `Reaction` objects.
"""
function build_reactions(
        df_reactions::AbstractDataFrame,
        species_list = [],
        list_symbols = nothing;
        verbose = false,
    )
    local_df_reactions = if isnothing(list_symbols)
        df_reactions
    else
        @view df_reactions[df_reactions.symbol .∈ Ref(list_symbols), :]
    end
    dict_species = Dict(symbol(s) => s for s in species_list)
    keylist = String[]
    reactions_list = Reaction[]
    _banner("Building reactions", :red)
    # A reactant by its symbol, or by that symbol with `_` read as `.`; then by
    # the other names `find_species` knows. A reactant found nowhere is an error:
    # no other species stands in for it.
    function choose_species(k)
        haskey(dict_species, k) && return dict_species[k]
        kdot = replace(k, "_" => ".")
        haskey(dict_species, kdot) && return dict_species[kdot]
        # `find_species` returns a new species when none matches.
        s = find_species(k, collect(values(dict_species)))
        any(x -> x === s, values(dict_species)) ||
            throw(ArgumentError("a reaction names $k, which the species given do not contain"))
        return s
    end
    progress = _progress(nrow(local_df_reactions))
    for row in eachrow(local_df_reactions)
        if verbose
            println(row[:symbol])
        end
        # The reactants by symbol, whatever the order of the keys of each entry.
        stoich = [(k, ν) for (k, ν) in _thermofun_reactants(row) if k != "e-"]
        whole = all(isinteger(last(x)) for x in stoich)
        reaction = Reaction(
            OrderedDict(choose_species(k) => (whole ? Int(ν) : ν) for (k, ν) in stoich);
            symbol = row.symbol,
        )
        complete_reaction_with_thermo_model!(reaction, row; verbose = verbose)
        key = row.symbol
        if key in keylist
            @warn("Symbol $key is used for multiple reactions")
        end
        push!(keylist, key)
        push!(reactions_list, reaction)
        _tick!(progress)
    end
    return reactions_list
end

function build_reactions(filename, species_list = [], list_symbols = nothing; verbose = false)
    _, _, df_reactions = read_thermofun_database(filename)
    return build_reactions(df_reactions, species_list, list_symbols; verbose = verbose)
end

"""
    get_compatible_species(df_substances::AbstractDataFrame, species_list; aggregate_states=[AS_AQUEOUS], exclude_species=[], union=false) -> DataFrame

Find species in the database compatible with a given list of species (sharing atoms).

# Arguments

  - `df_substances`: substance DataFrame.
  - `species_list`: list of target species symbols.
  - `aggregate_states`: filter for specific aggregate states (default: `[AS_AQUEOUS]`).
  - `exclude_species`: list of species symbols to exclude.
  - `union`: if true, includes the original `species_list` in the result (default: false).

# Returns

  - DataFrame of compatible substances.
"""
function get_compatible_species(
        df_substances::AbstractDataFrame,
        species_list;
        aggregate_states = [AS_AQUEOUS],
        exclude_species = [],
        union = false,
    )
    df_given_species = @view df_substances[df_substances.symbol .∈ Ref(species_list), :]
    involved_atoms = union_atoms(parse_formula.(df_given_species.formula))
    mask1 = last.(only.(df_substances.aggregate_state)) .∈ Ref(string.(aggregate_states))
    mask2 = issubset.(keys.(parse_formula.(df_substances.formula)), Ref(involved_atoms))
    mask3 = .!(df_substances.symbol .∈ Ref(exclude_species))
    df_compat = @view df_substances[mask1 .&& mask2 .&& mask3, :]
    if union
        return unique(vcat(df_given_species, df_compat))
    else
        return df_compat
    end
end

"""
    build_solid_solutions(toml_file, dict_species; skip_missing=true, instances=Dict()) -> Vector{SolidSolutionPhase}

Load solid solution phase definitions from a TOML file and assemble
[`SolidSolutionPhase`](@ref) objects from an existing species dictionary.

Each end-member species is automatically requalified to `SC_SSENDMEMBER` via
[`with_class`](@ref), regardless of the class stored in the database.

# Arguments

  - `toml_file`: path to a TOML file with `[[solid_solution]]` entries (see
    `data/solid_solutions.toml` for the format).
  - `dict_species`: `Dict{String, <:AbstractSpecies}` mapping symbol → species
    (typically built from `Dict(symbol(s) => s for s in build_species(...))`).
  - `skip_missing`: if `true` (default), silently skip phases whose end-members
    are not all present in `dict_species`; if `false`, throw an error.

# TOML format

```toml
[[solid_solution]]
name        = "CSHQ"
end_members = ["CSHQ-TobD", "CSHQ-TobH", "CSHQ-JenH", "CSHQ-JenD", "KSiOH", "NaSiOH"]
model       = "ideal"          # or "redlich_kister", "regular", "sublattice", "compound_energy"
```

A Redlich-Kister [RedlichKister1948](@cite) entry gives its parameters in J/mol
as `a0`, `a1`, `a2`, or names published dimensionless Guggenheim parameters with
`guggenheim = "<literature key>:<pair>"`, a row of the `guggenheim_parameters`
table of `data/literature/<key>.json`; they become `a = α R T` at 298.15 K, so
that no published value is copied into the file. The sign of `a1` follows the
order of `end_members`. A model that unmixes needs `instances = 2`, the number
of coexisting compositions [`SolidSolutionPhase`](@ref) may give the phase, or
`instances = "auto"`, which gives it the second only when a solve finds it
wanting to split; without either the phase is refused at construction, which is
what makes a gap visible. The keyword `instances`, a `Dict` from a phase name to
a number or `:auto`, replaces what the file declares for the phases it names.

```toml
[[solid_solution]]
name        = "AFm_SO4_OH"
end_members = ["C4AH13", "monosulphate12"]
model       = "redlich_kister"
guggenheim  = "Lothenbach2019:AFm SO4/OH"
instances   = 2
```

A sublattice entry names the published model it follows with
`sublattice = "<literature key>:<model>"`, read by [`sublattice_model`](@ref),
which also checks each end-member's formula against the one the source prints:

```toml
[[solid_solution]]
name        = "CNASH_ss"
end_members = ["T2C-CNASHss", "T5C-CNASHss", "TobH-CNASHss",
               "5CA", "5CNA", "INFCA", "INFCN", "INFCNA"]
model       = "sublattice"
sublattice  = "Myers2014:cnash"
```

An entry whose end members exist in one database only names it with
`database = "<file>"`. Where they are absent, that database was not loaded,
which is expected rather than an anomaly: the entry is skipped without a warning,
whatever `skip_missing` says. `CSHQ_Cl`, whose chloride end member is in the
database `cemdata18-chloride.json` alone, is declared so.

# Example

```julia
substances = build_species(datapath("cemdata18-thermofun.json"))
dict       = Dict(symbol(s) => s for s in substances)
ss_phases  = build_solid_solutions(datapath("solid_solutions.toml"), dict)
# The file holds alternatives (`CSHQ` or `CNASH_ss`, `Ettringite_ss` or
# `AFt_SO4_CO3`): a system takes the phases it needs by name.
chosen = [p for p in ss_phases if name(p) in ("CSHQ", "C3(AF)S0.84H", "Ettringite_ss")]
cs = ChemicalSystem(species, CEMDATA_PRIMARIES; solid_solutions = chosen)
```
"""
function build_solid_solutions(
        toml_file::AbstractString,
        dict_species::AbstractDict;
        skip_missing::Bool = true,
        instances::AbstractDict = Dict{String, Any}(),
    )
    data = TOML.parsefile(resolve_data_path(toml_file))
    entries = get(data, "solid_solution", [])
    phases = SolidSolutionPhase[]
    for entry in entries
        ss_name = entry["name"]
        em_symbols = entry["end_members"]

        # Check all end-members are available
        missing_syms = filter(sym -> !haskey(dict_species, sym), em_symbols)
        if !isempty(missing_syms) && haskey(entry, "database")
            @debug "build_solid_solutions: \"$ss_name\" needs $(entry["database"])"
            continue
        elseif !isempty(missing_syms)
            missing_str = join(missing_syms, ", ")
            if skip_missing
                @warn "build_solid_solutions: skipping \"$ss_name\" — " *
                    "end-members not found in dict_species: $missing_str"
                continue
            else
                error(
                    "build_solid_solutions: end-members not found for \"$ss_name\": " *
                        missing_str,
                )
            end
        end

        em_species = [dict_species[sym] for sym in em_symbols]

        # Build mixing model
        model_str = get(entry, "model", "ideal")
        mixing_model = if model_str == "redlich_kister" && haskey(entry, "guggenheim")
            _guggenheim_model(entry["guggenheim"], ss_name)
        elseif model_str == "redlich_kister"
            RedlichKisterModel(;
                a0 = get(entry, "a0", 0.0),
                a1 = get(entry, "a1", 0.0),
                a2 = get(entry, "a2", 0.0),
            )
        elseif model_str == "regular"
            # `W` as a list of rows (a full symmetric matrix), or `w` as the
            # single interaction parameter of a binary.
            if haskey(entry, "W")
                RegularSolutionModel(reduce(vcat, permutedims.(entry["W"])))
            else
                w = float(get(entry, "w", 0.0))
                RegularSolutionModel([0.0 w; w 0.0])
            end
        elseif model_str == "sublattice"
            haskey(entry, "sublattice") || error(
                "build_solid_solutions: \"$ss_name\" is a sublattice model and names no " *
                    "`sublattice = \"<literature key>:<model>\"` to read it from."
            )
            sublattice_model(entry["sublattice"], em_species)
        elseif model_str == "compound_energy"
            haskey(entry, "sublattice") || error(
                "build_solid_solutions: \"$ss_name\" is a compound-energy model and names no " *
                    "`sublattice = \"<literature key>:<model>\"` to read it from."
            )
            compound_energy_model(entry["sublattice"], em_species)
        elseif model_str == "ideal"
            IdealSolidSolutionModel()
        else
            @warn "build_solid_solutions: unknown model \"$model_str\" for " *
                "\"$ss_name\", defaulting to ideal"
            IdealSolidSolutionModel()
        end

        declared_instances = get(instances, ss_name, get(entry, "instances", 1))
        n_instances = declared_instances in ("auto", :auto) ? :auto : Int(declared_instances)
        acknowledge_degenerate = Bool(get(entry, "acknowledge_degenerate", false))
        push!(
            phases,
            SolidSolutionPhase(ss_name, em_species; model = mixing_model, instances = n_instances, acknowledge_degenerate),
        )
    end
    return phases
end

# Published dimensionless Guggenheim parameters, `"<key>:<pair>"`, as the
# Redlich-Kister model at 298.15 K.
function _guggenheim_model(ref::AbstractString, ss_name)
    parts = split(ref, ":"; limit = 2)
    length(parts) == 2 || error(
        "build_solid_solutions: \"$ss_name\" names guggenheim = \"$ref\"; " *
            "expected \"<literature key>:<pair>\", such as \"Lothenbach2019:AFm SO4/OH\"."
    )
    p = literature_row(String(parts[1]), "guggenheim_parameters", String(parts[2]))
    RT = R_GAS * T_STANDARD
    return RedlichKisterModel(; a0 = ustrip(p.alpha0) * RT, a1 = ustrip(p.alpha1) * RT)
end

# ── Holland and Powell (1998): the volume of a solid and its Landau transition ──
#
# Two methods of aq17 carry a mineral away from its reference state beyond its
# heat-capacity polynomial [HollandPowell1998; p. 312](@cite). Its volume follows
# temperature and pressure,
#
#     V(T, P) = V₁(T) (1 + k′P/k(T))^(−1/k′),
#     V₁(T) = V° [1 + a° (T − Tr) − 2c a° (√T − √Tr)],   k(T) = k₂₉₈ [1 − b (T − Tr)],
#
# the thermal expansion a° (1 − c/√T) rising to a° at high temperature, and the
# Gibbs energy gains ∫ V dP from the reference pressure (`mv_eos_murnaghan_hp98`).
# A mineral with an order-disorder transition gains the excess energy of
# tricritical Landau theory (`landau_holland_powell98`),
#
#     G_ex = h° − T s° + ∫ v dP + S_max [(T − T_c) Q² + T_c Q⁶ / 3],
#     h° = S_max T°_c (Q₀² − Q₀⁶/3),  s° = S_max Q₀²,  Q₀⁴ = 1 − Tr/T°_c,
#     Q⁴ = 1 − T/T_c below T_c (0 above),  T_c = T°_c + (V_max/S_max)(P − Pr),
#     v = V_max Q₀² [1 + a° (T − Tr) − 2c a° (√T − √Tr)],
#
# whose Gibbs energy, enthalpy and entropy vanish at the reference state, so that
# the record's values hold there; its volume there is V_max Q₀⁶/3. The article's
# pressures are absolute, its ∫ v dP starting from zero and its T_c moving with P;
# both are taken here from the reference pressure P°, which keeps the record's
# Gibbs energy at P° and differs from the article by the volume of disorder times
# one bar, about 0.1 J/mol. The constants b, c and k′ of the equations are read
# from `data/literature/HollandPowell1998.json`.
const _HP98 = let q(name) = ustrip(literature_value("HollandPowell1998", name))
    (; b = q("bulk_modulus_temperature_coefficient"), c = q("thermal_expansion_coefficient"), kprime = q("bulk_modulus_pressure_derivative"))
end

_hp98_v1(V0, a0, Tr, T) = V0 * (1 + a0 * (T - Tr) - 2 * _HP98.c * a0 * (sqrt(T) - sqrt(Tr)))
_hp98_dv1(V0, a0, T) = V0 * a0 * (1 - _HP98.c / sqrt(T))

# ∫_{Pr}^{P} v (1 + k′P/k)^(−1/k′) dP and its derivatives in T and in P, for `v`
# and `k` the values at T and `dv`, `dk` their derivatives in T.
function _murnaghan_integral(v, dv, k, dk, P, Pr)
    kp = _HP98.kprime
    e = (kp - 1) / kp
    u, ur = 1 + kp * P / k, 1 + kp * Pr / k
    I = v * k / (kp - 1) * (u^e - ur^e)
    # The derivative of k u^e in T, u depending on T through k.
    dkue(uu, PP) = dk * (uu^e - (kp - 1) * PP * uu^(e - 1) / k)
    dIdT = dv * k / (kp - 1) * (u^e - ur^e) + v / (kp - 1) * (dkue(u, P) - dkue(ur, Pr))
    dIdP = v * u^(-1 / kp)
    return I, dIdT, dIdP
end

# The Gibbs energy the two methods add to a mineral's, with its derivatives in T
# and P; `p` holds V° (m³/mol), a° (1/K), k₂₉₈ (Pa), Tr (K), Pr (Pa) and, for a
# transition, `landau = (T°_c, S_max, V_max)` in K, J/(mol K), m³/mol.
function _hp98_excess(p, T, P)
    k = p.k0 * (1 - _HP98.b * (T - p.Tr))
    dk = -p.k0 * _HP98.b
    G, dGdT, dGdP = _murnaghan_integral(_hp98_v1(p.V0, p.a0, p.Tr, T), _hp98_dv1(p.V0, p.a0, T), k, dk, P, p.Pr)
    if p.landau !== nothing
        Tc0, Smax, Vmax = p.landau
        Q0² = sqrt(1 - p.Tr / Tc0)
        h0 = Smax * Tc0 * (Q0² - Q0²^3 / 3)
        s0 = Smax * Q0²
        I, dI_T, dI_P = _murnaghan_integral(
            _hp98_v1(Vmax * Q0², p.a0, p.Tr, T), _hp98_dv1(Vmax * Q0², p.a0, T), k, dk, P, p.Pr,
        )
        Tc = Tc0 + Vmax / Smax * (P - p.Pr)
        Q² = _plain(T) < _plain(Tc) ? sqrt(1 - T / Tc) : zero(T / Tc)
        # Q is that of equilibrium, ∂G/∂Q = 0: the derivatives are the explicit ones.
        G += h0 - T * s0 + I + Smax * ((T - Tc) * Q² + Tc * Q²^3 / 3)
        dGdT += -s0 + dI_T + Smax * Q²
        dGdP += dI_P + Vmax * (Q²^3 / 3 - Q²)
    end
    return G, dGdT, dGdP
end

# The parameters of the two methods from a record of aq17, `nothing` when it
# gives no bulk modulus: the bulk modulus in kbar under `m_compressibility` and
# the thermal expansion in 1/K under `m_expansivity` (their unit labels are each
# other's in the file), T°_c in °C, S_max in J/(mol K) and V_max in J/bar under
# `m_landau_phase_trans_props`.
function _hp98_parameters(row, landau, Tr, Pr)
    value(key) = hasproperty(row, key) ? extract_value(row, key; with_units = false) : missing
    k0, a0, V0 = value(:m_compressibility), value(:m_expansivity), value(:sm_volume)
    (any(ismissing, (k0, a0, V0)) || iszero(k0)) && return nothing
    J_bar = u"J/(bar*mol)"
    return (;
        V0 = ustrip(us"m^3/mol", V0 * J_bar), a0 = Float64(a0), k0 = 1000 * k0 * ustrip(us"Pa", 1.0u"bar"), Tr, Pr,
        landau = landau === nothing ? nothing :
            (_celsius_to_kelvin(landau[1]), landau[2], ustrip(us"m^3/mol", landau[3] * J_bar)),
    )
end

# Adds the two methods to the functions of a mineral built from its heat
# capacity: G, H, S and C_p gain the excess and its derivatives, and V is that of
# the equation of state.
function _holland_powell98!(s, p)
    complete_thermo_functions!(s)
    G0, H0, S0, Cp0 = s[:ΔₐG⁰], s[:ΔₐH⁰], s[:S⁰], s[:Cp⁰]
    refs = (T = s.Tref, P = s.Pref)
    d2(T, P) = ForwardDiff.derivative(t -> _hp98_excess(p, t, P)[2], T)
    s[:ΔₐG⁰] = NumericFunc((T, P) -> G0(T = T, P = P) + first(_hp98_excess(p, T, P)), (:T, :P), refs, u"J/mol")
    s[:ΔₐH⁰] = NumericFunc((T, P) -> (x = _hp98_excess(p, T, P); H0(T = T, P = P) + x[1] - T * x[2]), (:T, :P), refs, u"J/mol")
    s[:S⁰] = NumericFunc((T, P) -> S0(T = T, P = P) - _hp98_excess(p, T, P)[2], (:T, :P), refs, u"J/(mol*K)")
    s[:Cp⁰] = NumericFunc((T, P) -> Cp0(T = T, P = P) - T * d2(T, P), (:T, :P), refs, u"J/(mol*K)")
    s[:V⁰] = NumericFunc((T, P) -> _hp98_excess(p, T, P)[3], (:T, :P), refs, u"m^3/mol")
    return s
end

# ── Akinfiev and Diamond (2003): a dissolved gas ─────────────────────────────
#
# The chemical potential of a dissolved nonelectrolyte, from that of its ideal
# gas, of pure water and of three parameters ξ, a, b [AkinfievDiamond2003; Eq. 18](@cite):
#
#     μ_aq = μ_g(T) − RT ln N_w + (1 − ξ) RT ln f°_w + ξ RT ln(R T ρ°_w / M_w)
#            + RT ρ°_w [a + b (10³/T)^½],
#
# N_w the moles of water in a kilogram, f°_w the fugacity of pure water (bar),
# ρ°_w its density (g/cm³) and R in cm³ bar/(mol K) inside the logarithm. The
# fugacity is that of the database's own water: RT ln f°_w is the Gibbs energy of
# its solvent less that of its water vapor as an ideal gas at 1 bar. The record
# gives the Gibbs energy, enthalpy and entropy of the solute at the reference
# state, and the heat capacity of its ideal gas: μ_g(T) is the function of that
# heat capacity through the values the equation gives the gas there, so that the
# solute keeps its own. The excess the equation adds to the record's functions is
# therefore the term above less its value and its slope in T at the reference
# state.
const _G_CM3_PER_KG_M3 = ustrip(us"g/cm^3", 1.0u"kg/m^3")

function _ad03_term(p, w, T, P)
    ρ = _hgk_density(T, P) * _G_CM3_PER_KG_M3
    Rv = ustrip(us"cm^3*bar/(mol*K)", R_GAS_Q)
    lnf = (w.water(T = T, P = P) - w.steam(T = T, P = w.Pr)) / (R_GAS * T)
    return R_GAS * T * (-log(w.Nw) + (1 - p.ξ) * lnf + p.ξ * log(Rv * T * ρ / w.Mw) + ρ * (p.a + p.b * sqrt(1000 / T)))
end

# The water and the water vapor of the database of `df`, as `_ad03_term` uses
# them, or `nothing` when it has not both.
function _ad03_water(df)
    df === nothing && return nothing
    wrow = findfirst(==("H2O@"), df.symbol)
    srow = findfirst(eachrow(df)) do r
        extract_classification(get(r, :aggregate_state, missing), AS_UNDEF) == AS_GAS &&
            replace(String(r.formula), r"\|.*$" => "") == "H2O"
    end
    (wrow === nothing || srow === nothing) && return nothing
    water, steam = (complete_thermo_functions!(_species_from_row(df[r, :])) for r in (wrow, srow))
    M = water[:M]
    return (;
        water = water[:ΔₐG⁰], steam = steam[:ΔₐG⁰], Pr = ustrip(us"Pa", steam.Pref),
        Mw = ustrip(us"g/mol", M), Nw = ustrip(us"mol/kg", 1 / M),
    )
end

# Applies the equation of Akinfiev and Diamond to a dissolved gas, with the water
# of the database `df`; without it, the gas stays usable at its reference state
# only.
function _akinfiev_diamond!(s, df)
    w = _ad03_water(df)
    # Refused until applied; with another method refused too, it stays so.
    (w === nothing || get(properties(s), :refused_method, "") != "solute_aknifiev_diamond03") && return s
    delete!(s.properties, :refused_method)
    ξ, a, b = s[:ad03]
    p = (; ξ, a, b)
    Tr, Pr = ustrip(us"K", s.Tref), ustrip(us"Pa", s.Pref)
    gr = _ad03_term(p, w, Tr, Pr)
    sr = ForwardDiff.derivative(t -> _ad03_term(p, w, t, Pr), Tr)
    return _add_excess!(s, (T, P) -> _ad03_term(p, w, T, P) - gr - sr * (T - Tr))
end

# Adds to the functions of `s` built from its record a Gibbs energy `x(T, P)`
# that vanishes at its reference state, with the enthalpy, entropy, heat
# capacity and volume that follow from it by differentiation.
function _add_excess!(s, x)
    complete_thermo_functions!(s)
    G0, H0, S0, Cp0 = s[:ΔₐG⁰], s[:ΔₐH⁰], s[:S⁰], s[:Cp⁰]
    dT(T, P) = ForwardDiff.derivative(t -> x(t, P), T)
    d2T(T, P) = ForwardDiff.derivative(t -> dT(t, P), T)
    refs = (T = s.Tref, P = s.Pref)
    s[:ΔₐG⁰] = NumericFunc((T, P) -> G0(T = T, P = P) + x(T, P), (:T, :P), refs, u"J/mol")
    s[:ΔₐH⁰] = NumericFunc((T, P) -> H0(T = T, P = P) + x(T, P) - T * dT(T, P), (:T, :P), refs, u"J/mol")
    s[:S⁰] = NumericFunc((T, P) -> S0(T = T, P = P) - dT(T, P), (:T, :P), refs, u"J/(mol*K)")
    s[:Cp⁰] = NumericFunc((T, P) -> Cp0(T = T, P = P) - T * d2T(T, P), (:T, :P), refs, u"J/(mol*K)")
    s[:V⁰] = NumericFunc((T, P) -> ForwardDiff.derivative(q -> x(T, q), P), (:T, :P), refs, u"m^3/mol")
    return s
end
