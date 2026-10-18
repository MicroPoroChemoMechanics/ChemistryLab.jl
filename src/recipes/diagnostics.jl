# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# The checks a system and a recipe deserve before a calculation: what a choice
# of species leaves out, or brings in twice, or uses outside its data. Each is
# something a calculation would otherwise do silently — set a constituent
# aside, leave a volume or a heat incomplete, read a species beyond the range
# its data cover — and each says what to do about it.

"""
    SystemDiagnostic(level, code, message, fix, subjects)

One finding of [`diagnose_system`](@ref): its `level` (`:blocking`, the
calculation cannot be trusted; `:warning`, a result will be incomplete or
doubtful; `:info`), a `code` naming the check, a `message` saying what was
found and what it costs, a `fix` saying what to do, and the `subjects`, the
species or constituents concerned.
"""
struct SystemDiagnostic
    level::Symbol
    code::Symbol
    message::String
    fix::String
    subjects::Vector{String}
end

function Base.show(io::IO, d::SystemDiagnostic)
    print(io, d.level, " (", d.code, "): ", d.message)
    isempty(d.fix) || print(io, " ", d.fix)
    return nothing
end

# Cemdata18 (Lothenbach et al. 2019, Table 1, note *): the Al-Si hydrogarnets
# "precipitate very slowly at 20 °C" and are "generally not included in
# calculations".
const _SLOW_HYDROGARNETS = ("C3AS0.41H5.18", "C3AS0.84H4.32")

# The temperature (K) below which the note on the hydrogarnets is reported as a
# warning rather than as information: the note concerns ambient conditions.
const _AMBIENT_LIMIT = 303.15

"""
    diagnose_system(cs; recipe = nothing, T = nothing, t = nothing, model = nothing,
                    reactants = String[], large = 250) -> Vector{SystemDiagnostic}

The findings worth reading before computing with system `cs`, and with
`recipe` (at age `t`, for laws of the degree of reaction that depend on it):

  - a constituent of the recipe that the system cannot hold (no species carries
    it, or one of its elements), set aside unreacted;
  - ettringite declared both as a pure phase and as a member of a solid
    solution, one phase counted twice;
  - a crystal used outside the temperature range its data declare (the
    temperature is `T`, or the recipe's), reported as information for the
    reactants, which only dissolve: the minerals of the recipe, and the species
    named in `reactants` (those of a phase list, `phase_list(name).reactants`);
  - the Al-Si hydrogarnets of Cemdata18, which it notes precipitate very slowly
    at 20 °C and are generally left out;
  - crystals without a standard molar volume (volumes and porosity incomplete)
    or without an enthalpy of formation (heat incomplete);
  - with a SIT `model`, the pairs of ions without an interaction coefficient;
  - a system of more than `large` species, whose first solve compiles longer.

An empty result means that none of these was found, not that the calculation
is right.

```julia
cs = phase_list_system("Portland paste (Lothenbach and Winnefeld 2006)", substances)
foreach(println, diagnose_system(cs; recipe))
```
"""
function diagnose_system(
        cs::ChemicalSystem; recipe = nothing, T = nothing, t = nothing, model = nothing,
        reactants = String[], large::Integer = 250,
    )
    out = SystemDiagnostic[]
    symbols = Set(symbol(s) for s in cs.species)
    crystals = [cs.species[i] for i in cs.idx_crystal]
    # the species that only dissolve: the minerals a recipe brings, and those named
    dissolving = Set(vcat(_recipe_minerals(recipe), String.(reactants)))

    aside = _set_aside(recipe, cs, t)
    isempty(aside) || push!(
        out, SystemDiagnostic(
            :warning, :not_in_system,
            "$(join(aside, ", ")) cannot be held by the system: set aside, unreacted.",
            "Add the species (or the species carrying their elements) to the system, or accept that they stay inert.",
            aside,
        )
    )

    if "ettringite" in symbols && "ettringite03_ss" in symbols
        push!(
            out, SystemDiagnostic(
                :warning, :declared_twice,
                "ettringite is declared both as a pure phase and, as ettringite03_ss (a third of it), in a solid solution: one phase counted twice.",
                "Keep one of the two.", ["ettringite", "ettringite03_ss"],
            )
        )
    end

    θ = _diagnosis_temperature(T, recipe)
    append!(out, _temperature_findings(θ, crystals, dissolving))

    slow = [h for h in _SLOW_HYDROGARNETS if h in symbols]
    if !isempty(slow)
        level = _hydrogarnet_level(θ)
        push!(
            out, SystemDiagnostic(
                level, :slow_hydrogarnet,
                "$(join(slow, " and ")): Cemdata18 notes that they precipitate very slowly at 20 °C and are generally not included in calculations (Lothenbach et al. 2019, Table 1).",
                "Remove them for a paste at ambient temperature.", slow,
            )
        )
    end

    no_volume = [symbol(s) for s in crystals if !_has_molar_volume(s)]
    isempty(no_volume) || push!(
        out, SystemDiagnostic(
            :warning, :missing_volume,
            "$(join(no_volume, ", ")) $(length(no_volume) == 1 ? "has" : "have") no standard molar volume: volumes and porosity leave them out.",
            "Use another database for these phases, or read the volumes as incomplete.", no_volume,
        )
    )

    no_enthalpy = [symbol(s) for s in crystals if _std_property(s, :ΔₐH⁰) === nothing]
    isempty(no_enthalpy) || push!(
        out, SystemDiagnostic(
            :warning, :missing_enthalpy,
            "$(join(no_enthalpy, ", ")) $(length(no_enthalpy) == 1 ? "has" : "have") no enthalpy of formation: a heat leaves them out.",
            "Read heats as incomplete, or use another database for these phases.", no_enthalpy,
        )
    )

    append!(out, _model_findings(model, cs))

    n = length(cs.species)
    n > large && push!(
        out, SystemDiagnostic(
            :info, :large_system, "The system holds $n species: the first equilibrium of a session takes longer to compile.",
            "Leave out the phases that cannot form in the material.", String[],
        )
    )
    return out
end

# The mineral species a recipe brings (none without a recipe).
_recipe_minerals(::Nothing) = String[]
_recipe_minerals(recipe::Recipe) =
    [symbol(c.species) for (m, _) in _material_masses(recipe) for c in m.constituents if c isa MineralConstituent]

# The constituents a recipe brings that the system cannot hold, at age `t`.
_set_aside(::Nothing, cs, t) = String[]
_set_aside(recipe::Recipe, cs, t) = unique(r.constituent for r in budget(recipe, cs; t).residual if r.reason === :not_in_system)

# The temperature the data are judged at: the one given, else the recipe's.
_diagnosis_temperature(T, recipe) = _kelvin(T)
_diagnosis_temperature(::Nothing, recipe::Recipe) = recipe.T
_diagnosis_temperature(::Nothing, ::Nothing) = nothing

# The crystals read outside the temperature range of their data, at `θ` (none
# without a temperature): a warning for those that may form, information for the
# reactants, which only dissolve. A comparison of values: numbers only.
_temperature_findings(::Nothing, crystals, dissolving) = SystemDiagnostic[]
function _temperature_findings(θ, crystals, dissolving)
    out = SystemDiagnostic[]
    outside = [symbol(s) for s in crystals if !(first(temperature_range(s)) <= θ <= last(temperature_range(s)))]
    formed = filter(∉(dissolving), outside)
    brought = filter(in(dissolving), outside)
    isempty(formed) || push!(
        out, SystemDiagnostic(
            :warning, :temperature_range,
            "$(join(formed, ", ")) $(length(formed) == 1 ? "is" : "are") used outside the temperature range $(length(formed) == 1 ? "its" : "their") data declare.",
            "Read the results for these phases with care, or leave them out at this temperature.", formed,
        )
    )
    isempty(brought) || push!(
        out, SystemDiagnostic(
            :info, :reactant_temperature_range,
            "The reactants $(join(brought, ", ")) are read outside the temperature range their data declare; they dissolve into the paste rather than form.",
            "", brought,
        )
    )
    return out
end

# The note on the Al-Si hydrogarnets concerns ambient conditions: a warning there
# (or when no temperature is known), information above.
_hydrogarnet_level(::Nothing) = :warning
_hydrogarnet_level(θ) = θ <= _AMBIENT_LIMIT ? :warning : :info

# What an activity model lacks for the ions of the system: the pairs of a SIT
# model without a coefficient; nothing to say of the others.
_model_findings(model, cs) = SystemDiagnostic[]
function _model_findings(model::SITActivityModel, cs)
    pairs = missing_epsilon_pairs(cs, model)
    isempty(pairs) && return SystemDiagnostic[]
    return [
        SystemDiagnostic(
            :warning, :missing_sit_pairs,
            "$(length(pairs)) pairs of ions have no SIT coefficient, such as $(join(("$(a)–$(b)" for (a, b) in first(pairs, 3)), ", ")): they are taken as not interacting.",
            "Use an activity model whose parameters cover these ions.", ["$(a)–$(b)" for (a, b) in pairs],
        ),
    ]
end
