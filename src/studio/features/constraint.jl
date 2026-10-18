# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2026 Jean-François Barthélémy (Cerema, UMR MCD)

# What the equilibrium holds fixed besides the composition: temperature and
# pressure by default, or a pH, an activity, saturated curing, an adiabatic
# cell. ChemistryLab's constraints (src/equilibrium/constraints.jl). A sealed
# volume is not offered: ChemistryLab refuses it on a condensed paste, whose
# volume barely depends on the pressure, and `porosity` already reports the
# volume change of sealed curing.

"""
    ConstraintChoice

The constraint of the equilibrium solved from the recipe.
"""
abstract type ConstraintChoice <: Choice end

"Temperature and pressure fixed: the usual equilibrium."
struct FixedTPChoice <: ConstraintChoice end
kind(::Type{FixedTPChoice}) = "fixed_tp"
fields(::Type{FixedTPChoice}) = Field[]

"A prescribed pH, held by adding acid or base."
struct FixedpHChoice <: ConstraintChoice
    pH::Float64
end
kind(::Type{FixedpHChoice}) = "fixed_ph"
fields(::Type{FixedpHChoice}) = [Field("pH", :number; default = 12.5, min = 0, max = 14)]

"A prescribed activity of one species, held by adding it (CO₂ of the air, say)."
struct FixedActivityChoice <: ConstraintChoice
    species::String
    activity::Float64
end
kind(::Type{FixedActivityChoice}) = "fixed_activity"
fields(::Type{FixedActivityChoice}) = [
    Field("species", :text; default = "CO2@"),
    Field("activity", :number; default = 1.0e-5, min = 0),
]

"Cured under water: the paste takes up the water its reactions shrink by."
struct SaturatedCuringChoice <: ConstraintChoice end
kind(::Type{SaturatedCuringChoice}) = "saturated_curing"
fields(::Type{SaturatedCuringChoice}) = Field[]

"An adiabatic cell: the enthalpy holds and the temperature follows."
struct AdiabaticChoice <: ConstraintChoice end
kind(::Type{AdiabaticChoice}) = "adiabatic"
fields(::Type{AdiabaticChoice}) = Field[]

members(::Type{ConstraintChoice}) =
    (FixedTPChoice, FixedpHChoice, FixedActivityChoice, SaturatedCuringChoice, AdiabaticChoice)

"""
    constraint_code(c) -> (lines, keyword)

The lines a constraint needs before the solve, and the keyword it adds to it
(empty for fixed temperature and pressure).
"""
constraint_code(::FixedTPChoice) = (String[], "")
constraint_code(c::FixedpHChoice) = (String[], "constraint = FixedpH($(julia_number(c.pH)))")
constraint_code(c::FixedActivityChoice) =
    (String[], "constraint = FixedActivity($(repr(c.species)), $(julia_number(c.activity)))")
constraint_code(::AdiabaticChoice) = (String[], "constraint = Adiabatic()")
constraint_code(::SaturatedCuringChoice) = (
    ["# Saturated curing refers to the fresh paste, before any reaction.", "fresh = budget(recipe, cs).state"],
    "constraint = SaturatedCuring(; reference = fresh)",
)

constraint_problems(c::ConstraintChoice) = bound_problems(c, "the constraint")
function constraint_problems(c::FixedActivityChoice)
    out = bound_problems(c, "the constraint")
    c.activity > 0 || push!(out, Problem("activity_not_positive", "The activity of $(c.species) must be positive.", Dict{String, Any}("species" => c.species)))
    isempty(strip(c.species)) && push!(out, Problem("no_species", "Name the species whose activity is held."))
    return out
end
