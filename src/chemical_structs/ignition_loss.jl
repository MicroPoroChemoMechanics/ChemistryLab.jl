# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using DynamicQuantities

# ── What a solid assemblage would lose on heating ────────────────────────────
#
# Thermogravimetry is the second observable `hydration_calibration.jl` asks for
# by name, and the reason it is worth having is that calorimetry constrains
# three combinations of six kinetic parameters and no more: a measurement that
# sees the PHASES rather than the heat breaks correlations heat cannot.
#
# Its total is computable from formulas alone, which is what is here. Its
# *curve* is not: turning a mass into a thermogram needs to know which phase
# releases what between which temperatures, and those windows are literature
# values, not consequences of the formulas. They are named as the missing input
# rather than guessed at.

"""
    bound_water(state::ChemicalState) -> Quantity

The water a solid assemblage would release on ignition, in kilograms, from the
hydrogen written into its phases.

Every hydrogen in a solid leaves as water — `Ca(OH)₂ → CaO + H₂O`,
`CaSO₄·2H₂O → CaSO₄ + 2H₂O` — so the amount is `H/2` per formula unit,
whatever form the hydrogen is written in. That is why this counts hydrogen and
not formula water: `Ca(OH)₂` has no `H₂O` in its formula and loses one on
ignition, and a rule that looked for `H₂O` would report zero.

**The aqueous phase is not included.** Pore solution is water, and evaporable
water is not what a hydrate bound — the distinction
[the water budget](@ref sec-theory-water-budget) is about. Only species in the
solid buckets are counted.

See also: [`ignition_loss`](@ref), [`bound_water_per_phase`](@ref).
"""
bound_water(state::ChemicalState) = ignition_loss(state).water

"""
    ignition_loss(state::ChemicalState) -> NamedTuple

What a solid assemblage would lose on ignition, as
`(; water, carbon_dioxide, total)` in kilograms.

Hydrogen leaves as water and carbon as carbon dioxide, which is the convention
a loss-on-ignition measurement reports and which is exact for the carbonates and
hydrates a cement or a clay assemblage is made of.

!!! warning "Carbon is assumed to be carbonate"
    `C` in a solid is counted as `CO₂`. For a cement, a clay or a soil mineral
    assemblage that is what it is. For organic carbon it is not — the mass is
    right and the species is not, and a sample with organic matter needs the
    split made before this is read.

# What this is and is not

It is the **total** a thermogram integrates to, computed from the formulas the
database gives. It is not a thermogram: which phase releases what between which
temperatures is a set of decomposition windows, and those are literature values
rather than consequences of a formula. This package does not carry them, and
inventing them would be fabricating the part of the measurement that does the
identifying.

# Example

```julia
loss = ignition_loss(state)
ustrip(us"g", loss.water), ustrip(us"g", loss.total)
```
"""
function ignition_loss(state::ChemicalState)
    system = state.system
    idx = _solid_indices(system)
    w = zero(0.0u"mol")
    c = zero(0.0u"mol")
    for i in idx
        sp = system.species[i]
        a = atoms(sp)
        n = state.n[i]
        h = get(a, :H, 0)
        cc = get(a, :C, 0)
        iszero(h) || (w += (h / 2) * n)
        iszero(cc) || (c += cc * n)
    end
    # The molar masses come from the species, never from a table written here.
    mw = _ignition_molar_mass(system, "H2O@", "H2O", _WATER_ATOMS)
    mc = _ignition_molar_mass(system, "CO2@", "CO2", _CO2_ATOMS)
    water = uconvert(us"kg", w * mw)
    co2 = uconvert(us"kg", c * mc)
    return (; water, carbon_dioxide = co2, total = water + co2)
end

"""
    _ignition_molar_mass(system, symbols..., atoms) -> Quantity

The molar mass of a released gas, taken from the system's own species when it
carries one under any of `symbols`, and computed from `atoms` otherwise.

Preferring the system's own is the rule this package keeps everywhere: a molar
mass written here would drift from the database the rest of the calculation
uses. The fallback exists because an assemblage may release a species it never
declared — a paste with no gas phase still loses water on ignition — and it is
weighed with the same atomic masses as every species.
"""
function _ignition_molar_mass(system::ChemicalSystem, a, b, atoms)
    for sp in system.species
        symbol(sp) in (a, b) && haskey(sp, :M) && return sp[:M]
    end
    return calculate_molar_mass(atoms)
end

const _WATER_ATOMS = OrderedDict(:H => 2, :O => 1)
const _CO2_ATOMS = OrderedDict(:C => 1, :O => 2)

"""
    bound_water_per_phase(state::ChemicalState) -> Vector{Pair{String,Quantity}}

The ignition water each solid species carries, largest first.

What a thermogram would resolve if the decomposition windows were known — so
this is the quantity a TGA observation operator needs per phase, and the place
those windows would attach.
"""
function bound_water_per_phase(state::ChemicalState)
    system = state.system
    mw = _ignition_molar_mass(system, "H2O@", "H2O", _WATER_ATOMS)
    out = Pair{String, typeof(uconvert(us"kg", 1.0u"mol" * mw))}[]
    for i in _solid_indices(system)
        sp = system.species[i]
        h = get(atoms(sp), :H, 0)
        iszero(h) && continue
        m = uconvert(us"kg", (h / 2) * state.n[i] * mw)
        ustrip(us"kg", m) > 0 || continue
        push!(out, symbol(sp) => m)
    end
    return sort(out; by = p -> -ustrip(us"kg", p.second))
end
