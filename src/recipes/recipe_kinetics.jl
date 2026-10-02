# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# ── A recipe as a kinetic problem ────────────────────────────────────────────
#
# The same recipe drives `hydrate`, where the extents are imposed, and
# `integrate`, where rate laws decide them: the constituents given a rate are the
# kinetic species, everything else is taken as the recipe says at the start.

"""
    KineticsProblem(recipe::Recipe, system::ChemicalSystem, rates, tspan; kwargs...)
        -> KineticsProblem

The hydration of `recipe` in `system` as a kinetic problem. `rates` maps the
name of a mineral constituent to its rate law (a [`KineticFunc`](@ref)): that
constituent enters whole and unreacted, and dissolves into the primaries of
`system` at that rate, by the reaction the conservation matrix of `system`
gives it. The other keywords go to the [`KineticsProblem`](@ref) of a list of
reactions (`activity_model`, `equilibrium_solver`, `calorimeter`).

Every other constituent is taken as [`budget`](@ref) takes it at `tspan[1]`: its
reacted part is in the equilibrium from the start, and its unreacted part stays
aside for the whole run, out of the problem. The reacted part of a constituent
known by its oxides enters as the primaries that carry its elements, with the
protons an oxide consumes written as hydroxide taken from the mixing water
(`CaO + H₂O → Ca²⁺ + 2 OH⁻`); its heat of dissolution is therefore not part of the
heat of the run.

A constituent given a rate must be a mineral constituent of the recipe whose
species is in `system`; a glass, known by its oxides only, has no formula to
dissolve and is refused.
"""
function KineticsProblem(
        recipe::Recipe, cs::ChemicalSystem, rates::AbstractDict, tspan::Tuple; kwargs...,
    )
    t0 = first(tspan)
    kinetic = Dict{String, Any}()        # constituent name => (material, constituent)
    for (m, _) in _material_masses(recipe), c in m.constituents
        haskey(rates, c.name) || continue
        c isa MineralConstituent || throw(
            ArgumentError(
                "KineticsProblem: $(c.name) of $(m.name) is known by its oxides only; it has no " *
                    "formula to dissolve, so it cannot be given a rate."
            )
        )
        haskey(kinetic, c.name) && throw(
            ArgumentError("KineticsProblem: $(c.name) is a constituent of more than one material; give each its own name.")
        )
        haskey(cs.dict_species, symbol(c.species)) || throw(
            ArgumentError("KineticsProblem: $(symbol(c.species)), the species of $(c.name), is not in the system.")
        )
        kinetic[c.name] = (m, c)
    end
    missing_ = setdiff(keys(rates), keys(kinetic))
    isempty(missing_) || throw(
        ArgumentError("KineticsProblem: no constituent of the recipe is named " * join(sort!(collect(missing_)), ", ") * ".")
    )

    bud = budget(recipe, cs; t = t0)
    A = Float64.(cs.SM.A)
    state = bud.state
    oxides = bud.b .- A * ustrip.(us"mol", state.n)
    # A kinetic constituent enters whole: its reacted part at `t0` is taken back
    # out of the equilibrium and the whole of it put in as the species.
    masses = Dict(m => x for (m, x) in _material_masses(recipe))
    for (_, (m, c)) in kinetic
        sym = symbol(c.species)
        mc = masses[m] * c.mass_fraction
        α = effective_extent(m, c, t0)
        M = ustrip(us"g/mol", c.species[:M])
        set_quantity!(state, sym, moles(state, sym) + (1 - α) * mc / M * u"mol")
    end
    _add_primaries!(state, oxides)

    reactions = KineticReaction[]
    for (name, (_, c)) in kinetic
        j = findfirst(s -> symbol(s) == symbol(c.species), cs.species)
        products = [cs.dict_species[symbol(prim)] for (q, prim) in enumerate(cs.SM.primaries) if !iszero(A[q, j])]
        rxn = Reaction([cs.species[j]], products; symbol = "$(symbol(c.species)) dissolution")
        rxn[:rate] = rates[name]
        push!(reactions, KineticReaction(cs, rxn))
    end
    return KineticsProblem(cs, reactions, state, tspan; kwargs...)
end

# The element budget `b` (over the primaries of the system) added to `state` as
# amounts of the primary species, a primary's column being a unit vector. Two
# negative amounts have a reading: protons, written as hydroxide taken from the
# water (a basic oxide, `K₂O + H₂O → 2 K⁺ + 2 OH⁻`), and water itself, taken from
# the mixing water (an acidic one, `SO₃ + H₂O → SO₄²⁻ + 2 H⁺`). Any other is
# refused, and so is a budget that would leave no water.
function _add_primaries!(state::ChemicalState, b::AbstractVector)
    cs = state.system
    prim = [symbol(p) for p in cs.SM.primaries]
    for (p, x) in zip(prim, b)
        abs(x) <= 1.0e-14 * max(1.0, maximum(abs, b)) && continue
        if x < 0 && p == "H+" && haskey(cs.dict_species, "OH-") && "H2O@" in prim
            set_quantity!(state, "OH-", moles(state, "OH-") - x * u"mol")
            set_quantity!(state, "H2O@", moles(state, "H2O@") + x * u"mol")
        elseif x >= 0 || p == "H2O@"
            set_quantity!(state, p, moles(state, p) + x * u"mol")
        else
            throw(ArgumentError("KineticsProblem: the oxides of the recipe need a negative amount of $p, which no species of the system can carry."))
        end
    end
    ustrip(us"mol", moles(state, "H2O@")) > 0 || throw(
        ArgumentError("KineticsProblem: the oxides of the recipe take more water than the recipe holds.")
    )
    return state
end
