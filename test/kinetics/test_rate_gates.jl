# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using ChemistryLab
using DynamicQuantities
using ForwardDiff
using OrderedCollections
using OrdinaryDiffEq
using Test

include("../reference_species.jl")

# Laws restricted to one direction, `dissolution_only` and `growth_only`: zero on
# the forbidden side of saturation, continuous with their slope at Ω = 1, the
# law itself away from it. The references are the closed form of the factor.

const _GATE_SUBS = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-thermofun.json"); verbose = false))

# Calcite in water, its dissolution, and the activities of a state whose ln Ω is
# shifted by `x` through the activity of Ca+2.
function _gate_calcite()
    sp = speciation(
        collect(values(_GATE_SUBS)), ["Cal", "H2O@"];
        aggregate_state = [AS_AQUEOUS], exclude_species = split("H2@ O2@ CH4@"),
    )
    cs = ChemicalSystem(sp, CEMDATA_PRIMARIES)
    rxn = Reaction(
        OrderedDict(cs["Cal"] => 1.0), OrderedDict(cs["Ca+2"] => 1.0, cs["CO3-2"] => 1.0);
        symbol = "calcite dissolution",
    )
    index = Dict(symbol(s) => i for (i, s) in enumerate(cs.species))
    sat = saturation_ratio(cs, rxn)
    lna0 = zeros(length(cs.species))
    # At saturation: ln a(Ca+2) chosen so that ln Ω = 0 with every other at zero.
    lna0[index["Ca+2"]] = -log(sat(298.15, 1.0e5, StateView(lna0, index)))
    lna(x) = (v = Vector{typeof(x)}(lna0); v[index["Ca+2"]] += x; StateView(v, index))
    return cs, rxn, index, lna, sat
end

@testsection "Laws that run one way" begin
    cs, rxn, index, lna, sat = _gate_calcite()
    n = StateView(fill(1.0, length(cs.species)), index)
    r0 = 1.0e-6
    law = (T, P, t, n, lna, n0) -> r0 * (1 - sat(T, P, lna))
    r(y) = law(298.15, 1.0e5, 0.0, n, lna(y), n)
    w = 1.0e-3
    diss = dissolution_only(law, cs, rxn; width = w)
    grow = growth_only(law, cs, rxn; width = w)
    @test diss isa KineticFunc && grow isa KineticFunc
    @test log(sat(298.15, 1.0e5, lna(0.0))) ≈ 0 atol = 1.0e-12
    s(x) = x > 0 ? x^2 / (x^2 + w^2) : 0.0          # the factor, in its driving force

    @testset "zero on the forbidden side, the law away from it" begin
        for x in (-0.5, -0.01, -2w)                 # undersaturated: ln Ω = x, to rounding
            lnΩ = log(sat(298.15, 1.0e5, lna(x)))
            @test diss(298.15, 1.0e5, 0.0, n, lna(x), n) ≈ law(298.15, 1.0e5, 0.0, n, lna(x), n) * s(-lnΩ) *
                (1 / (1 + ChemistryLab.SHRINK_FLOOR)) rtol = 1.0e-12
            @test grow(298.15, 1.0e5, 0.0, n, lna(x), n) == 0
        end
        for x in (0.0, 2w, 0.5)                     # saturated or supersaturated
            @test diss(298.15, 1.0e5, 0.0, n, lna(x), n) == 0
        end
        @test grow(298.15, 1.0e5, 0.0, n, lna(0.5), n) ≈ law(298.15, 1.0e5, 0.0, n, lna(0.5), n) * s(0.5) rtol = 1.0e-12
        # Far from saturation, the law itself to within (w/ln Ω)².
        @test diss(298.15, 1.0e5, 0.0, n, lna(-0.5), n) / law(298.15, 1.0e5, 0.0, n, lna(-0.5), n) ≈ 1 atol = 1.0e-5
    end

    @testset "continuous with its slope at saturation" begin
        d(x) = ForwardDiff.derivative(y -> diss(298.15, 1.0e5, 0.0, n, lna(y), n), x)
        @test d(0.0) == 0
        # The slope fades to nothing as the solution approaches saturation.
        @test abs(d(-1.0e-9)) < 1.0e-10 * abs(d(-w))
        # Against the closed form, away from zero: dr/dx = r′ s(−x) − r s′(−x).
        x = -w
        ds(y) = 2y * w^2 / (y^2 + w^2)^2
        g = 1 / (1 + ChemistryLab.SHRINK_FLOOR)
        @test d(x) ≈ g * (ForwardDiff.derivative(r, x) * s(-x) - r(x) * ds(-x)) rtol = 1.0e-10
    end

    @testset "a dissolving phase stops when it is gone" begin
        n0 = StateView(fill(1.0, length(cs.species)), index)
        empty_ = StateView([symbol(sp) == "Cal" ? 0.0 : 1.0 for sp in cs.species], index)
        @test diss(298.15, 1.0e5, 0.0, empty_, lna(-0.5), n0) == 0
        half = StateView([symbol(sp) == "Cal" ? 0.5 : 1.0 for sp in cs.species], index)
        @test diss(298.15, 1.0e5, 0.0, half, lna(-0.5), n0) ≈ r(-0.5) * s(0.5) * 0.5 / (0.5 + ChemistryLab.SHRINK_FLOOR) rtol = 1.0e-12
        # Even from nothing at all, no NaN.
        @test diss(298.15, 1.0e5, 0.0, empty_, lna(-0.5), empty_) == 0
        @test_throws ArgumentError dissolution_only(law, cs, rxn; width = 0.0)
    end

    @testset "gypsum in a supersaturated solution does not move" begin
        sp = speciation(
            collect(values(_GATE_SUBS)), ["Gp", "H2O@"];
            aggregate_state = [AS_AQUEOUS], exclude_species = split("H2@ O2@ H2S@ HS- S-2 HSO3- SO3-2 S2O3-2"),
        )
        cs2 = ChemicalSystem(sp, CEMDATA_PRIMARIES)
        rxn2 = Reaction(
            OrderedDict(cs2["Gp"] => 1.0), OrderedDict(cs2["Ca+2"] => 1.0, cs2["SO4-2"] => 1.0, cs2["H2O@"] => 2.0);
            symbol = "gypsum dissolution",
        )
        st = ChemicalState(cs2)
        set_quantity!(st, "H2O@", moles_of_water() * u"mol")
        set_quantity!(st, "Gp", 0.1u"mol")
        set_quantity!(st, "Ca+2", 0.04u"mol")
        set_quantity!(st, "SO4-2", 0.04u"mol")
        plain = first_order_rate(arrhenius_rate_constant(1.0e-6, 0.0), cs2, rxn2, FixedSurfaceArea(1.0))
        amounts = map((plain, dissolution_only(plain, cs2, rxn2))) do law
            kp = KineticsProblem(
                cs2, [KineticReaction(cs2, rxn2, law)], st, (0.0, 3600.0);
                equilibrium_solver = EquilibriumSolver(cs2, DiluteSolutionModel(), OptimaOptimizer()),
            )
            sol = integrate(kp, KineticsSolver(; ode_solver = Rodas5P(), reltol = 1.0e-10, abstol = 1.0e-16, saveat = [0.0, 3600.0]))
            ustrip(us"mol", last(speciated_states(sol, kp)).n[findfirst(==("Gp"), symbol.(cs2.species))])
        end
        @test amounts[1] > 0.1 + 1.0e-6               # the law itself grows it back
        @test amounts[2] == 0.1                       # restricted, it does not move
    end
end
