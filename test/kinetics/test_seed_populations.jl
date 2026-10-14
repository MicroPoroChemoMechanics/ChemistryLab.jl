# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using ChemistryLab
using DynamicQuantities
using ForwardDiff
using OrderedCollections
using OrdinaryDiffEq
using Test

include("../reference_species.jl")

# A phase that precipitates from its own seeds: the geometry of the seeds, the
# nucleation laws, and portlandite growing out of a supersaturated solution. The
# references are the formulas of the docstrings, written out again here.

const _SEED_SUBS = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-thermofun.json"); verbose = false))

# Portlandite and its reservoir of seeds in water, with the activities of a
# state whose ln Ω is shifted by `x` through Ca+2.
function _seed_system()
    h2o, hp, oh, ca = reference_species(("H2O@", "H+", "OH-", "Ca+2"))
    ch = _SEED_SUBS["Portlandite"]
    nuc = with_symbol(ch, "CH_nuc")
    cs = ChemicalSystem(AbstractSpecies[h2o, hp, oh, ca, ch, nuc], AbstractSpecies[h2o, hp, ca])
    rxn(s) = Reaction(OrderedDict(cs[s] => 1.0), OrderedDict(cs["Ca+2"] => 1.0, cs["OH-"] => 2.0); symbol = "$s dissolution")
    index = Dict(symbol(s) => i for (i, s) in enumerate(cs.species))
    sat = saturation_ratio(cs, rxn("Portlandite"))
    lna0 = zeros(length(cs.species))
    lna0[index["Ca+2"]] = -log(sat(298.15, 1.0e5, StateView(lna0, index)))
    lna(x) = (v = Vector{typeof(x)}(lna0); v[index["Ca+2"]] += x; StateView(v, index))
    return cs, rxn, index, lna
end

@testsection "Seeds, nucleation and growth" begin
    cs, rxn, index, lna = _seed_system()
    V_m = ustrip(us"m^3/mol", _SEED_SUBS["Portlandite"][:V⁰](T = 298.15u"K", P = 1.0e5u"Pa"; unit = true))
    amounts(nuc, ch) = StateView([symbol(s) == "CH_nuc" ? nuc : symbol(s) == "Portlandite" ? ch : 1.0 for s in cs.species], index)

    @testset "the area of spheres" begin
        r0 = 50.0e-9
        g = SphericalSeeds(50u"nm")
        @test seed_volume(g) ≈ 4π / 3 * r0^3 rtol = 1.0e-14
        area = SeededGrowthArea(cs["CH_nuc"], g)
        @test area.V_m ≈ V_m rtol = 1.0e-14
        A(nuc, ch) = ChemistryLab._reactive_area(area, 298.15, amounts(nuc, ch), amounts(nuc, ch), "Portlandite", 0.0)
        # The seeds alone: N spheres of radius r₀, N·4πr₀².
        N = 1.0e-6 * V_m / seed_volume(g)
        @test A(1.0e-6, 0.0) ≈ N * 4π * r0^2 rtol = 1.0e-7
        # Grown: N equal spheres holding the whole volume, (36πN)^{1/3} V^{2/3}.
        ref(nuc, ch) = cbrt(36π * nuc * V_m / seed_volume(g)) * ((nuc + ch) * V_m)^(2 / 3)
        @test A(1.0e-6, 3.0e-3) ≈ ref(1.0e-6, 3.0e-3) rtol = 1.0e-12
        @test A(1.0e-12, 3.0e-3) ≈ ref(1.0e-12, 3.0e-3) rtol = 1.0e-12   # few seeds, much grown
        @test A(0.0, 3.0e-3) == 0                       # no seeds, no growth
        @test A(0.0, 0.0) == 0
        # And its derivatives, against the formula.
        @test ForwardDiff.gradient(x -> A(x...), [1.0e-6, 3.0e-3]) ≈ ForwardDiff.gradient(x -> ref(x...), [1.0e-6, 3.0e-3]) rtol = 1.0e-12
        @test_throws "not of a site budget" total_area(area, 1.0, 1.0, 0.1)
    end

    @testset "the area of needles, up to their length" begin
        l0, a, lmax, w = 20.0e-9, 0.1, 400.0e-9, 1.0e-3
        g = NeedleSeeds(; length = 20u"nm", aspect = a, max_length = 400u"nm", width = w)
        @test seed_volume(g) ≈ π * a^2 * l0^3 rtol = 1.0e-14
        area = SeededGrowthArea(cs["CH_nuc"], g)
        A(nuc, ch) = ChemistryLab._reactive_area(area, 298.15, amounts(nuc, ch), amounts(nuc, ch), "Portlandite", 0.0)
        len(nuc, ch) = cbrt((nuc + ch) * V_m / (nuc * V_m / seed_volume(g) * π * a^2))
        cut(z) = z > 0 ? z^2 / (z + w) : 0.0
        ref(nuc, ch) = 3 * (nuc + ch) * V_m / len(nuc, ch) * cut(1 - len(nuc, ch) / lmax)
        for ch in (0.0, 1.0e-5, 1.0e-4)
            @test A(1.0e-7, ch) ≈ ref(1.0e-7, ch) rtol = 1.0e-12
        end
        # At their full length the needles stop growing.
        full = 1.0e-7 * ((lmax / l0)^3 - 1)
        @test len(1.0e-7, full) ≈ lmax rtol = 1.0e-12
        @test A(1.0e-7, full) < 1.0e-20 * A(1.0e-7, 0.0)   # zero but for the rounding of l
        @test A(1.0e-7, 1.1 * full) == 0
        @test_throws ArgumentError NeedleSeeds(; length = 20u"nm", aspect = a, max_length = 10u"nm")
    end

    @testset "the nucleation laws" begin
        cl = ClassicalNucleation(1.0e12, 2.0)
        J(lnΩ) = ChemistryLab._nucleation_frequency(cl, lnΩ)
        @test J(log(2.0)) ≈ 1.0e12 rtol = 1.0e-14
        @test J(log(3.0)) ≈ 1.0e12 * exp(2 / log(2)^2 - 2 / log(3)^2) rtol = 1.0e-14
        @test J(0.0) == 0 && J(-0.1) == 0 && J(1.0e-3) == 0
        # Every derivative vanishes at saturation: this one, at a hundredth.
        @test ForwardDiff.derivative(J, 0.05) ≈ 1.0e12 * exp(2 / log(2)^2 - 2 / 0.05^2) * 4 / 0.05^3 rtol = 1.0e-10
        th = ThresholdNucleation(1.0e10, 2.5; width = 1.0e-3)
        Jt(lnΩ) = ChemistryLab._nucleation_frequency(th, lnΩ)
        @test Jt(log(2.4)) == 0
        z = 3.0 / 2.5 - 1
        @test Jt(log(3.0)) ≈ 1.0e10 * z^2 / (z + 1.0e-3) rtol = 1.0e-12
        @test ForwardDiff.derivative(Jt, log(2.5)) == 0
        @test_throws ArgumentError ClassicalNucleation(1.0, 0.0)
        @test_throws ArgumentError ThresholdNucleation(1.0, -1.0)
    end

    @testset "the reservoir fills at the frequency times the seed volume" begin
        g = SphericalSeeds(50u"nm")
        law = ClassicalNucleation(1.0e12, 2.0)
        n = amounts(0.0, 1.0)                          # one mole of the substrate
        for (on, E) in ((FixedExtent(2.0), 2.0), (SubstrateAmount("Portlandite"), 1.0))
            r = nucleation_rate(law, cs, rxn("CH_nuc"), g; on)
            @test r isa KineticFunc
            @test r(298.15, 1.0e5, 0.0, n, lna(log(3.0)), n) ≈
                -1.0e12 * exp(2 / log(2)^2 - 2 / log(3)^2) * E * seed_volume(g) / V_m rtol = 1.0e-10
            @test r(298.15, 1.0e5, 0.0, n, lna(-0.5), n) == 0
        end
        # A saturation given by the caller, here twice the phase's.
        sat = saturation_ratio(cs, rxn("Portlandite"))
        r2 = nucleation_rate(law, cs, rxn("CH_nuc"), g; on = FixedExtent(1.0), saturation = (T, P, l) -> 2 * sat(T, P, l))
        @test r2(298.15, 1.0e5, 0.0, n, lna(log(1.5)), n) ≈ -1.0e12 * exp(2 / log(2)^2 - 2 / log(3)^2) * seed_volume(g) / V_m rtol = 1.0e-10
        @test_throws "not a species" nucleation_rate(law, cs, rxn("CH_nuc"), g; on = SubstrateAmount("C3S"))
    end

    @testset "portlandite out of a supersaturated solution" begin
        g = SphericalSeeds(50u"nm")
        k = arrhenius_rate_constant(1.0e-6, 0.0)
        grow = growth_only(
            transition_state([RateMechanism(k, 1.0, 1.0)], cs, rxn("Portlandite"), SeededGrowthArea(cs["CH_nuc"], g)),
            cs, rxn("Portlandite"),
        )
        nucleate = nucleation_rate(ClassicalNucleation(1.0e10, 2.0), cs, rxn("CH_nuc"), g; on = FixedExtent(1.0))
        krs = [KineticReaction(cs, rxn("CH_nuc"), nucleate), KineticReaction(cs, rxn("Portlandite"), grow)]
        function run(ca)
            st = ChemicalState(cs)
            set_quantity!(st, "H2O@", moles_of_water() * u"mol")
            set_quantity!(st, "Ca+2", ca * u"mol")
            set_quantity!(st, "OH-", 2ca * u"mol")
            kp = KineticsProblem(cs, krs, st, (0.0, 3600.0); equilibrium_solver = EquilibriumSolver(cs, DiluteSolutionModel(), OptimaOptimizer()))
            sol = integrate(kp, KineticsSolver(; ode_solver = Rodas5P(), reltol = 1.0e-8, abstol = 1.0e-18, saveat = [0.0, 600.0, 3600.0]))
            return st, speciated_states(sol, kp)
        end
        mol(s, x) = ustrip(us"mol", s.n[index[x]])
        # Undersaturated: nothing nucleates, nothing grows.
        st, states = run(0.005)
        @test all(s -> mol(s, "CH_nuc") + mol(s, "Portlandite") < 1.0e-20, states)
        # Supersaturated: seeds form, the phase grows on them, and the solution
        # falls back toward saturation; the calcium is all accounted for.
        st, states = run(0.04)
        sat = saturation_ratio(cs, rxn("Portlandite"))
        lnaof(s) = (l = log_activities(s, DiluteSolutionModel()); StateView([l[symbol(sp)] for sp in cs.species], index))
        Ω = [sat(298.15, 1.0e5, lnaof(s)) for s in states]
        @test Ω[1] > 2
        @test Ω[end] < Ω[2] < Ω[1]
        @test mol(states[end], "CH_nuc") > 0 && mol(states[end], "Portlandite") > mol(states[end], "CH_nuc")
        ca(s) = mol(s, "Ca+2") + mol(s, "CH_nuc") + mol(s, "Portlandite")
        @test maximum(abs(ca(s) - ca(st)) for s in states) < 1.0e-12
    end
end
