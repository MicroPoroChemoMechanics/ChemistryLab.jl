# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using ChemistryLab
using DynamicQuantities
using ForwardDiff
using Test

# Size classes of a ground mineral. The references are the closed forms of a
# lognormal distribution by volume: its specific surface as spheres,
# (6/ρ) e^{−ln d50 + σ²/2}, and the invariance of every class surface under a
# change of scale of the diameters.

@testsection "Size classes of a ground mineral" begin
    subs = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-thermofun.json"); verbose = false))
    c3s = subs["C3S"]
    M = ustrip(us"kg/mol", c3s[:M])
    ρ = ChemistryLab._species_density(c3s, 298.15, 1.0e5)
    @test 2500 < ρ < 3500                                # kg/m³, a clinker phase
    psd = LogNormalPSD(; d50 = 10u"μm", d90 = 30u"μm")
    # Φ⁻¹(0.9) = 1.2815515655446004, the 90th percentile of the standard normal.
    @test psd.σ ≈ log(3) / 1.2815515655446004 rtol = 1.0e-14
    whole = 6 / ρ * exp(-log(10.0e-6) + psd.σ^2 / 2)    # m²/kg

    @testset "the classes add up to the powder" begin
        for K in (1, 2, 5, 12)
            cl = size_classes(c3s, psd, K; mass = 1u"kg")
            @test length(cl) == K
            @test sum(c.amount for c in cl) ≈ 1 / M rtol = 1.0e-14
            @test all(c -> c.fraction ≈ 1 / K, cl)
            # Each class holds 1/K of the mass, so the surfaces add with weights 1/K.
            @test sum(ChemistryLab.specific_area(c.surface.initial) for c in cl) / K ≈ whole rtol = 1.0e-12
            # Finer classes first, each its own species.
            @test issorted([c.diameter for c in cl])
            @test allunique(symbol(c.species) for c in cl)
            @test cl[1].species === c3s
            @test all(c -> c.surface isa ShrinkingCoreArea, cl)
        end
        # One class is the mineral with the whole surface: the calculation of
        # one species, unchanged.
        one_ = only(size_classes(c3s, psd, 1; amount = 2.0u"mol"))
        @test one_.amount == 2.0
        @test ChemistryLab.specific_area(one_.surface.initial) ≈ whole rtol = 1.0e-14
        @test one_.diameter ≈ 6 / (ρ * whole) rtol = 1.0e-14
    end

    @testset "a change of scale of the diameters" begin
        # With σ held, every class surface scales as 1/d50: ∂s/∂d50 = −s/d50.
        s(d50) = [
            ChemistryLab.specific_area(c.surface.initial)
                for c in size_classes(c3s, LogNormalPSD(; d50, σ = psd.σ), 4; mass = 1.0, density = ρ)
        ]
        d = 10.0e-6
        @test ForwardDiff.derivative(s, d) ≈ -s(d) ./ d rtol = 1.0e-10
        # And in σ, against the derivative of the closed form for one class.
        w(σ) = ChemistryLab.specific_area(only(size_classes(c3s, LogNormalPSD(; d50 = d, σ), 1; mass = 1.0, density = ρ)).surface.initial)
        @test ForwardDiff.derivative(w, psd.σ) ≈ psd.σ * 6 / ρ * exp(-log(d) + psd.σ^2 / 2) rtol = 1.0e-10
    end

    @testset "a measured curve, sampled on the lognormal" begin
        # F linear in ln d between 400 points, from far below to far above:
        # the same classes, to the discretization.
        ds = exp.(range(log(10.0e-6) - 6psd.σ, log(10.0e-6) + 6psd.σ; length = 400))
        F = [ChemistryLab._normal_mass(-Inf, (log(x) - log(10.0e-6)) / psd.σ) for x in ds]
        F .-= F[1]
        tab = TabulatedPSD(ds .* u"m", 100 .* F)        # in percent
        for K in (1, 4)
            a = size_classes(c3s, psd, K; mass = 1.0, density = ρ)
            b = size_classes(c3s, tab, K; mass = 1.0, density = ρ)
            for (x, y) in zip(a, b)
                @test ChemistryLab.specific_area(y.surface.initial) ≈ ChemistryLab.specific_area(x.surface.initial) rtol = 2.0e-3
            end
        end
        # Exact for the curve itself: the classes add up to its own integral.
        b = size_classes(c3s, tab, 7; mass = 1.0, density = ρ)
        @test sum(ChemistryLab.specific_area(c.surface.initial) for c in b) / 7 ≈
            6 / ρ * inverse_diameter_moment(tab, 0.0, 1.0) rtol = 1.0e-12
    end

    @testset "what a distribution refuses" begin
        @test_throws ArgumentError LogNormalPSD(; d50 = 10u"μm")
        @test_throws ArgumentError LogNormalPSD(; d50 = 10u"μm", d90 = 5u"μm")
        @test_throws ArgumentError TabulatedPSD([1.0, 1.0] .* u"μm", [0.0, 1.0])
        @test_throws "start at zero" TabulatedPSD([1.0, 2.0] .* u"μm", [0.1, 1.0])
        @test_throws ArgumentError size_classes(c3s, psd, 0; mass = 1.0)
        @test_throws ArgumentError size_classes(c3s, psd, 2)
        @test_throws "molar volume" size_classes(Species("CaO"; aggregate_state = AS_CRYSTAL), psd, 2; mass = 1.0)
    end

    @testset "copies are found by their own symbol" begin
        # A copy shares its formula; the kinetic index takes symbols first, and a
        # formula only where it names one species.
        cal = subs["Cal"]
        a = with_symbol(cal, "CaCO3")       # a species whose symbol is its formula
        b = with_symbol(cal, "CaCO3_2")
        cs = ChemicalSystem([subs["H2O@"], a, b], ["H2O@", "CaCO3"])
        index = ChemistryLab._species_index(cs)
        @test index["CaCO3"] == 2
        @test index["CaCO3_2"] == 3
    end
end
