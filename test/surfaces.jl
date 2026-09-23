# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using ChemistryLab
using DynamicQuantities
using ForwardDiff
using Test

@testsection "Surface area models" begin

    # ── The numbers that were already shipped must not move ───────────────────
    #
    # `blaine_factor` is now `area_ratio` between two `BlaineSurfaceArea`, which
    # is the same division. These are the exact values the previous
    # implementation was pinned to, kept exact rather than approximate: a ratio
    # of two floats that divide evenly has no rounding to allow for.
    @testset "blaine_factor is unchanged, to the bit" begin
        @test blaine_factor(385u"m^2/kg") == 1.0
        @test blaine_factor(770u"m^2/kg") == 2.0
        @test blaine_factor(462u"m^2/kg") ≈ 1.2 rtol = 1.0e-12
        @test blaine_factor(2000u"m^2/kg"; blaine_ref = 400u"m^2/kg") == 5.0
        @test blaine_factor(385.0) == 1.0                 # plain Real, in m²/kg
        # and it is literally the ratio the new abstraction computes
        @test blaine_factor(462u"m^2/kg") ==
            area_ratio(BlaineSurfaceArea(462u"m^2/kg"), BlaineSurfaceArea(385.0u"m^2/kg"))
    end

    # ── The guard, in both directions ─────────────────────────────────────────
    @testset "two measurements of area do not compare" begin
        @test area_ratio(BlaineSurfaceArea(770.0), BlaineSurfaceArea(385.0)) == 2.0
        @test area_ratio(BETSurfaceArea(180.0), BETSurfaceArea(90.0)) == 2.0
        @test area_ratio(GeometricSurfaceArea(2.0), GeometricSurfaceArea(1.0)) == 2.0

        # The pair that has cost real hydration rates: silica fume's BET area
        # standing in for its effective Blaine fineness, a factor of ten.
        @test_throws ArgumentError area_ratio(
            BETSurfaceArea(20_000.0), BlaineSurfaceArea(385.0)
        )
        @test_throws ArgumentError area_ratio(
            BlaineSurfaceArea(385.0), BETSurfaceArea(20_000.0)
        )
        @test_throws ArgumentError area_ratio(
            GeometricSurfaceArea(1.0), BETSurfaceArea(90.0)
        )
        # and therefore through the front door as well
        @test_throws ArgumentError blaine_factor(BETSurfaceArea(20_000.0))

        # A typed Blaine area still goes through unchanged.
        @test blaine_factor(BlaineSurfaceArea(770.0)) == 2.0

        # The message names both measurements, which is what makes it useful.
        msg = try
            area_ratio(BETSurfaceArea(1.0), BlaineSurfaceArea(1.0))
            ""
        catch e
            sprint(showerror, e)
        end
        @test occursin("BET", msg) && occursin("Blaine", msg)
    end

    # ── Units are converted, or refused; never stripped in silence ────────────
    @testset "units convert or raise" begin
        @test FixedSurfaceArea(500.0u"cm^2").A ≈ 0.05 rtol = 1.0e-12
        @test BETSurfaceArea(0.09u"m^2/g").A_specific ≈ 90.0 rtol = 1.0e-12
        @test specific_area(BlaineSurfaceArea(380u"m^2/kg")) ≈ 380.0 rtol = 1.0e-12
        @test FixedSurfaceArea(0.5).A == 0.5           # plain Real is SI already

        # The case `safe_ustrip` would have swallowed: a specific area handed to
        # a constructor expecting a total one.
        @test_throws ArgumentError FixedSurfaceArea(20_000.0u"m^2/kg")
        @test_throws ArgumentError BETSurfaceArea(0.5u"m^2")
        @test_throws ArgumentError BlaineSurfaceArea(1.0u"K")

        # A total area has no mass to divide by, and says so.
        @test_throws MethodError specific_area(FixedSurfaceArea(0.5))

        @test area_method(BETSurfaceArea(1.0)) === :BET
        @test area_method(BlaineSurfaceArea(1.0)) === :Blaine
        @test area_method(GeometricSurfaceArea(1.0)) === :geometric
        @test area_method(FixedSurfaceArea(1.0)) === :fixed
    end

    # ── The element type is kept, so AD goes through ──────────────────────────
    @testset "constructors keep their element type" begin
        d = ForwardDiff.Dual(90.0, 1.0)
        @test BETSurfaceArea(d).A_specific === d
        @test FixedSurfaceArea(ForwardDiff.Dual(0.5, 1.0)).A isa ForwardDiff.Dual
        # The generated-constructor trap: Julia's automatic
        # `BETSurfaceArea(x::T) where T` is more specific than an untyped outer
        # constructor, so without an inner one this returned
        # `BETSurfaceArea{Int64}(90)` with the unit check never reached.
        @test BETSurfaceArea(90).A_specific === 90.0
        @test BETSurfaceArea(90) isa BETSurfaceArea{Float64}
        @test FixedSurfaceArea(1) isa FixedSurfaceArea{Float64}
        @test ShrinkingCoreArea(BETSurfaceArea(90.0), 1) isa
            ShrinkingCoreArea{BETSurfaceArea{Float64}, Float64}
    end

    # ── total_area, and the three-argument form kept for callers ──────────────
    @testset "total_area and surface_area agree where they must" begin
        fixed = FixedSurfaceArea(0.5)
        @test total_area(fixed, 0.01, 0.02, 0.1) == 0.5
        @test surface_area(fixed, 0.01, 0.1) == 0.5
        @test surface_area(fixed, 0.0, 0.1) == 0.5

        bet = BETSurfaceArea(90.0)
        @test total_area(bet, 0.01, 0.02, 0.1) ≈ 90.0 * 0.01 * 0.1 rtol = 1.0e-12
        @test surface_area(bet, 0.01, 0.1) ≈ 90.0 * 0.01 * 0.1 rtol = 1.0e-12
        @test surface_area(bet, 0.0, 0.1) == 0.0
        # An amount an ODE step pushed slightly negative gives no area, not a
        # negative one.
        @test surface_area(bet, -1.0e-12, 0.1) == 0.0
    end

    # ── A grain that shrinks ──────────────────────────────────────────────────
    @testset "ShrinkingCoreArea" begin
        n₀, M = 0.05, 0.1
        A₀ = 90.0 * n₀ * M

        # p = 1 is the mass-proportional law, i.e. exactly BET.
        linear = ShrinkingCoreArea(BETSurfaceArea(90.0); exponent = 1)
        for n in (0.0, 1.0e-6, 0.01, 0.05)
            @test total_area(linear, n, n₀, M) ≈ surface_area(BETSurfaceArea(90.0), n, M) rtol =
                1.0e-12
        end

        # p = 2/3 is the spherical geometry, recovered away from exhaustion. The
        # tolerance is the regularization's own bias, SHRINK_FLOOR/(3f), not a
        # convenience: at f = 1e-3 it is 3.3e-6.
        sphere = ShrinkingCoreArea(BETSurfaceArea(90.0))
        @test sphere.exponent ≈ 2 / 3
        for f in (1.0, 0.5, 0.1, 1.0e-3)
            expected = A₀ * f^(2 / 3)
            @test total_area(sphere, f * n₀, n₀, M) ≈ expected rtol =
                max(1.0e-12, SHRINK_FLOOR / (3 * f))
        end

        # The two fixed points are exact, not approximate.
        @test total_area(sphere, n₀, n₀, M) ≈ A₀ rtol = 1.0e-14
        @test total_area(sphere, 0.0, n₀, M) == 0.0
        @test total_area(sphere, -1.0e-12, n₀, M) == 0.0

        # It decreases, and it stays above the mass-proportional law: two thirds
        # of an exponent means area survives depletion longer than mass does.
        areas = [total_area(sphere, f * n₀, n₀, M) for f in 1.0:-0.1:0.1]
        @test issorted(areas; rev = true)
        @test all(
            total_area(sphere, f * n₀, n₀, M) > total_area(linear, f * n₀, n₀, M)
                for f in (0.9, 0.5, 0.1)
        )

        # Refused rather than wrong: no initial amount to be a fraction of, and
        # no way to answer the three-argument form honestly.
        @test_throws ArgumentError total_area(sphere, 0.0, 0.0, M)
        @test_throws ArgumentError surface_area(sphere, 0.01, M)

        @test area_method(sphere) === :BET
    end

    # ── Automatic differentiation, including where it used to be infinite ─────
    @testset "AD through every area model" begin
        n₀, M = 0.05, 0.1
        bet = BETSurfaceArea(90.0)
        sphere = ShrinkingCoreArea(bet)

        dbet = ForwardDiff.derivative(n -> total_area(bet, n, n₀, M), 0.01)
        @test dbet ≈ 90.0 * M rtol = 1.0e-12

        dsph = ForwardDiff.derivative(n -> total_area(sphere, n, n₀, M), 0.01)
        @test isfinite(dsph) && dsph > 0

        # The point of SHRINK_FLOOR: f^p has an infinite slope at exhaustion, and
        # an integrator does reach it.
        d0 = ForwardDiff.derivative(n -> total_area(sphere, n, n₀, M), 0.0)
        @test isfinite(d0)
        @test !isnan(d0)

        # Differentiating with respect to the area parameter, not only the amount.
        dA = ForwardDiff.derivative(a -> total_area(BETSurfaceArea(a), 0.01, n₀, M), 90.0)
        @test dA ≈ 0.01 * M rtol = 1.0e-12

        dfac = ForwardDiff.derivative(b -> blaine_factor(b), 385.0)
        @test dfac ≈ 1 / 385.0 rtol = 1.0e-12
    end

    # ── Surface: one support, one area ────────────────────────────────────────
    @testset "Surface names its host once" begin
        s = Surface("calcite", "Cal", BETSurfaceArea(90.0))
        @test s.name == "calcite"
        @test s.host == "Cal"
        @test area_method(s) === :BET

        prescribed = Surface("inert sorbent", FixedSurfaceArea(0.5))
        @test prescribed.host === nothing
        @test occursin("prescribed support", sprint(show, prescribed))
        @test occursin("calcite", sprint(show, s))
    end

end
