# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using ChemistryLab
using DynamicQuantities
using ForwardDiff
using LinearAlgebra
using Test

@testsection "sensitivities keep the traces and the conservation laws" begin

    # The regression cases of the sensitivity correction PoroMechanics.jl
    # carried as a local patch. A trace species keeps the conservation law of
    # the component it carries however small its amount, and its curvature,
    # `1/n`, up to 1e300 beside a conservation block of order one, must not make
    # a pseudo-inverse discard that law.
    @testset "a trace carries its component" begin
        for trace in (1.0e-13, 1.0e-100, 1.0e-300)
            A = Matrix{Float64}(I, 2, 2)
            n = [1.0e4, trace]
            db = [0.0, 1.0]
            dn = ChemistryLab._equilibrium_sensitivity(A, diagm(1 ./ n), zeros(2), db, n)
            @test dn ≈ db atol = 1.0e-10
            # Redundant conservation rows leave the multipliers free, not the
            # species' response.
            redundant = [1.0 0.0; 0.0 1.0; 0.0 2.0]
            dn = ChemistryLab._equilibrium_sensitivity(redundant, diagm(1 ./ n), zeros(2), redundant * db, n)
            @test redundant * dn ≈ redundant * db atol = 1.0e-10
        end
    end

    @testset "an absent pure phase is held, and an impossible perturbation raises" begin
        dn = ChemistryLab._equilibrium_sensitivity(
            [1.0 1.0], diagm([1.0, 0.0]), zeros(2), [1.0], [1.0, 0.0]; pinned = [false, true],
        )
        @test dn ≈ [1.0, 0.0]
        # The perturbation asks the held phase to move: no sensitivity solves it.
        @test_throws ErrorException ChemistryLab._equilibrium_sensitivity(
            Matrix{Float64}(I, 2, 2), diagm([1.0, 0.0]), zeros(2), [0.0, 1.0], [1.0, 0.0];
            pinned = [false, true],
        )
    end

    @testset "the derivative of a certified equilibrium through dual element amounts" begin
        # The call PoroMechanics.jl makes per cell: a state seeded with the dual
        # element type, the element amounts as dual numbers, one certified solve.
        # Calcite in water with CO2 present at trace level, differentiated with
        # respect to the calcium of the budget, against a centered difference.
        db = Dict(
            symbol(s) => s for s in build_species(
                    datapath("slop98-inorganic-thermofun.json"); verbose = false
                )
        )
        names = split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Cal")
        cs = ChemicalSystem([db[s] for s in names], ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"])
        st = ChemicalState(cs)
        set_quantity!(st, "H2O@", 1.0u"kg")
        set_quantity!(st, "Cal", 1.0e-3u"mol")
        A = Float64.(cs.SM.A)
        b0 = A * Float64[ustrip(us"mol", x) for x in st.n]
        ica = findfirst(p -> symbol(p) == "Ca+2", cs.SM.primaries)
        ical = findfirst(==("Cal"), symbol.(cs.species))
        e = zeros(length(b0)); e[ica] = 1.0
        function calcite(x)
            b = b0 .+ x .* e
            seed = [ustrip(us"mol", v) + zero(x) for v in st.n]
            eq, _ = equilibrate_certified(ChemicalState(cs, seed .* u"mol"); b)
            return ustrip(us"mol", eq.n[ical])
        end
        d_ad = ForwardDiff.derivative(calcite, 0.0)
        h = 1.0e-7
        d_fd = (calcite(h) - calcite(-h)) / 2h
        @test d_ad ≈ d_fd rtol = 1.0e-5
    end
end
