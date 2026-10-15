# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using ChemistryLab
using DynamicQuantities
using ForwardDiff
using OrderedCollections
using OrdinaryDiffEq
using Test

include("../reference_species.jl")

# A layer covering a dissolving mineral, counted in the sites it blocks. The
# references are the formulas of the docstrings and the logistic solution of the
# opening at a fixed saturation.

const _COVER_SUBS = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-thermofun.json"); verbose = false))

function _cover_system()
    h2o, hp, oh, ca = reference_species(("H2O@", "H+", "OH-", "Ca+2"))
    c3s = _COVER_SUBS["C3S"]
    layer = with_symbol(_COVER_SUBS["Portlandite"], "CH_layer")
    cs = ChemicalSystem(AbstractSpecies[h2o, hp, oh, ca, c3s, layer], AbstractSpecies[h2o, hp, ca, c3s])
    rxn = Reaction(OrderedDict(cs["CH_layer"] => 1.0), OrderedDict(cs["Ca+2"] => 1.0, cs["OH-"] => 2.0); symbol = "layer")
    index = Dict(symbol(s) => i for (i, s) in enumerate(cs.species))
    sat = saturation_ratio(cs, rxn)
    lna0 = zeros(length(cs.species))
    lna0[index["Ca+2"]] = -log(sat(298.15, 1.0e5, StateView(lna0, index)))
    lna(x) = (v = Vector{typeof(x)}(lna0); v[index["Ca+2"]] += x; StateView(v, index))
    return cs, rxn, index, lna
end

@testsection "A layer covering a surface" begin
    cs, rxn, index, lna = _cover_system()
    M = ustrip(us"kg/mol", cs["C3S"][:M])
    Γ, ssa = 1.0e-5, 400.0
    layer = CoveringLayer("CH_layer", "C3S", BETSurfaceArea(ssa * u"m^2/kg"), Γ * u"mol/m^2"; blockable = 0.99)
    N(n_s) = Γ * ssa * M * n_s                        # sites of the substrate, mol
    state(nℓ, n_s = 1.0) = StateView([symbol(s) == "CH_layer" ? nℓ : symbol(s) == "C3S" ? n_s : 1.0 for s in cs.species], index)

    @testset "the free fraction multiplies the substrate's rate" begin
        r0 = 1.0e-6
        law = covered((T, P, t, n, lna, n0) -> r0, cs, layer)
        for θ in (0.0, 0.3, 0.9)
            n = state(θ * N(1.0))
            @test law(298.15, 1.0e5, 0.0, n, lna(0.0), n) ≈ r0 * (1 - θ) atol = 2.0e-6 * r0
        end
        # Never below its floor, and continuous there.
        n = state(N(1.0))
        @test law(298.15, 1.0e5, 0.0, n, lna(0.0), n) ≈ r0 * 1.0e-6 rtol = 1.0e-12
        @test_throws "not a species" covered((T, P, t, n, lna, n0) -> r0, cs, CoveringLayer("X", "C3S", BETSurfaceArea(1.0), 1.0))
        @test_throws ArgumentError CoveringLayer("CH_layer", "C3S", BETSurfaceArea(1.0), 1.0; blockable = 1.5)
    end

    @testset "it forms on the blockable sites while its phase is stable" begin
        k, τ = 1.0e-3, 100.0
        r = covering_rate(layer, cs, rxn; k_form = k, τ_open = τ)
        x = log(3.0)
        n = state(0.2 * N(1.0))
        @test r(298.15, 1.0e5, 0.0, n, lna(x), n) ≈ -k * (0.99 - 0.2) * N(1.0) * (3 - 1) / (3 + 1) rtol = 1.0e-10
        # Saturation of the blockable sites: the amount tends to b·N and stays below.
        f(u, p, t) = [-r(298.15, 1.0e5, t, state(u[1]), lna(x), state(u[1]))]
        sol = solve(ODEProblem(f, [0.0], (0.0, 50 / (k / 2))), Tsit5(); reltol = 1.0e-10, abstol = 1.0e-20)
        @test sol.u[end][1] ≈ 0.99 * N(1.0) rtol = 1.0e-8
        @test all(u -> u[1] <= 0.99 * N(1.0) * (1 + 1.0e-12), sol.u)
    end

    @testset "it opens from the free sites once its phase is not stable" begin
        k, τ = 1.0e-3, 100.0
        r = covering_rate(layer, cs, rxn; k_form = k, τ_open = τ)
        x = log(0.4)                                    # Ω = 0.4
        a = (1 - 0.4) / τ
        θ0 = 0.9
        f(u, p, t) = [-r(298.15, 1.0e5, t, state(u[1]), lna(x), state(u[1]))]
        ts = range(0.0, 2000.0; length = 9)
        sol = solve(ODEProblem(f, [θ0 * N(1.0)], (0.0, 2000.0)), Tsit5(); reltol = 1.0e-11, abstol = 1.0e-22, saveat = ts)
        θ(t) = θ0 * exp(-a * t) / (1 - θ0 + θ0 * exp(-a * t))
        @test [u[1] / N(1.0) for u in sol.u] ≈ θ.(ts) rtol = 1.0e-7
        # Slowed by a factor of the run, here one half.
        half = covering_rate(layer, cs, rxn; k_form = k, τ_open = τ, opening = (T, P, t, n, lna, n0) -> 0.5)
        n = state(0.5 * N(1.0))
        @test half(298.15, 1.0e5, 0.0, n, lna(x), n) ≈ r(298.15, 1.0e5, 0.0, n, lna(x), n) / 2 rtol = 1.0e-14
    end

    @testset "differentiated, against the formula" begin
        r = covering_rate(layer, cs, rxn; k_form = 1.0e-3, τ_open = 100.0)
        # Opening, at Ω = 0.4: r = (n/τ)(1 − n/N)(1 − Ω), dr/dn = (1 − 2n/N)(1 − Ω)/τ.
        d = ForwardDiff.derivative(nℓ -> r(298.15, 1.0e5, 0.0, state(nℓ), lna(log(0.4)), state(1.0)), 0.3 * N(1.0))
        @test d ≈ (1 - 0.6) * 0.6 / 100.0 rtol = 1.0e-10
        # Forming, at Ω = 3: dr/dn = k (Ω − 1)/(Ω + 1).
        d = ForwardDiff.derivative(nℓ -> r(298.15, 1.0e5, 0.0, state(nℓ), lna(log(3.0)), state(1.0)), 0.3 * N(1.0))
        @test d ≈ 1.0e-3 * 0.5 rtol = 1.0e-10
        @test_throws ArgumentError covering_rate(layer, cs, rxn; k_form = 1.0, τ_open = 0.0)
    end
end
