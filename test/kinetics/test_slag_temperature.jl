# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# The clinker and the slag of Snellings et al. (2022) at 5, 20 and 40 °C
# (`scripts/snellings2022_kinetics.jl`): the isothermal integral of the Waller
# law the fits use checked against the rate of `waller` itself, and the numbers
# the page `examples/slag_temperature.md` states pinned.

isdefined(@__MODULE__, :sn22_slag_fit) || include(joinpath(pkgdir(ChemistryLab), "scripts", "snellings2022_kinetics.jl"))

@testsection "a slag-limestone cement at three temperatures" begin
    slag, clinker = sn22_degrees("slag"), sn22_degrees("clinker")
    # Fig. 6: six ages of three temperatures and three w/b, for each constituent.
    for d in (slag, clinker)
        @test length(d.age) == 54
        @test all(count((d.temperature_C .== T) .& (d.w_b .== w)) == 6 for T in (5, 20, 40), w in (0.4, 0.5, 0.6))
    end

    # The closed form is the integral of the package's law: at a constant
    # temperature the rate of `waller`, at the degree the closed form gives,
    # is the closed form's derivative, ceiling and Arrhenius factor included.
    θ = (τ = 5.0, n = 0.6, Ea = 67.0, α_max = 0.65)
    law = sn22_slag_law(θ)
    idx = Dict("slag" => 1)
    for T_C in (5.0, 20.0, 40.0), t in (0.3, 7.0, 90.0)
        α(s) = only(sn22_waller_degree(Tuple(θ), T_C, (s,))) / 100
        rate = law(ustrip(us"K", _sn22_kelvin(T_C)), 1.0e5, t * 86400, StateView([1 - α(t)], idx), StateView([0.0], idx), StateView([1.0], idx))
        @test rate ≈ ForwardDiff.derivative(α, t) / 86400 rtol = 1.0e-12
    end

    # The clinker under the published law: the water factor does not act at
    # these water/cement ratios, and the page's numbers at one day.
    @test sn22_clinker_degree(20, 0.4, (1, 7)) == sn22_clinker_degree(20, 0.6, (1, 7))
    one_day = [only(sn22_clinker_degree(T, 0.5, (1,))) for T in (5, 20, 40)]
    @test one_day ≈ [26.4956, 42.3543, 59.2364] atol = 1.0e-3

    # The slag: the fly-ash shape, then one ceiling per w/b.
    @test [sn22_slag_fit((w,); shared = (:τ,)).rms for w in (0.4, 0.5, 0.6)] ≈ [17.7728, 15.4653, 14.9238] atol = 1.0e-3
    @test sn22_slag_fit(; shared = (:τ,), per_wb = (:α_max,), τ₀ = 10.0).rms ≈ 3.51 atol = 0.01
    @test sn22_slag_fit(; shared = (:τ, :n), per_wb = (:α_max,), τ₀ = 10.0).rms ≈ 3.14 atol = 0.01
    fit = sn22_slag_fit(; per_wb = (:α_max,), τ₀ = 10.0)
    @test fit.rms ≈ 2.59 atol = 0.01
    @test fit.values ≈ [5.339, 0.5906, 67.4, 0.5412, 0.653, 0.7049] rtol = 1.0e-3
    @test fit.identifiability.rank == 6
    @test fit.identifiability.stderr[[1, 3]] ≈ [0.33, 0.17] atol = 0.005
    @test fit.θ[0.4].Ea == fit.θ[0.6].Ea

    # The law of w/b 0.4 on the round robin's slags.
    dz = sn22_durdzinski(fit.θ[0.4])
    largest(m, lab) = (k = (dz.material .== m) .& (dz.lab .== lab); maximum(abs, dz.model[k] .- dz.measured[k]))
    @test [largest("S1", "B"), largest("S2", "B"), largest("S1", "E"), largest("S2", "E")] ≈ [5.0985, 6.1574, 7.1574, 16.1574] atol = 1.0e-3
end
