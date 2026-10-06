# Alkali uptake by C-A-S-H, L'Hôpital et al. (2016), which the page
# `tutorials/validation_alkali_uptake.md` compares with CSHQ and CASH+NK. What
# is asserted: the transcription of Appendices A to C is consistent with
# itself, and the numbers the page states are pinned.

isdefined(@__MODULE__, :lh16_compute) || include(joinpath(pkgdir(ChemistryLab), "scripts", "lhopital2016_alkali.jl"))

@testsection "validation: alkali uptake by C-A-S-H (L'Hôpital et al. 2016)" begin
    # Appendix A: 2 g of solids in every batch, and the targets they make, from
    # the molar masses of the library.
    mix = lh16_table("mixing_proportions")
    M(s) = ustrip(us"g/mol", LH16_DB[s][:M])
    g(q) = ustrip(us"g", q)
    for i in eachindex(mix.CaO)
        @test g(mix.CaO[i]) + g(mix.SiO2[i]) + g(mix.CaO_Al2O3[i]) ≈ 2.0 atol = 0.0015
        si = g(mix.SiO2[i]) / M("Amor-Sl")
        ca = g(mix.CaO[i]) / M("Lim") + g(mix.CaO_Al2O3[i]) / M("CA")
        @test ca / si ≈ mix.Ca_Si_target[i] atol = 0.002
        @test 2 * g(mix.CaO_Al2O3[i]) / M("CA") / si ≈ mix.Al_Si_target[i] atol = 0.001
    end
    # Appendices B and C: every row read, and every sample of C a batch of B.
    @test length(lh16_table("pore_solution").element) == 682
    @test length(lh16_table("pH").pH) == 127
    sol = lh16_table("solid_composition")
    @test length(sol.Ca_Si) == 125
    batch(t, i) = (t.Al_Si_target[i], t.Ca_Si_target[i], t.alkali[i], ustrip(us"mol/L", t.alkali_concentration[i]))
    ph = lh16_table("pH")
    @test Set(batch(sol, i) for i in eachindex(sol.Ca_Si)) ⊆ Set(batch(ph, i) for i in eachindex(ph.pH))

    results = Dict{String, Any}()
    for gel in ("CSHQ", "CASH+NK")
        results[gel] = lh16_compute(gel)
    end
    for gel in ("CSHQ", "CASH+NK")
        @test length(results[gel]) == 49
        @test all(v.certified for v in values(results[gel]))
    end

    # Section 3: the alkali of the gel, computed over measured.
    ext(x) = (minimum(x), maximum(x))
    for (ca, q, n) in (
            (0.6, (0.12, 1.25), (0.82, 1.55)), (0.8, (0.24, 0.46), (1.0, 1.36)),
            (1.0, (0.32, 0.98), (0.9, 1.5)), (1.2, (0.41, 1.19), (0.62, 1.16)),
        )
        @test all(isapprox.(ext(lh16_uptake_ratios(results["CSHQ"], ca)), q; atol = 0.005))
        @test all(isapprox.(ext(lh16_uptake_ratios(results["CASH+NK"], ca)), n; atol = 0.005))
    end
    # Section 4: the solution, as medians of computed over measured.
    med(x) = sort(x)[cld(length(x), 2)]
    q, n = lh16_solution_ratios(results["CSHQ"]; low = true), lh16_solution_ratios(results["CASH+NK"]; low = true)
    # Each to the digit the page prints it to.
    near(x, y) = all(isapprox.(x, y; atol = 0.005))
    @test near([med(q.alkali), med(q.Si), med(q.Ca)], [1.34, 2.52, 2.41])
    @test near([med(n.alkali), med(n.Si), med(n.Ca)], [1.04, 0.96, 1.67])
    @test maximum(q.Si) ≈ 26.37 atol = 0.005
    Q, N = lh16_solution_ratios(results["CSHQ"]; low = false), lh16_solution_ratios(results["CASH+NK"]; low = false)
    @test all(isapprox.(ext(Q.alkali), (0.88, 1.14); atol = 0.005))
    @test near([med(Q.Si), med(N.alkali), med(N.Si), med(Q.Ca), med(N.Ca)], [3.61, 1.0, 1.39, 0.65, 0.79])
    @test maximum(Q.Si) ≈ 29.78 atol = 0.005
    # Section 5: the pH, where CSHQ leaves too much alkali, and the batch
    # without alkali neither model reaches.
    at(gel, b) = results[gel][b]
    @test at("CSHQ", (0.05, 0.6, "KOH", 0.05)).pH - _lh16_mean(lh16_measured_solution((0.05, 0.6, "KOH", 0.05), "pH")) ≈ 0.59 atol = 0.005
    pH08 = _lh16_mean(lh16_measured_solution((0.05, 0.8, "none", 0.0), "pH"))
    @test pH08 ≈ 10.27 atol = 0.005
    @test near([at("CSHQ", (0.05, 0.8, "none", 0.0)).pH, at("CASH+NK", (0.05, 0.8, "none", 0.0)).pH] .- pH08, [0.97, 0.55])
    @test maximum(q.ΔpH) ≈ 0.97 atol = 0.005
end
