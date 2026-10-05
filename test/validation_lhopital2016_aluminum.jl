# Aluminum uptake by C-S-H, L'Hôpital et al. (2016a), which the page
# `tutorials/validation_aluminum_uptake.md` compares with CNASH_ss and CSHQ.
# What is asserted: the transcription of Appendices A, B and D is consistent
# with itself, every synthesis certifies with both models, and the numbers the
# page states are pinned.

isdefined(@__MODULE__, :lh16a_compute) || include(joinpath(pkgdir(ChemistryLab), "scripts", "lhopital2016_aluminum.jl"))

@testsection "validation: aluminum uptake by C-S-H (L'Hôpital et al. 2016a)" begin
    # Appendix A: 2 g of solids in every synthesis, and the targets they make,
    # from the molar masses of the library.
    mix = lh16a_table("mixing_proportions")
    @test length(mix.CaO) == 34
    M(s) = ustrip(us"g/mol", LH16_DB[s][:M])
    g(q) = ustrip(us"g", q)
    for i in eachindex(mix.CaO)
        @test g(mix.CaO[i]) + g(mix.SiO2[i]) + g(mix.CaO_Al2O3[i]) ≈ 2.0 atol = 0.0015
        si = g(mix.SiO2[i]) / M("Amor-Sl")
        @test (g(mix.CaO[i]) / M("Lim") + g(mix.CaO_Al2O3[i]) / M("CA")) / si ≈ mix.Ca_Si_target[i] atol = 0.002
        @test 2 * g(mix.CaO_Al2O3[i]) / M("CA") / si ≈ mix.Al_Si_target[i] atol = 0.001
    end
    # Appendices B and D: every row read, each a synthesis of Appendix A, and
    # every synthesis has a measured gel.
    gel = lh16a_table("gel_composition")
    @test length(gel.Ca_Si) == 42
    @test length(lh16a_table("pore_solution").element) == 240
    @test length(lh16a_table("pH").pH) == 63
    synth = Set(zip(mix.Ca_Si_target, mix.Al_Si_target))
    @test Set(zip(gel.Ca_Si_target, gel.Al_Si_target)) == synth
    @test all(!ismissing(lh16a_measured(k..., :Ca_Si)) for k in synth)

    results = Dict{String, Any}()
    for model in ("CNASH_ss", "CSHQ")
        results[model] = lh16a_compute(model)
        @test length(results[model]) == 34
        @test all(v.certified for v in values(results[model]))
    end
    cn, cq = results["CNASH_ss"], results["CSHQ"]

    # Section 2: CNASH_ss takes all the aluminum up to 0.05 from Ca/Si 1.0 up,
    # and holds 0.10 to 0.12 beyond, whatever the Ca/Si.
    for ca in (1.0, 1.2, 1.4, 1.6), al in (0.03, 0.05)
        @test cn[(ca, al)].Al_Si ≈ al atol = 5.0e-4
    end
    @test all(isapprox.([cn[(0.8, 0.03)].Al_Si, cn[(0.8, 0.05)].Al_Si, cn[(0.6, 0.03)].Al_Si, cn[(0.6, 0.05)].Al_Si], [0.032, 0.053, 0.042, 0.071]; atol = 5.0e-4))
    beyond = [v.Al_Si for (k, v) in cn if k[2] >= 0.1]
    @test 0.0995 <= minimum(beyond) && maximum(beyond) <= 0.1225
    @test all(cq[k].Al_Si == 0 for k in synth)
    @test all(cn[(ca, 0.2)].straetlingite < 1.0e-9 for ca in (0.6, 0.8, 1.0, 1.2, 1.4, 1.6))

    # Section 3: where the rest goes, at Al/Si 0.2.
    near(x, y; atol = 0.05) = all(isapprox.(x, y; atol))
    @test near([cn[(ca, 0.2)].katoite for ca in (1.2, 1.4, 1.6)], [9.1, 11.1, 10.5])
    s06, s10, s16 = (lh16a_aluminum_split(cn[(ca, 0.2)].state, "CNASH_ss") for ca in (0.6, 1.0, 1.6))
    @test near([s06.gibbsite, s10.gibbsite, s16.katoite], [61.1, 38.5, 50.5])
    q06, q10, q16 = (lh16a_aluminum_split(cq[(ca, 0.2)].state, "CSHQ") for ca in (0.6, 1.0, 1.6))
    @test near([q06.gibbsite, q10.straetlingite, q16.katoite], [100.0, 99.6, 99.9])
    @test all(cq[(ca, 0.03)].straetlingite > 0 for ca in (0.8, 1.0))
    @test all(cq[(ca, 0.03)].katoite > 0 for ca in (1.2, 1.4, 1.6))

    # Section 4: the calcium of the gel at Al/Si 0.05.
    @test near([cn[(ca, 0.05)].Ca_Si for ca in (1.4, 1.6)], [1.19, 1.19]; atol = 0.005)
    @test near([cn[(ca, 0.05)].portlandite for ca in (1.4, 1.6)], [3.4, 10.5])
    gap = [abs(cq[(ca, 0.05)].Ca_Si - lh16a_measured(ca, 0.05, :Ca_Si)) for ca in (0.8, 1.0, 1.2, 1.4, 1.6)]
    @test argmax(gap) == 3
    @test maximum(gap) ≈ 0.052 atol = 5.0e-4
    @test all(v.portlandite < 1.0e-9 for v in values(cq))
end
