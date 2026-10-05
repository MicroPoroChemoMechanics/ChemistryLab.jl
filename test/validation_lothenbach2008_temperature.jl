# The two cements of Lothenbach, Matschei, Möschner and Glasser (2008) from 0
# to 60 °C, whose hydrates the page `examples/hydrates_temperature.md` compares
# with the article's calculation. What is asserted: the transcription is
# consistent, the constants of Parrott and Killoh the article prints are those
# the package ships, and the numbers the page states are pinned.

isdefined(@__MODULE__, :l08t_material) || include(joinpath(pkgdir(ChemistryLab), "scripts", "lothenbach2008_temperature.jl"))

@testsection "validation: the hydrates of SRPC and PLC from 0 to 60 °C (Lothenbach et al. 2008)" begin
    # Table 2 of the article is the Parrott-Killoh set of Lavergne et al. (2018).
    pk = l08t_table("parrott_killoh")
    for (phase, record) in L08T_CLINKER
        i = findfirst(==(phase), pk.phase)
        ref = getproperty(ChemistryLab, Symbol("PK84_PARAMS_", record))
        @test ustrip(pk.K1[i]) == ustrip(us"1/d", ref.k₁)
        @test ustrip(pk.K2[i]) == ustrip(us"1/d", ref.k₂)
        @test ustrip(pk.K3[i]) == ustrip(us"1/d", ref.k₃)
        @test (ustrip(pk.N1[i]), ustrip(pk.N3[i])) == (ref.n₁, ref.n₃)
    end
    # The volumes read from Figs. 5 and 6: seven layers on each side of the
    # transition, and each figure's top, the sum of its layers, lower above it.
    fig = l08t_table("hydrate_volumes")
    @test length(fig.phase) == 28
    top(c, T) = sum(ustrip(us"cm^3", fig.volume[i]) for i in eachindex(fig.phase) if fig.cement[i] == c && ustrip(fig.temperature_C[i]) == T)
    @test top("SRPC", 5) > top("SRPC", 58)
    @test top("PLC", 5) > top("PLC", 58)

    # The degrees of the clinker after 150 days, and the transition: no
    # monosulfate just below it, monosulfate just above.
    @test l08t_overall_degree("SRPC") ≈ 0.7111 atol = 1.0e-4
    @test l08t_overall_degree("PLC") ≈ 0.827 atol = 1.0e-4
    cs = l08t_system()
    srpc = l08t_material("SRPC")
    @test !l08t_has_monosulfate(l08t_state("SRPC", 52.9; cs, material = srpc))
    @test l08t_has_monosulfate(l08t_state("SRPC", 53.3; cs, material = srpc))
    plc = l08t_material("PLC")
    @test !l08t_has_monosulfate(l08t_state("PLC", 52.4; cs, material = plc))
    @test l08t_has_monosulfate(l08t_state("PLC", 52.8; cs, material = plc))
    # Inside the range the article gives for its own data.
    @test l08t_value("transition_temperature_low") < 52.4 && 53.3 < l08t_value("transition_temperature_high")

    # The SRPC on either side, as the page prints it.
    below = l08t_grouped(l08t_volumes(l08t_state("SRPC", 5.0; cs, material = srpc)))
    @test [below[g] for g in ("ettringite", "monocarbonate", "monosulfate", "portlandite", "C-S-H")] ≈
        [6.59, 1.53, 0.0, 10.57, 23.93] atol = 0.01
    above = l08t_grouped(l08t_volumes(l08t_state("SRPC", 58.0; cs, material = srpc)))
    @test [above[g] for g in ("ettringite", "monocarbonate", "monosulfate", "portlandite", "C-S-H")] ≈
        [4.61, 0.0, 2.67, 10.84, 23.55] atol = 0.01
    @test !haskey(below, "other") && !haskey(above, "other")
end
