# The pore solutions of Deschner et al. (2013), cured at 7 to 80 °C, speciated at
# their measured hydroxide and at their curing temperature, against the
# effective saturation indices the authors computed from them. What is asserted:
# the transcription (every solution of Tables A.1 and A.2, every index of Table
# B.1); every solution certifies; the hydroxide held is the one measured; and the
# agreements `docs/src/tutorials/validation_temperature_pore_solutions.md`
# states, pinned.

isdefined(@__MODULE__, :d13_solution) || include(joinpath(pkgdir(ChemistryLab), "scripts", "deschner_2013.jl"))

@testsection "validation: the pore solutions of Deschner et al. (2013), 7 to 80 °C" begin
    ps, b1 = d13_table("pore_solution"), d13_table("effective_saturation_indices")
    # Nine analytes, two pastes, five temperatures; the chloride of OPC-Qz at
    # 23 °C and at 50 °C after 7 days is not available.
    @test length(ps.system) == 488
    @test length(unique(zip(ps.system, ps.temperature, ps.age, ps.analyte))) == length(ps.system)
    @test all(>(0), ustrip.(us"mol/m^3", ps.concentration))
    @test sort(unique(ustrip.(us"K", ps.temperature))) ≈ [280.15, 296.15, 313.15, 323.15, 353.15]
    # Four indices for each of the 27 solutions of Table B.1, in both pastes.
    @test length(b1.system) == 4 * 27 * 2
    @test length(d13_solutions()) == 54
    mixes = d13_table("mixes")
    @test all(mixes.OPC .+ mixes.fly_ash .+ mixes.quartz .≈ 1)

    # Every solution of Table B.1, speciated at its hydroxide and its
    # temperature, against the indices of the paper. PINNED, not bounded: the
    # page states these figures, and a bound would let them go stale.
    systems = d13_systems()
    model = cemdata18_activity_model(:KOH)
    labels = ("CH", "Ettr", "Ms", "Strätl")
    diffs = Dict{Tuple{String, Float64, Float64}, Dict{String, Float64}}()
    for k in d13_solutions()
        r = d13_solution(systems.aqueous, k...; model)
        @test r.certificate.optimal
        @test temperature(r.state) ≈ k[2] * u"K"
        @test molalities(r.state)["OH-"] ≈ d13_concentrations(k...)["OH-"] rtol = 1.0e-6
        ix = d13_indices(r.state, systems, k[2]; model)
        diffs[k] = Dict(l => ix[l] - d13_published(k..., l) for l in labels)
    end
    worst(sel, ls = labels) = maximum(abs(diffs[k][l]) for k in keys(diffs) if sel(k) for l in ls)
    # From 7 to 50 °C: 0.08 at most, and a mean within 0.04 at each temperature.
    @test worst(k -> k[2] < 330) ≈ 0.08 atol = 2.5e-3
    for T in (280.15, 296.15, 313.15, 323.15), l in labels
        d = [diffs[k][l] for k in keys(diffs) if k[2] ≈ T]
        @test abs(sum(d) / length(d)) < 0.04
    end
    # At 80 °C: OPC-Qz to 0.12; OPC-FA up to 0.41 on portlandite, 0.23 on the others.
    @test worst(k -> k[2] > 350 && k[1] == "OPC-Qz") ≈ 0.121 atol = 2.5e-3
    @test worst(k -> k[2] > 350 && k[1] == "OPC-FA", ("CH",)) ≈ 0.405 atol = 2.5e-3
    @test worst(k -> k[2] > 350 && k[1] == "OPC-FA", ("Ettr", "Ms", "Strätl")) ≈ 0.227 atol = 2.5e-3
end
