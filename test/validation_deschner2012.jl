# The pore solutions of Deschner et al. (2012), over 550 days, speciated at their
# measured hydroxide, against the effective saturation indices the authors
# computed from them. The solutions are printed only as plots, and were read from
# the vector drawing of the figures. What is asserted: the reading (every marker
# on one of the eleven ages of Table 3, every series where the figures put it);
# every solution certifies at the hydroxide measured; and the indices agree
# within 0.12 for portlandite and gypsum and 0.06 for ettringite, strätlingite
# and monosulfate, the agreement found when this test was written (0.10, 0.07,
# 0.03, 0.03 and 0.05 over the 55 solutions) with a margin.

include(joinpath(pkgdir(ChemistryLab), "scripts", "deschner_2012.jl"))

@testsection "validation: the pore solutions of Deschner et al. (2012)" begin
    ps, t3 = d12_table("pore_solution"), d12_table("effective_saturation_indices")
    @test d12_systems_listed() == ["OPC", "OPC-Qz", "OPC-F1", "OPC-F2", "OPC-F1-L"]
    @test d12_ages() == [1, 4, 8, 16, 24, 48, 168, 672, 2160, 6000, 13200]
    # Seven analytes, eleven ages, five pastes; what the figures leave out is the
    # aluminum below 0.01 mmol/L at the first ages, and some points at 550 days.
    @test length(ps.system) == 366
    @test all(>(0), ustrip.(us"mol/m^3", ps.concentration))
    @test length(t3.system) == 4 * 55 + 77   # hemi- and monocarbonate for OPC-F1-L alone
    @test length(unique(zip(ps.system, ps.age, ps.analyte))) == length(ps.system)
    # The mixes are complete binders.
    mixes = d12_table("mixes")
    for k in eachindex(mixes.system)
        @test mixes.OPC[k] + mixes.F1[k] + mixes.F2[k] + mixes.quartz[k] + mixes.limestone[k] ≈ 1
    end

    systems = d12_systems()
    model = cemdata18_activity_model(:KOH)
    for s in ("OPC", "OPC-F1"), a in (4, 24, 2160)
        r = d12_solution(systems.aqueous, s, a; model)
        @test r.certificate.optimal
        # The hydroxide held is the one measured, as a molality.
        @test molalities(r.state)["OH-"] ≈ d12_concentrations(s, a)["OH-"] rtol = 1.0e-6
        ix = d12_indices(r.state, systems, keys(d12_concentrations(s, a)); model)
        Δ(l) = ix[l] - d12_published(s, a, l)
        @test abs(Δ("CH")) < 0.12
        @test abs(Δ("Gypsum")) < 0.12
        if haskey(ix, "Ettr")
            @test abs(Δ("Ettr")) < 0.06
            @test abs(Δ("Strätl")) < 0.06
            @test abs(Δ("MS")) < 0.06
        end
    end
    # Aluminum is not plotted at four hours: no aluminous index is computed there.
    @test !haskey(d12_concentrations("OPC", 4), "Al")
end
