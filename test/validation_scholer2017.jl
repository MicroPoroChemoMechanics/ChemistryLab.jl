# The early pore solutions of Schöler et al. (2017), speciated at their measured
# pH, against the saturation indices the authors computed from them. What is
# asserted: every solution certifies; the indices of portlandite and gypsum agree
# within 0.1 and that of the C-S-H within 0.25, the agreement found when this
# test was written (0.07, 0.08 and 0.21) with a margin; and the differences obey
# the stoichiometry of ettringite and monosulfate (ΔE − ΔMs = 2 ΔGp), which holds
# when the two calculations differ by a factor on each ion's activity.

include(joinpath(pkgdir(ChemistryLab), "scripts", "scholer_2017.jl"))

@testsection "validation: the early pore solutions of Schöler et al. (2017)" begin
    t6, t7 = s17_table("pore_solution"), s17_table("saturation_indices")
    @test length(s17_systems_listed()) == 8 && length(s17_ages()) == 6
    @test length(t6.system) == 48 * 6 - 1         # one silicon not detected
    @test length(t7.system) == 48 * 8 - 3         # three indices not calculated
    # The last column of Table 7 is the natural logarithm of the ion activity
    # product of alite, and the column before it that product against log K = -22.
    for s in s17_systems_listed(), a in s17_ages()
        ln = s17_published(s, a, "ln_IAP_C3S_exp")
        ln === nothing || @test ln / log(10) + 22.0 ≈ s17_published(s, a, "C3S_exp") atol = 0.03
    end

    systems = s17_systems()
    model = cemdata18_activity_model(:KOH)
    for s in ("C", "C-FA-\$", "C-Q"), a in (0.083, 6.0)
        r = s17_solution(systems.aqueous, s, a; model)
        @test r.certificate.optimal
        @test pH(r.state, model) ≈ s17_pH(s, a) atol = 0.01
        ix = s17_indices(r.state, systems; model)
        Δ(l) = ix[l] - s17_published(s, a, l)
        @test abs(Δ("CH")) < 0.1
        @test abs(Δ("Gp")) < 0.1
        s17_published(s, a, "CSH") === nothing || @test abs(Δ("CSH")) < 0.25
        @test Δ("E") - Δ("Ms") ≈ 2Δ("Gp") atol = 0.03
    end
end
