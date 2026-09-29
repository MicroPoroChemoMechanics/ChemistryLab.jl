# The Portland cement with 4 % limestone of Lothenbach et al. (2008), which the
# CASH+ page computes with CSHQ and with CASH+NK against its pore solution. What
# is asserted: the two transcriptions of the cement and of its kinetics agree, the
# recipe is the paper's, and the paste certifies at the first and last ages with
# both models of the C-S-H. The comparison with the analysis is on the page: a
# measurement is reported against, not asserted.

include(joinpath(pkgdir(ChemistryLab), "scripts", "lothenbach_2008.jl"))

@testsection "validation: the PC4 of Lothenbach et al. (2008), with CSHQ and CASH+NK" begin
    # Miron et al. (2022b) reprint the PC4 column and the kinetic constants; the
    # two transcriptions, made separately, are the same numbers.
    ours = l08_table("normative_phases")
    theirs = literature_table("Miron2022b", "pc4_normative_phases")
    for (p, x) in zip(ours.phase, ours.PC4)
        @test ustrip(theirs.PC4[findfirst(==(p), theirs.phase)]) == ustrip(x)
    end
    minors = l08_table("clinker_minor_oxides")
    for (ox, x) in zip(minors.oxide, minors.PC4)
        @test ustrip(theirs.PC4[findfirst(==("$ox in the clinker"), theirs.phase)]) == ustrip(x)
    end
    pk, pk22 = l08_table("parrott_killoh"), literature_table("Miron2022b", "parrott_killoh")
    for c in (:K1, :N1, :K2, :K3, :N3, :H)
        @test ustrip.(getproperty(pk, c)) == ustrip.(getproperty(pk22, c))
    end
    @test l08_value("blaine_PC4") == ustrip(literature_value("Miron2022b", "blaine_PC4"))
    # The normative phases sum to 100 g within the rounding of the table.
    @test sum(ustrip.(ours.PC)) ≈ 100 atol = 0.2
    @test sum(ustrip.(ours.PC4)) ≈ 100 atol = 0.2

    # The recipe: the minor oxides shared among the clinker phases are the
    # paper's totals, and the budget holds them.
    shares = l08_minor_distribution("PC4")
    for (ox, x) in zip(minors.oxide, minors.PC4)
        @test sum(s[ox] for s in values(shares)) ≈ ustrip(x) rtol = 1.0e-12
    end
    recipe = l08_recipe("PC4")
    @test sum(c.mass_fraction for c in recipe.binder[1].first.constituents) ≈ sum(ustrip.(ours.PC4)) / 100 rtol = 1.0e-12

    # Both models of the C-S-H certify the paste at the first and last ages.
    model = cemdata18_activity_model(:KOH)
    for gel in (:CSHQ, :CASHNK)
        cs = l08_system(gel)
        run = hydrate(recipe, cs, [1.0, 400.0]; model)
        @test all(rs -> rs.certificate.optimal, run.states)
        ps = pore_solution_mmol(run.states[end].state)
        @test ps["K"] > 0 && ps["Na"] > 0
    end
end
