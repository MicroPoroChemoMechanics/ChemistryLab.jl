# The silica fume shotcrete paste ESDRED of Lothenbach et al. (2014), which the
# page validation_silica_fume_paste.md computes with CSHQ and CASH+NK, with and
# without the formate of its set accelerator. What is asserted: the recipe is the
# paper's, the formate is the organic carbon of the accelerator and stays formate,
# and the paste certifies. The comparison with the analysis is on the page.

include(joinpath(pkgdir(ChemistryLab), "scripts", "esdred_2014.jl"))

@testsection "validation: the ESDRED paste of Lothenbach et al. (2014)" begin
    # The binder: 60 % CEM I and 40 % silica fume, the normative phases of the
    # CEM I scaled from the 100.31 g Table 1 prints to 100 g.
    recipe = es14_recipe()
    @test [last(p) for p in recipe.binder] == [0.6, 0.4]
    t = es14_table("normative_phases")
    @test sum(Float64.(ustrip.(t.percent))) ≈ 100.31 atol = 1.0e-9
    @test sum(es14_normative(p) for p in t.phase) ≈ 100 rtol = 1.0e-12
    # The minor oxides shared among the clinker phases are Table 1's totals.
    shares = es14_minor_distribution()
    minors = es14_table("alkalis_in_clinker")
    for (ox, x) in zip(minors.oxide, minors.percent)
        @test sum(s[ox] for s in values(shares)) ≈ ustrip(x) rtol = 1.0e-12
    end
    # The degrees of reaction of the NMR start at zero and grow.
    nmr = es14_nmr_extents()
    @test extent(nmr.silicates, 0) == 0 && extent(nmr.silica_fume, 0) == 0
    @test extent(nmr.silica_fume, 1310) ≈ 0.75 rtol = 1.0e-12

    # The formate is the organic carbon of the accelerator: 202 mM in the mixing
    # water, as Table 3 prints it, within the rounding of the analysis.
    acc, grams = es14_accelerator(; formate = true)
    oxides = only(acc.constituents).oxides
    n_formate = grams * oxides["CO"] / ustrip(us"g/mol", Species("CO")[:M])
    water_L = (es14_value("water_binder_ratio") * 100 + grams * oxides["H2O"]) / 1000
    @test 1000 * n_formate / water_L ≈ es14_value("initial_formate") rtol = 0.06

    # With the reduced species of sulfur and iron left out, nothing oxidizes it:
    # all the formate of the budget is formate at equilibrium, and the paste
    # certifies, with either model of the C-S-H.
    model = cemdata18_activity_model(:KOH)
    for gel in (:CSHQ, :CASHNK)
        cs = es14_system(gel; formate = true)
        @test !any(s -> symbol(s) in ES14_REDUCED, cs.species)
        rs = hydrate(es14_recipe(; formate = true), cs, [28.0]; model)
        @test rs.states[1].certificate.optimal
        n = Dict(symbol(s) => ustrip(us"mol", x) for (s, x) in zip(cs.species, rs.states[1].state.n))
        @test n["For-"] + n["ForH@"] ≈ n_formate rtol = 1.0e-8
    end
    # Without the formate, as Miron et al. (2022b) computed it, at the last age.
    rs = hydrate(es14_recipe(), es14_system(:CASHNK), [1310.0]; model)
    @test rs.states[1].certificate.optimal
end
