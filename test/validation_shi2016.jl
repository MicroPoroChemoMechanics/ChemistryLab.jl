# The carbonated mortars of Shi et al. (2016), against GEMS3K on the same budgets
# and the same phases (`test/reference/xgems_shi2016.json`, written by
# `test/reference/xgems_replay.py shi2016`). The comparison with their
# measurements and with their own calculation is on the carbonation page: a
# measurement is reported against, not asserted.

using JSON

include(joinpath(pkgdir(ChemistryLab), "scripts", "shi_2016.jl"))

@testsection "validation: the carbonated mortars of Shi et al. (2016) against GEMS3K" begin
    # The transcription of Table 5 holds together: the total CO2 binding capacity
    # is the CaO of portlandite, C-S-H, monocarbonate, straetlingite and
    # ettringite, and the effective one the total less the decalcified C-S-H, to
    # the rounding of the table.
    for mix in SHI16_MIXES
        parts = sum(shi16_computed(mix, p) for p in ("portlandite", "C-S-H", "monocarbonate", "straetlingite", "ettringite"))
        @test parts ≈ shi16_computed(mix, "total") atol = 0.0025
        @test shi16_computed(mix, "total") - shi16_computed(mix, "decalcified C-S-H") ≈ shi16_computed(mix, "effective") atol = 0.0015
    end
    # Table 2 sums to the binder.
    m = shi16_table("mixes")
    @test all(ustrip.(m.wpc .+ m.mk .+ m.ls) .≈ 100)

    fixture = JSON.parsefile(joinpath(@__DIR__, "reference", "xgems_shi2016.json"))
    cs = shi16_system()
    model = cemdata18_activity_model(:KOH)
    t = SHI16_AGE

    # The replay ran on the phases of this system, and on the budgets the recipes
    # and the carbon dioxide give.
    @test sort(String.(symbol.(cs.species))) == sort(fixture["species"])
    for row in fixture["rows"]
        paste = budget(shi16_recipe(row["mix"]), cs; t)
        elements = budget_elements(cs, shi16_carbonated_budget(cs, paste.b, row["co2_g"]))
        for (el, x) in row["elements"]
            @test elements[el] ≈ x rtol = 1.0e-12 atol = 1.0e-15
        end
    end
    @test all(row -> row["converged"], fixture["rows"])

    # The certified equilibria agree with GEMS3K, before carbonation and on the
    # way down, each paste carbonated as the page carbonates it: step by step,
    # each equilibrium started from the one before. (Started cold from the paste,
    # the cement alone at 20 g stops at a stationarity of 5.6e-6, uncertified,
    # its pH within 0.001 of GEMS3K.) When the fixture was written the largest
    # differences over the 44 budgets were 0.011 in pH and 7 % on an element above
    # 0.01 mmol/kg (aluminum); the tolerances are those, with a margin. An element
    # at a trace is left out: a relative difference means nothing there.
    M_CO2 = ustrip(us"g/mol", Species("CO2")[:M])
    grams = filter(<=(40.0), shi16_co2_grams())
    for mix in SHI16_MIXES
        rs, _ = equilibrate_certified(shi16_recipe(mix), cs; t, model)
        sweep = carbonate(rs, grams ./ M_CO2)
        @test all(s -> s.certificate.optimal, sweep)
        for row in fixture["rows"]
            (row["mix"] == mix && row["co2_g"] in (0.0, 20.0, 40.0)) || continue
            eq = sweep[findfirst(==(row["co2_g"]), grams)].state
            @test pH(eq, model) ≈ row["pH"] atol = 0.02
            ours = pore_solution_mmol(eq)
            for (el, x) in row["mmol_per_kg_water"]
                x > 0.01 && @test ours[el] ≈ x rtol = 0.08
            end
        end
    end
end
