# The ternary cements of De Weerdt et al. (2011), against GEMS3K on the same
# budgets and the same phases (`test/reference/xgems_deweerdt2011.json`, written
# by `test/reference/xgems_replay.py deweerdt2011`). The comparison with their
# measurements is on the validation page: a measurement is reported against, not
# asserted.

using JSON

include(joinpath(pkgdir(ChemistryLab), "scripts", "de_weerdt_2011.jl"))

@testsection "validation: the ternary cements of De Weerdt et al. (2011) against GEMS3K" begin
    # The transcription holds together: the column "%OPC reacted" of Table 7 is one
    # minus the four clinker phases over their sum at 0 days, to its rounding.
    for mix in DW11_MIXES
        c0 = sum(something(dw11_phase_content(mix, 0.0, p), 0.0) for p in ("C3S", "C2S", "C3A", "C4AF"))
        for d in (1.0, 7.0, 28.0, 90.0, 180.0)
            left = sum(something(dw11_phase_content(mix, d, p), 0.0) for p in ("C3S", "C2S", "C3A", "C4AF"))
            @test 100 * (1 - left / c0) ≈ dw11_phase_content(mix, d, "opc_reacted") atol = 0.35
        end
    end

    fixture = JSON.parsefile(joinpath(@__DIR__, "reference", "xgems_deweerdt2011.json"))
    cs = dw11_system()
    model = cemdata18_activity_model(:KOH)

    # The replay ran on the phases of this system, and on the budgets the recipes
    # give at each age of Table 8.
    @test sort(String.(symbol.(cs.species))) == sort(fixture["species"])
    @test all(row -> row["converged"], fixture["rows"])
    for row in fixture["rows"]
        elements = budget_elements(cs, budget(dw11_recipe(row["mix"]), cs; t = row["time_d"]).b)
        for (el, x) in row["elements"]
            @test elements[el] ≈ x rtol = 1.0e-12 atol = 1.0e-15
        end
    end

    # The certified equilibria agree with GEMS3K, the AFm binary declared as the
    # export declares it. When the fixture was written the largest differences
    # over the twenty budgets were 0.0008 in pH and 1.6 % on an element
    # (silicon); the tolerances are those of the Lothenbach and Winnefeld test.
    for row in fixture["rows"]
        row["time_d"] in (7.0, 90.0) || continue
        paste = budget(dw11_recipe(row["mix"]), cs; t = row["time_d"])
        eq, cert = equilibrate_certified(paste.state; model, b = paste.b)
        @test cert.optimal
        @test pH(eq, model) ≈ row["pH"] atol = 0.005
        ours = pore_solution_mmol(eq)
        for (el, x) in row["mmol_per_kg_water"]
            @test ours[el] ≈ x rtol = 0.03
        end
    end
end
