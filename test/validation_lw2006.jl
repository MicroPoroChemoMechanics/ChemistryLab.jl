# The CEM I 42.5 N of Lothenbach & Winnefeld (2006), against GEMS3K on the same
# budgets and the same phases (`test/reference/xgems_lw2006.json`, written by
# `test/reference/xgems_lw2006.py`). The comparison with their measured pore
# solution is on the validation page: a measurement is reported against, not
# asserted.

using JSON

include(joinpath(pkgdir(ChemistryLab), "scripts", "lothenbach_winnefeld_2006.jl"))

@testsection "validation: the CEM I 42.5 N of Lothenbach and Winnefeld (2006) against GEMS3K" begin
    fixture = JSON.parsefile(joinpath(@__DIR__, "reference", "xgems_lw2006.json"))
    cs = lw06_system()
    recipe = lw06_recipe()
    model = cemdata18_activity_model(:KOH)

    # The replay ran on the phases of this system, and on the budgets its recipe
    # gives at each age of Table 3.
    @test sort(String.(symbol.(cs.species))) == sort(fixture["species"])
    @test [r["time_h"] for r in fixture["rows"]] == lw06_hours()
    for row in fixture["rows"]
        elements = lw06_elements(cs, budget(recipe, cs; t = row["time_h"] / 24).b)
        for (el, x) in row["elements"]
            @test elements[el] ≈ x rtol = 1.0e-12 atol = 1.0e-15
        end
    end

    # The certified equilibria agree with GEMS3K. The two codes read different
    # releases of Cemdata18 (the export's, and ThermoHub v1.1.1), the export
    # keeps aqueous species of elements this budget lacks, and each computes the
    # Debye-Hueckel A and B from its own water model. When the fixture was written
    # the largest differences over the sixteen ages were 0.002 in pH and 1.6 % on
    # an element (silicon); the tolerances below are those, with a margin.
    for h in (1.0, 16.0, 26.0, 696.0)
        row = only(r for r in fixture["rows"] if r["time_h"] == h)
        @test row["converged"]
        paste = budget(recipe, cs; t = h / 24)
        eq, cert = equilibrate_certified(paste.state; model, b = paste.b)
        @test cert.optimal
        @test pH(eq, model) ≈ row["pH"] atol = 0.005
        ours = lw06_pore_solution(eq)
        for (el, x) in row["mmol_per_kg_water"]
            @test ours[el] ≈ x rtol = 0.03
        end
    end
end
