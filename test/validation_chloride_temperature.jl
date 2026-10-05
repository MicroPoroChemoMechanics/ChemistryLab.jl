# The chloride AFm phases of Balonis (2019) from 0 to 99 °C, which the page
# `examples/chloride_temperature.md` compares with the article's calculation.
# What is asserted: the records of Friedel's and Kuzel's salts are those of the
# article's Table 1, and the numbers the page states are pinned.

isdefined(@__MODULE__, :b19_scan) || include(joinpath(pkgdir(ChemistryLab), "scripts", "balonis2019_temperature.jl"))

@testsection "validation: Friedel's and Kuzel's salts from 0 to 99 °C (Balonis 2019)" begin
    # Table 1 is the Cemdata18 record of each salt at 298 K, to its last digit.
    t1 = literature_table(B19, "chloride_afm_thermodynamics")
    for (row, sym) in (("C4ACl2H10", "C4AClH10"), ("C4As0.5ClH12", "C4AsClH12"))
        i = findfirst(==(row), t1.phase)
        s = B19_DB[sym]
        at(q, u) = ustrip(u, s[q](T = 298.15u"K", P = 1.0e5u"Pa"; unit = true))
        @test at(:ΔₐG⁰, us"kJ/mol") ≈ ustrip(us"kJ/mol", t1.dfG[i]) atol = 0.11
        @test at(:ΔₐH⁰, us"kJ/mol") ≈ ustrip(us"kJ/mol", t1.dfH[i]) atol = 0.05
        @test at(:S⁰, us"J/(mol*K)") ≈ ustrip(us"J/(mol*K)", t1.S[i]) atol = 0.05
        @test at(:Cp⁰, us"J/(mol*K)") ≈ ustrip(us"J/(mol*K)", t1.Cp[i]) atol = 0.5
    end

    cs = b19_system()
    # Kuzel's salt against half monosulfate, half Friedel's salt and one water.
    g(s, T) = ustrip(us"kJ/mol", B19_DB[s][:ΔₐG⁰](T = (T + 273.15) * u"K", P = 1.0e5u"Pa"; unit = true))
    ΔrG(T) = g("monosulphate12", T) / 2 + g("C4AClH10", T) / 2 + g("H2O@", T) - g("C4AsClH12", T)
    @test [ΔrG(25), ΔrG(99)] ≈ [1.59, 0.34] atol = 0.005

    # The mixture with the most chloride and calcite, as the page scans it.
    scan = b19_scan(cs, "Fig. 12")
    @test b19_transition(scan, "Friedel's salt", "last") == (90, 95)
    @test b19_transition(scan, "ettringite", "last") == (90, 95)
    @test b19_transition(scan, "monosulfate", "first") === nothing
    @test 1 - b19_volumes(scan, 99)["total"] / b19_volumes(scan, 0)["total"] ≈ 0.26 atol = 0.001
    # The others, on either side of the changes the page reports.
    kuzel = b19_solids(cs, "Fig. 5", 25.0)
    @test kuzel["Kuzel's salt"] ≈ 2.309 atol = 0.001
    @test kuzel["Friedel's salt"] < 1.0e-3
    @test b19_solids(cs, "Fig. 5", 99.0)["Kuzel's salt"] ≈ 0.305 atol = 0.001
    @test b19_solids(cs, "Fig. 11", 50.0)["monocarbonate"] > 1.0e-3
    @test b19_solids(cs, "Fig. 11", 55.0)["monocarbonate"] < 1.0e-3
    @test !haskey(kuzel, "other")
end
