# The heat of the slag-limestone cement of Snellings et al. (2022) at its
# measured degrees of reaction (scripts/snellings2022_heat.jl), the glass of the
# slag at the enthalpy glass_enthalpy builds; the page examples/slag_heat.md
# compares every age and temperature.

include(joinpath(pkgdir(ChemistryLab), "scripts", "snellings2022_heat.jl"))

@testsection "validation: the heat of a slag-limestone cement (Snellings et al. 2022)" begin
    cs = sn22h_system()
    model = cemdata18_activity_model(:KOH)
    g = sn22h_glass(20)
    h = g.enthalpy
    αc, αs = sn22h_degree("clinker", 20, 0.5, 28), sn22h_degree("slag", 20, 0.5, 28)
    rs0 = sn22h_state(cs, 20, 0.5, 0.0, 0.0; h_glass = h, model)
    rs = sn22h_state(cs, 20, 0.5, αc, αs; h_glass = h, model)
    @test rs0.certificate.optimal && rs.certificate.optimal

    # The glass sets aside the oxides the system has no primary for as it
    # reacts; without their enthalpy the heat is not a number.
    @test isnan(heat_release(rs0, rs))
    aside = sn22h_set_aside()
    pc_g = sn22_value("pc_percent")
    q = heat_release(rs0, rs; set_aside = aside) / pc_g
    # 669 J per gram of Portland cement at 28 days, against 650 read off Fig. 3.
    @test isapprox(q, sn22h_measured(20, 0.5, 28).heat; rtol = 0.1)

    # The glass at the enthalpy of the crystals of its composition instead: the
    # heat differs by exactly the enthalpy of vitrification of the glass that
    # reacted.
    hc = h - g.vitrification
    q0 = heat_release(
        sn22h_state(cs, 20, 0.5, 0.0, 0.0; h_glass = hc, model),
        sn22h_state(cs, 20, 0.5, αc, αs; h_glass = hc, model); set_aside = aside,
    ) / pc_g
    glass = only(c for c in material_template("slag (Snellings 2022)", SN22_DB).constituents if c.name == "glass")
    m_glass = sn22_value("slag_percent") * glass.mass_fraction
    @test q - q0 ≈ αs * m_glass * ustrip(us"J/g", g.vitrification) / pc_g rtol = 1.0e-6

    # A constituent with a formula keeps the enthalpy of its record.
    pc = material_template("PC (Snellings 2022)", SN22_DB)
    @test_throws ArgumentError with_enthalpy(pc, Dict("Alite" => -12000.0); source = "x")
    @test_throws ArgumentError with_enthalpy(pc, Dict("no such" => -12000.0); source = "x")
end
