# The heat of the slag-limestone cement of Snellings et al. (2022) at its
# measured degrees of reaction (scripts/snellings2022_heat.jl), the glass of the
# slag at the enthalpy glass_enthalpy builds, and at the progress its bound
# water gives; the page examples/slag_heat.md compares every age and
# temperature.

include(joinpath(pkgdir(ChemistryLab), "scripts", "snellings2022_heat.jl"))

@testsection "validation: the heat of a slag-limestone cement (Snellings et al. 2022)" begin
    cs = sn22h_system()
    model = cemdata18_activity_model(:KOH)
    g = sn22h_glass(20)
    h = g.enthalpy
    αc, αs = sn22h_phase_degrees(20, 0.5, 28), sn22h_degree("slag", 20, 0.5, 28)
    rs0 = sn22h_state(cs, 20, 0.5, 0.0, 0.0; h_glass = h, model)
    rs = sn22h_state(cs, 20, 0.5, αc, αs; h_glass = h, model)
    @test rs0.certificate.optimal && rs.certificate.optimal

    # The glass sets aside the oxides the system has no primary for as it
    # reacts; without their enthalpy the heat is not a number.
    @test isnan(heat_release(rs0, rs))
    aside = sn22h_set_aside()
    pc_g = sn22_value("pc_percent")
    q = heat_release(rs0, rs; set_aside = aside) / pc_g
    # 696 J per gram of Portland cement at 28 days, each clinker phase at its own
    # degree (Fig. 5), against 650 read off Fig. 3.
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

    # The degrees of the four clinker phases (Fig. 5), weighted by their contents
    # at the mixing, give the degree of the clinker the paper plots (Fig. 6a).
    init = literature_table(SN22, "clinker_phase_initial")
    c0 = Dict(String(p) => ustrip(c) for (p, c) in zip(init.phase, init.percent))
    for T in (5, 20, 40), wb in (0.4, 0.5, 0.6), a in (1, 28)
        α = sn22h_phase_degrees(T, wb, a)
        weighted = sum(α[n] * c0[p] for (n, p) in SN22H_PHASES) / sum(values(c0))
        @test weighted ≈ sn22h_degree("clinker", T, wb, a) atol = 0.02
    end

    # One factor on the degrees, fitted on the bound water of the
    # thermogravimetry, at 20 C and one day: the bound water is the measured one,
    # the diffraction is ahead of it, and the heat at that progress is the
    # calorimetry's, which nothing was fitted on.
    o = sn22h_origin(cs, 20, 0.5; model)
    f = sn22h_progress(cs, 20, 1; model, origin = o)
    @test f.tga.bound_water ≈ sn22h_measured_tga("bound water", 20, 1) rtol = 1.0e-3
    @test 0.5 < f.λ < 0.8
    @test isapprox(o.heat(f.state), sn22h_measured(20, 0.5, 1).heat; rtol = 0.1)

    # A constituent with a formula keeps the enthalpy of its record.
    pc = material_template("PC (Snellings 2022)", SN22_DB)
    @test_throws ArgumentError with_enthalpy(pc, Dict("Alite" => -12000.0); source = "x")
    @test_throws ArgumentError with_enthalpy(pc, Dict("no such" => -12000.0); source = "x")
end
