# The enthalpy and heat capacity of a silicate glass from its oxides
# (glass_thermochemistry.jl), against the measurements it is built from.

@testsection "the enthalpy and heat capacity of a glass" begin
    M(f) = ustrip(us"g/mol", Species(f)[:M])
    Jg(x) = ustrip(us"J/g", x)
    aq = ChemistryLab._aq17_crystals()
    H(s) = ustrip(us"J/mol", aq[s][:ΔₐH⁰](T = 298.15u"K", P = 1.0e5u"Pa"; unit = true))
    # The mass fractions of a composition given in moles of oxides.
    function mass_fractions(mol)
        m = Dict(ox => v * M(ox) for (ox, v) in mol)
        tot = sum(values(m))
        return Dict(ox => v / tot for (ox, v) in m)
    end

    # A glass of the composition of a measured one is that glass: the enthalpy of
    # its crystal in aq17 plus its enthalpy of vitrification, read here from the
    # transcriptions independently of the code.
    dHv(key, table; kw...) = ustrip(us"J/mol", only(literature_table(key, table; kw...).dH_v))
    cases = (
        ("Gehlenite", Dict("CaO" => 2, "Al2O3" => 1, "SiO2" => 1), dHv("RichetBottinga1986", "vitrification_enthalpies"; formula = "Ca2Al2SiO7", T_s = 298.0u"K")),
        ("Akermanite", Dict("CaO" => 2, "MgO" => 1, "SiO2" => 2), dHv("RichetBottinga1986", "vitrification_enthalpies"; formula = "Ca2MgSi2O7", T_s = 298.0u"K")),
        ("Pseudowoll", Dict("CaO" => 1, "SiO2" => 1), dHv("RichetBottinga1986", "vitrification_enthalpies"; formula = "CaSiO3", T_s = 298.0u"K")),
        ("Anorthite", Dict("CaO" => 1, "Al2O3" => 1, "SiO2" => 2), dHv("RichetBottinga1984", "vitrification_enthalpies"; crystal = "anorthite", T_s = 298.0u"K")),
        ("Diopside", Dict("CaO" => 1, "MgO" => 1, "SiO2" => 2), dHv("RichetBottinga1984", "vitrification_enthalpies"; crystal = "diopside", T_s = 298.0u"K")),
        ("Quartz", Dict("SiO2" => 1), dHv("Navrotsky1980", "vitrification_enthalpies"; substance = "quartz", T = 298.0u"K", reference = "Kracek (1953)")),
    )
    for (crystal, mol, h) in cases
        Mf = sum(v * M(ox) for (ox, v) in mol)
        g = glass_enthalpy(mass_fractions(mol))
        @test Jg(g.formation_298) ≈ (H(crystal) + h) / Mf rtol = 1.0e-12
        @test isempty(g.unassigned) && Jg(g.span) == 0
        @test length(g.norm) == 1
        @test only(values(g.norm)) ≈ 1 / Mf rtol = 1.0e-12
    end

    # Anorthite plus akermanite, Ca3MgAl2Si4O15, is also gehlenite plus diopside
    # plus silica, and half gehlenite, half anorthite, diopside and half
    # pseudowollastonite: every combination with nothing left over is
    # Ak = 1 − Di, An = 1 − Ge, Pwo = Di − Ge, SiO2 = 2Ge − Di, for
    # Ge ≤ Di ≤ min(2Ge, 1), whose vertices are these three. Their enthalpies
    # of vitrification differ by what ideal mixing of the measured glasses
    # cannot reconcile; the result is the midpoint of the extremes and the span
    # their half-difference.
    v(crystal) = only(c[3] for c in cases if c[1] == crystal)
    mol = Dict("CaO" => 3, "MgO" => 1, "Al2O3" => 1, "SiO2" => 4)
    Mf = sum(n * M(ox) for (ox, n) in mol)
    g = glass_enthalpy(mass_fractions(mol))
    vertices = (
        v("Anorthite") + v("Akermanite"),
        v("Gehlenite") + v("Diopside") + v("Quartz"),
        (v("Gehlenite") + v("Anorthite") + v("Pseudowoll")) / 2 + v("Diopside"),
    ) ./ Mf
    lo, hi = extrema(vertices)
    @test Jg(g.vitrification) ≈ (lo + hi) / 2 rtol = 1.0e-10
    @test Jg(g.span) ≈ (hi - lo) / 2 rtol = 1.0e-10
    @test isempty(g.unassigned)

    # A slag holds more lime than the measured glasses can take: the norm is
    # unique, the excess counted as crystalline lime, and the enthalpy of
    # vitrification and its uncertainty are those of the norm.
    rows = literature_table("Snellings2022", "chemical_composition"; material = "GGBFS")
    slag = Dict(String(o) => ustrip(p) / 100 for (o, p) in zip(rows.oxide, rows.percent) if o != "Total")
    @test_throws ArgumentError glass_enthalpy(slag)            # P2O5 has no crystal in the reference
    s = glass_enthalpy(slag; ignore = ("P2O5",))
    @test s.ignored == Dict("P2O5" => slag["P2O5"])
    @test Set(keys(s.norm)) == Set(["gehlenite", "akermanite", "pseudowollastonite"])
    @test haskey(s.unassigned, "CaO") && Jg(s.span) == 0
    @test Jg(s.vitrification) ≈ s.norm["gehlenite"] * v("Gehlenite") + s.norm["akermanite"] * v("Akermanite") +
        s.norm["pseudowollastonite"] * v("Pseudowoll") rtol = 1.0e-12
    @test 140 < Jg(s.vitrification) < 175
    @test Set(keys(s.others)) == Set(["Na2O", "K2O", "Fe2O3", "TiO2", "MnO"])

    # The oxides on the scale of another database: the difference is the moles of
    # each oxide times the difference of its enthalpies of formation, exactly.
    cem = Dict(symbol(x) => x for x in build_species(datapath("cemdata18-thermofun.json"); verbose = false))
    hc(sym) = ustrip(us"J/mol", cem[sym][:ΔₐH⁰](T = 298.15u"K", P = 1.0e5u"Pa"; unit = true))
    s_cem = glass_enthalpy(slag; ignore = ("P2O5",), reference = ([cem["Lim"], cem["Qtz"]], ChemistryLab._default_glass_reference()...))
    @test s_cem.oxide_sources["CaO"] == 1 && s_cem.oxide_sources["SiO2"] == 1 && s_cem.oxide_sources["MgO"] == 2
    nCaO, nSiO2 = slag["CaO"] / M("CaO"), slag["SiO2"] / M("SiO2")
    @test Jg(s_cem.formation_298) - Jg(s.formation_298) ≈
        nCaO * (hc("Lim") - H("Lime")) + nSiO2 * (hc("Qtz") - H("Quartz")) rtol = 1.0e-9

    # Away from 298.15 K the heat capacity is integrated; Richet (1987)
    # reproduces the relative enthalpies of his Table A-I within about 1 %.
    kept = Dict(ox => f for (ox, f) in slag if ox != "P2O5")
    s40 = glass_enthalpy(slag; ignore = ("P2O5",), T = 313.15u"K")
    xs = range(298.15, 313.15; length = 301)
    cps = [ustrip(us"J/(g*K)", glass_heat_capacity(kept; T = x * u"K")) for x in xs]
    @test Jg(s40.enthalpy) - Jg(s.enthalpy) ≈ step(xs) * (sum(cps) - (cps[1] + cps[end]) / 2) rtol = 1.0e-6

    # SO3 has no crystal in the databases: it counts as anhydrite less lime, both
    # from the first collection that holds them, in the enthalpy as in the heat
    # capacity, read here from slop98 directly.
    slop = Dict(symbol(x) => x for x in ChemistryLab._default_glass_reference()[2])
    hs(sym) = ustrip(us"J/mol", slop[sym][:ΔₐH⁰](T = 298.15u"K", P = 1.0e5u"Pa"; unit = true))
    cs(sym) = ustrip(us"J/(mol*K)", slop[sym][:Cp⁰](T = 350.0u"K", P = 1.0e5u"Pa"; unit = true))
    ref = ([slop["Anh"], slop["Lim"]], ChemistryLab._default_glass_reference()...)
    sulfated = merge(slag, Dict("SO3" => 0.02))
    nSO3 = 0.02 / M("SO3")
    s0 = glass_enthalpy(slag; ignore = ("P2O5",), reference = ref)
    s1 = glass_enthalpy(sulfated; ignore = ("P2O5",), reference = ref)
    @test s1.oxide_sources["SO3"] == 1 && Jg(s1.vitrification) == Jg(s0.vitrification)
    @test Jg(s1.formation_298) - Jg(s0.formation_298) ≈ nSO3 * (hs("Anh") - hs("Lim")) rtol = 1.0e-9
    cpg(x) = ustrip(us"J/(g*K)", glass_heat_capacity(x; T = 350.0u"K", reference = ref))
    @test cpg(merge(kept, Dict("SO3" => 0.02))) - cpg(kept) ≈ nSO3 * (cs("Anh") - cs("Lim")) rtol = 1.0e-9
    # a collection that holds lime alone is passed over: both come from one database
    @test glass_enthalpy(sulfated; ignore = ("P2O5",), reference = ([slop["Lim"]], ref...)).oxide_sources["SO3"] == 2
    # and with no collection holding both, SO3 is refused like any oxide
    @test ChemistryLab._oxide_crystals("SO3", ([slop["Lim"]], [slop["Anh"]])) === nothing
    @test_throws ArgumentError glass_heat_capacity(Dict("SO3" => 0.02); reference = ([slop["Lim"]],))
    t = literature_table("Richet1987", "relative_enthalpies")
    for i in eachindex(t.glass)
        mol = Dict(ox => ustrip(getfield(t, Symbol(ox * "_percent"))[i]) / 100 for ox in ("Na2O", "K2O", "MgO", "CaO", "Al2O3", "SiO2"))
        Mmix = sum(n * M(ox) for (ox, n) in mol)
        frac = Dict(ox => n * M(ox) / Mmix for (ox, n) in mol if n > 0)
        ys = range(273.0, 1000.0; length = 401)
        c = [ustrip(us"J/(g*K)", glass_heat_capacity(frac; T = y * u"K")) for y in ys]
        Hrel = step(ys) * (sum(c) - (c[1] + c[end]) / 2) * Mmix
        @test Hrel ≈ ustrip(us"J/mol", t.H1000_minus_H273[i]) rtol = 0.01
    end

    @test_throws ArgumentError glass_enthalpy(Dict("CaO" => -0.1, "SiO2" => 1.1))
    @test_throws ArgumentError glass_enthalpy(Dict("CaO" => 0.0))
    @test_throws ArgumentError glass_heat_capacity(Dict("CaO" => -0.1))

    # A composition being calibrated. At 298.15 K the enthalpy per gram is
    # homogeneous of degree one in the fractions (every step, the norm
    # included, is linear in them at a fixed vertex), so its derivative along
    # the analysis scaled by x is the enthalpy itself. The moles per gram were
    # stored as Float64 and a dual fraction raised.
    slag2 = Dict("CaO" => 0.416, "SiO2" => 0.352, "Al2O3" => 0.113, "MgO" => 0.06)
    H0 = Jg(glass_enthalpy(slag2).enthalpy)
    dH = ForwardDiff.derivative(x -> Jg(glass_enthalpy(Dict(k => x * v for (k, v) in slag2)).enthalpy), 1.0)
    @test dH ≈ H0 rtol = 1.0e-12
end

@testsection "a glass species with an enthalpy, for a calorimeter" begin
    slag = Dict("CaO" => 0.416, "SiO2" => 0.352, "Al2O3" => 0.113, "MgO" => 0.06)
    sp = glass_species(slag; symbol = "BFS", M = 100.0u"g/mol", enthalpy = true)
    H(T) = sp[:ΔₐH⁰](; T, unit = false)
    Cp(T) = sp[:Cp⁰](; T, unit = false)
    # per formula unit of 100 g: glass_enthalpy per gram, at 298.15 K and away
    # from it, where glass_enthalpy integrates the same heat capacity
    @test H(298.15) ≈ 100 * ustrip(us"J/g", glass_enthalpy(slag).formation_298) rtol = 1.0e-12
    @test H(350.0) ≈ 100 * ustrip(us"J/g", glass_enthalpy(slag; T = 350.0u"K").enthalpy) rtol = 1.0e-9
    @test Cp(320.0) ≈ 100 * ustrip(us"J/(g*K)", glass_heat_capacity(slag; T = 320.0u"K")) rtol = 1.0e-12
    # the temperature of a semi-adiabatic cell is differentiated: dH/dT is Cp
    for T in (298.15, 320.0, 360.0)
        @test ForwardDiff.derivative(H, T) ≈ Cp(T) rtol = 1.0e-10
    end
    # and so is a composition being calibrated: at 298.15 K the enthalpy is
    # homogeneous of degree one in the fractions
    h(x) = glass_species(Dict(k => x * v for (k, v) in slag); symbol = "B", enthalpy = true)[:ΔₐH⁰](; T = 298.15, unit = false)
    @test ForwardDiff.derivative(h, 1.0) ≈ H(298.15) rtol = 1.0e-10
    # an oxide outside Richet's table takes its crystal's heat capacity, here too
    fe = Dict("CaO" => 0.4, "SiO2" => 0.4, "MnO" => 0.05)
    @test glass_species(fe; symbol = "G", enthalpy = true)[:Cp⁰](; T = 310.0, unit = false) ≈
        100 * ustrip(us"J/(g*K)", glass_heat_capacity(fe; T = 310.0u"K")) rtol = 1.0e-12
    # none by default: a glass has no enthalpy in a database
    @test !haskey(glass_species(slag; symbol = "BFS").properties, :ΔₐH⁰)
end
