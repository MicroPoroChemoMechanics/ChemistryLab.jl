# The recipe and process layer: extents, materials, recipes and what they put
# into an equilibrium, and the processes built on them. Every bookkeeping claim is
# checked against an identity (the oxides a decomposition reproduces, the mass a
# budget conserves) or against the hand-built recipe of `scripts/gruyaert2010.jl`.

using ChemistryLab, DynamicQuantities, ForwardDiff, OrderedCollections, Test

@testsection "Recipes and processes" begin
    substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
    db = Dict(symbol(s) => s for s in substances)

    @testset "extents" begin
        @test extent(ConstantExtent(0.4), 28) == 0.4
        @test_throws ArgumentError ConstantExtent(1.2)
        tab = TabulatedExtent([1, 7, 28], [0.1, 0.4, 0.7])
        @test extent(tab, 4) ≈ 0.25
        @test extent(tab, 0.5) == 0.1 && extent(tab, 90) == 0.7
        @test extent(tab, 7u"d") ≈ 0.4
        lg = LogisticExtent(; final = 0.8, half_time = 10, slope = 1.5)
        @test extent(lg, 10) ≈ 0.4
        @test extent(CappedExtent(ConstantExtent(0.9), powers_alpha_max(0.3)), 1) ≈ 0.3 / 0.42
        pk = ParrottKillohExtent("C3S")
        a = [extent(pk, t) for t in (0.1, 1, 7, 28, 365)]
        @test issorted(a) && 0 < first(a) && last(a) <= 1
        # The tabulated curve is the integral of the law: its slope is the rate,
        # to the linear interpolation between grid points 0.4 % apart in time.
        rate = ChemistryLab.parrott_killoh_avrami(ChemistryLab._pk84_params("C3S"), "C3S")
        t = 7.0
        slope = ForwardDiff.derivative(τ -> extent(pk, τ), t) / 86400
        α = extent(pk, t)
        @test slope ≈ rate(293.15, 1.0e5, 0.0, Dict("C3S" => 1 - α), nothing, Dict("C3S" => 1.0)) rtol = 1.0e-2
        # The w/c factor of Parrott and Killoh: no effect below 1.333 w/c, and a
        # hydration that stops at (1 + 4.444 w/c)/3.333.
        pk3 = ParrottKillohExtent("C3S"; w_c = 0.3)
        @test extent(pk3, 0.5) == extent(pk, 0.5)
        @test extent(pk3, 3650) < extent(pk, 3650)
        @test extent(pk3, 3650) <= (1 + 4.444 * 0.3) / 3.333
        # The constants and the critical degree a paper fits. Given the default
        # constants, the law is the default one, to the last bit; a larger
        # critical degree lets the phase hydrate further at the same w/c; a faster
        # diffusion stage of belite (Lothenbach et al. 2008) hydrates it further.
        same = ParrottKillohExtent("C3S"; parameters = (k₁ = 1.5, n₁ = 0.7, k₂ = 0.05, k₃ = 1.1, n₃ = 3.3))
        @test [extent(same, t) for t in (1, 28, 365)] ≈ [extent(pk, t) for t in (1, 28, 365)] rtol = 1.0e-12
        # The rate constants are per day, as printed (read as per second they
        # hydrated the phase within the first day).
        @test extent(same, 1) < 0.5
        @test extent(ParrottKillohExtent("C3S"; parameters = (k₁ = 1.5u"1/d",)), 1) ≈ extent(pk, 1) rtol = 1.0e-12
        @test extent(ParrottKillohExtent("C3S"; w_c = 0.3, H = 1.8), 3650) > extent(pk3, 3650)
        @test extent(ParrottKillohExtent("C3S"; w_c = 0.3, H = 1.8), 0.5) == extent(pk, 0.5)
        c2s = ParrottKillohExtent("C2S")
        @test extent(c2s, 365) < extent(ParrottKillohExtent("C2S"; parameters = (k₂ = 0.02, k₃ = 0.7)), 365) < 1
        @test_throws ArgumentError ParrottKillohExtent("C3S"; parameters = (k4 = 1.0,))
    end

    @testset "oxides of a phase, Bogue, and decompositions" begin
        c = oxide_content(db["C3S"], ("CaO", "SiO2"))
        @test sum(values(c)) ≈ 1 rtol = 1.0e-12
        ox = literature_oxides("Durdzinski2017", "chemical_composition", "PC")
        @test ox["CaO"] ≈ 0.637
        wide = literature_oxides("Gruyaert2010", "oxides", "OPC-CAL")
        @test wide["CaO"] ≈ 0.622 && !haskey(wide, "blaine")
        bg = bogue(ox, db; sulfate = "Anh")
        # The phases reproduce the six oxides they are computed from, exactly.
        for o in ("CaO", "SiO2", "Al2O3", "Fe2O3", "SO3")
            @test sum(f * oxide_content(db[p], (o,))[o] for (p, f) in bg.phases) ≈ ox[o] rtol = 1.0e-10
        end
        # Gypsum would bring water the analysis counts in its loss on ignition.
        @test sum(values(bogue(ox, db).phases)) > sum(values(bg.phases))
        @test_throws ArgumentError bogue(ox, db; sulfate = "Cal")
        @test haskey(bg.remainder, "MgO") && !haskey(bg.remainder, "CaO")
        # A decomposition recovers the fractions an analysis was made from.
        truth = OrderedDict("C3S" => 0.6, "C2S" => 0.2, "C3A" => 0.1, "Cal" => 0.1)
        made = OrderedDict(o => sum(f * oxide_content(db[p], (o,))[o] for (p, f) in truth) for o in ("CaO", "SiO2", "Al2O3", "CO2"))
        d = decompose(made, collect(keys(truth)), db)
        @test all(isapprox(d.fractions[p], truth[p]; atol = 1.0e-10) for p in keys(truth))
        @test d.misfit < 1.0e-12
        # And when no non-negative combination fits, the nearest one, with its misfit.
        d2 = decompose(OrderedDict("CaO" => 0.2, "SiO2" => 0.8), ["C3S", "Qtz"], db)
        @test all(>=(0), values(d2.fractions))
        # The glass by difference and its crystals give back the analysis.
        sfa = literature_oxides("Durdzinski2017", "chemical_composition", "SFA")
        cr = OrderedDict("Qtz" => 0.149)
        g = reactive_part(sfa, cr, db)
        @test g.mass_fraction ≈ 0.851
        @test g.mass_fraction * g.oxides["SiO2"] + 0.149 * oxide_content(db["Qtz"], ("SiO2",))["SiO2"] ≈ sfa["SiO2"]
    end

    @testset "materials refuse what they cannot be" begin
        @test_throws ArgumentError OxideConstituent("g", Dict("CaO" => 40.0))     # percent, not fractions
        @test_throws ArgumentError OxideConstituent("g", Dict("NaCl" => 0.1))     # not an oxide
        @test_throws ArgumentError Material(
            "x", :cement; constituents = [
                MineralConstituent(db["C3S"]; mass_fraction = 0.7), MineralConstituent(db["C2S"]; mass_fraction = 0.4),
            ]
        )
        @test_throws ArgumentError Material(
            "salt", :salt; constituents = [
                MineralConstituent(db["Gp"]; mass_fraction = 1.0, extent = 0.5),
            ]
        )
        @test_throws ArgumentError Material("x", :nonsense; constituents = [MineralConstituent(db["C3S"]; mass_fraction = 1.0)])
        @test "S1 slag (Durdzinski 2017)" in material_templates()
        @test_throws KeyError material_template("no such material", db)
        @test_throws ErrorException ChemistryLab._material_from_entry(Dict("name" => "incomplete"), db)
        # The fly ashes: the glass by difference beside its crystals, which are
        # inert in the siliceous one and reactive in the Rietveld route.
        sfa = material_template("siliceous fly ash (Durdzinski 2017)", db)
        glass = only(c for c in sfa.constituents if c.name == "glass")
        @test glass.mass_fraction ≈ 1 - 0.149 - 0.193
        @test all(c -> extent(c.extent, 28) == 0, (c for c in sfa.constituents if c.name != "glass"))
        cfa = material_template("calcareous fly ash (Durdzinski 2017)", db)
        @test only(c for c in cfa.constituents if c.name == "glass").mass_fraction ≈ 0.897
        @test sum(c.mass_fraction for c in cfa.constituents) ≈ 1 rtol = 1.0e-9
        # A clinker by its phases and the rest of its analysis; a fly ash with a
        # hematite; a limestone whose calcite (by TGA) holds more lime than its
        # analysis (by XRF) does, the oxides left still fractions of their mass.
        clinker = material_template("clinker (De Weerdt 2011)", db)
        @test [c.name for c in clinker.constituents] == ["C2S", "C3S", "C3A", "C4AF", "minor oxides"]
        @test only(c for c in clinker.constituents if c.name == "minor oxides").mass_fraction ≈ 0.08
        ash = material_template("siliceous fly ash (De Weerdt 2011)", db)
        @test "Hematite" in [c.name for c in ash.constituents]
        @test only(c for c in ash.constituents if c.name == "glass").mass_fraction ≈ 0.68
        rest = only(c for c in material_template("limestone (De Weerdt 2011)", db).constituents if c.name == "minor oxides")
        @test sum(values(rest.oxides)) ≈ 1 && rest.mass_fraction ≈ 0.19
        @test_throws ErrorException ChemistryLab._material_from_entry(
            Dict("name" => "x", "phases" => "DeWeerdt2011:mineral_composition:clinker", "remainder" => true), db,
        )
        # A clinker whose phases, counted with their pure formulas, leave less mass
        # than its minor oxides: `remainder = "analysis"` keeps the oxides no phase
        # holds at the amounts of the analysis, the phases scaled to make room.
        # The white cement of Shi et al. (2016) leaves 1.3 % for 3.9 %.
        wpc = material_template("white Portland cement (Shi 2016)", db)
        @test sum(c.mass_fraction for c in wpc.constituents) ≈ 1 rtol = 1.0e-12
        minor = only(c for c in wpc.constituents if c.name == "minor oxides")
        xrf = literature_oxides("Shi2016", "chemical_composition", "wPc")
        for ox in ("K2O", "Na2O", "MgO", "TiO2")
            @test minor.mass_fraction * minor.oxides[ox] ≈ xrf[ox] rtol = 1.0e-12
        end
        @test sum(values(minor.oxides)) ≈ 1 rtol = 1.0e-12
        c3s = only(c for c in wpc.constituents if c.name == "C3S").mass_fraction
        @test 0.95 * 0.649 < c3s < 0.649
        # The materials of Schöler et al. (2015). The glass of the slag, found by
        # difference, is the authors' own estimate of it (their Table 3) to within
        # a percent of each major oxide, once both are brought to 100 %: the
        # authors normalize theirs, and the template does not, so as not to
        # invent the part of the analysis it does not report.
        bfs = material_template("blast-furnace slag (Schöler 2015)", db)
        g15 = only(c for c in bfs.constituents if c.name == "glass")
        @test g15.mass_fraction ≈ 0.985 rtol = 1.0e-12
        t3 = literature_row("Scholer2015", "glass_composition", "BFS")
        total = sum(values(g15.oxides))
        @test 0.95 < total < 1
        for ox in ("SiO2", "CaO", "Al2O3", "MgO")
            @test g15.oxides[ox] / total ≈ ustrip(getproperty(t3, Symbol(ox))) / 100 atol = 0.01
        end
        opc15 = material_template("OPC (Schöler 2015)", db)
        @test sum(c.mass_fraction for c in opc15.constituents) ≈ 1 rtol = 1.0e-12
        # The two polymorphs of C2S and of C3A of the Rietveld analysis are one
        # constituent each, named by the database symbol.
        names15 = [c.name for c in opc15.constituents]
        @test issubset(["C3S", "C2S", "C3A", "C4AF", "Bassanite", "Syngenite"], names15)
        @test allunique(names15) && !any(n -> occursin("C2S", n) && n != "C2S", names15)
        p15 = ChemistryLab._literature_phases("Scholer2015:phases:OPC")
        c2s = only(c for c in opc15.constituents if c.name == "C2S").mass_fraction
        c3s = only(c for c in opc15.constituents if c.name == "C3S").mass_fraction
        @test c2s / c3s ≈ (p15["alpha' C2S"] + p15["beta C2S"]) / p15["C3S"] rtol = 1.0e-12
        fa15 = material_template("siliceous fly ash (Schöler 2015)", db)
        @test only(c for c in fa15.constituents if c.name == "glass").mass_fraction ≈ 0.687 rtol = 1.0e-12
        @test material_template("limestone (Schöler 2015)", db).constituents[1] isa ChemistryLab.OxideConstituent
        @test material_template("CEM I 52.5 N (Gruyaert 2010), Bogue", db).kind === :cement
        # With `remainder = true` the same oxides are squeezed into what the phases
        # leave, and the potassium falls to a third of the analysis.
        crystals = ChemistryLab._literature_phases("Shi2016:phase_composition:wPc")
        squeezed = ChemistryLab._glass("Shi2016:chemical_composition:wPc", crystals, db)
        @test squeezed.mass_fraction * squeezed.oxides["K2O"] < 0.4 * xrf["K2O"]
        # Where the phases leave room (the clinker of De Weerdt et al.), the two
        # rules give the same material.
        dw = ChemistryLab._material_from_entry(
            merge(ChemistryLab._template_entry("clinker (De Weerdt 2011)"), Dict("remainder" => "analysis")), db,
        )
        @test [(c.name, c.mass_fraction) for c in dw.constituents] == [(c.name, c.mass_fraction) for c in clinker.constituents]
        dwm = only(c for c in dw.constituents if c.name == "minor oxides")
        @test all(dwm.oxides[k] ≈ v for (k, v) in only(c for c in clinker.constituents if c.name == "minor oxides").oxides)
    end

    @testset "the non-negative least squares meets its optimality conditions" begin
        # x ≥ 0, the gradient Aᵀ(Ax − b) zero where x > 0 and non-negative where
        # x = 0: the conditions of the minimum, whatever path found it.
        seed = UInt64(0x2026_0927_0000_0003)
        rnd() = (seed ⊻= seed << 13; seed ⊻= seed >> 7; seed ⊻= seed << 17; (seed >> 11) / Float64(1 << 53))
        for _ in 1:50
            A = [rnd() for _ in 1:4, _ in 1:6]
            b = [rnd() - 0.3 for _ in 1:4]
            x = ChemistryLab._nnls(A, b)
            g = transpose(A) * (A * x - b)
            @test all(>=(0), x)
            @test all(abs(g[j]) < 1.0e-9 for j in eachindex(x) if x[j] > 0)
            @test all(g[j] > -1.0e-9 for j in eachindex(x) if x[j] == 0)
        end
    end

    # A small paste: clinker, gypsum, slag, water.
    pure = split(
        "C3S C2S C3A C4AF Gp Anh Cal Portlandite ettringite monosulphate12 " *
            "monocarbonate hemicarbonate C4AH13 C3AH6 C3FH6 straetlingite " *
            "hydrotalcite Mg2AlC0.5OH Brc FeOOHmic AlOHmic Amor-Sl K2SO4 Qtz"
    )
    csh = ["CSHQ-JenD", "CSHQ-JenH", "CSHQ-TobD", "CSHQ-TobH", "KSiOH", "NaSiOH"]
    sp = speciation(substances, vcat(pure, csh); aggregate_state = [AS_AQUEOUS], exclude_species = split("H2@ O2@ CH4@"))
    cs = ChemicalSystem(sp, CEMDATA_PRIMARIES; solid_solutions = [SolidSolutionPhase("CSHQ", [db[m] for m in csh])])
    A = Float64.(cs.SM.A)
    model = cemdata18_activity_model(:KOH)

    pc = with_extents(
        material_template("PC (Durdzinski 2017), Bogue", cs),
        Dict("C3S" => 0.8, "C2S" => 0.5, "C3A" => 0.9, "C4AF" => 0.6),
    )
    slag = material_template("S1 slag (Durdzinski 2017)", cs)
    slag = with_extents(slag, Dict{String, Any}(); material_extent = 0.4)

    @testset "a recipe's budget conserves what it is given" begin
        @test_throws ArgumentError Recipe(pc => 0.7, slag => 0.2; w_b = 0.4)
        r = Recipe(pc => 0.0, slag => 0.4; w_b = 0.4, balance = pc.name)
        @test last(r.binder[1]) ≈ 0.6
        bud = budget(r, cs)
        # Every oxygen-free element the materials hold: in the budget as reacted,
        # in the residue as not.
        reacted_Ca = bud.b[findfirst(==("Ca+2"), [symbol(p) for p in cs.SM.primaries])]
        total_Ca = 0.0
        for (m, mass) in ChemistryLab._material_masses(r), c in m.constituents
            mc = mass * c.mass_fraction
            α = effective_extent(m, c, nothing)
            if c isa MineralConstituent
                total_Ca += α * mc * get(atoms(c.species), :Ca, 0) / ustrip(us"g/mol", c.species[:M])
            else
                total_Ca += α * mc * get(c.oxides, "CaO", 0.0) / ChemistryLab._oxide_molar_mass("CaO")
            end
        end
        @test reacted_Ca ≈ total_Ca rtol = 1.0e-12
        unreacted = [x for x in bud.residual if x.reason === :unreacted]
        @test sum(x.mass for x in unreacted) ≈
            sum(mass * c.mass_fraction * (1 - effective_extent(m, c, nothing)) for (m, mass) in ChemistryLab._material_masses(r) for c in m.constituents) rtol = 1.0e-12
        # The titanium of the cement has no primary in this system: its reacted
        # part is kept aside, and said so, rather than refused or dropped.
        aside = [x for x in bud.residual if x.reason === :not_in_system]
        @test any(x -> endswith(x.constituent, "TiO2"), aside)
        # A mineral of the residue has its volume from its molar volume.
        @test all(x -> x.volume !== nothing, (x for x in unreacted if x.constituent == "C3S"))
        # A glass has no density unless its source gives one: its volume is missing.
        @test any(x -> x.volume === nothing, bud.residual)
    end

    @testset "the layer poses the problem the hand-built recipe posed" begin
        # `scripts/gruyaert2010.jl` builds a CEM I 52.5 N and slag paste by hand;
        # the same materials through the layer give the same element budget, and
        # so the same certified answer.
        include(joinpath(pkgdir(ChemistryLab), "scripts", "gruyaert2010.jl"))
        gcs = gruyaert_system()
        c = gruyaert_bogue("OPC-CAL")
        α = 0.74
        opc = Material(
            "OPC-CAL", :cement;
            constituents = vcat(
                [MineralConstituent(G10_DB[p]; mass_fraction = f, extent = α) for (p, f) in c.clinker],
                [MineralConstituent(G10_DB["Gp"]; mass_fraction = c.gypsum), MineralConstituent(G10_DB["Cal"]; mass_fraction = c.calcite)],
            ),
        )
        bfs = oxide_material("BFS-CAL", filter(p -> first(p) in ("CaO", "SiO2", "Al2O3", "Fe2O3", "MgO", "SO3"), literature_oxides("Gruyaert2010", "oxides", "BFS-CAL")); extent = 0.3)
        r = Recipe(opc => 0.5, bfs => 0.5; w_b = 0.5, T = 293.15u"K")
        bud = budget(r, gcs)
        _, hand = gruyaert_budget(gcs; slag = 0.5, alpha_cement = α, alpha_slag = 0.3)
        @test bud.b ≈ hand rtol = 1.0e-12
    end

    @testset "a recipe at equilibrium, and what is read off it" begin
        r = Recipe(pc => 0.7, slag => 0.3; w_b = 0.45)
        rs, cert = equilibrate_certified(r, cs; model)
        @test cert.optimal
        @test A * ustrip.(us"mol", rs.state.n) ≈ rs.b rtol = 1.0e-10 atol = 1.0e-12
        pm = phase_masses(rs)
        @test pm["Portlandite"] > 0
        @test any(startswith("unreacted "), keys(pm))
        # Started from that answer, the same paste at 5 °C is solved at 5 °C:
        # a start carries its amounts, not its temperature.
        cold = Recipe(pc => 0.7, slag => 0.3; w_b = 0.45, T = 278.15u"K")
        rs5, cert5 = equilibrate_certified(cold, cs; model, start = rs.state)
        @test cert5.optimal
        @test temperature(rs5.state) == 278.15u"K"
        # The slag glass has no sourced density: the volume of its residue, and
        # with it the porosity, are missing rather than estimated.
        @test isnan(porosity(rs).total)
        # Given a density (a test value, not a published one), both are defined.
        slag_d = oxide_material(
            "S1 with a density", literature_oxides("Durdzinski2017", "chemical_composition", "S1");
            density = 2.9, extent = 0.4,
        )
        # The cement by its phases, every element of which the system holds.
        pcp = with_extents(
            material_template("PC (Durdzinski 2017), phases", cs),
            Dict("C3S" => 0.8, "C2S" => 0.45, "C3A" => 0.9, "C4AF" => 0.6),
        )
        rs_d, _ = equilibrate_certified(Recipe(pcp => 0.7, slag_d => 0.3; w_b = 0.45), cs; model)
        @test isempty(volume(rs_d).missing)
        p = porosity(rs_d)
        @test 0 < p.total < 1 && p.void >= 0
        @test volume(rs_d).residual ≈ sum(x.mass for x in rs_d.residual if x.constituent == "S1 with a density") / 2.9 +
            sum(x.volume for x in rs_d.residual if x.constituent != "S1 with a density" && x.volume !== nothing) rtol = 1.0e-12
        # The fractions of the paste: relative to its initial volume, residue
        # included, and closed by the void of `porosity`.
        vf = volume_fractions(rs_d)
        @test sum(values(vf)) ≈ 1 rtol = 1.0e-12
        @test vf["void"] ≈ p.void rtol = 1.0e-10 atol = 1.0e-14
        V0 = volume(rs_d.initial).total
        @test sum(v for (k, v) in vf if startswith(k, "unreacted ")) ≈ volume(rs_d).residual / (ChemistryLab._in_unit(us"cm^3", V0) + volume(rs_d).residual) rtol = 1.0e-12
        @test vf["Portlandite"] > 0
        @test_throws ArgumentError volume_fractions(rs)
        @test 0 < bound_water(rs) < 0.4
        # Over a temperature window, from decomposition windows (test values,
        # not published ones): sharp steps, so a range holding every midpoint
        # recovers the whole, and one holding some of them recovers theirs.
        held = vcat(
            [p.first => ChemistryLab._in_unit(us"g", p.second) for p in bound_water_per_phase(rs.state)],
            ChemistryLab._unreacted_water(rs),
        )
        mids = Dict(ph => 400.0 + 10k for (k, (ph, _)) in enumerate(held))
        steps = [DecompositionWindow(ph, T, 1.0e-3) for (ph, T) in mids]
        @test bound_water(rs; window = (300.0, 2000.0), windows = steps) ≈ bound_water(rs) rtol = 1.0e-12
        cut = 400.0 + 10 * (length(held) ÷ 2) + 5
        @test bound_water(rs; window = (300.0u"K", cut * u"K"), windows = steps) ≈
            sum(g for (ph, g) in held if mids[ph] < cut) / rs.recipe.binder_mass rtol = 1.0e-12
        big = first(held).first
        @test_throws ArgumentError bound_water(rs; window = (300.0, 2000.0), windows = filter(w -> w.phase != big, steps))
        @test_throws ArgumentError bound_water(rs; window = (300.0, 2000.0))
        @test_throws ArgumentError bound_water(rs; windows = steps)
        @test_throws ArgumentError bound_water(rs; window = (500.0, 400.0), windows = steps)
        # An unreacted mineral that holds water counts with it: gypsum half
        # reacted (a test extent), whose unreacted half keeps its two waters.
        gypsum = Material("gypsum", :other; constituents = [MineralConstituent(db["Gp"]; mass_fraction = 1.0, extent = 0.5)])
        rg, _ = equilibrate_certified(Recipe(pc => 0.95, gypsum => 0.05; w_b = 0.45), cs; model)
        Mw = ustrip(us"g/mol", Species("H2O")[:M])
        unreacted = 0.5 * 0.05 * rg.recipe.binder_mass / ustrip(us"g/mol", db["Gp"][:M]) * 2 * Mw
        @test bound_water(rg) ≈ (ChemistryLab._in_unit(us"g", ignition_loss(rg.state).water) + unreacted) / rg.recipe.binder_mass rtol = 1.0e-12
        ps = pore_solution(rs)
        @test ps.pH > 12 && ps.elements[:K] > 0
        # The residue has no sourced enthalpy (the glass): the heat is not complete.
        @test isnan(enthalpy(rs))
        @test isnan(volume(rs).residual) && "S1 slag (Durdzinski 2017)" in volume(rs).missing

        # The heat between two states of one paste counts the residue by what
        # changed. The cement alone: its titanium is set aside, unsourced, with
        # the same mass in both states, so it cancels; the unreacted clinker
        # changes, with the enthalpy of its records, and counts.
        early = with_extents(pc, Dict("C3S" => 0.4, "C2S" => 0.1, "C3A" => 0.5, "C4AF" => 0.3))
        rs1, _ = equilibrate_certified(Recipe(early => 1.0; w_b = 0.45), cs; model)
        rs2, _ = equilibrate_certified(Recipe(pc => 1.0; w_b = 0.45), cs; model)
        @test isnan(enthalpy(rs1)) && any(x -> x.enthalpy === nothing && x.reason === :not_in_system, rs1.residual)
        q = heat_release(rs1, rs2)
        sourced(r) = sum(x.enthalpy for x in r.residual if x.enthalpy !== nothing)
        H(r) = ustrip(us"J", enthalpy(r.state))
        @test q ≈ -((H(rs2) + sourced(rs2)) - (H(rs1) + sourced(rs1))) rtol = 1.0e-12
        @test 0 < q / 100 < 600   # J per g of cement: a heat of hydration's order
        @test heat_release(rs2, rs2) == 0
        # With the slag glass reacting further between the two, its unreacted
        # mass changes without an enthalpy: no heat can be given.
        slag6 = with_extents(slag, Dict{String, Any}(); material_extent = 0.6)
        rs6, _ = equilibrate_certified(Recipe(pc => 0.7, slag6 => 0.3; w_b = 0.45), cs; model)
        @test isnan(heat_release(rs, rs6))
        hot, _ = equilibrate_certified(Recipe(pc => 1.0; w_b = 0.45, T = 313.15u"K"), cs; model)
        @test_throws ArgumentError heat_release(rs2, hot)
    end

    @testset "processes" begin
        r = Recipe(pc => 1.0; w_b = 0.45)
        pc_t = with_extents(
            material_template("PC (Durdzinski 2017), Bogue", cs),
            Dict("C3S" => TabulatedExtent([1, 28], [0.3, 0.8]), "C2S" => TabulatedExtent([1, 28], [0.05, 0.4])),
        )
        h = hydrate(Recipe(pc_t => 1.0; w_b = 0.45), cs, [1, 7, 28]; model)
        @test all(rs -> rs.certificate.optimal, h)
        ch = [phase_masses(rs)["Portlandite"] for rs in h]
        @test issorted(ch)
        tbl = process_table(h; phases = ["Portlandite", "ettringite"])
        @test size(tbl, 1) == 3 && tbl.Portlandite == ch
        @test h[end] === h[3] && h[begin] === h[1]
        @test findfirst(rs -> phase_masses(rs)["Portlandite"] > ch[1], h) == 2

        bl = blend(r, slag, [0.0, 0.3], cs; model)
        @test all(rs -> rs.certificate.optimal, bl)
        @test phase_masses(bl[2])["Portlandite"] < phase_masses(bl[1])["Portlandite"]

        rs, _ = equilibrate_certified(r, cs; model)
        co = carbonate(rs, [0.0, 0.2, 0.6])
        @test all(s -> s.certificate.optimal, co)
        @test pH(co[3].state, model) < pH(co[1].state, model)
        for s in co
            @test A * ustrip.(us"mol", s.state.n) ≈ s.b rtol = 1.0e-10 atol = 1.0e-12
        end

        # A salt enters by its formula, which need not be a species of the system:
        # calcium carbonate as the column of calcite. One the primaries cannot
        # express (no chlorine here), or no formula at all, is refused.
        ical = findfirst(s -> symbol(s) == "Cal", cs.species)
        sa = add_salt(rs, "CaCO3", [0.05])
        @test only(sa.states).certificate.optimal
        @test only(sa.states).b ≈ rs.b .+ 0.05 .* A[:, ical] rtol = 1.0e-14
        @test_throws ArgumentError add_salt(rs, "NaCl", [0.1])
        @test_throws ArgumentError titrate(rs, "not a formula", [0.1])

        le = leach(rs, 2)
        @test all(s -> s.certificate.optimal, le)
        @test pore_solution(le[2]).elements[:K] < pore_solution(rs).elements[:K]
    end

    @testset "derivatives through a recipe are exact" begin
        iw = findfirst(s -> symbol(s) == "H2O@", cs.species)
        Mw = ustrip(us"g/mol", cs.species[iw][:M])
        # Water enters the budget as water, with the H⁺ and OH⁻ of neutral pH a
        # state seeds in proportion to it; the reacted part of a mineral as that
        # mineral, and its unreacted part leaves the residue as it reacts.
        bw(w) = budget(Recipe(pc => 0.7, slag => 0.3; w_b = w), cs).b
        st45 = budget(Recipe(pc => 0.7, slag => 0.3; w_b = 0.45), cs).state
        iH, iOH = (findfirst(s -> symbol(s) == x, cs.species) for x in ("H+", "OH-"))
        seed = ustrip(us"mol", moles(st45, "H+")) / ustrip(us"mol", moles(st45, "H2O@"))
        @test seed > 0 && moles(st45, "OH-") == moles(st45, "H+")
        dbw = ForwardDiff.derivative(bw, 0.45)
        @test dbw ≈ (A[:, iw] .+ seed .* (A[:, iH] .+ A[:, iOH])) * 100 / Mw rtol = 1.0e-12
        ic = findfirst(s -> symbol(s) == "C3S", cs.species)
        f = only(c.mass_fraction for c in pc.constituents if c.name == "C3S")
        Mc = ustrip(us"g/mol", cs.species[ic][:M])
        at(α) = budget(Recipe(with_extents(pc, Dict("C3S" => α)) => 1.0; w_b = 0.45), cs)
        @test ForwardDiff.derivative(α -> at(α).b, 0.8) ≈ A[:, ic] * 100 * f / Mc rtol = 1.0e-12
        @test ForwardDiff.derivative(α -> sum(x.mass for x in at(α).residual if x.constituent == "C3S"), 0.8) ≈ -100 * f rtol = 1.0e-12

        # The law of an extent: with every rate constant scaled by λ, the law is
        # the same in a time scaled by λ, and so is it at a temperature `T` in
        # the time scaled by its Arrhenius factor. At λ = 1 and T = T_ref the
        # derivatives are the time times the rate, to the tabulation.
        p0 = ChemistryLab._pk84_params("C3S")
        t = 7.0
        α = extent(ParrottKillohExtent("C3S"), t)
        rate = ChemistryLab.parrott_killoh_avrami(p0, "C3S")
        α̇(T, a) = rate(T, 1.0e5, 0.0, Dict("C3S" => 1 - a), nothing, Dict("C3S" => 1.0))
        scaled(λ) = extent(ParrottKillohExtent("C3S"; parameters = (k₁ = λ * p0.k₁, k₂ = λ * p0.k₂, k₃ = λ * p0.k₃)), t)
        @test ForwardDiff.derivative(scaled, 1.0) ≈ t * 86400 * α̇(293.15, α) rtol = 1.0e-3
        Tr = ChemistryLab.safe_ustrip(us"K", p0.T_ref)
        Ea = ChemistryLab.safe_ustrip(us"J/mol", p0.Ea)
        αr = extent(ParrottKillohExtent("C3S"; T = Tr), t)
        dT = ForwardDiff.derivative(T -> extent(ParrottKillohExtent("C3S"; T), t), Tr)
        @test dT ≈ t * 86400 * Ea / (ChemistryLab.R_GAS * Tr^2) * α̇(Tr, αr) rtol = 1.0e-3

        # Through the equilibrium: the answer keeps the mass balance in its
        # derivatives, and is the derivative of the same equilibrium in its
        # budget, by the chain rule. The readers carry the derivatives along.
        r(w) = Recipe(pc => 0.7, slag => 0.3; w_b = w)
        rs0, _ = equilibrate_certified(r(0.45), cs; model)
        ich = findfirst(s -> symbol(s) == "Portlandite", cs.species)
        Mch = ustrip(us"g/mol", cs.species[ich][:M])
        read_off(rs) = vcat(
            ustrip.(us"mol", rs.state.n), phase_masses(rs)["Portlandite"], bound_water(rs),
            pore_solution(rs).elements[:K], residual_mass(rs),
        )
        out = w -> read_off(first(equilibrate_certified(r(w), cs; model)))
        Tg = typeof(ForwardDiff.Tag(out, Float64))
        lifted = out(ForwardDiff.Dual{Tg}(0.45, 1.0))
        # The values are those of the plain recipe: the solve ran on them.
        @test ForwardDiff.value.(lifted) ≈ read_off(rs0) rtol = 1.0e-12
        d = ForwardDiff.partials.(lifted, 1)
        ns = length(cs.species)
        dn = d[1:ns]
        @test A * dn ≈ dbw rtol = 1.0e-8 atol = 1.0e-12
        dn_b = ForwardDiff.derivative(
            s -> ustrip.(us"mol", first(equilibrate_certified(rs0.initial; model, b = rs0.b .+ s .* dbw)).n), 0.0,
        )
        @test dn ≈ dn_b rtol = 1.0e-10 atol = 1.0e-14
        @test d[ns + 1] ≈ dn[ich] * Mch rtol = 1.0e-12
        @test d[ns + 4] == 0   # the residue does not depend on the water
    end

    @testset "a recipe as a kinetic problem" begin
        # The clinker silicates given their rates, everything else as the recipe
        # has it at the start: C3S and C2S whole, the other constituents reacted
        # as `budget` takes them, the oxides of the cement as the primaries that
        # carry them.
        rates = Dict("C3S" => parrott_killoh_avrami(PK84_PARAMS_C3S, "C3S"), "C2S" => parrott_killoh_avrami(PK84_PARAMS_C2S, "C2S"))
        r = Recipe(pc => 1.0; w_b = 0.45)
        kp = KineticsProblem(r, cs, rates, (0.0, 86400.0))
        whole = with_extents(pc, Dict("C3S" => 1.0, "C2S" => 1.0))
        n0 = ustrip.(us"mol", kp.initial_state.n)
        @test A * n0 ≈ budget(Recipe(whole => 1.0; w_b = 0.45), cs).b rtol = 1.0e-12 atol = 1.0e-12
        @test all(>=(0), n0)
        c3s = only(c for c in pc.constituents if c.name == "C3S")
        @test ustrip(us"mol", moles(kp.initial_state, "C3S")) ≈
            r.binder_mass * c3s.mass_fraction / ustrip(us"g/mol", c3s.species[:M]) rtol = 1.0e-12
        # One dissolution per rate, each conserving the elements.
        @test length(kp.kinetic_reactions) == 2
        @test maximum(abs, A * transpose(kp.ν)) < 1.0e-12
        # The kinetics alone (no equilibrium solver): the silicates dissolve.
        sol = integrate(kp, KineticsSolver())
        ic3s = findfirst(s -> symbol(s) == "C3S", cs.species)
        k3 = findfirst(==(ic3s), kp.idx_kinetic)
        @test sol.u[end][k3] < sol.u[1][k3]
        # Dissolved into the primaries without a partition, they consume `H+`
        # below zero, which is hydroxide: the total of a component, not a
        # co-reactant run out, and the run is a success.
        p = sol.prob.p
        iH = findfirst(s -> symbol(s) == "H+", cs.species)
        kH = findfirst(==(iH), p.idx_equilibrium)
        @test p.n_initial_full[iH] + sum(p.νe[j, kH] * sol.u[end][p.n_nk + j] for j in 1:(p.n_rxn_state)) < 0
        @test SciMLBase.successful_retcode(sol)
        # The reacted oxides enter as the primaries: a basic oxide's protons as
        # hydroxide, an acidic oxide's water taken from the mixing water, and
        # anything else refused.
        st = deepcopy(kp.initial_state)
        ox(f) = ChemistryLab.oxide_budget(OrderedDict(f => 1.0), cs.SM.primaries; mass = 1.0u"g")
        oh0, w0 = moles(st, "OH-"), moles(st, "H2O@")
        ChemistryLab._add_primaries!(st, ox("K2O"))
        @test moles(st, "OH-") > oh0 && moles(st, "H2O@") < w0
        w1 = moles(st, "H2O@")
        ChemistryLab._add_primaries!(st, ox("SO3"))
        @test moles(st, "H2O@") < w1
        @test A * ustrip.(us"mol", st.n) ≈ A * n0 .+ ox("K2O") .+ ox("SO3") rtol = 1.0e-12 atol = 1.0e-14
        neg = zeros(length(cs.SM.primaries)); neg[findfirst(p -> symbol(p) == "Ca+2", cs.SM.primaries)] = -1.0
        @test_throws ArgumentError ChemistryLab._add_primaries!(deepcopy(kp.initial_state), neg)
        dry = zeros(length(cs.SM.primaries)); dry[findfirst(p -> symbol(p) == "H2O@", cs.SM.primaries)] = -1.0e6
        @test_throws ArgumentError ChemistryLab._add_primaries!(deepcopy(kp.initial_state), dry)
        # What cannot be given a rate is refused by name.
        @test_throws ArgumentError KineticsProblem(r, cs, Dict("no such" => rates["C3S"]), (0.0, 1.0))
        glass = first(c.name for c in slag.constituents if c isa ChemistryLab.OxideConstituent)
        err = try
            KineticsProblem(Recipe(pc => 0.7, slag => 0.3; w_b = 0.45), cs, Dict(glass => rates["C3S"]), (0.0, 1.0))
        catch e
            e
        end
        @test err isa ArgumentError
        # The refusal names the way round it.
        @test occursin("glass_species", sprint(showerror, err))
    end

    @testset "phase lists" begin
        # Every symbol a list names is a record of the database it is written for,
        # and every solid solution a phase of data/solid_solutions.toml.
        declared = [e["name"] for e in ChemistryLab.TOML.parsefile(datapath("solid_solutions.toml"))["solid_solution"]]
        for list in phase_lists()
            pl = phase_list(list)
            names = Set(symbol.(build_species(datapath(pl.database); verbose = false)))
            @test isempty(setdiff(vcat(pl.reactants, pl.products, pl.exclude_aqueous), names))
            @test isempty(setdiff(pl.solid_solutions, declared))
            @test isempty(setdiff(keys(pl.instances), pl.solid_solutions))
        end
        @test_throws KeyError phase_list("no such paste")

        list = "Portland paste (Lothenbach and Winnefeld 2006)"
        pl = phase_list(list)
        pcs = phase_list_system(list, substances)
        syms = String.(symbol.(pcs.species))
        @test issubset(vcat(pl.reactants, pl.products), syms)
        @test isempty(intersect(pl.exclude_aqueous, syms))
        # The list declares the AFm binary with one composition and room for a
        # second, where data/solid_solutions.toml declares two.
        afm = only(p for p in pcs.solid_solutions if name(p) == "AFm_SO4_OH")
        @test afm.instances == 1 && afm.max_instances == 2
        filed = only(p for p in build_solid_solutions(datapath("solid_solutions.toml"), db) if name(p) == "AFm_SO4_OH")
        @test filed.instances == 2
        @test issubset(String.(symbol.(end_members(afm))), syms)
        @test Set(name.(pcs.solid_solutions)) == Set(pl.solid_solutions)

        # Departures from the list are named: a phase added, a phase removed, and
        # a removal of what the list does not hold refused.
        friedel = phase_list_system(list, substances; add = ["C4AClH10"], remove = ["hemicarbonate"])
        fsyms = String.(symbol.(friedel.species))
        @test "C4AClH10" in fsyms && "Cl-" in fsyms && !("hemicarbonate" in fsyms)
        @test_throws ArgumentError phase_list_system(list, substances; remove = ["C4AClH10"])
        @test_throws ArgumentError phase_list_system(list, substances; add = ["no such phase"])
        # A database without a member of the list's solid solutions is refused,
        # after the warning of `build_solid_solutions` that names the phase.
        err = @test_logs (:warn, r"CSHQ") match_mode = :any try
            phase_list_system(list, [s for s in substances if symbol(s) != "KSiOH"])
        catch e
            e
        end
        @test err isa ArgumentError && occursin("CSHQ", sprint(showerror, err))
        # Another model of the gel in place of the list's, and aqueous species left
        # out besides the list's.
        cashdb = build_species(datapath("cemdata18-cashplus.json"); verbose = false)
        swapped = phase_list_system(list, cashdb; replace = Dict("CSHQ" => "CASH+NK"), exclude_aqueous = ["NaOH@", "KOH@"])
        @test Set(name.(swapped.solid_solutions)) == Set(["CASH+NK", "AFm_SO4_OH"])
        ssyms = String.(symbol.(swapped.species))
        @test "TCNh" in ssyms && !("KSiOH" in ssyms) && !("NaOH@" in ssyms) && !("KOH@" in ssyms)
        @test_throws ArgumentError phase_list_system(list, cashdb; replace = Dict("CNASH_ss" => "CASH+NK"))

        # A leaching step reads the system of the last answer. Once a step has
        # given the AFm its second instance, that system has two more species than
        # the paste's, and the renewal is built in it.
        big = with_instances(pcs, "AFm_SO4_OH" => 2; T = 293.15)
        st0 = ChemicalState(pcs; T = 293.15u"K", n = fill(1.0e-3, length(pcs.species)) .* u"mol")
        stbig = with_instances(st0, big)
        @test length(big.species) == length(pcs.species) + 2
        st1, bren = ChemistryLab._renewal(stbig, 50.0)
        @test length(st1.n) == length(big.species)
        # The budget is the solids of the last answer and 50 g of pure water.
        nexp = ustrip.(us"mol", stbig.n)
        iw = findfirst(s -> symbol(s) == "H2O@", big.species)
        nexp[big.idx_aqueous] .= 0.0
        nexp[iw] = 50.0 / ustrip(us"g/mol", big.species[iw][:M])
        @test bren ≈ Float64.(big.SM.A) * nexp rtol = 1.0e-14
    end
end
