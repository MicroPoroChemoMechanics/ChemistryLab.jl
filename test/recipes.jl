# The recipe and process layer: extents, materials, recipes and what they put
# into an equilibrium, and the processes built on them. Every bookkeeping claim is
# checked against an identity (the oxides a decomposition reproduces, the mass a
# budget conserves) or against the hand-built recipe of `scripts/gruyaert2010.jl`.

using ChemistryLab, DynamicQuantities, OrderedCollections, Test

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
        t, h = 7.0, 1.0e-3
        slope = (extent(pk, t + h) - extent(pk, t - h)) / (2h) / 86400
        α = extent(pk, t)
        @test slope ≈ rate(293.15, 1.0e5, 0.0, Dict("C3S" => 1 - α), nothing, Dict("C3S" => 1.0)) rtol = 1.0e-2
        # The w/c factor of Parrott and Killoh: no effect below 1.333 w/c, and a
        # hydration that stops at (1 + 4.444 w/c)/3.333.
        pk3 = ParrottKillohExtent("C3S"; w_c = 0.3)
        @test extent(pk3, 0.5) == extent(pk, 0.5)
        @test extent(pk3, 3650) < extent(pk, 3650)
        @test extent(pk3, 3650) <= (1 + 4.444 * 0.3) / 3.333
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
        @test 0 < bound_water(rs) < 0.4
        ps = pore_solution(rs)
        @test ps.pH > 12 && ps.elements[:K] > 0
        # The residue has no sourced enthalpy (the glass): the heat is not complete.
        @test isnan(enthalpy(rs))
        @test isnan(volume(rs).residual) && "S1 slag (Durdzinski 2017)" in volume(rs).missing
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

        le = leach(rs, 2)
        @test all(s -> s.certificate.optimal, le)
        @test pore_solution(le[2]).elements[:K] < pore_solution(rs).elements[:K]
    end
end
