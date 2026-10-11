# An oxide analysis as an element budget — the entry route for a glass.

@testsection "oxide budget" begin
    substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
    sp = speciation(
        substances, ["Portlandite", "Cal", "Gp", "Amor-Sl", "Brc", "AlOHmic"];
        aggregate_state = [AS_AQUEOUS],
    )
    cs = ChemicalSystem(sp, CEMDATA_PRIMARIES)
    prim = cs.SM.primaries

    # A ground granulated blast-furnace slag, in the shape a datasheet reports.
    slag = Dict("CaO" => 0.41, "SiO2" => 0.36, "Al2O3" => 0.11, "MgO" => 0.08)

    @testset "every element arrives, and only what was weighed in" begin
        b = oxide_budget(slag, prim; mass = 60.0u"g")
        # Calcium enters only through CaO, so the budget must equal the moles of
        # CaO weighed in -- to machine precision, since it is the same division.
        for (oxide, frac, comp) in (("CaO", 0.41, "Ca+2"), ("MgO", 0.08, "Mg+2"))
            i = findfirst(p -> symbol(p) == comp, prim)
            n = 60.0 * frac / ustrip(us"g/mol", Species(oxide)[:M])
            @test b[i] ≈ n rtol = 1.0e-12
        end
        # Silicon likewise, through SiO2 alone.
        i = findfirst(p -> symbol(p) == "SiO2@", prim)
        @test b[i] ≈ 60.0 * 0.36 / ustrip(us"g/mol", Species("SiO2")[:M]) rtol = 1.0e-12
    end

    @testset "the analysis is not renormalized" begin
        # It sums to 0.96: the missing 4 % is loss on ignition and minor oxides
        # the sheet does not report. Scaling it to 1 would invent material.
        @test sum(values(slag)) ≈ 0.96 atol = 1.0e-12
        b1 = oxide_budget(slag, prim; mass = 100.0u"g")
        b2 = oxide_budget(slag, prim; mass = 200.0u"g")
        @test b2 ≈ 2 .* b1 rtol = 1.0e-12          # linear in the mass, as it must be
        # and NOT equal to what a renormalized analysis would give
        scaled = Dict(k => v / 0.96 for (k, v) in slag)
        @test !isapprox(oxide_budget(scaled, prim; mass = 100.0u"g"), b1; rtol = 1.0e-3)
    end

    @testset "additivity" begin
        # Two oxides together contribute what each contributes alone: the budget
        # is a sum over the analysis, and nothing couples the entries.
        a = Dict("CaO" => 0.41)
        c = Dict("MgO" => 0.08)
        @test oxide_budget(merge(a, c), prim; mass = 50.0u"g") ≈
            oxide_budget(a, prim; mass = 50.0u"g") .+
            oxide_budget(c, prim; mass = 50.0u"g") rtol = 1.0e-12
        # A zero fraction contributes nothing rather than failing.
        @test oxide_budget(Dict("CaO" => 0.41, "TiO2" => 0.0), prim; mass = 50.0u"g") ≈
            oxide_budget(a, prim; mass = 50.0u"g") rtol = 1.0e-12
    end

    @testset "refusals" begin
        # An oxide the primaries cannot express has no decomposition, and a
        # least-squares approximation of one would put elements into the budget
        # that the oxide does not contain.
        @test_throws ArgumentError oxide_budget(
            Dict("TiO2" => 0.01), prim; mass = 10.0u"g"
        )
        @test_throws ArgumentError oxide_budget(
            Dict("CaO" => -0.1), prim; mass = 10.0u"g"
        )
    end

    @testset "primary_decomposition itself" begin
        # Portlandite over the primaries: one calcium, and the protons that the
        # oxide convention implies.
        x = primary_decomposition(Species("CaO"), prim)
        i = findfirst(p -> symbol(p) == "Ca+2", prim)
        @test x[i] ≈ 1.0 rtol = 1.0e-10
        @test_throws ArgumentError primary_decomposition(Species("TiO2"), prim)
    end

    @testset "an element below the valence of its primary goes through the unit charge" begin
        # Carbon (II) over the primaries of carbon (IV): CO is
        # CO3-2 + 4 H+ − 2 H2O − 2 Zz, its column of the conservation matrix.
        withz = [Species(s) for s in ("H2O@", "H+", "CO3-2", "Ca+2", "FeO2-", "Zz")]
        nozz = withz[1:(end - 1)]
        coef(x, s) = x[findfirst(p -> symbol(p) == s, withz)]
        x = primary_decomposition(Species("CO"), withz)
        @test coef(x, "CO3-2") ≈ 1 && coef(x, "H+") ≈ 4 && coef(x, "H2O@") ≈ -2
        @test coef(x, "Zz") ≈ -2
        # Written over the primaries, it carries the species' own charge.
        @test sum(x .* charge.(withz)) ≈ 0 atol = 1.0e-12
        @test coef(primary_decomposition(Species("FeO"), withz), "Zz") ≈ -1
        # An element at the valence of its primary leaves Zz untouched, exactly.
        @test coef(primary_decomposition(Species("CaO"), withz), "Zz") === 0.0
        @test coef(primary_decomposition(Species("Fe2O3"), withz), "Zz") === 0.0
        # Without Zz the system cannot hold carbon (II): refused, rather than
        # given the charge of carbon (IV) in silence.
        @test_throws ArgumentError primary_decomposition(Species("CO"), nozz)
        # A recipe does not drop such an oxide as if its element were absent:
        # the refusal reaches the caller.
        @test ChemistryLab._representable("FeO", nozz)
        @test !ChemistryLab._representable("TiO2", nozz)
    end

    @testset "carbon (II) of an oxide budget becomes formate" begin
        # 0.01 mol of CO and 0.005 mol of Na2O in a kilogram of water, over the
        # formate and the carbonate of SUPCRT: the carbon of valence two stays
        # formate, nothing in the system being able to oxidize it.
        org = Dict(symbol(s) => s for s in build_species(datapath("slop98-organic-thermofun.json"); verbose = false))
        inorg = Dict(symbol(s) => s for s in build_species(datapath("slop98-inorganic-thermofun.json"); verbose = false))
        sp = vcat([inorg[s] for s in ("H2O@", "H+", "OH-", "Na+", "CO3-2", "HCO3-", "CO2@")], [org["For-"], org["ForH@"]])
        cs = ChemicalSystem(sp)
        prim = cs.SM.primaries
        M(f) = ustrip(us"g/mol", Species(f)[:M])
        m_CO, m_Na2O = 0.01 * M("CO"), 0.005 * M("Na2O")
        m = m_CO + m_Na2O
        st = ChemicalState(cs)
        set_quantity!(st, "H2O@", 1.0u"kg")
        b = Float64.(cs.SM.A) * ustrip.(us"mol", st.n) .+
            oxide_budget(Dict("CO" => m_CO / m, "Na2O" => m_Na2O / m), prim; mass = m * u"g")
        eq, cert = equilibrate_certified(st; b)
        @test cert.optimal
        n = Dict(symbol(s) => ustrip(us"mol", x) for (s, x) in zip(cs.species, eq.n))
        @test n["For-"] + n["ForH@"] ≈ 0.01 rtol = 1.0e-9
        @test sum(n[symbol(s)] * charge(s) for s in cs.species) ≈ 0 atol = 1.0e-12
    end
end

@testsection "glass_species — a material with no formula unit" begin
    slag = Dict(
        "CaO" => 0.41, "SiO2" => 0.36, "Al2O3" => 0.11,
        "MgO" => 0.08, "SO3" => 0.02,
    )
    sp = glass_species(slag; symbol = "GGBS", M = 95.0u"g/mol")
    a = atoms(sp)

    # Every element the analysis reports is carried, MAGNESIUM INCLUDED, which
    # is the whole reason this exists: a slag written as a representative
    # mineral -- anorthite, `CaAl2Si2O8` -- has no magnesium, so no hydrotalcite
    # can form from it, and hydrotalcite is the one phase a slag is certain to
    # make.
    for el in (:Ca, :Si, :Al, :Mg, :S, :O)
        @test haskey(a, el) && a[el] > 0
    end
    @test symbol(sp) == "GGBS"
    @test ustrip(us"g/mol", sp[:M]) ≈ 95.0
    # The analysis is not renormalized; what it does not report is recorded.
    @test sp[:modeled_mass_fraction] ≈ 0.98

    # The element RATIOS are the analysis's, whatever formula-unit size is asked
    # for: `M` scales the unit, it does not change the material.
    sp2 = glass_species(slag; symbol = "GGBS2", M = 190.0u"g/mol")
    a2 = atoms(sp2)
    @test a2[:Ca] / a[:Ca] ≈ 2.0
    @test a2[:Mg] / a[:Mg] ≈ 2.0

    # And the ratio is the one the analysis implies, computed from molar masses
    # the database supplies rather than from a table written here.
    @test a[:Ca] / a[:Si] ≈
        (0.41 / ustrip(us"g/mol", Species("CaO")[:M])) /
        (0.36 / ustrip(us"g/mol", Species("SiO2")[:M])) rtol = 1.0e-10

    @test_throws ArgumentError glass_species(Dict("CaO" => -0.1); symbol = "X")
    @test_throws ArgumentError glass_species(Dict("CaO" => 0.0); symbol = "X")
    @test_throws ArgumentError glass_species(slag; symbol = "X", M = -1.0u"g/mol")

    # A composition being calibrated: the element counts carry the derivative of
    # the mass fraction, one mole of CaO per M(CaO) grams. They were stored as
    # Float64 and a dual fraction raised.
    M_CaO = ustrip(us"g/mol", Species("CaO")[:M])
    dCa = ForwardDiff.derivative(
        x -> atoms(glass_species(Dict("CaO" => x, "SiO2" => 0.36); symbol = "G", M = 95.0u"g/mol"))[:Ca],
        0.41,
    )
    @test dCa ≈ 95.0 / M_CaO rtol = 1.0e-12
end
