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
end
