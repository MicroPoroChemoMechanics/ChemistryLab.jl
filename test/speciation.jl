@testsection "Speciation" begin
    @testsection "union_atoms from dicts" begin
        d1 = Dict(:H => 2, :O => 1)
        d2 = Dict(:C => 1, :O => 2)
        atoms = union_atoms([d1, d2])

        @test :H in atoms
        @test :C in atoms
        @test :O in atoms
        @test length(atoms) == 3   # no duplicates

        # Order follows ATOMIC_ORDER (C before H)
        idx_C = findfirst(==(:C), atoms)
        idx_H = findfirst(==(:H), atoms)
        @test !isnothing(idx_C) && !isnothing(idx_H)
        @test idx_C < idx_H
    end

    @testsection "union_atoms from species list" begin
        sp_list = [
            Species("H2O"),
            Species("CaCO3"),
            Species("Na+"),
        ]
        atoms = union_atoms(sp_list)
        @test :H in atoms
        @test :O in atoms
        @test :Ca in atoms
        @test :C in atoms
        @test :Na in atoms
        @test length(atoms) == length(unique(atoms))   # no duplicates
    end

    @testsection "union_atoms from atom vectors" begin
        v1 = [:H, :O]
        v2 = [:Ca, :O, :C]
        atoms = union_atoms([v1, v2])
        @test :H in atoms && :Ca in atoms && :C in atoms && :O in atoms
        @test length(atoms) == 4
    end

    @testsection "same_components dispatch" begin
        # Regular Species → atoms_charge (exported)
        sp = [Species("CaCO3"), Species("H2O")]
        @test ChemistryLab.same_components(sp) === atoms_charge

        # CemSpecies → oxides_charge (exported)
        cem = [CemSpecies("C3S"), CemSpecies("CH")]
        @test ChemistryLab.same_components(cem) === oxides_charge
    end

    @testsection "item_order dispatch" begin
        sp = [Species("H2O")]
        @test ChemistryLab.item_order(sp) === ATOMIC_ORDER

        cem = [CemSpecies("C3S")]
        @test ChemistryLab.item_order(cem) === OXIDE_ORDER
    end

    @testsection "speciation by atom list" begin
        all_sp = [
            Species("H2O"; aggregate_state = AS_AQUEOUS),
            Species("NaCl"; aggregate_state = AS_CRYSTAL),
            Species("CO2"; aggregate_state = AS_GAS),
            Species("CaCO3"; aggregate_state = AS_CRYSTAL),
        ]

        # Filter to species containing only H and O atoms
        filtered = speciation(all_sp, [:H, :O])
        syms = symbol.(filtered)
        @test "H2O" in syms
        @test "NaCl" ∉ syms
        @test "CO2" ∉ syms   # CO2 has C which is not in [:H, :O]

        # Filter to species containing Ca, C, O
        filtered2 = speciation(all_sp, [:Ca, :C, :O])
        syms2 = symbol.(filtered2)
        @test "CaCO3" in syms2
        @test "H2O" ∉ syms2
    end

    @testsection "speciation by aggregate_state filter" begin
        all_sp = [
            Species("H2O"; aggregate_state = AS_AQUEOUS),
            Species("NaCl"; aggregate_state = AS_CRYSTAL),
            Species("CO2"; aggregate_state = AS_GAS),
        ]
        # Include only aqueous species (H2O has only H and O)
        filtered = speciation(all_sp, [:H, :O, :Na, :Cl, :C]; aggregate_state = [AS_AQUEOUS])
        syms = symbol.(filtered)
        @test "H2O" in syms
        @test "NaCl" ∉ syms
        @test "CO2" ∉ syms
    end

    @testsection "speciation by seed species list" begin
        all_sp = [
            Species("H2O"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT),
            Species("Na+"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE),
            Species("CaCO3"; aggregate_state = AS_CRYSTAL),
            Species("Mg+2"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE),
        ]
        # Seed with Na+ — should include species made of Na (and/or Zz)
        seed = [Species("Na+"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)]
        filtered = speciation(all_sp, seed)
        # Na+ must be in the filtered set (it's a seed)
        @test any(s -> symbol(s) == "Na+", filtered)
    end

    # The three traps of 0.38: a view, the seeds kept whatever the filters, and
    # duplicates when two results are concatenated.
    @testsection "speciation returns a vector, and says what it keeps" begin
        all_sp = [
            Species("H2O"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT),
            Species("H+"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE),
            Species("OH-"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE),
            Species("Ca+2"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE),
            Species("CaOH+"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE),
            Species("Ca(OH)2"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT),
            Species("CaCO3"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT),
        ]
        before = copy(all_sp)

        # A vector of its own: filtering it leaves the list as it was.
        sp = speciation(all_sp, ["Ca(OH)2", "H2O"]; aggregate_state = [AS_AQUEOUS])
        @test sp isa Vector
        filter!(s -> aggregate_state(s) == AS_AQUEOUS, sp)
        @test all_sp == before

        # The seeds are kept whatever the filters, a solid one included...
        sp = speciation(all_sp, ["Ca(OH)2", "H2O"]; aggregate_state = [AS_AQUEOUS])
        @test "Ca(OH)2" in symbol.(sp)
        @test "CaCO3" ∉ symbol.(sp)
        # ...unless asked to go through them like every other species.
        sp = speciation(all_sp, ["Ca(OH)2", "H2O"]; aggregate_state = [AS_AQUEOUS], keep_seeds = false)
        @test "Ca(OH)2" ∉ symbol.(sp)
        @test issubset(["H2O", "H+", "OH-", "Ca+2", "CaOH+"], symbol.(sp))

        # The caller's `include_species` is added to the seeds, not substituted.
        sp = speciation(all_sp, ["Ca(OH)2"]; aggregate_state = [AS_AQUEOUS], include_species = ["CaCO3"])
        @test issubset(["Ca(OH)2", "CaCO3"], symbol.(sp))

        # A seed that is not in the list is an error, not a smaller space.
        @test_throws "Portlandit" speciation(all_sp, ["Portlandit", "H2O"])

        # Two results join without duplicates through `union`; through `vcat`
        # they carry each shared species twice, which `ChemicalSystem` reports
        # and drops.
        a = speciation(all_sp, ["Ca+2", "H2O"]; aggregate_state = [AS_AQUEOUS])
        b = speciation(all_sp, ["Ca(OH)2", "H2O"]; aggregate_state = [AS_AQUEOUS])
        u = union(a, b)
        @test length(u) == length(unique(symbol.(u)))
        cs_u = ChemicalSystem(u)
        cs_v = @test_logs_own (:warn, r"listed more than once") ChemicalSystem(vcat(a, b))
        @test symbol.(cs_v.species) == symbol.(cs_u.species)
        @test size(cs_v.SM.A) == size(cs_u.SM.A)

        # Two different species under one symbol are refused.
        @test_throws "two different species" ChemicalSystem(
            [all_sp[1], Species("CaCO3"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE), all_sp[7]],
        )
    end
end
