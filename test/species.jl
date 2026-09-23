@testsection "Species" begin
    # Test basic Species construction
    @test_nowarn Species("H2O")
    water = Species("H2O", name = "Water", symbol = "H2O", aggregate_state = AS_AQUEOUS)
    @test name(water) == "Water"
    @test symbol(water) == "H2O"
    @test aggregate_state(water) == AS_AQUEOUS

    # Test atoms and composition
    @test atoms(water) == Dict(:H => 2, :O => 1)
    @test charge(water) == 0

    # Test properties
    @test water.M ≈ 18.015u"g/mol"  # Molar mass of water

    # Test property manipulation
    water[:custom_prop] = 42
    @test water[:custom_prop] == 42
    @test haskey(water, :custom_prop)
    @test !haskey(water, :nonexistent)

    # Test species equality (based on formula + aggregate_state + class)
    water2 = Species("H2O", name = "Water2", aggregate_state = AS_AQUEOUS)
    @test water == water2

    # Test species inequality: same formula, different aggregate_state
    vapour = Species("H₂O"; name = "Vapour", symbol = "H₂O⤴", aggregate_state = AS_GAS, class = SC_GASFLUID)
    @test vapour != water

    # Test ionic species
    nacl = Species("Na+")
    @test charge(nacl) == 1
    @test atoms(nacl) == Dict(:Na => 1)

    # Test complex formula
    calcium_carbonate = Species("CaCO3", aggregate_state = AS_CRYSTAL)
    @test atoms(calcium_carbonate) == Dict(:Ca => 1, :C => 1, :O => 3)
    @test aggregate_state(calcium_carbonate) == AS_CRYSTAL

    # Test aqueous solvent class
    h2o_solvent = Species("H2O"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
    @test class(h2o_solvent) == SC_AQSOLVENT

    # Test aqueous solute class
    na_plus = Species("Na+"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    @test class(na_plus) == SC_AQSOLUTE

    # Test component class
    sio2 = Species("SiO2"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)
    @test class(sio2) == SC_COMPONENT

    # Test mendeleev_filter (not exported — use ChemistryLab.mendeleev_filter)
    valid_species = Species("CaCO3"; aggregate_state = AS_CRYSTAL)
    @test !isnothing(ChemistryLab.mendeleev_filter(valid_species))
    invalid_species = Species("Xx"; aggregate_state = AS_UNDEF)
    @test isnothing(ChemistryLab.mendeleev_filter(invalid_species))

    # Test apply: double all stoichiometric coefficients
    doubled = apply(x -> x * 2, water)
    @test atoms(doubled) == Dict(:H => 4, :O => 2)

    # Test aqueous dissolved species with @ suffix
    @test_nowarn Species("CO2@"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    co2_aq = Species("CO2@"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    @test aggregate_state(co2_aq) == AS_AQUEOUS

    # Test CemSpecies construction
    c3s = CemSpecies("C3S")
    @test !isempty(atoms(c3s))
    @test :Ca in keys(atoms(c3s))
    @test :Si in keys(atoms(c3s))

    # CemSpecies oxide composition
    ch = CemSpecies("CH")
    ch_oxides = oxides(ch)
    @test :C in keys(ch_oxides)   # CaO
    @test :H in keys(ch_oxides)   # H2O

    # CemSpecies components returns oxides
    @test oxides(c3s) == components(c3s)

    # CemSpecies hash consistency
    c3s_copy = CemSpecies("C3S")
    @test c3s == c3s_copy
    @test hash(c3s) == hash(c3s_copy)
end

@testsection "with_class" begin
    # Base case: SC_COMPONENT → SC_SSENDMEMBER
    s = Species("CaCO3"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)
    s2 = with_class(s, SC_SSENDMEMBER)

    @test class(s2) == SC_SSENDMEMBER
    @test class(s) == SC_COMPONENT          # original unchanged

    # All other fields preserved
    @test name(s2) == name(s)
    @test symbol(s2) == symbol(s)
    @test formula(s2) == formula(s)
    @test aggregate_state(s2) == aggregate_state(s)
    @test properties(s2) === properties(s)  # same dict object (shared by reference)

    # Round-trip: any class can be set
    s3 = with_class(s2, SC_AQSOLUTE)
    @test class(s3) == SC_AQSOLUTE
    @test class(s2) == SC_SSENDMEMBER       # s2 still unchanged

    # Works with a species carrying properties
    sp = Species("Ca+2"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    sp[:custom] = 42
    sp2 = with_class(sp, SC_UNDEF)
    @test class(sp2) == SC_UNDEF
    @test sp2[:custom] == 42

    # with_class result can be used in a SolidSolutionPhase
    em1 = with_class(
        Species("Em1"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT),
        SC_SSENDMEMBER,
    )
    em2 = with_class(
        Species("Em2"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT),
        SC_SSENDMEMBER,
    )
    @test_nowarn SolidSolutionPhase("SS", [em1, em2])

    # SolidSolutionPhase also accepts SC_COMPONENT directly (auto-requalifies)
    raw1 = Species("Raw1"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)
    raw2 = Species("Raw2"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)
    ss_auto = @test_nowarn SolidSolutionPhase("SSAuto", [raw1, raw2])
    @test class(end_members(ss_auto)[1]) == SC_SSENDMEMBER
end

@testsection "what makes two species the same species" begin
    # The contract: formula, aggregate state, class, and whatever the SYMBOL adds
    # to that. The symbol is in it because the first three do not separate
    # polymorphs, and a polymorph is a different substance.

    cal = Species("CaCO3"; symbol = "Cal", aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)
    arg = Species("CaCO3"; symbol = "Arg", aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)

    # SAME formula, state and class -- which is exactly why the symbol has to
    # count. In CEMDATA18 their standard Gibbs energies differ by 821 J/mol.
    @test formula(cal) == formula(arg)
    @test aggregate_state(cal) == aggregate_state(arg)
    @test class(cal) == class(arg)

    @test cal != arg
    @test hash(cal) != hash(arg)
    @test !haskey(Dict(cal => 1), arg)
    @test length(unique([cal, arg])) == 2
    @test length(Set([cal, arg])) == 2
    # A generic CaCO3 names no polymorph, so it is neither of them.
    @test Species("CaCO3"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT) != cal

    # THE OTHER DIRECTION, which is what `_identity_symbol` exists for: a symbol
    # that merely spells the species' own formula adds nothing, so two spellings
    # are one species -- with equal hashes, or `Dict` would disagree with `==`.
    #
    # `unicode(formula(s))` would NOT do as the canonical form: it is not
    # constant on the equality classes of `Formula`, which are composition and
    # charge. `Formula("e") == Formula("e-")` while their unicode spellings are
    # "e" and "e⁻", and `Formula("Ca+2") == Formula("Ca⁺²")` while theirs are
    # "Ca²⁺" and "Ca⁺²". Canonicalizing would have split those pairs; dropping
    # the symbol entirely does not, because the formula is compared anyway.
    for (a, b) in (("H2O", "H₂O"), ("e", "e-"), ("Ca+2", "Ca⁺²"), ("CO3-2", "CO₃²⁻"))
        s1, s2 = Species(a), Species(b)
        @test s1 == s2
        @test hash(s1) == hash(s2)
        @test Dict(s1 => 1)[s2] == 1
        @test length(unique([s1, s2])) == 1
        # and the stored symbol is untouched -- it is a lookup key, not an
        # identity, and `cs["H2O"]` depends on it being what the caller typed.
        @test symbol(s1) == a
        @test symbol(s2) == b
    end
    @test ChemistryLab._identity_symbol(Species("H2O")) === nothing
    @test ChemistryLab._identity_symbol(cal) == "Cal"

    # The lookup key, asserted here because canonicalizing the stored symbol --
    # the obvious way to do all of the above -- breaks it.
    cs = ChemicalSystem([Species("H2O"; aggregate_state = AS_AQUEOUS)])
    @test cs["H2O"] == Species("H2O"; aggregate_state = AS_AQUEOUS)

    # The invariant `isequal ⟹ hash`, over a whole shipped database rather than
    # over a handful of cases chosen to pass.
    subs = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
    violations = [
        (symbol(subs[i]), symbol(subs[j]))
            for i in eachindex(subs) for j in eachindex(subs)
            if i < j && isequal(subs[i], subs[j]) && hash(subs[i]) != hash(subs[j])
    ]
    @test isempty(violations)
end
