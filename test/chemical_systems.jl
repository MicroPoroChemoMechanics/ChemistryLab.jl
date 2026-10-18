@testsection "ChemicalSystem" begin
    # ── Shared fixtures ──────────────────────────────────────────────────────
    h2o = Species("H2O"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
    hplus = Species("H+"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    oh = Species("OH-"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    nacl = Species("NaCl"; aggregate_state = AS_CRYSTAL)
    co2g = Species("CO2"; aggregate_state = AS_GAS, class = SC_GASFLUID)
    sio2 = Species("SiO2"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)

    @testsection "Basic construction" begin
        cs = ChemicalSystem([h2o, hplus, oh])
        @test length(cs) == 3
        @test size(cs) == (3,)
    end

    @testsection "Indexing" begin
        cs = ChemicalSystem([h2o, hplus, oh])
        @test cs[1] == h2o
        @test cs["H2O"] == h2o
        @test cs["H+"] == hplus
        @test_throws KeyError cs["nonexistent"]
    end

    @testsection "Iteration" begin
        cs = ChemicalSystem([h2o, hplus, oh])
        collected = collect(cs)
        @test length(collected) == 3
        @test h2o in collected
    end

    @testsection "Views by aggregate state" begin
        cs = ChemicalSystem([h2o, hplus, nacl, co2g])
        @test length(aqueous(cs)) == 2   # h2o + hplus
        @test length(crystal(cs)) == 1   # nacl
        @test length(gas(cs)) == 1       # co2g
        @test all(s -> aggregate_state(s) == AS_AQUEOUS, aqueous(cs))
        @test all(s -> aggregate_state(s) == AS_CRYSTAL, crystal(cs))
        @test all(s -> aggregate_state(s) == AS_GAS, gas(cs))
    end

    @testsection "Views by class" begin
        cs = ChemicalSystem([h2o, hplus, oh, co2g, sio2])
        @test length(solutes(cs)) == 2      # hplus + oh
        @test class(solvent(cs)) == SC_AQSOLVENT
        @test length(gasfluid(cs)) == 1     # co2g
        @test length(components(cs)) == 1   # sio2
    end

    @testsection "Construction with explicit primaries" begin
        # WATER ALONE CANNOT DESCRIBE THIS SYSTEM, and until 0.18.0 nothing said
        # so. The three species span a two-dimensional space — over (H, O,
        # charge) the determinant of their columns is exactly zero, so the rank
        # is 2 — and a single component cannot express the other two. The
        # decomposition runs through `pinv`, which projects rather than refuses,
        # so the matrix came out reading `H+ = 0.4 H2O` and `OH- = 0.6 H2O`.
        # That is not an approximation: it conserves a quantity that is not
        # conserved, and it lets a solver make `H+` out of water with no `OH-`
        # to go with it, charge and all.
        sp = [h2o, hplus, oh]
        err = try
            ChemicalSystem(sp, [h2o])
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("cannot be written over", sprint(showerror, err))

        # Two components, and the same species list is expressible: this is the
        # smallest basis that describes an acid-base system, and it is what the
        # refusal points the caller at.
        cs = ChemicalSystem(sp, [h2o, hplus])
        @test length(cs.SM.primaries) == 2
        @test Set(String.(symbol.(cs.SM.primaries))) == Set(["H2O", "H+"])
    end

    @testsection "Construction from string primaries" begin
        sp = [h2o, hplus, oh]
        cs = ChemicalSystem(sp, ["H2O", "H+"])
        primaries_syms = symbol.(cs.SM.primaries)
        @test "H2O" in primaries_syms
        @test "H+" in primaries_syms
        @test "OH-" ∉ primaries_syms
    end

    @testsection "Pre-computed matrices" begin
        cs = ChemicalSystem([h2o, hplus, oh])
        @test !isnothing(cs.CSM)   # canonical matrix built
        @test !isnothing(cs.SM)    # primary-relative matrix built
        @test size(cs.CSM.A, 2) == 3   # 3 species → 3 columns
    end

    @testsection "Reactions via kinetic_species" begin
        # With no kinetic_species → no reactions
        cs_no_kin = ChemicalSystem([h2o, hplus, oh])
        @test isempty(cs_no_kin.reactions)
        @test isempty(cs_no_kin.idx_kinetic)

        # Given as a vector of pairs rather than a dictionary, over components
        # that include one which is not a species (dissolved CO2). The null
        # space is then computed by rational elimination, and calcite enters two
        # of the vectors it returns (with Ca²⁺ and HCO₃⁻, and with Ca²⁺ and
        # CO₃²⁻): isolating it in one column takes a column operation.
        calcite = Species("CaCO3"; symbol = "Calcite", aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)
        ca = Species("Ca+2"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
        hco3 = Species("HCO3-"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
        co3 = Species("CO3-2"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
        co2 = Species("CO2"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
        sp = [calcite, h2o, hplus, oh, ca, hco3, co3]
        comps = [h2o, hplus, ca, co2]
        calcite_rate(T, P, t, n, lna, n0) = 1.0e-6

        plain = StoichMatrix(sp, comps)
        @test co2 ∉ plain.species                     # a component that is not a species
        @test iszero(plain.A * plain.N)
        @test count(!iszero, plain.N[1, :]) == 2      # calcite in two vectors

        cs = ChemicalSystem(sp, comps; kinetic_species = ["Calcite" => calcite_rate])
        @test cs.idx_kinetic == [1]
        @test symbol.(kinetic_species(cs)) == ["Calcite"]
        A, N = cs.SM.A, cs.SM.N
        # Still a basis of the null space, with calcite in one vector only.
        @test size(N, 2) == length(sp) - size(A, 1) == 3
        @test iszero(A * N)
        @test count(!iszero, N[1, :]) == 1

        # That vector is the kinetic reaction. It conserves every element and
        # the charge, counted on the formulas of its species rather than on the
        # matrix it came from.
        rxn = only(cs.reactions)
        ν = rxn[calcite]
        @test !iszero(ν)
        balance = Dict{Symbol, Any}()
        for (s, c) in rxn, (el, k) in atoms(s)
            balance[el] = get(balance, el, 0) + c * k
        end
        @test all(iszero, values(balance))
        @test iszero(sum(c * charge(s) for (s, c) in rxn))
        # Its rate is the one given, per mole of calcite.
        @test rxn[:rate](298.15, 1.0e5, 0.0, nothing, nothing, nothing) ==
            calcite_rate(298.15, 1.0e5, 0.0, nothing, nothing, nothing) / abs(ν)
    end

    @testsection "merge two systems" begin
        cs1 = ChemicalSystem([h2o, hplus])
        cs2 = ChemicalSystem([oh, h2o])   # h2o duplicate — cs1 wins
        cs = merge(cs1, cs2)
        @test length(cs) == 3             # H2O, H+, OH- (no duplicate H2O)
        @test cs["H2O"] === cs1["H2O"]    # cs1 wins on conflict
    end

    @testsection "merge variadic" begin
        cs1 = ChemicalSystem([h2o])
        cs2 = ChemicalSystem([hplus])
        cs3 = ChemicalSystem([oh])
        cs = merge(cs1, cs2, cs3)
        @test length(cs) == 3
    end

    @testsection "merge keeps the solid solutions" begin
        # Until 0.37.0 the merged system was built from the species alone: the
        # solid solutions of both systems were dropped.
        em1 = Species("AFm1"; aggregate_state = AS_CRYSTAL, class = SC_SSENDMEMBER)
        em2 = Species("AFm2"; aggregate_state = AS_CRYSTAL, class = SC_SSENDMEMBER)
        cs1 = ChemicalSystem([em1, em2]; solid_solutions = [SolidSolutionPhase("AFm", [em1, em2])])
        cs = merge(ChemicalSystem([h2o, hplus, oh]), cs1)
        @test length(cs) == 5
        @test [name(ss) for ss in cs.solid_solutions] == ["AFm"]
        @test cs.ss_groups == [[4, 5]]
        # A phase declared twice is declared again, its copies built once.
        two = ChemicalSystem(
            [em1, em2];
            solid_solutions = [SolidSolutionPhase("AFm", [em1, em2]; model = RedlichKisterModel(a0 = 20_000.0), instances = 2)],
        )
        cs2 = merge(two, ChemicalSystem([h2o]))
        @test sort([name(ss) for ss in cs2.solid_solutions]) == ["AFm", "AFm#2"]
        @test count(s -> symbol(s) == "AFm1#2", cs2.species) == 1
        @test length(cs2) == length(two) + 1
    end
end
