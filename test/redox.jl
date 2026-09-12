using Logging

# Redox: the electron as a component, and the potentials read off a couple.
#
# Nothing in the package solved a multi-valence system before this file existed,
# so these tests are the first exercise of the branch of `StoichMatrix` that
# keeps the charge row when an element carries several oxidation states.

@testsection "redox" begin
    substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
    byname = Dict(symbol(s) => s for s in substances)

    redox_species = ["SO4-2", "HS-", "Fe+2", "Fe+3", "Ca+2", "O2@", "H2@"]
    sp = speciation(substances, redox_species; aggregate_state = [AS_AQUEOUS])
    cs = ChemicalSystem(sp, CEMDATA_PRIMARIES)

    st = ChemicalState(cs)
    set_quantity!(st, "H2O@", 1.0u"kg")
    for (s, n) in (
            "SO4-2" => 1.0e-3, "HS-" => 1.0e-3, "Fe+2" => 1.0e-6,
            "Fe+3" => 1.0e-6, "Ca+2" => 1.0e-3,
        )
        set_quantity!(st, s, n * u"mol")
    end
    model = DiluteSolutionModel()

    @testset "the electron carries the conventional standard state" begin
        # A convention, not a measurement -- and the whole of pe rests on it.
        @test ELECTRON.ΔₐG⁰ == 0.0
        @test ELECTRON.ΔₐH⁰ == 0.0
        @test ELECTRON.S⁰ == 0.0
        @test symbol(ELECTRON) == "e-"
        @test charge(ELECTRON) == -1
    end

    @testset "a half-reaction balances itself, electron included" begin
        r = half_reaction(st, "SO4-2", "HS-")
        # The textbook sulfate/sulfide half-reaction, obtained from the element
        # and charge balance alone -- no coefficient is transcribed here.
        @test occursin("8e⁻", r.equation)
        @test occursin("9H⁺", r.equation)
        # Published log K at 25 C is 33.66 (PHREEQC llnl/wateq4f); this one is
        # computed from CEMDATA18's own Gibbs energies, so agreement is a check
        # that the two datasets share a reference state.
        @test r.logK⁰(T = 298.15) ≈ 33.66 atol = 0.1

        rf = half_reaction(st, "Fe+3", "Fe+2")
        @test occursin("e⁻", rf.equation)
        # E0 = 0.771 V for Fe(III)/Fe(II), i.e. log K = 0.771 / 0.05916 = 13.03
        @test rf.logK⁰(T = 298.15) ≈ 13.03 atol = 0.1
    end

    @testset "pe and Eh" begin
        # With the two iron activities equal, every activity term cancels and
        # pe must fall exactly on log K. That is an identity, so it tests the
        # assembly of the sum rather than the data.
        pe_fe = pe(st, model; couple = "Fe+3" => "Fe+2")
        @test pe_fe ≈ half_reaction(st, "Fe+3", "Fe+2").logK⁰(T = 298.15) atol = 1.0e-6

        # Sulfate/sulfide at equal molality and this pH is strongly reducing.
        pe_s = pe(st, model)
        @test pe_s < 0

        # Nernst: 0.05916 V per pe unit at 25 C.
        @test Eh(st, model) ≈ 0.05916 * pe_s atol = 1.0e-3
        @test Eh(st, model; couple = "Fe+3" => "Fe+2") ≈ 0.05916 * pe_fe atol = 1.0e-3

        # The two couples disagree by a wide margin, which is the physical
        # point: a paste has no single redox state unless the couples are at
        # mutual equilibrium.
        @test abs(pe_fe - pe_s) > 10
    end

    @testset "refusals" begin
        @test_throws ArgumentError half_reaction(st, "SO4-2", "NotASpecies")
        @test_throws ArgumentError half_reaction(st, "NotASpecies", "HS-")
        # A species cannot be both members. Asked for it, the balance closes by
        # creating matter -- `∅ = Ca²⁺ + 2e⁻` -- and returns a plausible-looking
        # potential (-49.9) from no chemistry at all, which is worse than an
        # error.
        @test_throws ArgumentError pe(st, model; couple = "Ca+2" => "Ca+2")
        @test_throws ArgumentError half_reaction(st, "Ca+2", "Ca+2")
    end

    @testset "the charge row survives as a conservation component" begin
        # This is what makes a redox calculation possible at all: with sulfur at
        # -II and +VI in the same list, charge is no longer a combination of the
        # element rows, so `StoichMatrix` keeps it.
        comps = String.(symbol.(cs.SM.primaries))
        @test "Zz" in comps

        # And with one valence per element it does not: charge is then a fixed
        # combination of the element rows, and keeping it would add a redundant
        # one. The list has to be given explicitly to show this -- `speciation`
        # derives species from the ATOMS, so asking for `SO4-2` alone pulls in
        # HS-, SO3-2 and S2O3-2 with it, and the system has redox freedom
        # whether or not the caller wanted any.
        cs1 = ChemicalSystem(
            [
                Species("H2O"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT),
                Species("Na+"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE),
                Species("Cl-"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE),
            ],
        )
        @test !("Zz" in String.(symbol.(cs1.SM.primaries)))
    end
end

@testsection "redox as a prescribed condition" begin
    substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
    sp = speciation(
        substances, ["SO4-2", "HS-", "Ca+2", "O2@", "H2@"];
        aggregate_state = [AS_AQUEOUS],
    )
    cs = ChemicalSystem(sp, CEMDATA_PRIMARIES)
    st = ChemicalState(cs)
    set_quantity!(st, "H2O@", 1.0u"kg")
    set_quantity!(st, "SO4-2", 1.0e-3u"mol")
    set_quantity!(st, "HS-", 1.0e-3u"mol")
    set_quantity!(st, "Ca+2", 1.0e-3u"mol")
    model = DiluteSolutionModel()
    b = Float64.(cs.SM.A) * ustrip.(us"mol", st.n)

    @testset "the prescribed potential is the one obtained" begin
        for target in (-6.0, -3.0, 0.0)
            eq, cert = equilibrate_certified(
                st; model = model, b = b, constraint = FixedpE(target)
            )
            @test cert.optimal
            # Read back through `pe`, which knows nothing of the constraint:
            # agreement is a check on the whole path, not on one formula.
            @test pe(eq, model) ≈ target atol = 1.0e-3
        end
    end

    @testset "raising the potential oxidizes the sulfur" begin
        ratios = map((-6.0, -3.0, 0.0)) do target
            eq, _ = equilibrate_certified(
                st; model = model, b = b, constraint = FixedpE(target)
            )
            ustrip(us"mol", moles(eq, "SO4-2")) / ustrip(us"mol", moles(eq, "HS-"))
        end
        # S(VI)/S(-II) must increase with pe, and by orders of magnitude.
        @test issorted(ratios)
        @test ratios[end] / ratios[1] > 100
    end

    @testset "FixedEh is FixedpE through Nernst" begin
        # 0.05916 V per pe unit at 25 C.
        eq_e, cert_e = equilibrate_certified(
            st; model = model, b = b, constraint = FixedEh(-0.1775u"V")
        )
        @test cert_e.optimal
        @test pe(eq_e, model) ≈ -3.0 atol = 5.0e-3
    end

    @testset "refusals" begin
        # A couple at one oxidation state prescribes nothing.
        @test_throws ArgumentError equilibrate_certified(
            st; model = model, b = b,
            constraint = FixedpE(0.0; couple = "Ca+2" => "Ca+2"),
        )
        # The titrant has to be a species of the system.
        @test_throws ArgumentError equilibrate_certified(
            st; model = model, b = b, constraint = FixedpE(0.0; titrant = "NotHere"),
        )
    end
end

@testsection "a couple that does not buffer is said so" begin
    # A slag paste puts all of its sulfur into an AFm phase, leaving the aqueous
    # sulfide at the solver floor. The `pe` computed from such a couple is set by
    # `ϵ` and not by the chemistry, and it looks like an ordinary answer -- which
    # is exactly why it has to announce itself.
    substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
    sp = speciation(
        substances, ["SO4-2", "HS-", "Ca+2", "O2@"]; aggregate_state = [AS_AQUEOUS]
    )
    cs = ChemicalSystem(sp, CEMDATA_PRIMARIES)
    model = DiluteSolutionModel()

    both = ChemicalState(cs)
    set_quantity!(both, "H2O@", 1.0u"kg")
    set_quantity!(both, "SO4-2", 1.0e-3u"mol")
    set_quantity!(both, "HS-", 1.0e-3u"mol")
    set_quantity!(both, "Ca+2", 1.0e-3u"mol")
    @test_logs min_level = Logging.Warn pe(both, model)      # both present: silent

    starved = ChemicalState(cs)
    set_quantity!(starved, "H2O@", 1.0u"kg")
    set_quantity!(starved, "SO4-2", 1.0e-3u"mol")
    set_quantity!(starved, "Ca+2", 1.0e-3u"mol")             # no sulfide at all
    @test_logs (:warn,) match_mode = :any pe(starved, model)
end
