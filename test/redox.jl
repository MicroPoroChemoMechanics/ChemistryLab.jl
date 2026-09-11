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
        # Ca(II) has one valence here, so the "couple" balances with no electron.
        @test_throws ArgumentError pe(st, model; couple = "Ca+2" => "Ca+2")
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
