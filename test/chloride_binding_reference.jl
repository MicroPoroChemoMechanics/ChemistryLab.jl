# Chloride binding, against a paper that published its own inventory.
#
# Guo, Hong, Qiao, Ou & Li, Journal of Colloid and Interface Science 531 (2018)
# 56-63, https://doi.org/10.1016/j.jcis.2018.07.005, Fig. 1(b).
#
# Most cement papers report a paste and leave the modeler to reconstruct what
# was in it. This one states the inventory outright — per liter of concrete,
# C-S-H 225 g, CH 90 g, AFm 9 g, AFt 22.5 g, porosity 14.6 % — and then sweeps
# the chloride. That makes it reproducible without a single assumed number, and
# what it reports is a stoichiometric statement worth checking:
#
#     chloride replaces the sulfate of the AFm; the sulfate released converts
#     more AFm into ettringite; by about 1 % NaCl no AFm is left.
#
# Guo's thermodynamic data are CEMDATA18's, cited as such, so the phases can be
# matched one to one — and the molar masses of their Table 3 confirm it:
# AFm 622.5, AFt 1255.1, Friedel's salt 561.3, CH 74.1 all reproduce from the
# CEMDATA18 formulas to better than 0.03 g/mol.

@testsection "Chloride binding in a hydrated paste" begin

    subs = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
    byname = Dict(symbol(s) => s for s in subs)
    molar(s) = ustrip(us"g/mol", byname[s][:M])
    cshq = ["CSHQ-TobD", "CSHQ-TobH", "CSHQ-JenH", "CSHQ-JenD", "NaSiOH", "KSiOH"]

    @testset "Guo's Table 3 is CEMDATA18, phase for phase" begin
        for (phase, published) in (
                "monosulphate12" => 622.5, "ettringite" => 1255.1,
                "C4AClH10" => 561.3, "Portlandite" => 74.1,
            )
            @test molar(phase) ≈ published atol = 0.03
        end
        # Their CaCO3 is printed as 100.9; calcite is 100.09. A typo, and it
        # touches nothing in this figure, which carries no carbonate.
        @test molar("Cal") ≈ 100.09 atol = 0.01
    end

    function build(pure)
        sp = speciation(
            subs, vcat(pure, cshq);
            aggregate_state = [AS_AQUEOUS], exclude_species = split("H2@ O2@ CH4@")
        )
        return ChemicalSystem(
            sp, CEMDATA_PRIMARIES;
            solid_solutions = [SolidSolutionPhase("CSHQ", [byname[m] for m in cshq])],
        )
    end

    # 146 g of pore solution per liter of concrete, of which `nacl_frac` is salt.
    # The C-S-H enters as its bulk oxides at Guo's Ca/Si = 1.67; only the
    # element vector reaches the solver, so lime and silica here are
    # bookkeeping.
    pore_g = 146.0
    function charged(cs, nacl_frac)
        st = ChemicalState(cs)
        n_csh = 225 / 191.4                       # Guo's C-S-H, M = 191.4
        set_quantity!(st, "Lim", (1.67 * n_csh) * u"mol")
        set_quantity!(st, "Amor-Sl", n_csh * u"mol")
        set_quantity!(st, "Portlandite", (90 / molar("Portlandite")) * u"mol")
        set_quantity!(st, "monosulphate12", (9 / molar("monosulphate12")) * u"mol")
        set_quantity!(st, "ettringite", (22.5 / molar("ettringite")) * u"mol")
        set_quantity!(st, "H2O@", (pore_g * (1 - nacl_frac) / 18.015) * u"mol")
        salt = pore_g * nacl_frac / 58.44
        if salt > 0
            set_quantity!(st, "Na+", salt * u"mol")
            set_quantity!(st, "Cl-", salt * u"mol")
        end
        return st
    end

    # Descending, for the reason recorded in `limestone_blending_reference.jl`:
    # a warm start survives only while the assemblage holds.
    fractions = [0.02, 0.01, 0.005, 0.001, 0.0]
    function sweep(pure)
        cs = build(pure)
        idx = Dict(symbol(s) => i for (i, s) in enumerate(cs.species))
        budgets = [
            Float64.(cs.SM.A) * ustrip.(us"mol", charged(cs, f).n) for f in fractions
        ]
        states, certs = equilibrate_path(
            charged(cs, first(fractions)), budgets; model = HKFActivityModel(
                å = 0.0, Ḃ = 0.097637, Kₙ = 0.0
            )
        )
        @test all(c.optimal for c in certs)
        rows = map(states) do eq
            n = ustrip.(us"mol", eq.n)
            a(s) = haskey(idx, s) ? n[idx[s]] : 0.0
            (
                AFm = a("monosulphate12") + a("monosulphate14"),
                AFm12 = a("monosulphate12"), AFm14 = a("monosulphate14"),
                AFt = a("ettringite") + a("ettringite30"),
                FS = a("C4AClH10"), Kuzel = a("C4AsClH12"),
            )
        end
        return Dict(fractions[i] => rows[i] for i in eachindex(fractions))
    end

    # Guo's own phase list — their Tables 1 and 2 carry C-S-H, CH, AFm, AFt and
    # Friedel's salt, and nothing else.
    guo = sweep(split("Portlandite monosulphate12 ettringite C4AClH10 Lim Amor-Sl"))

    @testset "the figure, on the phase list that drew it" begin
        # The inventory is preserved where there is no chloride to disturb it:
        # 9 g of AFm and 22.5 g of AFt are still 9 g and 22.5 g.
        @test guo[0.0].AFm ≈ 9 / molar("monosulphate12") rtol = 1.0e-3
        @test guo[0.0].AFt ≈ 22.5 / molar("ettringite") rtol = 1.0e-3
        @test guo[0.0].FS == 0

        # Consumed by 1 %, which is where Guo places it and where the XRD of
        # Hirao et al. that they cite loses the AFm reflection.
        @test guo[0.01].AFm == 0
        @test 0 < guo[0.005].AFm < guo[0.001].AFm < guo[0.0].AFm

        # The plateau. Guo read 0.023 and 0.010 mol/L off Fig. 1(b).
        @test guo[0.02].AFt ≈ 0.023 rtol = 0.02
        @test guo[0.02].FS ≈ 0.01 rtol = 0.05
        @test guo[0.02].AFt ≈ guo[0.01].AFt rtol = 1.0e-3      # flat once AFm is gone

        # And the two conservation statements behind those numbers, which is
        # what makes the agreement more than a coincidence of two curves.
        # Sulfate: the AFm gives up one, ettringite takes three.
        consumed = guo[0.0].AFm - guo[0.02].AFm
        @test guo[0.02].AFt - guo[0.0].AFt ≈ consumed / 3 rtol = 1.0e-2
        # Aluminum: what the AFm held is split between Friedel's salt and the
        # ettringite that grew, both of which carry two aluminums per formula.
        @test guo[0.02].FS + (guo[0.02].AFt - guo[0.0].AFt) ≈ consumed rtol = 1.0e-2
    end

    full = sweep(
        split(
            "Portlandite monosulphate12 monosulphate14 ettringite ettringite30 " *
                "C4AClH10 C4AsClH12 C4AH13 C4AH19 C3AH6 straetlingite AlOHmic " *
                "Gp Anh Amor-Sl Lim syngenite Na2SO4"
        )
    )

    @testset "what a fuller phase list adds, and where it stops mattering" begin
        # Two things Guo's five phases cannot express.
        #
        # FIRST, the hydration state. CEMDATA18 carries monosulfate at 9, 10.5,
        # 12, 14 and 16 waters, and at these conditions the stable one is the
        # 14-hydrate — not the 12 Guo used. Their own constant says as much
        # without their noticing: Guo's Table 2 gives log K = -29.2628 for a
        # phase they write as Ca4Al2(SO4)(OH)12·6H2O, M = 622.5, which is the
        # 12-hydrate. Cemdata18's value for the 12-hydrate is -29.23; -29.26 is
        # the 14-hydrate. They paired one hydrate's constant with another's
        # formula, and it made no difference because the two are 0.03 log units
        # apart.
        @test full[0.0].AFm14 > 0
        @test full[0.0].AFm12 == 0
        @test full[0.0].AFm ≈ guo[0.0].AFm rtol = 1.0e-3     # same amount either way

        # SECOND, Kuzel's salt — half a chloride and half a sulfate per AFm
        # layer — which is the phase the transition actually goes through. It
        # holds the whole low-chloride range and peaks where Guo has Friedel's
        # salt just starting.
        @test full[0.005].Kuzel > 0.01
        @test full[0.005].FS == 0
        @test full[0.001].Kuzel > 0

        # And where it stops mattering: once the chloride is high enough to take
        # the last sulfate out of the AFm layer, Kuzel's salt is gone and the
        # two phase lists agree to three digits. Guo's plateau is right even
        # though the path to it is not theirs.
        @test full[0.02].Kuzel == 0
        @test full[0.02].FS ≈ guo[0.02].FS rtol = 1.0e-3
        @test full[0.02].AFt ≈ guo[0.02].AFt rtol = 1.0e-3
    end

    # NOT CHECKED, and worth saying why.
    #
    #  - pH. Guo reports 13.213 falling to 13.128 across the sweep; this system
    #    sits at 12.5, which is portlandite in alkali-free water. The difference
    #    is their pore solution's alkalis, whose concentrations the paper takes
    #    from a reference it does not reproduce. Nothing here can supply them.
    #
    #  - The 5 % end of their abscissa. At that salt loading the ionic strength
    #    reaches 1.4 mol/L, past the range the B-dot term in `HKFActivityModel`
    #    was fitted for; the solve stops certifying and returns pH 15.3. The
    #    sweep above stops at 2 %, where I is still about 0.4. A Pitzer model
    #    would be the instrument for the rest.
    #
    #  - Everything in their §2 that is surface complexation. `SC_SURFCOMPLEX`
    #    and the site families arrived in v0.20.0 and that half of the paper —
    #    the five `≡SiOH` reactions, the diffuse layer, the Ca > Cl > Na > K
    #    ordering — is now expressible and is not attempted here.
end
