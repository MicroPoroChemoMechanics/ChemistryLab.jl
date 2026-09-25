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
# GUO'S CONSTANTS ARE NOT CEMDATA18'S. Their §2 says so plainly -- "Cemdata2007
# gives Kp and DrGT0 for nearly all phases in cement hydrate" -- and the note
# under their dissolution table points at Lothenbach, Matschei, Moschner &
# Glasser (2008), which is Cemdata07. So this is a comparison ACROSS database
# versions, and the interesting question is where the two versions still agree.
#
# The phases themselves do map one to one, and the molar masses of their Table 3
# settle it: AFm 622.5, AFt 1255.1, Friedel's salt 561.3, CH 74.1 all reproduce
# from the CEMDATA18 formulas to better than 0.03 g/mol.

include("reference_species.jl")

@testsection "Chloride binding in a hydrated paste" begin

    subs = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
    byname = Dict(symbol(s) => s for s in subs)
    molar(s) = ustrip(us"g/mol", byname[s][:M])
    cshq = ["CSHQ-TobD", "CSHQ-TobH", "CSHQ-JenH", "CSHQ-JenD", "NaSiOH", "KSiOH"]

    # Guo's inventory and tables, from data/literature/Guo2018.json.
    guo_g(q) = ustrip(u"g", literature_value("Guo2018", q))       # g per liter of concrete
    guo_M(phase) = ustrip(u"g/mol", literature_row("Guo2018", "hydrates", phase).molar_mass)

    @testset "Guo's Table 3 is CEMDATA18, phase for phase" begin
        for (phase, printed) in (
                "monosulphate12" => "AFm", "ettringite" => "AFt",
                "C4AClH10" => "Friedel's salt", "Portlandite" => "CH",
            )
            @test molar(phase) ≈ guo_M(printed) atol = 0.03
        end
        # Their CaCO3 is printed as 100.9; calcite is 100.09. A typo, and it
        # touches nothing in this figure, which carries no carbonate.
        @test abs(molar("Cal") - guo_M("CaCO3")) > 0.8
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
    # The pore solution fills the porosity, at a density of 1 g/cm³.
    pore_g = 10 * literature_value("Guo2018", "porosity_percent")   # g per liter
    # Their C-S-H is 5(CaO)·3(SiO2)·6.3(H2O) (Table 2), and its molar mass (Table
    # 3) is that of one third of it: moles of silicon, then lime and water.
    csh = let t = literature_table("Guo2018", "csh_formula")
        Dict(zip(t.oxide, t.coefficient ./ t.coefficient[findfirst(==("SiO2"), t.oxide)]))
    end
    M_water = molar("H2O@")
    M_salt = ustrip(us"g/mol", Species("NaCl")[:M])
    function charged(cs, nacl_frac)
        st = ChemicalState(cs)
        n_csh = guo_g("csh_per_liter") / guo_M("CSH")
        set_quantity!(st, "Lim", (csh["CaO"] * n_csh) * u"mol")
        set_quantity!(st, "Amor-Sl", n_csh * u"mol")
        # AND ITS STRUCTURAL WATER. Guo's dissolution reaction is written for
        # (CaO)5(SiO2)3(H2O)6.3, M = 574.1; their Table 3's 191.4 is one third
        # of that, so 225 g is 1.1755 mol of (CaO)1.667(SiO2)(H2O)2.1 and
        # carries 2.469 mol -- 44.5 g -- of water in the SOLID. Leaving it out
        # makes the system 44.5 g short: the solver then draws that water out
        # of the pore solution to hydrate the C-S-H, and the pore volume comes
        # out at 84 mL against the 146 mL the paper's porosity states.
        set_quantity!(st, "Portlandite", (guo_g("ch_per_liter") / molar("Portlandite")) * u"mol")
        set_quantity!(st, "monosulphate12", (guo_g("afm_per_liter") / molar("monosulphate12")) * u"mol")
        set_quantity!(st, "ettringite", (guo_g("aft_per_liter") / molar("ettringite")) * u"mol")
        set_quantity!(
            st, "H2O@",
            (pore_g * (1 - nacl_frac) / M_water + csh["H2O"] * n_csh) * u"mol"
        )
        salt = pore_g * nacl_frac / M_salt
        if salt > 0
            set_quantity!(st, "Na+", salt * u"mol")
            set_quantity!(st, "Cl-", salt * u"mol")
        end
        return st
    end

    # Descending, for the reason recorded in `limestone_blending_reference.jl`:
    # a warm start survives only while the assemblage holds.
    fractions = [0.02, 0.015, 0.01, 0.0075, 0.005, 0.0025, 0.001, 0.0]
    function sweep(pure)
        cs = build(pure)
        idx = Dict(symbol(s) => i for (i, s) in enumerate(cs.species))
        budgets = [
            Float64.(cs.SM.A) * ustrip.(us"mol", charged(cs, f).n) for f in fractions
        ]
        states, certs = equilibrate_path(
            charged(cs, first(fractions)), budgets; model = HKFActivityModel(
                å = 0.0, Ḃ = gems_bdot(), Kₙ = 0.0
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

    # Guo's own phase list. Their dissolution table carries C-S-H, CH, AFm,
    # AFt and Friedel's salt, and nothing else. (It is their Table 2 by the
    # captions and their Table 1 by the body text, which disagree.)
    guo = sweep(split("Portlandite monosulphate12 ettringite C4AClH10 Lim Amor-Sl"))

    @testset "the figure, on the phase list that drew it" begin
        # The inventory is preserved where there is no chloride to disturb it:
        # 9 g of AFm and 22.5 g of AFt are still 9 g and 22.5 g.
        @test guo[0.0].AFm ≈ guo_g("afm_per_liter") / molar("monosulphate12") rtol = 1.0e-3
        @test guo[0.0].AFt ≈ guo_g("aft_per_liter") / molar("ettringite") rtol = 1.0e-3
        @test guo[0.0].FS == 0

        # Consumed by 1 %, which is where Guo places it and where the XRD of
        # Hirao et al. that they cite loses the AFm reflection.
        @test guo[0.01].AFm == 0
        @test 0 < guo[0.005].AFm < guo[0.0].AFm

        # And nothing at all has happened at 0.1 %, which is Guo's own claim:
        # below about 0.5 % "the concentration of chloride is so low that it
        # cannot form Friedel's salt". The inventory is still untouched there.
        @test guo[0.001].FS == 0
        @test guo[0.001].AFm ≈ guo[0.0].AFm rtol = 1.0e-4

        # The plateau. Guo read 0.023 and 0.010 mol/L off Fig. 1(b).
        @test guo[0.02].AFt ≈ ustrip(u"mol", literature_value("Guo2018", "aft_plateau")) rtol = 0.02
        @test guo[0.02].FS ≈ ustrip(u"mol", literature_value("Guo2018", "friedel_plateau")) rtol = 0.05
        # Those two are THEIR numbers, read off a figure by eye, so the
        # agreement is a tolerance against a reading. The values this package
        # computes are printed on the page beside them, and are pinned here.
        @info "chloride: the plateau this package computes" AFt = guo[0.02].AFt FS = guo[0.02].FS
        @test guo[0.02].AFt ≈ 0.02274 atol = 5.0e-5
        @test guo[0.02].FS ≈ 0.00965 atol = 5.0e-6
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
        # 14-hydrate -- not the 12 Guo wrote. That is a difference between the
        # two database versions and not a slip of theirs: Guo's log K of
        # -29.2628 is Cemdata07's value for the 12-hydrate, and Cemdata18
        # RECALCULATED it to -29.23 (its Table 2 marks that entry ***,
        # "recalculated in this paper from DfG values") while carrying -29.26
        # for the 14-hydrate. Their number therefore lands within 0.003 of
        # CEMDATA18's 14-hydrate by arithmetic coincidence, and the amount is
        # unaffected either way.
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

        # AND THE TABLE, at the three chloride loadings this sweep computes. The
        # page prints six columns; three of them are rungs no solve here visits,
        # and they are marked there rather than presented as results. What is
        # computed is pinned at half the last printed digit — a bound of the
        # form `Kuzel > 0.01` leaves `0.0107` free to become `0.0140` with
        # nothing going red.
        @info "chloride: the computed columns of the Kuzel table" columns = [
            (f, full[f].AFm14, full[f].Kuzel, full[f].FS) for f in (0.001, 0.005, 0.01)
        ]
        for (frac, ms, kuzel, fs) in (
                (0.001, 0.0131, 0.0011, 0.0),
                (0.0025, 0.0085, 0.0047, 0.0),
                (0.0075, 0.0, 0.0116, 0.0),
                (0.015, 0.0, 0.0, 0.0096),
                (0.005, 0.001, 0.0107, 0.0),
                (0.01, 0.0, 0.0058, 0.0048),
            )
            @test full[frac].AFm14 ≈ ms atol = 5.0e-5
            @test full[frac].Kuzel ≈ kuzel atol = 5.0e-5
            @test full[frac].FS ≈ fs atol = 5.0e-5
        end
        # The 14-hydrate is the stable one at these conditions, not the 12 Guo
        # wrote: asserted rather than left to the prose above.
        @test full[0.001].AFm12 == 0
    end

    # NOT CHECKED, and worth saying why.
    #
    #  - The pore volume, which comes out at 129 mL against the 146 mL their
    #    14.6 % porosity states. The 17 mL is the difference between Guo's
    #    C-S-H, which carries 2.1 H2O per (CaO)1.667(SiO2), and CEMDATA18's
    #    CSHQ end members, which carry more and take the rest out of the pore.
    #    Nothing can close that gap without replacing the C-S-H model.
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
