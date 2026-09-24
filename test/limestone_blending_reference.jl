# Blending a Portland cement with limestone, against two published figures.
#
# Kulik, Winnefeld, Kulik, Miron & Lothenbach, CemGEMS, RILEM Technical Letters
# 6 (2021) 36-52, https://doi.org/10.21809/rilemtechlett.2021.140, Fig. 7A; and
# [Lothenbach2019] Figs. 13-14.
#
# CemGEMS is the closest oracle there is for this package: same CEMDATA18, same
# Gibbs-energy minimization, and the paper states that it and GEM-Selektor agree
# exactly. What it does NOT give is the cement it ran on — Fig. 7A cites a
# composition published elsewhere — so what is checked here is the part that
# does not depend on it: the ORDER the phases appear in, the point at which
# added limestone stops reacting, and where the iron goes.
#
# The sequence both papers report is
#
#     monosulfate  →  hemicarbonate  →  monocarbonate + ettringite  →  free calcite
#
# and reproducing it turned out to depend on something neither paper states in
# those terms. See the second testset.

@testsection "Limestone blending of a CEM I" begin

    subs = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
    byname = Dict(symbol(s) => s for s in subs)

    clinker = CemSpecies.(split("C3S C2S C3A C4AF"))
    bogue = inv(mass_matrix(CanonicalStoichMatrix(clinker)).A)
    oxide_mass(ox) = ustrip(us"g/mol", Species(ox)[:M])

    pure = split(
        "C3S C2S C3A C4AF AlOHmic C3AH6 C3FH6 C4FH13 CAH10 C2AH7.5 " *
            "straetlingite monocarbonate hemicarbonate Femonocarbonate " *
            "Fe-hemicarbonate Fe-monosulphate C3FS1.34H3.32 Cal Portlandite " *
            "Anh Gp hydrotalcite Brc Mgs K2SO4 syngenite Na2SO4 K2O Na2O " *
            "Amor-Sl FeOOHmic Lim C3AFS0.84H4.32 C3FS0.84H4.32"
    )
    solutions = [
        "CSHQ" => ["CSHQ-JenD", "CSHQ-JenH", "CSHQ-TobD", "CSHQ-TobH", "KSiOH", "NaSiOH"],
        "C3(AF)S0.84H" => ["C3AFS0.84H4.32", "C3FS0.84H4.32"],
        "AFt_SO4" => ["ettringite", "ettringite30"],
        "AFm_SO4_OH" => ["C4AH13", "monosulphate12"],
        "hydrotalc" => ["Mg3AlC0.5OH", "Mg3FeC0.5OH"],
    ]
    species = speciation(
        subs, vcat(pure, reduce(vcat, last.(solutions)), ["CO2@"]);
        aggregate_state = [AS_AQUEOUS]
    )
    cs = ChemicalSystem(
        species, CEMDATA_PRIMARIES;
        solid_solutions = [
            SolidSolutionPhase(n, [byname[m] for m in mem]) for (n, mem) in solutions
        ],
    )
    idx = Dict(symbol(s) => i for (i, s) in enumerate(cs.species))
    model = HKFActivityModel(å = 0.0, Ḃ = 0.097637, Kₙ = 0.0)

    # `limestone_g` grams of calcite in exchange for clinker, the binder held at
    # 100 g, w/b = 0.5, complete hydration — the convention of Fig. 7A, whose
    # abscissa is "the addition of dry SCMs in grams in exchange to dry PC so
    # that the mass of the binder remains constant at 100 g".
    function charged(limestone_g; Fe2O3 = 4.49)
        oxides = Dict(
            "CaO" => 65.03 + (4.49 - Fe2O3), "SiO2" => 21.4, "Al2O3" => 3.84,
            "Fe2O3" => Fe2O3, "MgO" => 1.0, "K2O" => 0.46, "Na2O" => 0.13,
            "SO3" => 2.3,
        )
        f = (100 - limestone_g) / sum(values(oxides))
        n_ox = Dict(k => v * f / oxide_mass(k) for (k, v) in oxides)
        # Lime combined with the sulfate is not available to the silicates.
        cao_free = oxides["CaO"] - oxides["SO3"] * oxide_mass("CaO") / oxide_mass("SO3")
        g = bogue * [cao_free, oxides["SiO2"], oxides["Al2O3"], oxides["Fe2O3"]]

        st = ChemicalState(cs)
        for (i, ph) in enumerate(clinker)
            n = g[i] * f / ustrip(us"g/mol", byname[symbol(ph)][:M])
            n > 0 && set_quantity!(st, symbol(ph), n * u"mol")
        end
        set_quantity!(st, "Gp", n_ox["SO3"] * u"mol")
        set_quantity!(st, "Brc", n_ox["MgO"] * u"mol")
        set_quantity!(st, "K2O", n_ox["K2O"] * u"mol")
        set_quantity!(st, "Na2O", n_ox["Na2O"] * u"mol")
        limestone_g > 0 &&
            set_quantity!(st, "Cal", (limestone_g / oxide_mass("CaCO3")) * u"mol")
        free_water = 0.5 * 100 / oxide_mass("H2O") - 2n_ox["SO3"] - n_ox["MgO"]
        set_quantity!(st, "H2O@", free_water * u"mol")
        return st
    end

    function read_off(eq)
        n = ustrip.(us"mol", eq.n)
        amount(s) = n[idx[s]]
        fe(s) = n[idx[s]] * Float64(get(atoms_charge(byname[s]), :Fe, 0))
        return (
            Ms = amount("monosulphate12"), Hc = amount("hemicarbonate"),
            Mc = amount("monocarbonate"),
            AFt = amount("ettringite") + amount("ettringite30"),
            calcite = amount("Cal"),
            fe_in_hg = (fe("C3AFS0.84H4.32") + fe("C3FS0.84H4.32")) /
                sum(
                n[i] * Float64(get(atoms_charge(cs.species[i]), :Fe, 0))
                    for i in eachindex(n)
            ),
        )
    end

    # ONE path, walked DOWNHILL, and the order is not cosmetic.
    #
    # A cold solve of this system costs about 26 s; a warm one, started from a
    # neighboring answer, 0.2-0.7 s. But a warm start only helps while the phase
    # ASSEMBLAGE holds: at a point where a phase appears or vanishes the
    # certificate refuses the warm answer and `equilibrate_certified` falls back
    # to its full multi-start cascade.
    #
    # Ascending from 0 g, the monosulfate-to-carboaluminate switch at 0.5 g cost
    # 125 s on its own. Descending from 4 g — where the assemblage is simple and
    # stable — the same sweep costs 52 s in total, with only the first point and
    # the carbonate-free end paying a cascade. Same answers, three times faster.
    #
    # The intermediate rungs at 0.8, 0.6 and 0.3 g are there to keep each step
    # inside one assemblage; they cost 0.2-0.7 s each and are not optional.
    # Dropping the 0.3 g rung alone — so that the last step crosses the
    # monosulfate boundary in one jump — put the run back up to 2.5 minutes.
    ladder = [
        (4.0, 4.49), (0.0, 4.49),                        # the ferriferous clinker
        (4.0, 2.5), (2.0, 2.5), (1.0, 2.5),              # then the ordinary one,
        (0.8, 2.5), (0.6, 2.5), (0.5, 2.5), (0.3, 2.5), (0.0, 2.5),
    ]
    budgets = [
        Float64.(cs.SM.A) * ustrip.(us"mol", charged(g; Fe2O3 = fe).n)
            for (g, fe) in ladder
    ]
    states, certs = equilibrate_path(
        charged(first(ladder)[1]; Fe2O3 = first(ladder)[2]), budgets; model = model
    )
    @test all(c.optimal for c in certs)
    at = Dict(ladder[i] => read_off(states[i]) for i in eachindex(ladder))

    @testset "the iron goes where Cemdata18 says it goes" begin
        # Fig. 14D of [Lothenbach2019]: with Cemdata18 "close to 100 % of the
        # iron is bound by the siliceous hydrogarnet solid solution", against
        # about 80 % in hemi-/monocarbonate with Cemdata07. That is the sharpest
        # quantitative claim the figure makes, and it does not depend on the
        # cement — the phase is simply the only stable home for Fe(III) here.
        for key in keys(at)
            @test at[key].fe_in_hg ≈ 1.0 atol = 1.0e-3
        end
    end

    @testset "the carboaluminate sequence, and what gates it" begin
        # THE SEQUENCE DOES NOT APPEAR ON EVERY CEM I, and finding out why was
        # the work. On the oxide analysis this repository uses elsewhere —
        # Fe2O3 = 4.49 % — no carbonate AFm forms at ANY limestone content, the
        # added calcite stays inert, and every aluminum ends in the mixed Al-Fe
        # siliceous hydrogarnet. That is the behavior [Lothenbach2019] §3.1
        # attributes to Cemdata07 and says Cemdata18 cures.
        #
        # It is not a defect. `C3AFS0.84H4.32` takes ONE aluminum per iron, so
        # the phase is capped by the iron available. Measured here at 4 g of
        # limestone, sweeping Fe2O3 with the balance going to CaO:
        #
        #     Fe2O3 %   4.49     3.50     2.50     1.50     1.00     0.50
        #     hydrogarnet 0.0487   0.0396   0.0283   0.0170   0.0113   0.0057
        #     monocarbonate 0.0000 0.0031   0.0083   0.0135   0.0161   0.0187
        #
        # The hydrogarnet tracks the iron one for one. Where there is enough
        # iron to pair with every aluminum, nothing is left to make a
        # carboaluminate; below about 3.5 % there is, and the published sequence
        # returns. Fig. 7A's cement is evidently on that side.
        @test at[(0.0, 4.49)].Mc == 0
        @test at[(4.0, 4.49)].Mc == 0
        @test at[(4.0, 4.49)].calcite ≈ 4.0 / oxide_mass("CaCO3") rtol = 1.0e-3

        # At 2.5 % Fe2O3 the figure is reproduced step by step.
        a0, a05, a1, a4 = at[(0.0, 2.5)], at[(0.5, 2.5)], at[(1.0, 2.5)], at[(4.0, 2.5)]

        # Without carbonate: monosulfate, and ettringite held down.
        @test a0.Ms > 1.0e-3
        @test a0.Mc == 0 && a0.Hc == 0
        @test a0.calcite == 0

        # Half a gram is enough to destroy the monosulfate outright, and
        # hemicarbonate appears as the intermediate both papers describe.
        @test a05.Ms == 0
        @test a05.Hc > 1.0e-3
        # Ettringite is stabilized by the same step — the mechanism the whole
        # figure is about. It gains more than 40 %.
        @test a05.AFt > 1.4 * a0.AFt

        # By one gram the hemicarbonate is spent, the monocarbonate is at its
        # plateau, and calcite starts to survive undissolved.
        @test a1.Hc == 0
        @test a1.Mc > 1.0e-3
        @test a1.calcite > 0

        # Past saturation the limestone is a filler: what is added is what is
        # left over, to within the little the plateau keeps drifting.
        @test a4.calcite - a1.calcite ≈ 3.0 / oxide_mass("CaCO3") rtol = 0.05
        @test a4.Mc ≈ a1.Mc rtol = 0.05
    end
end
