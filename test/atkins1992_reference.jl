# Agreement with a measured solution, not with another model.
#
# Atkins, Bennett, Dawes, Glasser, Kindness & Read, "A thermodynamic model for
# blended cements", Cem. Concr. Res. 22 (1992) 497-502.
#
# Nearly every cement-chemistry reference in this suite is a calculation:
# Cemdata18's tables are derived quantities, and a cross-check against GEMS or
# Reaktoro compares two codes reading the same database. Atkins et al. is
# different — their Table 2 reports ANALYZED solution compositions for slurries
# of synthetic hydrates in CO2-free water at 25 °C, two to four solids at a
# time. That is an experiment, and it can disagree with the database.
#
# Only one of their mixtures is usable quantitatively, and the reason the others
# are not is recorded at the end of this file: it is worth knowing before anyone
# spends an afternoon on them.

include("reference_species.jl")

@testsection "Atkins et al. (1992) hydrate slurries" begin

    subs = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
    byname = Dict(symbol(s) => s for s in subs)

    # The Ca-Al-Si-S-H-O corner of CEMDATA18. Gibbsite is left out on the
    # database's own instruction (Cemdata18 §2.1: its precipitation "should be
    # suppressed for calculations at ambient temperatures, where
    # microcrystalline Al(OH)3 will form instead"), and leaving it in does
    # change the answer — it takes over from AH3 in the sulfate-poor mixtures.
    pure = split(
        "ettringite C3AH6 C4AH13 C4AH19 C2AH7.5 CAH10 monosulphate12 " *
            "monosulphate14 straetlingite AlOHmic Portlandite Gp Anh " *
            "hemihydrate Amor-Sl C3AS0.41H5.18 C3AS0.84H4.32 Lim"
    )
    cshq = ["CSHQ-TobD", "CSHQ-TobH", "CSHQ-JenH", "CSHQ-JenD"]

    species = speciation(
        subs, vcat(pure, cshq);
        aggregate_state = [AS_AQUEOUS], exclude_species = split("H2@ O2@ CH4@")
    )
    cs = ChemicalSystem(
        species, CEMDATA_PRIMARIES;
        solid_solutions = [SolidSolutionPhase("CSHQ", [byname[m] for m in cshq])],
    )
    idx = Dict(symbol(s) => i for (i, s) in enumerate(cs.species))

    # CEMDATA18 carries no ion-size parameter, so a GEM-Selektor run of it
    # starts from å = 0; this is the same model the CEM I cross-check uses.
    model = HKFActivityModel(å = 0.0, Ḃ = gems_bdot(), Kₙ = 0.0)

    # Atkins' Table 2, from data/literature/Atkins1992.json: the analyzed
    # solution of mixture 32, in mmol/L.
    measured = let t = literature_table("Atkins1992", "mixture_32_solution")
        Dict(zip(t.element, ustrip.(u"mol/m^3", t.measured)))   # mol/m³ = mmol/L
    end
    Ca_Si = literature_value("Atkins1992", "mixture_32_Ca_Si")

    # Experiment 32: ettringite + a C-S-H of nominal Ca/Si = 0.9, in water.
    # The C-S-H enters as its bulk oxides — only the element vector `b` reaches
    # the solver, so lime and amorphous silica here are bookkeeping and not a
    # claim about which solids exist.
    state = ChemicalState(cs)
    set_quantity!(state, "ettringite", 0.002u"mol")
    set_quantity!(state, "Lim", 0.018u"mol")        # Ca
    set_quantity!(state, "Amor-Sl", 0.018u"mol" / Ca_Si)   # Si, so Ca/Si = 0.9
    set_quantity!(state, "H2O@", 1.0u"kg")
    b = Float64.(cs.SM.A) * ustrip.(us"mol", state.n)

    eq, cert = equilibrate_certified(state; model = model, b = b)
    @test cert.optimal

    n = ustrip.(us"mol", eq.n)
    # Dissolved totals, in mmol/L, of any equilibrium state.
    function totals(e)
        m = ustrip.(us"mol", e.n)
        vol = ustrip(uconvert(us"L", volume(e).liquid))
        return el -> 1000 * sum(
            m[i] * Float64(get(atoms_charge(cs.species[i]), el, 0)) for i in cs.idx_aqueous
        ) / vol
    end
    total = totals(eq)

    @testset "the assemblage Atkins reports" begin
        # AFt and C-S-H both survive, and no portlandite: that much the
        # experiment and the calculation agree on.
        @test n[idx["ettringite"]] > 1.0e-4
        @test sum(n[idx[m]] for m in cshq) > 1.0e-3
        @test n[idx["Portlandite"]] < 1.0e-9

        # Not quite "no new phases", though. A little strätlingite comes out,
        # which Atkins did not see. It is ~2 % of the solid and does not carry
        # the disagreements below, but it is there.
        @test 0 < n[idx["straetlingite"]] < 0.02 * sum(n[idx[m]] for m in cshq)
    end

    @testset "aluminium and pH: agreement" begin
        # Measured at 6 months: Al 0.136 mmol/L, pH 11.0.
        #
        # Al is the strong result. Across a sweep of the C-S-H Ca/Si from 0.75
        # to 1.0 and a fourfold change in solid loading it never leaves
        # 0.138-0.159 mmol/L, so this is the database answering, not a fit.
        @test total(:Al) ≈ measured["Al"] rtol = 0.15
        @test pH(eq, model) ≈ literature_value("Atkins1992", "mixture_32_pH") atol = 0.4
    end

    @testset "silicon: over-predicted, for a reason the paper gives" begin
        # Measured 0.076 mmol/L; calculated ~0.40, a factor of five.
        #
        # This is not a transcription slip and not a ChemistryLab artifact:
        # Atkins' own model gave 0.448 mmol/L on the same mixture, and they
        # explain it. Electron microscopy of the solid showed the ettringite had
        # taken up silicon, "substituting on average for 40 % of the available
        # SO4 sites", which "would tend to lower aqueous Si concentrations, and
        # increase SO4 levels as observed". CEMDATA18 carries no Si-bearing AFt
        # end-member either, so the same silicon has nowhere to go and stays in
        # solution. A database that gained one would move this number, which is
        # why it is pinned rather than skipped.
        @test total(:Si) / measured["Si"] > 3.5
        @test total(:Si) / measured["Si"] < 7.0
    end

    @testset "calcium: a disagreement the model cannot absorb" begin
        # Measured 1.95 mmol/L; calculated 2.72. Worth stating plainly, because
        # the obvious escape does not work: Atkins report that their C-S-H
        # dissolved incongruently over the test, which would have lowered its
        # Ca/Si below the nominal 0.9, and a lower Ca/Si ought to mean less
        # calcium in solution. Swept from Ca/Si 0.75 to 1.0 the calculated
        # calcium goes through a MINIMUM of about 2.7 mmol/L near 0.85-0.90 and
        # rises to 3.4 at 0.75 — it never approaches 1.95. The floor is a
        # property of the CSHQ model, not a free parameter, so the gap is real.
        @test total(:Ca) ≈ 2.72 rtol = 0.1
        @test total(:Ca) > 1.3 * measured["Ca"]
        # The page prints a "calculated" column beside Atkins' measured one.
        # Those five numbers are this package's answer, and until now only the
        # calcium had an assertion anchored on it — at 10 %, which is wider than
        # the three digits printed. Pinned here at the printed precision. The
        # assertions above stay: they say something different, namely how far
        # the answer is from the MEASUREMENT, which is what the page is about.
        @test total(:Ca) ≈ 2.72 atol = 5.0e-3
        @test total(:Al) ≈ 0.149 atol = 5.0e-4
        @test total(:Si) ≈ 0.396 atol = 5.0e-4
        @test total(:S) ≈ 1.18 atol = 5.0e-3
        @test pH(eq, model) ≈ 11.33 atol = 5.0e-3

        # Sulfate lands close (1.18 against 1.08 measured), but it is the one
        # number here that moves freely with the C-S-H Ca/Si — 0.55 at Ca/Si 1.0,
        # 3.2 at 0.75 — so the agreement is not evidence of much. Atkins' own
        # 1992 model gave 0.597, so CEMDATA18 is nearer; that is the honest
        # claim, and it is bounded loosely on purpose.
        @test 0.5 < total(:S) < 2.5
    end

    @testset "the calcium floor, across the C-S-H Ca/Si" begin
        # Atkins report that their C-S-H dissolved incongruently, which lowers
        # its Ca/Si below the nominal 0.9. The lime is held and the silica
        # varied, walked from the nominal mixture already solved.
        ratios = [1.0, 0.85, 0.8, 0.75]
        function budget(r)
            st = ChemicalState(cs)
            set_quantity!(st, "ettringite", 0.002u"mol")
            set_quantity!(st, "Lim", 0.018u"mol")
            set_quantity!(st, "Amor-Sl", 0.018u"mol" / r)
            set_quantity!(st, "H2O@", 1.0u"kg")
            return Float64.(cs.SM.A) * ustrip.(us"mol", st.n)
        end
        sweep, certs = equilibrate_path(eq, budget.(ratios); model = model)
        @test all(c.optimal for c in certs)
        # Pinned at the precision the page prints, half the last digit.
        printed = Dict(
            1.0 => (2.72, 0.142, 0.339, 0.76, 11.47),
            0.85 => (2.83, 0.152, 0.429, 1.5, 11.24),
            0.8 => (3.06, 0.154, 0.465, 1.92, 11.13),
            0.75 => (3.46, 0.157, 0.504, 2.51, 11.0),
        )
        for (r, e) in zip(ratios, sweep)
            t = totals(e)
            Ca, Al, Si, S, pHr = printed[r]
            @test t(:Ca) ≈ Ca atol = 5.0e-3
            @test t(:Al) ≈ Al atol = 5.0e-4
            @test t(:Si) ≈ Si atol = 5.0e-4
            @test t(:S) ≈ S atol = 5.0e-3
            @test pH(e, model) ≈ pHr atol = 5.0e-3
            # The floor: nowhere on the sweep does the calcium reach the measurement.
            @test t(:Ca) > 1.3 * measured["Ca"]
        end
        # The calcium rises as the Ca/Si falls below the nominal ratio, which
        # is the direction incongruent dissolution would take it, and is flat
        # between 0.9 and 1.0.
        Ca_at = [totals(e)(:Ca) for e in sweep]                      # 1.0, 0.85, 0.8, 0.75
        @test Ca_at[4] > Ca_at[3] > Ca_at[2] > total(:Ca)
        @test abs(Ca_at[1] - total(:Ca)) < 0.01
    end

    # ── Why the rest of Table 2 is not here ──────────────────────────────────
    #
    # Ten mixtures are tabulated; nine are not testable, and the reasons divide
    # in three. Recorded so that the next reader does not rediscover them:
    #
    #  - The solid-to-water ratio is not reported, and these assemblages are not
    #    invariant. Wherever a C-S-H buffers the calcium, a fourfold change in
    #    solid loading moves the answer by a factor of two: mixture 1
    #    (AFt + CSH 1.7) gives Ca = 8.0 mmol/L at one loading and 16.0 at four
    #    times that, against 2.67 measured. Nothing can be concluded from a
    #    number that depends on an unreported quantity.
    #
    #  - Four mixtures grew a phase Atkins identifies as metastable, so the
    #    measured solution is not an equilibrium one and should not match: AFm
    #    in mixture 3, C4AH13 in 8 and 21, siliceous hydrogarnet in 30. The
    #    paper says so itself and uses them to argue the point.
    #
    #  - Mixtures 10 and 21 are in 0.4 M NaOH. At that ionic strength the B-dot
    #    term in `HKFActivityModel` is past the range it was fitted for
    #    (`I ≲ 1 mol/kg` is the documented limit, and less for the water
    #    activity); a Pitzer model would be the honest instrument there.
    #
    # Their Table 3 — pore fluids of 5-year-old OPC, 30 % BFS and 30 % FA pastes
    # — is a better target for a future test, since alkali and pH are reported
    # and the binder compositions are given; it needs the alkali uptake of the
    # CSHQ Na/K end-members and belongs with the CEM I/II examples.
end
