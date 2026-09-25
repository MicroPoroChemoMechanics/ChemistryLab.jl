# Alkali uptake by C-S-H, against measured solutions.
#
# Hong & Glasser, "Alkali binding in cement pastes: Part I. The C-S-H phase",
# Cem. Concr. Res. 29 (1999) 1893-1903. Synthetic C-S-H at four Ca/Si ratios,
# 0.6 g of soft-dried gel in 9 mL of NaOH or KOH solution, 1 to 300 mM, at
# 20 °C; 48 analyzed solutions (their Table 1) and the distribution ratios Rd
# computed from them (their Table 2), from data/literature/HongGlasser1999.json.
#
# Read the circularity first. Cemdata18 models alkali uptake with the Na and K
# end members of its CSHQ solid solution, and its paper (§2.7, Fig. 10) says
# the Gibbs energies of its alkali end members were fine-tuned on these very
# isotherms. Reproducing the alkali checks that the model is applied as it was
# fitted; the calcium and the pH were not fitted, and they are the part of the
# comparison that can disagree.

include("reference_species.jl")

@testsection "Hong & Glasser (1999) alkali uptake by C-S-H" begin

    subs = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
    byname = Dict(symbol(s) => s for s in subs)
    cshq = ["CSHQ-TobD", "CSHQ-TobH", "CSHQ-JenH", "CSHQ-JenD"]
    model = HKFActivityModel(å = 0.0, Ḃ = gems_bdot(), Kₙ = 0.0)

    hg(q) = literature_value("HongGlasser1999", q)
    T = (273.15 + ustrip(hg("temperature_C"))) * u"K"
    solid_g = ustrip(us"g", hg("solid_mass"))
    volume_mL = ustrip(u"mL", hg("solution_volume"))
    water_fraction(r) = r < 1.35 ? ustrip(hg("water_content_low_CaSi")) :
        ustrip(hg("water_content_high_CaSi"))
    M(s) = ustrip(us"g/mol", Species(s)[:M])

    # Concentrations in mmol/L (= mol/m³), Rd in mL/g.
    mM(v) = ustrip.(u"mol/m^3", v)
    measured = let t = literature_table("HongGlasser1999", "aqueous_phase")
        (
            system = t.system, Ca_Si = t.Ca_Si, target = mM(t.target), alkali = mM(t.alkali),
            Ca = mM(t.Ca), cation_charge = mM(t.cation_charge), OH = mM(t.OH), pH = t.pH,
        )
    end
    binding = let t = literature_table("HongGlasser1999", "binding")
        (
            system = t.system, Ca_Si = t.Ca_Si, target = mM(t.target), initial = mM(t.initial),
            final = mM(t.final), Rd = ustrip.(u"mL/g", t.Rd),
        )
    end

    @testset "the published rows are what they say" begin
        # The printed cation sum is Na (or K) + 2 Ca, to the rounding of the
        # print, in every row but one: K at 15 mM on the Ca/Si 1.8 gel prints
        # 45.0 for 14.6 + 2 × 15.5 = 45.6.
        off = [
            (s, r, round(t)) for (s, r, t, a, ca, c) in zip(
                    measured.system, measured.Ca_Si, measured.target, measured.alkali,
                    measured.Ca, measured.cation_charge,
                ) if abs(c - (a + 2ca)) > 0.051 + 0.005 * c
        ]
        @test off == [("K", 1.8, 15.0)]
        # And its gap to the measured hydroxide is the silicate nobody measured:
        # the authors call ±5 % satisfactory at high ionic strength. Two of the
        # sixteen rows at 100 and 300 mM fall outside it, one at 9.1 %.
        gap = [
            abs(c - oh) / c for (c, oh, t) in
                zip(measured.cation_charge, measured.OH, measured.target) if t > 99
        ]
        @test length(gap) == 16
        @test count(>(ustrip(hg("charge_balance_tolerance"))), gap) == 2
        @test maximum(gap) ≈ 0.091 atol = 5.0e-4
        # In the dilute solutions, 1 and 5 mM, the same gap reaches 55 % and is
        # 16 % on average: nothing compared there carries much.
        dilute = [
            abs(c - oh) / c for (c, oh, t) in
                zip(measured.cation_charge, measured.OH, measured.target) if t < 6
        ]
        @test maximum(dilute) ≈ 0.55 atol = 5.0e-3
        @test sum(dilute) / length(dilute) ≈ 0.16 atol = 5.0e-3
    end

    # One system per alkali: a budget with no potassium would leave KSiOH, and
    # the K+ species, pinned at the floor with nothing to say.
    function system(alkali)
        member = alkali == "Na" ? "NaSiOH" : "KSiOH"
        sp = speciation(
            subs, vcat(["Portlandite", "Amor-Sl", "Lim"], cshq, [member]);
            aggregate_state = [AS_AQUEOUS], exclude_species = split("H2@ O2@ CH4@")
        )
        return ChemicalSystem(
            sp, CEMDATA_PRIMARIES;
            solid_solutions = [SolidSolutionPhase("CSHQ", [byname[m] for m in vcat(cshq, [member])])],
        )
    end

    # One experiment, brought to one kilogram of solution water. The gel enters
    # as lime, silica and its own water (14 % of the soft-dried mass at Ca/Si
    # 0.85 and 1.2, 18 % above), from the molar masses of the library; the 9 mL
    # of solution are taken as 9 g of water.
    function charged(cs, alkali, r, c0_mM)
        scale = 1000 / volume_mL
        w = water_fraction(r)
        n_Si = solid_g * (1 - w) / (r * M("CaO") + M("SiO2")) * scale
        st = ChemicalState(cs; T = T)
        set_quantity!(st, "Lim", r * n_Si * u"mol")
        set_quantity!(st, "Amor-Sl", n_Si * u"mol")
        set_quantity!(st, "H2O@", (1000 + solid_g * w * scale) / M("H2O") * u"mol")
        set_quantity!(st, alkali == "Na" ? "Na+" : "K+", c0_mM * 1.0e-3 * u"mol")
        set_quantity!(st, "OH-", c0_mM * 1.0e-3 * u"mol")
        return st
    end

    results = Dict{Tuple{String, Float64, Float64}, NamedTuple}()
    for alkali in ("Na", "K")
        cs = system(alkali)
        idx = Dict(symbol(s) => i for (i, s) in enumerate(cs.species))
        A = Float64.(cs.SM.A)
        el = Symbol(alkali)
        member = alkali == "Na" ? "NaSiOH" : "KSiOH"
        per_member = Float64(atoms(cs.species[idx[member]])[el])
        for r in sort(unique(measured.Ca_Si))
            rows = [
                i for i in eachindex(binding.system)
                    if binding.system[i] == alkali && binding.Ca_Si[i] == r
            ]
            # Walked down from the most concentrated solution: the assemblage is
            # the same C-S-H throughout, so each solve starts from the last.
            sort!(rows; by = i -> -binding.initial[i])
            states = [charged(cs, alkali, r, binding.initial[i]) for i in rows]
            budgets = [A * ustrip.(us"mol", st.n) for st in states]
            eqs, certs = equilibrate_path(first(states), budgets; model)
            @test all(c.optimal for c in certs)
            for (i, eq) in zip(rows, eqs)
                n = ustrip.(us"mol", eq.n)
                vol_L = ustrip(uconvert(us"L", volume(eq).liquid))
                dissolved(e) = 1000 * sum(
                    n[k] * Float64(get(atoms_charge(cs.species[k]), e, 0)) for k in cs.idx_aqueous
                ) / vol_L
                c_alk = dissolved(el)                           # mmol/L
                # Rd as the paper defines it: alkali in the solid per gram of the
                # soft-dried gel, over alkali per mL of solution.
                in_solid = n[idx[member]] * per_member * 1000 * volume_mL / 1000   # mmol per experiment
                Rd = (in_solid / solid_g) / (c_alk / 1000)
                results[(alkali, r, round(binding.target[i]))] = (
                    alkali = c_alk, Ca = dissolved(:Ca), pH = pH(eq, model), Rd = Rd,
                    portlandite = n[idx["Portlandite"]],
                )
            end
        end
    end

    # Keyed by the nominal concentration: 100 mM reads back as 99.999... mol/m³.
    at(s, r, t) = results[(s, r, Float64(t))]
    pts = [
        (s, r, round(t), a, ca, ph) for (s, r, t, a, ca, ph) in
            zip(measured.system, measured.Ca_Si, measured.target, measured.alkali, measured.Ca, measured.pH)
    ]
    measured_Rd = Dict(
        (s, r, round(t)) => rd for (s, r, t, rd) in
            zip(binding.system, binding.Ca_Si, binding.target, binding.Rd)
    )

    @testset "the 100 mM rung, as the chapter prints it" begin
        # (alkali, Ca) in mmol/L, pH, Rd in mL/g; pinned at the printed digit.
        printed = Dict(
            ("Na", 0.85) => (69.7, 0.15, 12.805, 6.41), ("Na", 1.2) => (83.1, 0.43, 12.94, 2.99),
            ("Na", 1.5) => (91.5, 3.81, 12.999, 1.3), ("Na", 1.8) => (93.7, 8.33, 13.033, 0.93),
            ("K", 0.85) => (76.1, 0.15, 12.835, 4.94), ("K", 1.2) => (87.6, 0.41, 12.964, 2.36),
            ("K", 1.5) => (94.6, 3.62, 13.015, 1.03), ("K", 1.8) => (96.4, 8.06, 13.046, 0.74),
        )
        for ((s, r), (a, ca, ph, rd)) in printed
            c = at(s, r, 100)
            @test c.alkali ≈ a atol = 0.05
            @test c.Ca ≈ ca atol = 0.005
            @test c.pH ≈ ph atol = 5.0e-4
            @test c.Rd ≈ rd atol = 0.005
        end
    end

    @testset "where the model and the measurements part" begin
        # The alkali: less of it stays in solution than was measured at 46 of the
        # 48 points, so the model binds more than the gel did.
        @test count(p -> at(p[1], p[2], p[3]).alkali < p[4], pts) == 46
        # Rd at 100 mM, computed over measured: from 0.95 to 2.7.
        ratio = [at(s, r, 100).Rd / measured_Rd[(s, r, 100.0)] for (s, r) in unique([(p[1], p[2]) for p in pts])]
        @test minimum(ratio) ≈ 0.95 atol = 0.005
        @test maximum(ratio) ≈ 2.74 atol = 0.005

        # The pH, which nothing was fitted to: within 0.072 of the electrode from
        # 15 to 100 mM wherever the gel is below Ca/Si 1.8, and low at 300 mM,
        # by 0.08 to 0.21 at every Ca/Si.
        mid = [abs(at(p[1], p[2], p[3]).pH - p[6]) for p in pts if p[2] <= 1.5 && p[3] in (15, 50, 100)]
        @test maximum(mid) ≈ 0.072 atol = 5.0e-4
        low = [p[6] - at(p[1], p[2], p[3]).pH for p in pts if p[3] == 300]
        @test all(>(0), low)
        @test minimum(low) ≈ 0.078 atol = 5.0e-4
        @test maximum(low) ≈ 0.21 atol = 5.0e-4

        # The calcium: low at Ca/Si 1.2 and 1.5 from 15 mM on, high at 1.8. At
        # 1.8 the CSHQ model cannot hold all of the calcium and portlandite
        # precipitates, 3.7 to 6.8 % of the solid, where the gel carried 0.2 %.
        @test all(at(p[1], p[2], p[3]).Ca < p[5] for p in pts if p[2] in (1.2, 1.5) && 15 <= p[3] <= 100)
        @test all(at(p[1], p[2], p[3]).Ca > p[5] for p in pts if p[2] == 1.8)
        ch = [
            100 * at(p[1], p[2], p[3]).portlandite * M("Ca(OH)2") / (solid_g * 1000 / volume_mL)
                for p in pts if p[2] == 1.8
        ]
        @test minimum(ch) ≈ 3.7 atol = 0.05
        @test maximum(ch) ≈ 6.8 atol = 0.05
        @test all(at(p[1], p[2], p[3]).portlandite == 0 for p in pts if p[2] < 1.8)
    end
end
