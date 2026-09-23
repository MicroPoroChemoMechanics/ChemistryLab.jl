# Agreement with the Cemdata18 paper's own tables.
#
# Lothenbach, Kulik, Matschei, Balonis, Baquerizo, Dilnesa, Miron & Myers,
# "Cemdata18: A chemical thermodynamic database for hydrated Portland cements
# and alkali-activated materials", Cem. Concr. Res. 115 (2019) 472-506,
# https://doi.org/10.1016/j.cemconres.2018.04.018
#
# `data/cemdata18-thermofun.json` is a vendored copy of a file maintained
# elsewhere. Nothing in the rest of the suite would notice if one Gibbs energy
# in it drifted by a few kJ/mol: every equilibrium would still converge, every
# figure would still be drawn, and the numbers would simply be wrong. The paper
# publishes BOTH the standard formation properties (Table 1, Table 3, Appendix
# D) and the solubility products derived from them (Table 2, Table 3), which
# makes the file checkable against a source rather than against itself.
#
# The same round trip is already run on the zeolite extension in
# `test/zeolites.jl`. This file does it for the base.

using JSON

@testsection "Cemdata18 reference tables" begin

    sp = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-thermofun.json")))
    raw = JSON.parsefile(datapath("cemdata18-thermofun.json"); dicttype = Dict{String, Any})
    rec = Dict(String(s["symbol"]) => s for s in raw["substances"])

    # ── Table 2 and Table 3: log Ks0 against ΔfG° ────────────────────────────
    #
    # Table 2 writes its dissolution reactions over Al(OH)4⁻, Fe(OH)4⁻ and
    # SiO(OH)3⁻ (H3SiO4⁻ for thaumasite). CEMDATA18's GEMS primaries are AlO2⁻,
    # FeO2⁻ and HSiO3⁻, which differ from those by 2, 2 and 1 H2O respectively.
    # Rather than transcribe the substitution AND the paper's water
    # coefficients — two chances to get it wrong, one of which would be silently
    # absorbed by the other — only the non-water products are transcribed here,
    # and the water is recovered from the hydrogen balance. The oxygen and
    # charge balances are then free checks: they can only close if the
    # transcribed products are right.
    #
    # Coefficients are the paper's, phase by phase, in its own order.
    table2 = [
        # AFt
        ("ettringite", -44.9, Dict("Ca+2" => 6, "AlO2-" => 2, "SO4-2" => 3, "OH-" => 4)),
        ("tricarboalu", -46.5, Dict("Ca+2" => 6, "AlO2-" => 2, "CO3-2" => 3, "OH-" => 4)),
        ("Fe-ettringite", -44.0, Dict("Ca+2" => 6, "FeO2-" => 2, "SO4-2" => 3, "OH-" => 4)),
        ("thaumasite", -24.75, Dict("Ca+2" => 3, "HSiO3-" => 1, "SO4-2" => 1, "CO3-2" => 1, "OH-" => 1)),
        # Hydrogarnet
        ("C3AH6", -20.5, Dict("Ca+2" => 3, "AlO2-" => 2, "OH-" => 4)),
        ("C3AS0.41H5.18", -25.35, Dict("Ca+2" => 3, "AlO2-" => 2, "HSiO3-" => 0.41, "OH-" => 3.59)),
        ("C3AS0.84H4.32", -26.7, Dict("Ca+2" => 3, "AlO2-" => 2, "HSiO3-" => 0.84, "OH-" => 3.16)),
        ("C3FH6", -26.3, Dict("Ca+2" => 3, "FeO2-" => 2, "OH-" => 4)),
        ("C3FS0.84H4.32", -32.5, Dict("Ca+2" => 3, "FeO2-" => 2, "HSiO3-" => 0.84, "OH-" => 3.16)),
        ("C3AFS0.84H4.32", -30.2, Dict("Ca+2" => 3, "AlO2-" => 1, "FeO2-" => 1, "HSiO3-" => 0.84, "OH-" => 3.16)),
        ("C3FS1.34H3.32", -34.2, Dict("Ca+2" => 3, "FeO2-" => 2, "HSiO3-" => 1.34, "OH-" => 2.66)),
        # AFm
        ("C4AH19", -25.45, Dict("Ca+2" => 4, "AlO2-" => 2, "OH-" => 6)),
        ("C4AH13", -25.25, Dict("Ca+2" => 4, "AlO2-" => 2, "OH-" => 6)),
        ("C2AH7.5", -13.8, Dict("Ca+2" => 2, "AlO2-" => 2, "OH-" => 2)),
        ("CAH10", -7.6, Dict("Ca+2" => 1, "AlO2-" => 2)),
        ("hemicarbonate", -29.13, Dict("Ca+2" => 4, "AlO2-" => 2, "CO3-2" => 0.5, "OH-" => 5)),
        ("monocarbonate", -31.47, Dict("Ca+2" => 4, "AlO2-" => 2, "CO3-2" => 1, "OH-" => 4)),
        ("monosulphate14", -29.26, Dict("Ca+2" => 4, "AlO2-" => 2, "SO4-2" => 1, "OH-" => 4)),
        ("monosulphate12", -29.23, Dict("Ca+2" => 4, "AlO2-" => 2, "SO4-2" => 1, "OH-" => 4)),
        ("straetlingite", -19.7, Dict("Ca+2" => 2, "AlO2-" => 2, "HSiO3-" => 1, "OH-" => 1)),
        ("C4AClH10", -27.27, Dict("Ca+2" => 4, "AlO2-" => 2, "Cl-" => 2, "OH-" => 4)),      # Friedel's salt
        ("C4AsClH12", -28.53, Dict("Ca+2" => 4, "AlO2-" => 2, "Cl-" => 1, "SO4-2" => 0.5, "OH-" => 4)),  # Kuzel's salt
        ("mononitrate", -28.67, Dict("Ca+2" => 4, "AlO2-" => 2, "NO3-" => 2, "OH-" => 4)),
        # Fe-AFm
        ("C4FH13", -30.75, Dict("Ca+2" => 4, "FeO2-" => 2, "OH-" => 6)),
        ("Fe-hemicarbonate", -30.83, Dict("Ca+2" => 4, "FeO2-" => 2, "CO3-2" => 0.5, "OH-" => 5)),
        ("Femonocarbonate", -34.59, Dict("Ca+2" => 4, "FeO2-" => 2, "CO3-2" => 1, "OH-" => 4)),
        ("Fe-monosulphate", -31.57, Dict("Ca+2" => 4, "FeO2-" => 2, "SO4-2" => 1, "OH-" => 4)),
        # Sulfates
        ("Anh", -4.357, Dict("Ca+2" => 1, "SO4-2" => 1)),
        ("Gp", -4.581, Dict("Ca+2" => 1, "SO4-2" => 1)),
        ("hemihydrate", -3.59, Dict("Ca+2" => 1, "SO4-2" => 1)),
        ("syngenite", -7.2, Dict("K+" => 2, "Ca+2" => 1, "SO4-2" => 2)),
        # (Hydr)oxides
        ("AlOHam", 0.24, Dict("AlO2-" => 1, "OH-" => -1)),
        ("AlOHmic", -0.67, Dict("AlO2-" => 1, "OH-" => -1)),
        ("Gbs", -1.12, Dict("AlO2-" => 1, "OH-" => -1)),
        ("FeOOHmic", -5.6, Dict("FeO2-" => 1, "OH-" => -1)),
        ("Gt", -8.6, Dict("FeO2-" => 1, "OH-" => -1)),
        ("Portlandite", -5.2, Dict("Ca+2" => 1, "OH-" => 2)),
        ("Amor-Sl", -2.714, Dict("SiO2@" => 1)),
        ("Qtz", -3.746, Dict("SiO2@" => 1)),
        # Hydrotalcite-pyroaurite, M-S-H
        ("Mg3AlC0.5OH", -33.29, Dict("Mg+2" => 3, "AlO2-" => 1, "CO3-2" => 0.5, "OH-" => 4)),
        ("Mg3FeC0.5OH", -33.64, Dict("Mg+2" => 3, "FeO2-" => 1, "CO3-2" => 0.5, "OH-" => 4)),
        ("M075SH", -28.8, Dict("Mg+2" => 1.5, "SiO2@" => 2, "OH-" => 3)),
        ("M15SH", -23.57, Dict("Mg+2" => 1.5, "SiO2@" => 1, "OH-" => 3)),
        # Zeolites
        ("zeoliteP_Ca", -20.3, Dict("Ca+2" => 1, "AlO2-" => 2, "SiO2@" => 2)),
        ("natrolite", -30.2, Dict("Na+" => 2, "AlO2-" => 2, "SiO2@" => 3)),
        ("chabazite", -25.8, Dict("Ca+2" => 1, "AlO2-" => 2, "SiO2@" => 4)),
        ("zeoliteX", -20.1, Dict("Na+" => 2, "AlO2-" => 2, "SiO2@" => 2.5)),
        ("zeoliteY", -25.0, Dict("Na+" => 2, "AlO2-" => 2, "SiO2@" => 4)),
        # Table 3 — the two alternative hydrotalcite modules. `hydrotalcite` is
        # the single phase recommended for PC; the `M*A-OH-LDH` trio is the
        # ideal solid solution recommended for alkali-activated materials.
        ("hydrotalcite", -56.02, Dict("Mg+2" => 4, "AlO2-" => 2, "OH-" => 6)),
        ("M4A-OH-LDH", -49.7, Dict("Mg+2" => 4, "AlO2-" => 2, "OH-" => 6)),
        ("M6A-OH-LDH", -72.0, Dict("Mg+2" => 6, "AlO2-" => 2, "OH-" => 10)),
        ("M8A-OH-LDH", -94.3, Dict("Mg+2" => 8, "AlO2-" => 2, "OH-" => 14)),
    ]

    # The two M-S-H end-members do not close. Their tabulated ΔfG° and their
    # published log Ks0 disagree by 2.76 and 2.26 kJ/mol when read through
    # CEMDATA18's own aqueous energies — about half a log unit. Every other
    # phase in the table closes to better than 0.04, so this is a property of
    # the source and not of the transcription.
    #
    # It is one of two disagreements these two phases carry, and the second is
    # in the testset below: their (ΔfG°, ΔfH°, S°) triplet is not self
    # consistent either. Cemdata18's Table 1 footnote r says why — S° and Cp°
    # for M-S-H were "estimated from Cp and S of talc, chrysotile and H2O"
    # rather than measured ([Nied2016]) — and an estimated entropy that was
    # never reconciled with the tabulated Gibbs energy is exactly what this
    # looks like.
    #
    # Pinned at the observed offset rather than hidden behind a loose tolerance,
    # so that a corrected upstream file shows up as a failing test, not silence.
    msh_offset = Dict("M075SH" => 0.4828, "M15SH" => 0.3951)

    R = 8.31446261815324                     # J/(mol·K), CODATA
    RTln10 = R * 298.15 * log(10)
    # The TABULATED ΔfG°, read from the file, not `ΔₐG⁰(T = 298.15)`. The two
    # are the same number for 220 of the 228 substances and the testset after
    # this one pins which eight they are not; using the tabulated value here
    # keeps this check on the question it is asking — does the vendored file
    # agree with the paper — instead of mixing it with how the package rebuilds
    # a Gibbs energy.
    Gf(k) = Float64(rec[k]["sm_gibbs_energy"]["values"][1])
    nel(k, e) = Float64(get(atoms_charge(sp[k]), e, 0))

    @testset "Table 2/3: log Ks0 closes against ΔfG° (298.15 K, 1 bar)" begin
        # 48 rows of Table 2 plus the 4 of Table 3.
        @test length(table2) == 52

        worst_ordinary = 0.0
        for (s, logK_published, products) in table2
            @test haskey(sp, s)

            # Water from the hydrogen balance, then oxygen and charge as checks.
            nH2O = (nel(s, :H) - sum(ν * nel(k, :H) for (k, ν) in products)) / 2
            @test sum(ν * nel(k, :O) for (k, ν) in products) + nH2O ≈ nel(s, :O) atol = 1.0e-9
            @test sum(ν * nel(k, :Zz) for (k, ν) in products) ≈ nel(s, :Zz) atol = 1.0e-9

            ΔrG = sum(ν * Gf(k) for (k, ν) in products) + nH2O * Gf("H2O@") - Gf(s)
            logK = -ΔrG / RTln10

            if haskey(msh_offset, s)
                @test logK - logK_published ≈ msh_offset[s] atol = 0.005
            else
                @test logK ≈ logK_published atol = 0.05
                worst_ordinary = max(worst_ordinary, abs(logK - logK_published))
            end
        end

        # Not merely "each within tolerance": the whole table agrees to better
        # than the rounding of the published values.
        @test worst_ordinary < 0.05
    end

    # ── Where the rebuilt Gibbs energy meets the tabulated one ───────────────
    #
    # `ΔₐG⁰(T)` is the apparent Gibbs energy of formation, and the package forms
    # it from ΔfH° and S° rather than reading ΔfG° off the file. At the
    # reference point T = 298.15 K the two must therefore agree — and they do,
    # to the last bit, for 220 of the 228 substances. For eight they do not.
    #
    # The eight are not a random selection, and they are not the ones missing a
    # heat-capacity block: six of them have one. What they share is that
    # Cemdata18 says their entropy and heat capacity were ESTIMATED rather than
    # measured — the six alkali C-S-H end members by the linear Ca/Si relations
    # of Table 4 (the paper's Eqs 2a and 2b), the two M-S-H end members from
    # talc, chrysotile and water (Table 1, footnote r, after [Nied2016]). An
    # estimated S° that was never reconciled with the tabulated ΔfG° leaves the
    # triplet (ΔfG°, ΔfH°, S°) inconsistent, and rebuilding the third from the
    # other two is what makes that visible: 0.70 J/K/mol of entropy for the
    # alkali members, 4.6 for M075SH.
    #
    # The control is the five zeolites, which carry no heat-capacity block at
    # all and still land on their tabulated value exactly. So this is the data,
    # not the code path.
    #
    # Asserted as an exhaustive list. A ninth substance joining it is a change
    # worth being told about.
    @testset "ΔₐG⁰(298.15 K) is the tabulated ΔfG°, with eight exceptions" begin
        rebuilt = Dict(                                      # J/mol
            "ECSH1-KSH" => -243.7, "ECSH2-KSH" => -243.7,
            "ECSH1-NaSH" => -207.6, "ECSH2-NaSH" => -207.6,
            "KSiOH" => -208.6, "NaSiOH" => -208.6,
            "M075SH" => -1364.8, "M15SH" => -1088.6,
        )
        found = String[]
        for (k, entry) in rec
            haskey(entry, "sm_gibbs_energy") || continue
            gap = sp[k].ΔₐG⁰(T = 298.15) - Float64(entry["sm_gibbs_energy"]["values"][1])
            abs(gap) > 1.0 || continue
            push!(found, k)
            @test haskey(rebuilt, k)
            @test gap ≈ get(rebuilt, k, NaN) atol = 1.0
        end
        @test sort(found) == sort(collect(keys(rebuilt)))

        # Six of the eight carry a heat-capacity block, so a missing one is not
        # the explanation.
        @test count(k -> haskey(rec[k], "TPMethods"), keys(rebuilt)) == 6

        # And the control.
        for k in ("chabazite", "natrolite", "zeoliteP_Ca", "zeoliteX", "zeoliteY")
            @test !haskey(rec[k], "TPMethods")
            @test sp[k].ΔₐG⁰(T = 298.15) ≈
                Float64(rec[k]["sm_gibbs_energy"]["values"][1]) atol = 1.0
        end
    end

    # ── Appendix D: standard properties and HKF parameters ───────────────────
    #
    # Table 2 exercises ΔfG° of the aqueous primaries hard — a drift in any one
    # of them would break dozens of reactions at once — but says nothing about
    # S°, Cp°, V° or the HKF equation-of-state coefficients, which are what
    # carry the database away from 25 °C and 1 bar. Those are checked here
    # directly against Table D.1 and Table D.2.
    #
    # Two conventions of the printed tables have to be undone:
    #
    #  - Table D.1 scales its HKF columns as a1·10, a2·10⁻², a4·10⁻⁴, c2·10⁻⁴
    #    and ω0·10⁻⁵, in the calorimetric units the HKF papers use.
    #  - Table D.1's volume column is headed "V⁰ (J/bar)" but carries cm³/mol:
    #    Ca²⁺ is listed at -18.44, and -18.44 J/bar would be -184.4 cm³/mol.
    #    The file stores J/bar (-1.8439), which is the same -18.44 cm³/mol.
    #    Table D.2's gas column really is J/bar (2479 J/bar = 24.79 L/mol, the
    #    ideal-gas molar volume at 298.15 K and 1 bar), so the two tables share
    #    a header and not a unit.
    #
    # ΔG⁰ and ΔH⁰ below are in kJ/mol as printed; V⁰ in cm³/mol.
    tableD1 = [
        #  symbol          ΔG⁰       ΔH⁰       S⁰       Cp⁰      V⁰      a1·10    a2·10⁻²   a3       a4·10⁻⁴  c1        c2·10⁻⁴   ω0·10⁻⁵
        ("Al+3", -483.71, -530.63, -325.1, -128.7, -45.24, -3.3802, -17.0071, 14.5185, -2.0758, 10.7, -8.06, 2.753),
        ("AlO2-", -827.48, -925.57, -30.21, -49.04, 9.47, 3.7221, 3.9954, -1.5879, -2.9441, 15.2391, -5.4585, 1.7418),
        ("AlOH+2", -692.6, -767.27, -184.93, 55.97, -2.73, 2.0469, -2.7813, 6.8376, -2.6639, 29.7923, -0.3457, 1.7247),
        ("Ca+2", -552.79, -543.07, -56.48, -30.92, -18.44, -0.1947, -7.252, 5.2966, -2.4792, 9.0, -2.522, 1.2366),
        ("Ca(SO4)@", -1310.38, -1448.43, 20.92, -104.6, 4.7, 2.4079, -1.8992, 6.4895, -2.7004, -8.4942, -8.1271, -0.001),
        ("CaOH+", -717.02, -751.65, 28.03, 6.05, 5.76, 2.7243, -1.1303, 6.1958, -2.7322, 11.1286, -2.7493, 0.4496),
        ("Cl-", -131.29, -167.11, 56.74, -122.49, 17.34, 4.032, 4.801, 5.563, -2.847, -4.4, -5.714, 1.456),
        ("CO2@", -386.02, -413.84, 117.57, 243.08, 32.81, 6.2466, 7.4711, 2.8136, -3.0879, 40.0325, 8.8004, -0.02),
        ("CO3-2", -527.98, -675.31, -50.0, -289.33, -6.06, 2.8524, -3.9844, 6.4142, -2.6143, -3.3206, -17.1917, 3.3914),
        ("Fe+2", -91.5, -92.24, -105.86, -32.44, -22.64, -0.7867, -9.6969, 9.5479, -2.378, 14.786, -4.6437, 1.4382),
        ("FeO2-", -368.26, -443.82, 44.35, -234.93, 0.45, 2.3837, -1.9602, 6.5182, -2.6979, -13.3207, -14.5028, 1.4662),
        ("HCO3-", -586.94, -690.01, 98.45, -34.85, 24.21, 7.5621, 1.1505, 1.2346, -2.8266, 12.9395, -4.7579, 1.2733),
        ("HSiO3-", -1014.6, -1144.68, 20.92, -87.2, 4.53, 2.9735, -0.5181, 5.9467, -2.7575, 8.1489, -7.3123, 1.5511),
        ("K+", -282.46, -252.14, 101.04, 8.39, 9.01, 3.559, -1.473, 5.435, -2.712, 7.4, -1.791, 0.1927),
        ("Mg+2", -453.99, -465.93, -138.07, -21.66, -22.01, -0.8217, -8.599, 8.39, -2.39, 20.8, -5.892, 1.5372),
        ("Na+", -261.88, -240.28, 58.41, 38.12, -1.21, 1.839, -2.285, 3.256, -2.726, 18.18, -2.981, 0.3306),
        ("NO3-", -110.91, -206.89, 146.94, -66.8, 28.66, 7.3161, 6.7824, -4.6838, -3.0594, 7.7, -6.725, 1.0977),
        ("OH-", -157.27, -230.01, -10.71, -136.34, -4.71, 1.2527, 0.0738, 1.8423, -2.7821, 4.15, -10.346, 1.7246),
        ("SO4-2", -744.46, -909.7, 18.83, -266.09, 12.92, 8.3014, -1.9846, -6.2122, -2.697, 1.64, -17.998, 3.1463),
    ]

    prop(s, key) = Float64(rec[s][key]["values"][1])

    @testset "Table D.1: aqueous standard properties" begin
        for (s, G, H, S, Cp, V) in ((r[1], r[2], r[3], r[4], r[5], r[6]) for r in tableD1)
            @test haskey(rec, s)
            @test prop(s, "sm_gibbs_energy") / 1000 ≈ G atol = 0.01
            @test prop(s, "sm_enthalpy") / 1000 ≈ H atol = 0.01
            @test prop(s, "sm_entropy_abs") ≈ S atol = 0.01
            @test prop(s, "sm_heat_capacity_p") ≈ Cp atol = 0.01
            # J/bar in the file, cm³/mol in the printed table.
            @test 10 * prop(s, "sm_volume") ≈ V atol = 0.01
        end
    end

    @testset "Table D.1: HKF equation-of-state coefficients" begin
        # Undo the printed scaling; see the note above.
        scale = (10.0, 1.0e-2, 1.0, 1.0e-4, 1.0, 1.0e-4, 1.0e-5)
        for r in tableD1
            s = r[1]
            printed = r[7:13]
            coeffs = rec[s]["TPMethods"][1]["eos_hkf_coeffs"]["values"]
            @test String(first(keys(rec[s]["TPMethods"][1]["method"]))) == "3" ||
                haskey(rec[s]["TPMethods"][1], "eos_hkf_coeffs")
            for i in 1:7
                expected = printed[i] / scale[i]
                @test Float64(coeffs[i]) ≈ expected rtol = 1.0e-4 atol = 1.0e-6
            end
        end
    end

    # Table D.1's last block leaves HKF behind: these species are extrapolated
    # by Cp(T) integration or by a log K(T) fit, and the paper gives their heat
    # capacity polynomial instead of EOS coefficients. Footnote ** records that
    # SiO2@ borrows S° and C°p from quartz, which is why a0, a1 and a2 below are
    # quartz's.
    @testset "Table D.1: the non-HKF tail" begin
        @test prop("SiO2@", "sm_gibbs_energy") / 1000 ≈ -833.41 atol = 0.01
        @test prop("SiO2@", "sm_entropy_abs") ≈ 41.34 atol = 0.01
        @test prop("SiO2@", "sm_heat_capacity_p") ≈ 44.47 atol = 0.01
        cp = rec["SiO2@"]["TPMethods"][1]["m_heat_capacity_ft_coeffs"]["values"]
        @test Float64(cp[1]) ≈ 46.94 atol = 0.01
        @test Float64(cp[2]) ≈ 0.034 atol = 0.001
        @test Float64(cp[3]) ≈ -1.13e6 rtol = 3.0e-3

        @test prop("SiO3-2", "sm_gibbs_energy") / 1000 ≈ -938.51 atol = 0.01
        @test prop("SiO3-2", "sm_entropy_abs") ≈ -80.2 atol = 0.01
        @test prop("SiO3-2", "sm_heat_capacity_p") ≈ 119.83 atol = 0.01

        @test prop("CaSiO3@", "sm_gibbs_energy") / 1000 ≈ -1517.56 atol = 0.01
        @test prop("CaSiO3@", "sm_entropy_abs") ≈ -136.68 atol = 0.01
        @test prop("CaSiO3@", "sm_heat_capacity_p") ≈ 88.9 atol = 0.01
    end

    # Table D.2. The gas volume column really is J/bar here: 2479 J/bar is
    # 24.79 L/mol, the ideal-gas molar volume at 298.15 K and 1 bar.
    tableD2 = [
        #  symbol   ΔG⁰      ΔH⁰      S⁰       Cp⁰     V⁰(J/bar)  a0      a1        a2
        ("CH4", -50.66, -74.81, 186.26, 35.75, 2479, 23.64, 0.0479, -192464),
        ("CO2", -394.39, -393.51, 213.74, 37.15, 2479, 44.22, 0.0088, -861904),
        ("H2", 0.0, 0.0, 130.68, 28.82, 2479, 27.28, 0.0033, 50208),
        ("H2O", -228.68, -242.4, 187.25, 40.07, 2479, 52.99, -0.0435, 5472),
        ("H2S", -33.75, -20.63, 205.79, 34.2, 2479, 32.68, 0.0124, -192464),
        ("N2", 0.0, 0.0, 191.61, 29.13, 2479, 28.58, 0.0038, -50208),
        ("O2", 0.0, 0.0, 205.14, 29.32, 2479, 29.96, 0.0042, -167360),
    ]

    @testset "Table D.2: gaseous standard properties" begin
        for (s, G, H, S, Cp, V, a0, a1, a2) in tableD2
            @test haskey(rec, s)
            @test prop(s, "sm_gibbs_energy") / 1000 ≈ G atol = 0.01
            @test prop(s, "sm_enthalpy") / 1000 ≈ H atol = 0.01
            @test prop(s, "sm_entropy_abs") ≈ S atol = 0.01
            @test prop(s, "sm_heat_capacity_p") ≈ Cp atol = 0.01
            @test prop(s, "sm_volume") ≈ V atol = 1.0
            cp = rec[s]["TPMethods"][1]["m_heat_capacity_ft_coeffs"]["values"]
            @test Float64(cp[1]) ≈ a0 atol = 0.01
            @test Float64(cp[2]) ≈ a1 atol = 1.0e-4
            @test Float64(cp[3]) ≈ a2 atol = 1.0
        end
    end

    # ── What the shipped file does not carry ─────────────────────────────────
    #
    # Two rows of Table 2 cannot be checked, and the reason is worth recording
    # where it will be seen: the phases are in the printed table but not in the
    # vendored file. Asserting their absence keeps this note honest — if a later
    # database update adds them, this test fails and the note gets updated along
    # with the coverage.
    @testset "Table 2 rows the vendored file does not cover" begin
        # Nitrite-AFm: the solid is present, but NO2⁻ is not, so its dissolution
        # reaction cannot be written over the file's own primaries.
        @test haskey(sp, "mononitrite")
        @test !haskey(sp, "NO2-")

        # Fe-Friedel's salt (C4FCl2H10, log Ks0 = -28.62) is absent altogether,
        # as are amorphous and microcrystalline Fe(OH)3.
        @test !any(startswith(k, "C4FCl") for k in keys(sp))
        @test !haskey(sp, "FeOHam")
        @test !haskey(sp, "FeOHmic")
    end
end
