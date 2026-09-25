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
    # Coefficients are the paper's, phase by phase, in its own order, as
    # transcribed in data/literature/Lothenbach2019.json. Its rows with no
    # package symbol are the four phases the vendored file does not carry.
    solubility = literature_table("Lothenbach2019", "solubility_products")
    products_of(phase) = let r = literature_table(
            "Lothenbach2019", "dissolution_products"; phase
        )
        Dict(zip(r.species, r.coefficient))
    end
    table2 = [
        (s, logK, products_of(s)) for (s, logK) in zip(solubility.phase, solubility.log_Ks0)
            if !isempty(s)
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

    RTln10 = R_GAS * 298.15 * log(10)
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
        worst_row = ""
        ordinary = Float64[]
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
                dev = abs(logK - logK_published)
                push!(ordinary, dev)
                dev > worst_ordinary && (worst_ordinary = dev; worst_row = s)
            end
        end

        # PINNED, NOT BOUNDED. `worst_ordinary < 0.05` is what this asserted
        # first, and the page quotes the number it displays — 0.041 — from a run
        # rather than from the assertion. A threshold BOUNDS a disagreement; it
        # does not PIN it. If the worst row drifted to 0.047 the suite would stay
        # green while the page printed a stale 0.041, which is the one failure
        # mode a comparison page cannot afford.
        @test worst_ordinary ≈ 0.0406 atol = 5.0e-4
        @test worst_row == "M8A-OH-LDH"
        # 48 of the 50 are an order of magnitude better again, and the two that
        # are not are the two layered double hydroxides whose published values
        # are quoted to one decimal.
        @test count(<(0.005), ordinary) == 48
        @test sort(ordinary)[end - 1] ≈ 0.0204 atol = 5.0e-4
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
        # "220 of the 228" on the page is this subtraction, so the 228 is pinned
        # here rather than recalled: a database update that adds a substance
        # changes the sentence, and this is what says so.
        @test count(k -> haskey(rec[k], "sm_gibbs_energy"), keys(rec)) == 228
        # The entropy the gap implies, which is the column the page prints and
        # nothing asserted. `ΔfG° = ΔfH° − T S°`, so a gap on `ΔfG°` at the
        # reference temperature is `−T` times an inconsistency in `S°`.
        # The page prints this column to two decimals, so that is the precision
        # it is pinned at — half the last displayed digit. Writing more digits
        # here than the page shows would assert something the page does not say.
        implied = Dict(
            "ECSH1-KSH" => 0.82, "ECSH2-KSH" => 0.82,
            "ECSH1-NaSH" => 0.7, "ECSH2-NaSH" => 0.7,
            "KSiOH" => 0.7, "NaSiOH" => 0.7,
            "M075SH" => 4.58, "M15SH" => 3.65,
        )
        for (k, gap) in rebuilt
            @test -gap / 298.15 ≈ implied[k] atol = 5.0e-3
        end

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
    # Read back in the printed units: kJ/mol, and cm³/mol for V⁰.
    d1 = literature_table("Lothenbach2019", "aqueous_species")
    tableD1 = [
        (
            d1.species[i], ustrip(u"kJ/mol", d1.dfG[i]), ustrip(u"kJ/mol", d1.dfH[i]),
            ustrip(u"J/(mol*K)", d1.S[i]), ustrip(u"J/(mol*K)", d1.Cp[i]),
            ustrip(u"cm^3/mol", d1.V[i]), d1.a1_e1[i], d1.a2_em2[i], d1.a3[i],
            d1.a4_em4[i], d1.c1[i], d1.c2_em4[i], d1.omega_em5[i],
        ) for i in eachindex(d1.species)
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
        tail = literature_table("Lothenbach2019", "aqueous_species_cp")
        for (k, sym) in enumerate(tail.species)
            @test prop(sym, "sm_gibbs_energy") / 1000 ≈ ustrip(u"kJ/mol", tail.dfG[k]) atol = 0.01
            @test prop(sym, "sm_entropy_abs") ≈ ustrip(u"J/(mol*K)", tail.S[k]) atol = 0.01
            @test prop(sym, "sm_heat_capacity_p") ≈ ustrip(u"J/(mol*K)", tail.Cp[k]) atol = 0.01
        end
        quartz = literature_row("Lothenbach2019", "aqueous_cp_polynomial", "SiO2@")
        cp = rec["SiO2@"]["TPMethods"][1]["m_heat_capacity_ft_coeffs"]["values"]
        @test Float64(cp[1]) ≈ ustrip(u"J/(mol*K)", quartz.a0) atol = 0.01
        @test Float64(cp[2]) ≈ ustrip(u"J/(mol*K^2)", quartz.a1) atol = 0.001
        @test Float64(cp[3]) ≈ ustrip(u"J*K/mol", quartz.a2) rtol = 3.0e-3
    end

    # Table D.2. The gas volume column really is J/bar here: 2479 J/bar is
    # 24.79 L/mol, the ideal-gas molar volume at 298.15 K and 1 bar.
    d2 = literature_table("Lothenbach2019", "gases")
    tableD2 = [
        (
            d2.species[i], ustrip(u"kJ/mol", d2.dfG[i]), ustrip(u"kJ/mol", d2.dfH[i]),
            ustrip(u"J/(mol*K)", d2.S[i]), ustrip(u"J/(mol*K)", d2.Cp[i]),
            ustrip(u"J/bar", d2.V[i]), ustrip(u"J/(mol*K)", d2.a0[i]),
            ustrip(u"J/(mol*K^2)", d2.a1[i]), ustrip(u"J*K/mol", d2.a2[i]),
        ) for i in eachindex(d2.species)
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
