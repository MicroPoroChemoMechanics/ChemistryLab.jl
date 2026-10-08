# Agreement with the Cemdata18 paper's own tables.
#
# Lothenbach, Kulik, Matschei, Balonis, Baquerizo, Dilnesa, Miron & Myers,
# "Cemdata18: A chemical thermodynamic database for hydrated Portland cements
# and alkali-activated materials", Cem. Concr. Res. 115 (2019) 472-506,
# https://doi.org/10.1016/j.cemconres.2018.04.018
#
# `cemdata18-thermofun.json` is a file maintained elsewhere, obtained from its
# publisher. Nothing in the rest of the suite would notice if one Gibbs energy
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
    # package symbol are the two phases the database does not carry.
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

    # The two M-S-H end-members do not close. The file defines them by their
    # reactions, whose log K coefficients give -28.634 and -23.460 at 25 °C, the
    # temperature of Table 2, where the records state the published -28.80 and
    # -23.57: computed from their reactions, as ThermoFun computes them, they
    # disagree with Table 2 by 0.166 and 0.110, 0.95 and 0.63 kJ/mol. Every other
    # phase in the table closes to better than 0.04, so this is a property of the
    # source and not of the transcription. GEMS, run on the CemGEMS export of the
    # database, does not close them either (0.17 and 0.09). Their tabulated
    # energies, which the package used before it computed such a substance from
    # its reaction, disagree by 0.244 and 0.205.
    #
    # Pinned at the observed offset rather than hidden behind a loose tolerance,
    # so that a corrected upstream file shows up as a failing test, not silence.
    msh_offset = Dict("M075SH" => 0.1659, "M15SH" => 0.1095)

    RTln10 = R_GAS * 298.15 * log(10)
    # Each energy at 25 °C, as the package forms it from the file. A record is
    # tabulated at its own reference temperature `Tst`, which is 298.15 K for 230
    # of the 238 and 293.15 K for eight of them, two of which are the M-S-H end
    # members of this table: reading their tabulated ΔfG° as a 25 °C value, as
    # this test once did, put a 20 °C energy into a 25 °C constant and doubled
    # their offset (0.48 and 0.40). The testset after this one pins that
    # `ΔₐG⁰(Tst)` is the tabulated value for every record.
    G25(k) = sp[k].ΔₐG⁰(T = 298.15)
    nel(k, e) = Float64(get(atoms_charge(sp[k]), e, 0))

    @testset "Table 2/3: log Ks0 closes against the energies at 298.15 K, 1 bar" begin
        # 50 rows of Table 2 plus the 4 of Table 3.
        @test length(table2) == 54

        worst_ordinary = 0.0
        worst_row = ""
        ordinary = Float64[]
        for (s, logK_published, products) in table2
            @test haskey(sp, s)

            # Water from the hydrogen balance, then oxygen and charge as checks.
            nH2O = (nel(s, :H) - sum(ν * nel(k, :H) for (k, ν) in products)) / 2
            @test sum(ν * nel(k, :O) for (k, ν) in products) + nH2O ≈ nel(s, :O) atol = 1.0e-9
            @test sum(ν * nel(k, :Zz) for (k, ν) in products) ≈ nel(s, :Zz) atol = 1.0e-9

            ΔrG = sum(ν * G25(k) for (k, ν) in products) + nH2O * G25("H2O@") - G25(s)
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
        # first, and the page quotes the number it displays — 0.040 — from a run
        # rather than from the assertion. A threshold BOUNDS a disagreement; it
        # does not PIN it. If the worst row drifted to 0.047 the suite would stay
        # green while the page printed a stale 0.040, which is the one failure
        # mode a comparison page cannot afford.
        @test worst_ordinary ≈ 0.04 atol = 5.0e-4
        @test worst_row == "M8A-OH-LDH"
        # 50 of the 52 are an order of magnitude better again, and the two that
        # are not are the two layered double hydroxides whose published values
        # are quoted to one decimal.
        @test count(<(0.005), ordinary) == 50
        @test sort(ordinary)[end - 1] ≈ 0.0199 atol = 5.0e-4
    end

    # ── Each record at its own reference temperature ─────────────────────────
    #
    # `ΔₐG⁰(T)` is anchored to the tabulated ΔfG° at the record's own `Tst`, so at
    # `Tst` the two agree by construction, for every record. That is a check of
    # the code path, not of the data. Eight records are tabulated at 293.15 K
    # rather than 298.15 K: the six alkali C-S-H end members and the two M-S-H
    # end members. GEMS anchors them the same way (its standard energy of KSiOH
    # at 293.15 K is the tabulated −440800 J/mol), so their energies at 25 °C
    # differ from the tabulated numbers by the temperature step alone, about
    # `−S° × 5 K`. Comparing their tabulated ΔfG° with `ΔₐG⁰(298.15)`, as an
    # earlier version of this test did, measured that step and attributed it to
    # an inconsistent estimated entropy.
    @testset "ΔₐG⁰ at each record's Tst is the tabulated ΔfG°, except by a reaction" begin
        Tst(k) = Float64(get(rec[k], "Tst", 298.15))
        tabulated(k) = Float64(rec[k]["sm_gibbs_energy"]["values"][1])
        withG = [k for k in keys(rec) if haskey(rec[k], "sm_gibbs_energy")]
        # "238" on the page is this count.
        @test length(withG) == 238
        @test sort([k for k in withG if Tst(k) != 298.15]) == sort(
            [
                "ECSH1-KSH", "ECSH2-KSH", "ECSH1-NaSH", "ECSH2-NaSH",
                "KSiOH", "NaSiOH", "M075SH", "M15SH",
            ]
        )
        @test all(k -> Tst(k) == 293.15 || Tst(k) == 298.15, withG)
        # To 0.3 J/mol: the evaluation of the HKF equations of state at their
        # reference point rounds at that level; the other records agree to 1e-9.
        # The 17 records the file defines by a reaction follow their reaction,
        # as ThermoFun computes them, and their tabulated energies are not used.
        byreaction = [k for k in withG if haskey(rec[k], "reaction")]
        @test length(byreaction) == 17
        @test maximum(k -> abs(sp[k].ΔₐG⁰(T = Tst(k)) - tabulated(k)), setdiff(withG, byreaction)) < 0.5
    end

    # ── Whether a record's three formation properties agree ──────────────────
    #
    # ΔfG°, ΔfH° and S° are related by ΔfG° = ΔfH° − Tst (S° − Σ S°(elements)),
    # with the element entropies the file itself carries. The package forms
    # ΔₐG⁰ from ΔfG°, so this relation is never used; it is the data's own
    # consistency that it measures. The six alkali C-S-H end members close it
    # exactly at their 293.15 K, and not at 298.15 K, which is what says their
    # numbers are 20 °C numbers. The two M-S-H end members close it at neither,
    # by kilojoules: that disagreement belongs to the source.
    @testset "the triplet (ΔfG°, ΔfH°, S°) of each record" begin
        Sel = Dict(String(e["symbol"]) => Float64(e["entropy"]["values"][1]) for e in raw["elements"])
        Tst(k) = Float64(get(rec[k], "Tst", 298.15))
        prop1(k, key) = Float64(rec[k][key]["values"][1])
        function residual(k, T)
            ΣS = sum(Float64(ν) * Sel[String(e)] for (e, ν) in atoms_charge(sp[k]))
            return prop1(k, "sm_enthalpy") - T * (prop1(k, "sm_entropy_abs") - ΣS) - prop1(k, "sm_gibbs_energy")
        end
        for k in ("KSiOH", "NaSiOH", "ECSH1-KSH", "ECSH2-KSH", "ECSH1-NaSH", "ECSH2-NaSH")
            @test abs(residual(k, 293.15)) < 0.1
            @test residual(k, 298.15) > 790
        end
        # Printed on the page to the joule per mole.
        @test residual("M075SH", 293.15) ≈ -6611.3 atol = 0.5
        @test residual("M075SH", 298.15) ≈ -1793.2 atol = 0.5
        @test residual("M15SH", 293.15) ≈ -5698.2 atol = 0.5
        @test residual("M15SH", 298.15) ≈ -1726.5 atol = 0.5
        # The relation holds for a phase measured at 298.15 K.
        @test abs(residual("C3AH6", 298.15)) < 0.1
        # Over the crystalline records, the counts the page quotes.
        crystals = [
            k for k in keys(rec) if haskey(rec[k], "sm_gibbs_energy") &&
                haskey(rec[k], "sm_enthalpy") && haskey(rec[k], "sm_entropy_abs") &&
                aggregate_state(sp[k]) == AS_CRYSTAL
        ]
        r = Dict(k => abs(residual(k, Tst(k))) for k in crystals)
        @test length(crystals) == 143
        @test count(<(1), values(r)) == 78
        @test count(<(100), values(r)) == 126
        worst = first.(sort(collect(r); by = last, rev = true)[1:7])
        @test worst == ["M075SH", "M15SH", "C4AF", "C3A", "CA2", "C12A7", "CA"]
    end

    # The page states these numbers in prose; each one is asserted above, and
    # this says the page states the asserted ones.
    @testset "the validation page quotes the pinned numbers" begin
        page = read(joinpath(pkgdir(ChemistryLab), "docs", "src", "tutorials", "published_data_validation.md"), String)
        for quoted in (
                "**52 of the 54 rows close.**", "`0.040` on\n`M8A-OH-LDH`", "50 are inside\n`0.005`",
                "| **+0.166** |", "| **+0.110** |", "For 230 of the 238 records",
                "Over the 143 crystalline records", "`1 J/mol` for 78 and to `100 J/mol` for 126",
                "| −6.61 kJ/mol | −1.79 kJ/mol |", "| −5.70 kJ/mol | −1.73 kJ/mol |",
                "| `0.000` |", "| `0.011` |", "closes to `−0.02 kJ/mol`",
            )
            @test occursin(quoted, page)
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

    # ── What the database does not carry ────────────────────────────────────
    #
    # Two rows of Table 2 cannot be checked, and the reason is worth recording
    # where it will be seen: the phases are in the printed table but not in the
    # database. Asserting their absence keeps this note honest — if a later
    # release adds them, this test fails and the note gets updated along with
    # the coverage. (Amorphous and microcrystalline Fe(OH)3, absent from earlier
    # exports of the database, are in the release ChemistryLab reads, and are
    # checked with the other rows above.)
    @testset "Table 2 rows the database does not cover" begin
        # Nitrite-AFm: the solid is present, but NO2⁻ is not, so its dissolution
        # reaction cannot be written over the file's own primaries.
        @test haskey(sp, "mononitrite")
        @test !haskey(sp, "NO2-")

        # Fe-Friedel's salt (C4FCl2H10, log Ks0 = -28.62) is absent altogether.
        @test !any(startswith(k, "C4FCl") for k in keys(sp))
        @test haskey(sp, "Fe(OH)3(am)") && haskey(sp, "Fe(OH)3(mic)")
    end

    # ── The same two rows, checked where they can be ────────────────────────
    #
    # Nitrite-AFm through the NO2⁻ of slop98, whose aqueous ions share the
    # reference state of Cemdata18's: the Ca²⁺, NO3⁻ and water of the two files
    # are the same records. Their OH⁻ are not, 27 J/mol apart, and the reaction
    # takes Cemdata18's, with which the log Ks0 of Table 2 closes. Fe-Friedel's salt in the database cemdata18-chloride.json,
    # which ChemistryLab builds from Table 1 of the paper (src/databases/derived.jl).
    # Each log Ks0 is recomputed from the species at 25 °C, over the products
    # transcribed with the rest of Table 2.
    @testset "Table 2 rows checked outside the shipped file" begin
        slop = Dict(symbol(s) => s for s in build_species(datapath("slop98-inorganic-thermofun.json"); verbose = false))
        ext = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-chloride.json"); verbose = false))
        @test all(isapprox(slop[k].ΔₐG⁰(T = 298.15), G25(k); rtol = 1.0e-12) for k in ("Ca+2", "NO3-", "H2O@"))
        @test slop["OH-"].ΔₐG⁰(T = 298.15) - G25("OH-") ≈ -27.0 atol = 1.0e-6
        # Pinned, as the rows above: docs/src/tutorials/published_data_validation.md
        # quotes both offsets.
        for (phase, printed, species, offset) in (
                ("mononitrite", "Nitrite-AFm", merge(sp, Dict("NO2-" => slop["NO2-"])), 0.0),
                ("C4FCl2H10", "Fe-Friedel's salt", ext, 0.011),
            )
            products = products_of(phase)
            published = solubility.log_Ks0[only(findall(==(printed), solubility.printed))]
            n(k, e) = Float64(get(atoms_charge(species[k]), e, 0))
            nH2O = (n(phase, :H) - sum(ν * n(k, :H) for (k, ν) in products)) / 2
            @test sum(ν * n(k, :O) for (k, ν) in products) + nH2O ≈ n(phase, :O) atol = 1.0e-9
            @test sum(ν * n(k, :Zz) for (k, ν) in products) ≈ n(phase, :Zz) atol = 1.0e-9
            g(k) = species[k].ΔₐG⁰(T = 298.15)
            ΔrG = sum(ν * g(k) for (k, ν) in products) + nH2O * g("H2O@") - g(phase)
            @test -ΔrG / RTln10 - published ≈ offset atol = 1.0e-3
        end

        # The record of Fe-Friedel's salt is the row of Table 1: its energies,
        # entropy and volume at 298.15 K, and its heat capacity the printed
        # polynomial at any temperature, a T^-1/2 term included.
        r = literature_row("Lothenbach2019", "solid_standard_properties", "C4FCl2H10")
        ff = ext["C4FCl2H10"]
        @test ff[:ΔₐG⁰](T = 298.15) ≈ ustrip(us"J/mol", r.dfG) rtol = 1.0e-12
        @test ff[:ΔₐH⁰](T = 298.15) ≈ ustrip(us"J/mol", r.dfH) rtol = 1.0e-12
        @test ff[:S⁰](T = 298.15) ≈ ustrip(us"J/(mol*K)", r.S) rtol = 1.0e-12
        @test ff[:V⁰](T = 298.15, P = 1.0e5) ≈ ustrip(us"m^3/mol", r.V) rtol = 1.0e-12
        for T in (283.15, 298.15, 323.15)
            @test ff[:Cp⁰](T = T) ≈ ustrip(r.a0) + ustrip(r.a1) * T + ustrip(r.a2) / T^2 + ustrip(r.a3) / sqrt(T) rtol = 1.0e-12
        end
        # Its enthalpy is the one note l of the table says was recalculated from
        # its Gibbs energy and entropy: with the elements' entropies of the
        # Cemdata18 file, ΔfH − T ΔfS − ΔfG closes to 0.02 kJ/mol.
        S_el = Dict(Symbol(e["symbol"]) => Float64(e["entropy"]["values"][1]) for e in raw["elements"])
        ΔfS = ustrip(us"J/(mol*K)", r.S) - sum(c * S_el[el] for (el, c) in atoms_charge(ff) if el !== :Zz)
        @test ustrip(us"J/mol", r.dfH) - 298.15 * ΔfS - ustrip(us"J/mol", r.dfG) ≈ -24.4 atol = 0.5

        # Friedel's salt and Fe-Friedel's salt form the ideal solid solution of
        # note k, built from the database that has both.
        fr = only(filter(p -> name(p) == "Friedel_AlFe", build_solid_solutions(datapath("solid_solutions.toml"), ext)))
        @test model(fr) isa IdealSolidSolutionModel
        @test symbol.(end_members(fr)) == ["C4AClH10", "C4FCl2H10"]
        @test !any(p -> name(p) == "Friedel_AlFe", build_solid_solutions(datapath("solid_solutions.toml"), sp))
    end
end
