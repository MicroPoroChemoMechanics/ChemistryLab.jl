# Cemdata07's standard properties, against the database this package ships.
#
# Lothenbach, Matschei, Möschner & Glasser, Cem. Concr. Res. 38 (2008) 1-18,
# Table 4: cemdata2007, the generation before
# Cemdata18, with the Gibbs energy, enthalpy, entropy and molar volume of 37
# solids and water, from data/literature/Lothenbach2008.json.
#
# The energies of the two generations differ where Cemdata18 revised them, and
# this test separates the revisions from the constants carried over. The molar
# volumes are another matter: most were calculated from unit cells, and nothing
# in a revision of the solubility products moves them. Agreement there checks
# the volumes of the solids the package's volume balances rest on, which no
# other test covers.

include("reference_species.jl")

@testsection "Lothenbach et al. (2008) cemdata2007 against CEMDATA18" begin

    subs = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
    byname = Dict(symbol(s) => s for s in subs)
    T, P = 298.15u"K", 1.0e5u"Pa"
    G(s) = ustrip(us"J/mol", s[:ΔₐG⁰](T = T, P = P; unit = true)) / 1000        # kJ/mol
    V(s) = 1.0e6 * ustrip(us"m^3/mol", s[:V⁰](T = T, P = P; unit = true))      # cm³/mol
    t = literature_table("Lothenbach2008", "thermodynamic_data")
    kJ(x) = ustrip(us"J/mol", x) / 1000
    cm3(x) = 1.0e6 * ustrip(us"m^3/mol", x)

    # Paired on composition, never on the number. A formula two phases share
    # is resolved by the name: amorphous silica is not quartz, calcite is not
    # aragonite, the amorphous Al(OH)3 is not gibbsite; ettringite is the pure
    # phase, not the copy CEMDATA18 keeps for its solid solution; and M4AH10 is
    # the hydrotalcite CEMDATA18 recommends for Portland cement, not its
    # alkali-activated M4A-OH-LDH. The two C-S-H are paired by name too, since
    # the table abbreviates them to C1.67SH2.1 and C0.83SH1.3.
    by_name = Dict(
        "SiO2,am" => "Amor-Sl", "Al(OH)3(am)" => "AlOHam", "calcite" => "Cal",
        "(Al-)ettringite" => "ettringite", "M4AH10" => "hydrotalcite",
        "C1.67SH2.1 (jennite)" => "Jennite", "C0.83SH1.3 (tobermorite)" => "Tob-II",
        "H2O" => "H2O@",
    )
    solids = [s for s in subs if aggregate_state(s) == AS_CRYSTAL]
    count_of(a, e) = Float64(get(a, e, 0))
    function counterpart(phase, formula)
        haskey(by_name, phase) && return (byname[by_name[phase]], 1.0)
        a = atoms(Species(formula))
        e1 = first(keys(a))
        found = Tuple{Any, Float64}[]
        for s in solids
            b = atoms(s)
            Set(keys(a)) == Set(keys(b)) || continue
            k = count_of(b, e1) / count_of(a, e1)
            all(isapprox(count_of(b, e), k * count_of(a, e); rtol = 1.0e-4) for e in keys(a)) &&
                push!(found, (s, k))
        end
        isempty(found) && return nothing
        # The whole formula unit where CEMDATA18 carries it; two whole units of
        # one composition would be an unresolved polymorph, and `only` says so.
        whole = [f for f in found if f[2] ≈ 1]
        return isempty(whole) ? first(found) : only(whole)
    end

    rows = [
        (; phase, match = counterpart(phase, formula), dfG = kJ(g), V = cm3(v), tentative = string(tn) == "true")
            for (phase, formula, g, v, tn) in zip(t.phase, t.formula, t.dfG, t.V, t.tentative)
    ]
    matched = [r for r in rows if r.match !== nothing]
    dG(r) = G(r.match[1]) / r.match[2] - r.dfG
    dV(r) = V(r.match[1]) / r.match[2] - r.V
    formula_of = Dict(zip(t.phase, t.formula))

    @testset "the energies: carried over, revised, or another compound" begin
        # Twenty-two are Cemdata07's to the printed 0.01 kJ/mol: 21 solids and
        # water. Seven were revised, none by more than 10 kJ/mol. Nine have no
        # CEMDATA18 solid of the same composition.
        same = [r for r in matched if abs(dG(r)) <= 0.005 + 1.0e-6]
        @test length(same) == 22
        shift = Dict(symbol(r.match[1]) => dG(r) for r in matched if abs(dG(r)) > 0.005 + 1.0e-6)
        printed = Dict(
            "C3AH6" => 1.94, "C3FH6" => -6.53, "C4AH13" => 0.87, "monosulphate12" => 0.14,
            "C4FH13" => -7.71, "Fe-monosulphate" => 9.3, "Femonocarbonate" => 5.19,
        )
        @test Set(keys(shift)) == Set(keys(printed))
        for (k, v) in printed
            @test shift[k] ≈ v atol = 5.0e-3
        end
        @test sort([r.phase for r in rows if r.match === nothing]) == sort(
            [
                "C3AS0.8H4.4", "C2AH8", "C2FH8", "C4FC̄0.5H12", "C2FSH8", "M4AcH9", "M4FH10",
                "Al2O3", "Fe(OH)3(mic)",
            ]
        )
    end

    @testset "the molar volumes" begin
        # 24 of the 28 paired solids to half the printed digit, 0.5 cm³/mol;
        # the other four, three of them iron phases, by up to 1.7 cm³/mol.
        solids_m = [r for r in matched if r.phase != "H2O"]
        @test length(solids_m) == 28
        off = Dict(symbol(r.match[1]) => dV(r) for r in solids_m if abs(dV(r)) > 0.5 + 1.0e-6)
        @test Set(keys(off)) == Set(["Fe-ettringite", "monosulphate12", "Fe-monosulphate", "Femonocarbonate"])
        @test maximum(abs, values(off)) ≈ 1.67 atol = 5.0e-3
    end

    @testset "another hydration state is priced as water" begin
        # Three rows have a CEMDATA18 counterpart at another hydration state,
        # the difference counted on the hydrogen rather than read off a name.
        # The Gibbs energy between the two is the water's, to within 3 %.
        gw = G(byname["H2O@"])
        for (p08, s18, printed) in (
                ("C2AH8", "C2AH7.5", -234.4), ("C4FC̄0.5H12", "Fe-hemicarbonate", -243.7),
                ("Fe(OH)3(mic)", "FeOOHmic", -231.5),
            )
            r = only(x for x in rows if x.phase == p08)
            dn = (count_of(atoms(Species(formula_of[p08])), :H) - count_of(atoms(byname[s18]), :H)) / 2
            @test dn > 0
            per = (r.dfG - G(byname[s18])) / dn
            @test per ≈ printed atol = 0.05
            @test per / gw ≈ 1 atol = 0.03
        end
    end

    @testset "the same generation as Lothenbach (2010)" begin
        # Table 1 of Lothenbach (2010) and this Table 4 print the same Cemdata07
        # solubility products, phase for phase, under different names.
        s08 = literature_table("Lothenbach2008", "solubility_products")
        s10 = literature_table("Lothenbach2010", "solubility_products")
        k08 = Dict(zip(s08.phase, s08.log_Ks0))
        k10 = Dict(zip(s10.mineral, s10.log_Ks0))
        same = Dict(
            "(Al-)ettringite" => "Ettringite", "Tricarboaluminate" => "Tricarboaluminate",
            "Fe-ettringite" => "Fe-ettringite", "C3AH6" => "C3AH6",
            "C3AS0.8H4.4" => "Siliceous hydrogarnet", "C3FH6" => "C3FH6", "C4AH13" => "C4AH13",
            "C2AH8" => "C2AH8", "C4AS̄H12" => "Monosulfoaluminate",
            "C4AC̄H11" => "Monocarboaluminate", "C4AC̄0.5H12" => "Hemicarboaluminate",
            "C2ASH8" => "Stratlingite", "C4FH13" => "C4FH13", "C2FH8" => "C2FH8",
            "C4FS̄H12" => "Fe-monosulfate", "C4FC̄H12" => "Fe-monocarbonate",
            "C4FC̄0.5H12" => "Fe-hemicarbonate", "C2FSH8" => "Fe-stratlingite",
            "M4AH10" => "M4AH10", "M4AcH9" => "M4ACH9", "M4FH10" => "M4FH10",
            "C1.67SH2.1 (jennite)" => "Jennite-type C-S-H",
            "C0.83SH1.3 (tobermorite)" => "Tobermorite-type C-S-H", "SiO2,am" => "SiO2,am",
            "Al(OH)3(am)" => "Al(OH)3,am", "Fe(OH)3(mic)" => "Fe(OH)3,mic",
        )
        @test length(same) == 26
        for (a, b) in same
            @test k08[a] == k10[b]
        end
    end
end
