# Cemdata07's solubility products, and a measured pore solution.
#
# Lothenbach, "Thermodynamic equilibrium calculations in cementitious systems",
# Mater. Struct. 43 (2010) 1413-1433, from data/literature/Lothenbach2010.json.
# Its Table 1 is the Cemdata07 generation of the database this package ships as
# Cemdata18: 29 solubility products with their dissolution reactions, and no
# Gibbs energy or molar volume. Its Table 2 is the pore solution of an ordinary
# Portland cement paste after 69 days.
#
# Two uses. The generation gap: which constants Cemdata18 kept and which it
# revised, read from our own energies — and with it the identity of the
# constants Guo et al. (2018) used. And the claim, made on the same page as
# Table 2, that the pore solutions of old pastes are often found a little above
# saturation for portlandite and ettringite.

include("reference_species.jl")

@testsection "Lothenbach (2010) Cemdata07 and a 69-day pore solution" begin

    subs = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
    byname = Dict(symbol(s) => s for s in subs)
    G(s) = ustrip(us"J/mol", s[:ΔₐG⁰](T = 298.15u"K", P = 1.0e5u"Pa"; unit = true))
    RTln10 = R_GAS * 298.15 * log(10)

    table1 = literature_table("Lothenbach2010", "solubility_products")
    products_of(m) = let r = literature_table("Lothenbach2010", "dissolution_products"; mineral = m)
        Dict(zip(r.species, r.coefficient))
    end
    # The paper writes aluminate, ferrite and silicate hydrated; CEMDATA18's
    # species differ from those by water only, which the hydrogen balance
    # restores.
    ours = Dict("Al(OH)4-" => "AlO2-", "Fe(OH)4-" => "FeO2-", "SiO(OH)3-" => "HSiO3-", "H3SiO4-" => "HSiO3-")
    aq(sp) = byname[get(ours, sp, sp)]
    count_of(a, e) = Float64(get(a, e, 0))

    # Each mineral of Table 1 against the CEMDATA18 solid of the same
    # composition, matched on the atoms and never on the number: `k` is how
    # many of the paper's formula units the CEMDATA18 one holds. The paper
    # prints 5/3 as 1.6667, so compositions agree to four decimals. A formula
    # two polymorphs share is resolved by the name: the paper's SiO2 and
    # Al(OH)3 are the amorphous ones, not quartz or gibbsite.
    polymorph = Dict("SiO2,am" => "Amor-Sl", "Al(OH)3,am" => "AlOHam")
    solids = [s for s in subs if aggregate_state(s) == AS_CRYSTAL]
    function counterpart(mineral, formula)
        haskey(polymorph, mineral) && return (byname[polymorph[mineral]], 1.0)
        a = atoms_charge(Species(formula))
        found = []
        for s in solids
            b = atoms_charge(s)
            Set(keys(a)) == Set(keys(b)) || continue
            k = count_of(b, first(keys(a))) / count_of(a, first(keys(a)))
            all(isapprox(count_of(b, e), k * count_of(a, e); rtol = 1.0e-4) for e in keys(a)) &&
                push!(found, (s, k))
        end
        isempty(found) && return nothing
        # The whole formula unit when CEMDATA18 has it, not a half used in a
        # solid solution.
        i = findfirst(f -> f[2] ≈ 1, found)
        return i === nothing ? first(found) : found[i]
    end

    rows = []
    @testset "every reaction of Table 1 balances" begin
        for (m, logK, dec, tentative, formula) in zip(
                table1.mineral, table1.log_Ks0, table1.decimals, table1.tentative, table1.formula,
            )
            a = atoms_charge(Species(formula))
            prod = products_of(m)
            # Water from the hydrogen balance; oxygen and charge are then checks.
            nH2O = (count_of(a, :H) - sum(ν * count_of(atoms_charge(aq(s)), :H) for (s, ν) in prod)) / 2
            # 2e-4: the C-S-H rows print 5/3 and 7/3 to four decimals.
            @test sum(ν * count_of(atoms_charge(aq(s)), :O) for (s, ν) in prod) + nH2O ≈ count_of(a, :O) atol = 2.0e-4
            @test sum(ν * count_of(atoms_charge(aq(s)), :Zz) for (s, ν) in prod) ≈ 0 atol = 2.0e-4
            match = counterpart(m, formula)
            ours_logK = if match === nothing
                nothing
            else
                s, k = match
                ΔrG = sum(ν * G(aq(x)) for (x, ν) in prod) + nH2O * G(byname["H2O@"]) - G(s) / k
                -ΔrG / RTln10
            end
            push!(rows, (; m, logK, dec, tentative, match, ours = ours_logK))
        end
    end

    kept = [r for r in rows if r.ours !== nothing && abs(r.ours - r.logK) <= 0.5 * 10.0^(-r.dec) + 1.0e-3]
    moved = [r for r in rows if r.ours !== nothing && abs(r.ours - r.logK) > 0.5 * 10.0^(-r.dec) + 1.0e-3]
    absent = [r for r in rows if r.ours === nothing]
    @testset "what Cemdata18 kept, revised and dropped" begin
        # Eleven constants are Cemdata07's to the printed digit, from our own
        # energies; nine moved; nine have no CEMDATA18 solid of the same
        # composition (a different hydration state, or a phase not carried).
        @test sort([symbol(r.match[1]) for r in kept]) == sort(
            [
                "ettringite", "tricarboalu", "Fe-ettringite", "monocarbonate", "hemicarbonate",
                "straetlingite", "hydrotalcite", "Tob-II", "Amor-Sl", "syngenite", "AlOHam",
            ]
        )
        # The revisions, ours minus Cemdata07, to the digit the chapter prints.
        shift = Dict(symbol(r.match[1]) => r.ours - r.logK for r in moved)
        printed = Dict(
            "thaumasite" => -0.1, "C3AH6" => 0.34, "C3FH6" => -1.14, "C4AH13" => 0.15,
            "monosulphate12" => 0.026, "C4FH13" => -1.35, "Fe-monosulphate" => 1.63,
            "CAH10" => -0.1, "Jennite" => 0.007,
        )
        @test Set(keys(shift)) == Set(keys(printed))
        for (k, v) in printed
            @test shift[k] ≈ v atol = (k in ("monosulphate12", "Jennite") ? 5.0e-4 : 5.0e-3)
        end
        @test sort([r.m for r in absent]) == sort(
            [
                "Siliceous hydrogarnet", "C2AH8", "C2FH8", "Fe-monocarbonate", "Fe-hemicarbonate",
                "Fe-stratlingite", "M4ACH9", "M4FH10", "Fe(OH)3,mic",
            ]
        )
        # All eight tentative values but one are among the moved or the absent.
        @test count(r -> string(r.tentative) == "true", kept) == 1   # hydrotalcite, M4AH10
    end

    @testset "Guo's constants are Cemdata07's" begin
        # Guo et al. (2018) print ettringite, monosulfate and their C-S-H at
        # -44.9085, -29.2628 and -13.1659; Cemdata07 at -44.9, -29.26 and -13.17,
        # the same numbers to the digit Lothenbach prints.
        guo = literature_table("Guo2018", "dissolution")
        l10(m) = only(r.logK for r in rows if r.m == m)
        for (g, m, dec) in (
                ("AFt", "Ettringite", 1), ("AFm", "Monosulfoaluminate", 2),
                ("CSH", "Jennite-type C-S-H", 2),
            )
            @test only(k for (p, k) in zip(guo.phase, guo.log_K) if p == g) ≈ l10(m) atol = 0.5 * 10.0^(-dec)
        end
    end

    # ── Table 2: the saturation of a 69-day pore solution ─────────────────────
    #
    # The totals as analyzed (mmol/L, taken as mmol per kg of water), the
    # hydroxide from the charge balance, lithium and strontium left out
    # (CEMDATA18 carries neither; 0.7 and 0.05 mmol/L). No temperature is given
    # with the table, so the database's own, 25 °C.
    pore = literature_table("Lothenbach2010", "pore_solution")
    mM(v) = ustrip(u"mol/m^3", v)
    sp = speciation(
        subs, [byname[s] for s in ("H2O@", "Na+", "K+", "Ca+2", "AlO2-", "HSiO3-", "SO4-2")];
        aggregate_state = [AS_AQUEOUS], exclude_species = split("H2@ O2@ CH4@")
    )
    cs = ChemicalSystem(sp, CEMDATA_PRIMARIES)
    idx = Dict(symbol(s) => i for (i, s) in enumerate(cs.species))
    models = (gems = HKFActivityModel(å = 0.0, Ḃ = gems_bdot(), Kₙ = 0.0), davies = DaviesActivityModel())
    portlandite = Dict("Ca+2" => 1.0, "OH-" => 2.0)
    si(lna, s, prod, nions) = let
        a = atoms_charge(byname[s])
        nH2O = (count_of(a, :H) - sum(ν * count_of(atoms_charge(aq(x)), :H) for (x, ν) in prod)) / 2
        lnIAP = sum(ν * lna[symbol(aq(x))] for (x, ν) in prod) + nH2O * lna["H2O@"]
        lnK = -(sum(ν * G(aq(x)) for (x, ν) in prod) + nH2O * G(byname["H2O@"]) - G(byname[s])) / (R_GAS * 298.15)
        (lnIAP - lnK) / log(10) / nions
    end
    sat = []
    for i in eachindex(pore.pressure_MPa)
        c = Dict(
            "Na+" => mM(pore.Na[i]), "K+" => mM(pore.K[i]), "Ca+2" => mM(pore.Ca[i]),
            "AlO2-" => mM(pore.Al[i]), "HSiO3-" => mM(pore.Si[i]), "SO4-2" => mM(pore.S[i]),
        )
        st = ChemicalState(cs)
        for (s, v) in c
            set_quantity!(st, s, v * 1.0e-3 * u"mol")
        end
        oh = c["Na+"] + c["K+"] + 2c["Ca+2"] - c["AlO2-"] - c["HSiO3-"] - 2c["SO4-2"]
        set_quantity!(st, "OH-", oh * 1.0e-3 * u"mol")
        set_quantity!(st, "H2O@", 1.0u"kg")
        for (name, model) in pairs(models)
            eq, cert = equilibrate_certified(st; model)
            @test cert.optimal
            lna = log_activities(eq, model)
            push!(
                sat, (
                    model = name, p = pore.pressure_MPa[i],
                    free_OH = 1000 * ustrip(us"mol", eq.n[idx["OH-"]]), OH = mM(pore.OH[i]),
                    CH = si(lna, "Portlandite", portlandite, 3),
                    AFt = si(lna, "ettringite", products_of("Ettringite"), 15),
                    pH = pH(eq, model),
                ),
            )
        end
    end
    @testset "saturated or not depends on the activity model" begin
        # The page carrying Table 2 reports old pore solutions a little above
        # saturation for portlandite and ettringite. At
        # I ≈ 0.5 mol/kg that is a statement about the activity model. With
        # Davies it holds: effective indices of +0.14 to +0.16 and +0.18 to
        # +0.23. With the B-dot model identified from GEMS (å = 0) the same
        # analyses come out just UNDERsaturated, -0.10 to -0.11 and -0.15 to
        # -0.18. Either way within a quarter of a log unit of equilibrium.
        g = [x for x in sat if x.model == :gems]
        d = [x for x in sat if x.model == :davies]
        @test all(abs(a - b) < 5.0e-4 for (a, b) in zip(extrema(x.CH for x in g), (-0.113, -0.098)))
        @test all(abs(a - b) < 5.0e-4 for (a, b) in zip(extrema(x.AFt for x in g), (-0.176, -0.153)))
        @test all(abs(a - b) < 5.0e-4 for (a, b) in zip(extrema(x.CH for x in d), (0.135, 0.163)))
        @test all(abs(a - b) < 5.0e-4 for (a, b) in zip(extrema(x.AFt for x in d), (0.181, 0.227)))
        # The free hydroxide the analysis reports comes back within 6 % from the
        # charge balance of the totals, with the GEMS model.
        @test maximum(abs(x.free_OH / x.OH - 1) for x in g) < 0.06
    end
end
