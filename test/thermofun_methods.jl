using JSON
using Logging

# The methods of a ThermoFun database beyond the heat-capacity polynomial and the
# HKF equations, against ThermoFun itself (`test/reference/thermofun_methods.py`):
# the substances defined by a reaction, the reactions with their entropy, the
# entropy obtained by integrating the heat capacity, and the methods of aq17: the
# volume and the Landau transition of Holland and Powell's minerals, Akinfiev and
# Diamond's dissolved gases, and the fluids of Holland and Powell, whose departure
# from the ideal gas is the gas phase's.
#
# Two conventions separate the values, and bound the agreement below. ThermoFun
# takes the gas constant as 8.31451 J/(mol K), the package as the exact
# 8.314462618: a Gibbs energy of reaction of 240 kJ/mol differs by 1.4 J/mol from
# that alone. ThermoFun takes the solvent from its equation of state, -237182.28
# J/mol and -285832 J/mol at 25 °C, where the package anchors it on the record,
# -237183 and -285881: a substance whose reaction involves n waters carries n
# times these differences, 0.72 J/mol in its Gibbs energy (6 waters for
# chabazite) and 49 J/mol in its enthalpy, which is therefore compared as a change
# from the reference state.

@testsection "ThermoFun methods against ThermoFun" begin
    oracle = JSON.parsefile(joinpath(@__DIR__, "reference", "thermofun_methods.json"))["databases"]
    at(f, row) = f(T = row["T"], P = row["P"])
    at(f, row::NamedTuple) = f(T = row.T, P = row.P)

    for (file, substances) in sort!(collect(oracle); by = first)
        rows = substances["substances"]
        _, subs, reacs = read_thermofun_database(datapath(file))
        methods(s) = [only(values(m.method)) for m in coalesce(only(subs[subs.symbol .== s, :TPMethods]), [])]
        @testset "$file" begin
            sp = Dict(symbol(s) => s for s in build_species(subs, collect(keys(rows))))
            for (s, points) in rows
                ref = first(points)
                @test (ref["T"], ref["P"]) == (298.15, 1.0e5)
                @test !haskey(ChemistryLab.properties(sp[s]), :refused_method)
                record = only(eachrow(subs[subs.symbol .== s, :]))
                if "mv_eos_murnaghan_hp98" in methods(s) && !ismissing(record.m_compressibility)
                    # Holland and Powell's minerals. ThermoFun integrates their
                    # volume from zero pressure and leaves the integral out at the
                    # reference state: away from it, its energies exceed these by
                    # the integral up to the reference pressure P°, P° V and
                    # P° (V − T ∂V/∂T) at P° (10 J/mol for albite).
                    Pr = 1.0e5
                    landau = only(m for m in record.TPMethods if only(values(m.method)) == "landau_holland_powell98")
                    Tc, _, Vmax = Float64.(landau.m_landau_phase_trans_props.values)
                    Q0² = sqrt(1 - 298.15 / (Tc + T_ZERO_CELSIUS))
                    for row in points
                        T = row["T"]
                        V(t) = sp[s][:V⁰](T = t, P = Pr)
                        away = (T, row["P"]) != (298.15, Pr)
                        oG = away ? Pr * V(T) : 0.0
                        oH = away ? Pr * (V(T) - T * ForwardDiff.derivative(V, T)) : 0.0
                        @test at(sp[s][:ΔₐG⁰], row) + oG ≈ row["G"] atol = 0.05
                        @test at(sp[s][:ΔₐH⁰], row) + oH ≈ row["H"] atol = 0.05
                        @test at(sp[s][:S⁰], row) ≈ row["S"] atol = 1.0e-3
                        @test at(sp[s][:Cp⁰], row) ≈ row["Cp"] atol = 1.0e-3
                        if away
                            # ThermoFun reports a volume without part of the
                            # Landau excess its Gibbs energy carries, 0.2 % of
                            # albite's; the Gibbs energies agree at 1 and 100 bar.
                            @test at(sp[s][:V⁰], row) ≈ 1.0e-5 * row["V"] rtol = 3.0e-3
                        else
                            # At the reference state ThermoFun reports the
                            # record's volume, where the equations of the article
                            # add the excess volume of disorder V_max Q₀⁶/3.
                            @test at(sp[s][:V⁰], row) ≈ 1.0e-5 * (row["V"] + Vmax * Q0²^3 / 3) rtol = 1.0e-5
                        end
                    end
                    continue
                elseif "solute_aknifiev_diamond03" in methods(s) || s == "H2O@"
                    # Akinfiev and Diamond's dissolved gases, and the water they
                    # are computed with. ThermoFun takes the water of IAPWS-95,
                    # the package that of Haar, Gallagher and Kell anchored on the
                    # record: 1.3 J/mol apart in G and 50 J/mol in H at 25 °C,
                    # their densities by up to 7e-5 here, on a term
                    # RT ρ (a + b √(10³/T)) of some 29 kJ/mol for CO2. The gases
                    # agree to the same order: 8 J/mol in G up to 150 °C, against
                    # the 12 kJ/mol the equation adds to the ideal gas there.
                    gas = s != "H2O@"
                    for row in points
                        @test at(sp[s][:ΔₐG⁰], row) ≈ row["G"] atol = gas ? 10.0 : 1.5
                        @test at(sp[s][:ΔₐH⁰], row) - at(sp[s][:ΔₐH⁰], ref) ≈ row["H"] - ref["H"] atol = gas ? 40.0 : 5.0
                        @test at(sp[s][:S⁰], row) ≈ row["S"] atol = gas ? 0.15 : 0.01
                        @test at(sp[s][:Cp⁰], row) ≈ row["Cp"] atol = gas ? 1.0 : 0.1
                        @test at(sp[s][:V⁰], row) ≈ 1.0e-5 * row["V"] rtol = gas ? 2.0e-3 : 1.0e-4
                    end
                    continue
                end
                for row in points
                    @test at(sp[s][:ΔₐG⁰], row) ≈ row["G"] atol = 6.0
                    # A fluid of Holland and Powell carries in its ThermoFun
                    # enthalpy, entropy and volume the departure from the ideal
                    # gas, which here is the gas phase's.
                    aggregate_state(sp[s]) == AS_GAS && continue
                    @test at(sp[s][:ΔₐH⁰], row) - at(sp[s][:ΔₐH⁰], ref) ≈ row["H"] - ref["H"] atol = 2.0
                    @test at(sp[s][:S⁰], row) ≈ row["S"] atol = 1.0e-2
                    @test at(sp[s][:Cp⁰], row) ≈ row["Cp"] atol = 1.0e-2
                    # ThermoFun reports volumes in J/bar; the volume of water
                    # differs in the same way as its energies. A solid keeps the
                    # volume of its record, where ThermoFun gives that of its
                    # reaction.
                    aggregate_state(sp[s]) == AS_AQUEOUS || continue
                    @test at(sp[s][:V⁰], row) ≈ 1.0e-5 * row["V"] atol = 1.0e-9 rtol = 1.0e-5
                end
            end
            # The reactions: their log K, which ThermoFun takes from the
            # coefficients of `logk_fpt_function` (two records of Cemdata18,
            # M075SH and M15SH, state a value at 25 °C that their coefficients
            # do not give), at 1 bar, where it carries no volume of reaction;
            # their Gibbs energy from those of their species at every pressure;
            # and their entropy as the record states it (`drsm_entropy`).
            isempty(substances["reactions"]) && continue
            everything = build_species(subs)
            for (r, points) in substances["reactions"]
                # Cemdata18 gives some reactions twice, the second copy listing
                # each reactant twice: read by symbol, the two are one reaction.
                rxns = build_reactions(reacs, everything, [r])
                @test allequal(Dict(symbol(k) => v for (k, v) in x) for x in rxns)
                rxn = first(rxns)
                for row in points
                    row["P"] == 1.0e5 && @test at(rxn.logKr, row) ≈ row["logK"] rtol = 1.0e-6
                    @test at(rxn.ΔᵣG⁰, row) ≈ row["G"] atol = 6.0
                end
                record = first(eachrow(reacs[reacs.symbol .== r, :]))
                @test ustrip(rxn.ΔᵣS⁰_Tref) == only(record.drsm_entropy["values"])
            end
        end
    end
    # The substances PSI/Nagra defines by a reaction say which one.
    _, subs, _ = read_thermofun_database(datapath("psinagra-12-07-thermofun.json"))
    s = only(build_species(subs, ["CaSiO3@"]))
    @test s[:defining_reaction] == "CaSiO3@"
    # Every method of aq17 is computed: nothing is restricted to its reference
    # state.
    # Read through `import_database`, which takes the reader from the extension.
    _, aq17, _ = import_database(datapath("aq17-thermofun.json"))
    built = @test_logs min_level = Logging.Warn build_species(aq17)
    @test any(x -> symbol(x) == "Calcite", built)
    @test any(x -> symbol(x) == "H2O@", built)
    # A dissolved gas of Akinfiev and Diamond without the water of its database
    # is known at its reference state, and nowhere else.
    nowater = aq17[aq17.symbol .!= "H2O@", :]
    co2 = @test_logs (:warn, r"reference state only.*CO2@ \(`solute_aknifiev_diamond03`\)") match_mode = :any only(build_species(nowater, ["CO2@"]))
    @test co2[:refused_method] == "solute_aknifiev_diamond03"
    @test co2[:ΔₐG⁰](T = 298.15, P = 1.0e5) ≈ -386030 atol = 1.0e-6
    err = try
        co2[:ΔₐG⁰](T = 333.15, P = 1.0e5)
    catch e
        e
    end
    @test err isa ArgumentError && occursin("solute_aknifiev_diamond03", sprint(showerror, err))
end

@testsection "ThermoFun methods: log K forms and reaction records" begin
    # The seven terms of `logk_fpt_function`, the last one A₆/√T, and their
    # derivatives against the closed form.
    rec = Dict{String, Any}(
        "symbol" => "X",
        "TPMethods" => [Dict{String, Any}("logk_ft_coeffs" => Dict{String, Any}("values" => [1.0, 2.0e-3, -300.0, 0.5, 1.0e4, 1.0e-6, 3.0, 0.0]))],
    )
    f = ChemistryLab._thermofun_log10K(rec)
    L(T) = 1.0 + 2.0e-3 * T - 300.0 / T + 0.5 * log(T) + 1.0e4 / T^2 + 1.0e-6 * T^2 + 3.0 / sqrt(T)
    for T in (280.0, 350.0)
        @test f(T)[1] ≈ L(T)
        @test f(T)[2] ≈ ForwardDiff.derivative(L, T)
        @test f(T)[3] ≈ ForwardDiff.derivative(t -> ForwardDiff.derivative(L, t), T)
    end
    rec["TPMethods"][1]["logk_ft_coeffs"]["values"][8] = 1.0
    @test_throws ArgumentError ChemistryLab._thermofun_log10K(rec)
    # Without coefficients: the value at the reference, van 't Hoff's equation
    # with a constant heat capacity of reaction.
    bare = Dict{String, Any}(
        "symbol" => "Y", "TPMethods" => Any[],
        "logKr" => Dict("values" => [2.0]), "drsm_enthalpy" => Dict("values" => [-10000.0]),
        "drsm_heat_capacity_p" => Dict("values" => [50.0]),
    )
    g = ChemistryLab._thermofun_log10K(bare)
    k = ChemistryLab.R_GAS * log(10)
    M(T) = 2.0 + 10000.0 / k * (1 / T - 1 / 298.15) + 50.0 / k * (298.15 / T - 1 + log(T / 298.15))
    @test g(298.15)[1] ≈ 2.0
    @test g(340.0)[1] ≈ M(340.0)
    @test g(340.0)[2] ≈ ForwardDiff.derivative(M, 340.0)
    @test g(340.0)[3] ≈ ForwardDiff.derivative(t -> ForwardDiff.derivative(M, t), 340.0)
    # A reaction that lists a species twice counts it once; two coefficients refuse it.
    twice = Dict{String, Any}("symbol" => "Z", "reactants" => [Dict("symbol" => "A", "coefficient" => -1), Dict("symbol" => "A", "coefficient" => -1), Dict("symbol" => "Z", "coefficient" => 1)])
    @test ChemistryLab._thermofun_reactants(twice) == OrderedDict("A" => -1.0, "Z" => 1.0)
    twice["reactants"][2]["coefficient"] = -2
    @test_throws ArgumentError ChemistryLab._thermofun_reactants(twice)
end

@testsection "ThermoFun methods: refusals and the reactants of a reaction" begin
    aq17 = JSON.parsefile(datapath("aq17-thermofun.json"); dicttype = Dict{String, Any})
    cem = JSON.parsefile(datapath("cemdata18-thermofun.json"); dicttype = Dict{String, Any})
    record(db, s) = deepcopy(only(r for r in db["substances"] if r["symbol"] == s))
    # A method ThermoFun does not list, and a Landau transition without the bulk
    # modulus its volume needs: both refused by name, the substance known at its
    # reference state only.
    odd = record(aq17, "Gibbsite")
    odd["symbol"] = "Odd"
    push!(odd["TPMethods"], Dict{String, Any}("method" => Dict{String, Any}("99" => "an_unknown_method")))
    nobulk = record(aq17, "Quartz")
    nobulk["symbol"] = "NoBulk"
    delete!(nobulk, "m_compressibility")
    # Reactants found by symbol, by the symbol with `_` read as `.`, and by
    # formula; one found nowhere.
    lime = record(cem, "Lim")
    lime["symbol"] = "Lim.e"
    reaction(sym, reactants) = Dict{String, Any}(
        "symbol" => sym, "equation" => "", "Tst" => 298.15, "Pst" => 100000,
        "logKr" => Dict("values" => [-22.8]), "drsm_entropy" => Dict("values" => [0.0]),
        "drsm_enthalpy" => Dict("values" => [0.0]), "drsm_heat_capacity_p" => Dict("values" => [0.0]),
        "TPMethods" => Any[], "reactants" => [Dict{String, Any}("symbol" => s, "coefficient" => c) for (s, c) in reactants],
    )
    db = Dict{String, Any}(
        "elements" => cem["elements"],
        "substances" => [odd, nobulk, lime, record(cem, "Portlandite"), record(cem, "H2O@")],
        "reactions" => [
            reaction("slaking", [("Ca(OH)2", -1), ("Lim_e", 1), ("H2O@", 1)]),
            reaction("nowhere", [("Portlandite", -1), ("Nope", 1)]),
        ],
    )
    file = joinpath(mktempdir(), "tiny-thermofun.json")
    write(file, JSON.json(db))
    _, subs, reacs = read_thermofun_database(file)
    built = @test_logs (:warn, r"reference state only.*Odd \(`an_unknown_method`\).*NoBulk \(`landau_holland_powell98`\)") match_mode = :any build_species(subs)
    sp = Dict(symbol(s) => s for s in built)
    @test sp["Odd"][:refused_method] == "an_unknown_method"
    @test sp["NoBulk"][:refused_method] == "landau_holland_powell98"
    @test_throws ArgumentError sp["Odd"][:ΔₐG⁰](T = 333.15, P = 1.0e5)
    rxn = only(build_reactions(reacs, built, ["slaking"]))
    @test Set(symbol(k) for (k, _) in rxn) == Set(["Portlandite", "Lim.e", "H2O@"])
    @test_throws "a reaction names Nope" build_reactions(reacs, built, ["nowhere"])
end
