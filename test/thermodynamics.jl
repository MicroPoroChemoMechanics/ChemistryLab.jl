using JSON

@testsection "Thermodynamics" begin
    @testsection "THERMO_MODELS registry" begin
        # Known models must be present
        @test haskey(THERMO_MODELS, :cp_ft_equation)
        @test haskey(THERMO_MODELS, :logk_fpt_function)

        # cp_ft_equation must have Cp, H, S, G expressions
        model = THERMO_MODELS[:cp_ft_equation]
        @test haskey(model, :Cp)
        @test haskey(model, :H)
        @test haskey(model, :S)
        @test haskey(model, :G)
        @test haskey(model, :units)

        # logk_fpt_function must have logKr
        model_logk = THERMO_MODELS[:logk_fpt_function]
        @test haskey(model_logk, :logKr)
    end

    @testsection "THERMO_FACTORIES populated at __init__" begin
        @test !isempty(THERMO_FACTORIES)
        @test haskey(THERMO_FACTORIES, :cp_ft_equation)
        @test haskey(THERMO_FACTORIES, :logk_fpt_function)

        factories = THERMO_FACTORIES[:cp_ft_equation]
        @test haskey(factories, :Cp)
        @test haskey(factories, :H)
        @test haskey(factories, :S)
        @test factories[:Cp] isa ThermoFactory
    end

    @testsection "SymbolicFunc from constant Number" begin
        tf = SymbolicFunc(42.0)
        @test tf isa SymbolicFunc
        @test tf() ≈ 42.0
    end

    @testsection "SymbolicFunc from constant Quantity" begin
        tf = SymbolicFunc(100.0u"J/mol")
        @test tf isa SymbolicFunc
        val = tf(; unit = true)
        @test val isa AbstractQuantity
    end

    @testsection "SymbolicFunc from expression" begin
        # Simple linear function: f(T) = 2*T
        tf = SymbolicFunc(:(2 * T), [:T])
        @test tf isa SymbolicFunc
        @test tf(; T = 3.0) ≈ 6.0
        @test tf(; T = 10.0) ≈ 20.0
    end

    @testsection "SymbolicFunc arithmetic" begin
        tf1 = SymbolicFunc(:(2 * T), [:T])
        tf2 = SymbolicFunc(:(3 * T), [:T])

        # Addition
        tf_sum = tf1 + tf2
        @test tf_sum(; T = 1.0) ≈ 5.0

        # Subtraction
        tf_diff = tf2 - tf1
        @test tf_diff(; T = 1.0) ≈ 1.0

        # Scalar multiplication
        tf_scaled = tf1 * 2.0
        @test tf_scaled(; T = 1.0) ≈ 4.0

        # Unary negation
        tf_neg = -tf1
        @test tf_neg(; T = 1.0) ≈ -2.0
    end

    @testsection "ThermoFactory construction and call" begin
        factory = ThermoFactory(:(a * T + b), [:T])
        @test factory isa ThermoFactory
        @test haskey(factory.params, :a)
        @test haskey(factory.params, :b)
        @test haskey(factory.vars, :T)

        # Instantiate with parameters
        tf = factory(; a = 2.0, b = 5.0)
        @test tf isa SymbolicFunc
        @test tf(; T = 10.0) ≈ 25.0   # 2*10 + 5
    end

    @testsection "ThermoFactory called from several tasks at once" begin
        # Independent calculations build their species on threads, and every
        # species compiles its functions through the same global factories: the
        # memo of each factory then receives new keys from several threads at
        # once. Unguarded, those concurrent insertions corrupted the `Dict` and
        # crashed a documentation build. On one thread this checks that the
        # calls neither deadlock nor disagree with the serial ones; on several it
        # is the race itself.
        factory = ThermoFactory(:(a * T^2 + b * T + c), [:T])
        keys_ = [(a = 1.0e-3 * i, b = 0.5 * i, c = -2.0 * i) for i in 1:48]
        concurrent = fetch.(
            [Threads.@spawn [factory(; k...)(; T = 300.0) for k in keys_[j:6:end]] for j in 1:6]
        )
        serial = [[k.a * 300.0^2 + k.b * 300.0 + k.c for k in keys_[j:6:end]] for j in 1:6]
        @test concurrent ≈ serial
        @test length(factory.cache) == length(keys_)
    end

    if Base.get_extension(ChemistryLab, :SymbolicNumericIntegrationExt) !== nothing
        @testsection "add_thermo_model" begin
            model_name = :test_linear_model
            Cpexpr = :(c₀ + c₁ * T)

            # Register a new model from a Cp expression
            @test_nowarn add_thermo_model(model_name, Cpexpr)

            @test haskey(THERMO_MODELS, model_name)
            @test haskey(THERMO_FACTORIES, model_name)

            factories = THERMO_FACTORIES[model_name]
            @test haskey(factories, :Cp)
            @test haskey(factories, :H)   # integrated automatically

            # Clean up: remove the test model so it doesn't pollute other tests
            delete!(THERMO_MODELS, model_name)
            delete!(THERMO_FACTORIES, model_name)
        end
    end

    @testsection "build_thermo_functions" begin
        # Build thermodynamic functions from cp_ft_equation with minimal
        # parameters: a constant heat capacity and the reference properties of
        # liquid water, read from the shipped SLOP98 rather than recalled.
        water = only(
            s for s in JSON.parsefile(datapath("slop98-inorganic-thermofun.json"))["substances"]
                if s["symbol"] == "H2O@"
        )
        w(key) = float(only(water[key]["values"]))
        params = [
            :a₀ => w("sm_heat_capacity_p") * u"J/(mol*K)",
            :a₁ => 0.0u"J/(mol*K^2)",
            :a₂ => 0.0u"J*K/mol",
            :a₃ => 0.0u"J/(mol*K^0.5)",
            :a₄ => 0.0u"J/(mol*K^3)",
            :a₅ => 0.0u"J/(mol*K^4)",
            :a₆ => 0.0u"J/(mol*K^5)",
            :a₇ => 0.0u"J*K^2/mol",
            :a₈ => 0.0u"J/mol",
            :a₉ => 0.0u"J/(mol*K^1.5)",
            :a₁₀ => 0.0u"J/(mol*K)",
            :T => 298.15u"K",
            :S⁰ => w("sm_entropy_abs") * u"J/(mol*K)",
            :ΔfH⁰ => w("sm_enthalpy") * u"J/mol",
            :ΔfG⁰ => w("sm_gibbs_energy") * u"J/mol",
        ]
        thermo = build_thermo_functions(:cp_ft_equation, params)

        @test haskey(thermo, :Cp⁰)
        @test haskey(thermo, :ΔₐH⁰)
        @test haskey(thermo, :S⁰)
        @test haskey(thermo, :ΔₐG⁰)

        @test thermo[:Cp⁰] isa SymbolicFunc

        # Cp at Tref is the constant term a₀ alone
        cp_val = thermo[:Cp⁰](; T = 298.15)
        @test isapprox(cp_val, w("sm_heat_capacity_p"); rtol = 1.0e-6)
    end

    @testsection "entries without a heat-capacity model" begin
        # Made-up reference values: only the way they are extrapolated matters.
        Tr, T = 298.15, 350.0
        S0, H0, G0, cp0 = 100.0, -1.0e5, -1.3e5, 50.0
        ref = [:ΔₐH⁰ => H0 * u"J/mol", :ΔₐG⁰ => G0 * u"J/mol", :T => Tr * u"K", :P => 1.0e5u"Pa"]
        function entry(params)
            s = Species("CaO"; aggregate_state = AS_CRYSTAL)
            s[:thermo_params] = params
            return s
        end

        # A heat capacity given at Tref only is held constant.
        const_cp = entry([:Cp⁰ => cp0 * u"J/(mol*K)"; :S⁰ => S0 * u"J/(mol*K)"; ref])
        @test const_cp[:ΔₐH⁰](T = T, unit = false) ≈ H0 + cp0 * (T - Tr)
        @test const_cp[:ΔₐG⁰](T = T, unit = false) ≈ G0 - S0 * (T - Tr) + cp0 * (T - Tr - T * log(T / Tr))

        # No heat capacity but an entropy: a zero heat capacity, so that ΔₐG⁰
        # still falls as -S⁰ instead of staying at its tabulated value.
        no_cp = entry([:Cp⁰ => missing; :S⁰ => S0 * u"J/(mol*K)"; ref])
        @test no_cp[:Cp⁰](T = T, unit = false) == 0
        @test no_cp[:ΔₐH⁰](T = T, unit = false) ≈ H0
        @test no_cp[:ΔₐG⁰](T = T, unit = false) ≈ G0 - S0 * (T - Tr)

        # Without an entropy there is nothing to extrapolate with.
        bare = entry([:Cp⁰ => missing; :S⁰ => missing; ref])
        @test bare[:ΔₐG⁰](T = T, unit = false) ≈ G0
    end

    @testsection "a heat capacity given on several intervals" begin
        # Quartz lists one polynomial below its alpha-beta transition and one
        # above. The functions are anchored at Tref, on the interval containing it.
        qtz = only(
            s for s in JSON.parsefile(datapath("cemdata18-thermofun.json"))["substances"]
                if s["symbol"] == "Qtz"
        )
        coeffs(m) = float.(m["m_heat_capacity_ft_coeffs"]["values"])
        cps = [m for m in qtz["TPMethods"] if haskey(m, "m_heat_capacity_ft_coeffs")]
        @test length(cps) == 2
        low = only(m for m in cps if m["limitsTP"]["lowerT"] <= 298.15 <= m["limitsTP"]["upperT"])
        # the eleven terms of `:cp_ft_equation`
        powers(T) = [1, T, T^-2, T^-0.5, T^2, T^3, T^4, T^-3, T^-1, sqrt(T), log(T)]
        species = only(build_species(datapath("cemdata18-thermofun.json"), ["Qtz"]; verbose = false))
        T = 400.0
        @test species[:Cp⁰](T = T, unit = false) ≈ sum(coeffs(low) .* powers(T))
        @test all(!(species[:Cp⁰](T = T, unit = false) ≈ sum(coeffs(m) .* powers(T))) for m in cps if m !== low)

        # The interval is chosen by its limits, not by its rank in the list.
        methods = JSON.parse(
            """[{"method": {"0": "cp_ft_equation"}, "limitsTP": {"lowerT": 500.0, "upperT": 900.0},
             "m_heat_capacity_ft_coeffs": {"values": [2.0]}},
            {"method": {"0": "cp_ft_equation"}, "limitsTP": {"lowerT": 273.15, "upperT": 500.0},
             "m_heat_capacity_ft_coeffs": {"values": [1.0]}},
            {"method": {"0": "mv_constant"}}]"""
        )
        @test ChemistryLab._reference_cp_interval(methods, 298.15) === methods[2]
        @test ChemistryLab._reference_cp_interval(methods, 1000.0) === methods[1]
        @test ChemistryLab._reference_cp_interval(methods[3:3], 298.15) === nothing
    end

    @testsection "ForwardDiff — SymbolicFunc AD" begin
        using ForwardDiff

        sf = SymbolicFunc(:(a + b * T + c / T^2); a = 30.0, b = 5.0e-3, c = -1.5e5)

        # ∂f/∂T must be finite and equal to b - 2c/T³
        df_dT = ForwardDiff.derivative(T -> sf(; T = T), 298.15)
        expected = 5.0e-3 - 2 * (-1.5e5) / 298.15^3
        @test isapprox(df_dT, expected; rtol = 1.0e-6)
    end

    @testsection "ForwardDiff — SymbolicFunc arithmetic AD" begin
        using ForwardDiff

        f1 = SymbolicFunc(:(a₀ + a₁ * T); a₀ = 30.0, a₁ = 2.0e-3)
        f2 = SymbolicFunc(:(b₀ + b₁ * T); b₀ = 10.0, b₁ = 1.0e-3)

        f_sum = f1 + f2
        f_diff = f1 - f2
        f_scal = 2.5 * f1

        # derivatives of sum/diff/scalar-mult
        dsum = ForwardDiff.derivative(T -> f_sum(; T = T), 300.0)
        @test isapprox(dsum, 2.0e-3 + 1.0e-3; rtol = 1.0e-6)

        ddiff = ForwardDiff.derivative(T -> f_diff(; T = T), 300.0)
        @test isapprox(ddiff, 2.0e-3 - 1.0e-3; rtol = 1.0e-6)

        dscal = ForwardDiff.derivative(T -> f_scal(; T = T), 300.0)
        @test isapprox(dscal, 2.5 * 2.0e-3; rtol = 1.0e-6)
    end

    @testsection "ForwardDiff — cross-type (SymbolicFunc ± NumericFunc) AD" begin
        using ForwardDiff

        sf = SymbolicFunc(:(a + b * T); a = 100.0, b = 0.5)
        nf = NumericFunc(
            (T, P) -> 50.0 + 0.3 * T,
            (:T, :P),
            (T = 298.15u"K", P = 1.0e5u"Pa"),
            1.0u"1",   # dimensionless — must match sf
        )

        # Cross-type sum → NumericFunc
        cross_sum = sf + nf
        @test cross_sum isa NumericFunc

        d_cross = ForwardDiff.derivative(T -> cross_sum(; T = T, P = 1.0e5), 300.0)
        @test isapprox(d_cross, 0.5 + 0.3; rtol = 1.0e-6)

        # Cross-type diff → NumericFunc
        cross_diff = sf - nf
        d_diff = ForwardDiff.derivative(T -> cross_diff(; T = T, P = 1.0e5), 300.0)
        @test isapprox(d_diff, 0.5 - 0.3; rtol = 1.0e-6)

        # Division (e.g. G / (R·T) pattern used for logK)
        cross_div = sf / nf
        d_div = ForwardDiff.derivative(T -> cross_div(; T = T, P = 1.0e5), 300.0)
        @test isfinite(d_div)
    end
end

@testsection "a heat capacity followed across its intervals" begin
    # Hematite in slop98: three heat-capacity intervals, a transition at the top
    # of the first (950 K) with its enthalpy and entropy, none at the top of the
    # second. Only the first interval used to be read.
    raw = JSON.parsefile(datapath("slop98-inorganic-thermofun.json"); dicttype = Dict{String, Any})
    rec = only(r for r in raw["substances"] if r["symbol"] == "Hem")
    hem = only(s for s in build_species(datapath("slop98-inorganic-thermofun.json"); verbose = false) if symbol(s) == "Hem")
    G(T) = hem[:ΔₐG⁰](T = T)
    H(T) = hem[:ΔₐH⁰](T = T)
    S(T) = hem[:S⁰](T = T)
    Cp(T) = hem[:Cp⁰](T = T)
    tr = only(m["m_phase_trans_props"] for m in rec["TPMethods"] if haskey(m, "m_phase_trans_props") && m["m_phase_trans_props"]["values"][1] < 951)
    Tt, dS, dH = tr["values"][1], tr["values"][2], tr["values"][3]
    @test dH > 0 && dS > 0
    # The tabulated value at the reference temperature.
    @test G(Float64(rec["Tst"])) ≈ Float64(rec["sm_gibbs_energy"]["values"][1]) atol = 1.0e-6
    # In each interval, the functions are those of one heat capacity.
    for T in (500.0, 1000.0, 1300.0)
        @test -ForwardDiff.derivative(G, T) ≈ S(T) rtol = 1.0e-12
        @test ForwardDiff.derivative(H, T) ≈ Cp(T) rtol = 1.0e-12
    end
    # G − H + T S does not depend on T, across the intervals and the transition.
    c(T) = G(T) - H(T) + T * S(T)
    @test c(1000.0) ≈ c(500.0) rtol = 1.0e-12
    @test c(1300.0) ≈ c(500.0) rtol = 1.0e-12
    # The transition adds its enthalpy and its entropy, the Gibbs energy stays.
    δ = 1.0e-7
    @test H(Tt + δ) - H(Tt - δ) ≈ dH atol = 1.0e-3
    @test S(Tt + δ) - S(Tt - δ) ≈ dS atol = 1.0e-6
    @test abs(G(Tt + δ) - G(Tt - δ)) < 1.0e-3 + abs(dH - Tt * dS)
    # The top of the second interval carries no transition: H continuous.
    T2 = 1049.9999755859
    @test H(T2 + δ) - H(T2 - δ) ≈ 0 atol = 1.0e-3
    # Past 950 K the heat capacity is the second interval's polynomial.
    @test Cp(1000.0) ≈ 150.62399291992 rtol = 1.0e-9

    # Below the reference interval. No shipped record has an interval there,
    # so hematite is given one: a copy of its reference interval under 280 K,
    # with the transition of 950 K moved to its top. Going down across it takes
    # the transition's enthalpy and entropy away and changes nothing else:
    # H₂ = H − ΔH, S₂ = S − ΔS and G₂ = G − ΔH + T ΔS below 280 K, and the
    # reference interval is untouched.
    low = deepcopy(rec["TPMethods"][1])
    low["limitsTP"]["lowerT"], low["limitsTP"]["upperT"] = 200.0, 280.0
    low["m_phase_trans_props"]["values"][1] = 280.0
    rec2 = deepcopy(rec)
    rec2["TPMethods"][1]["limitsTP"]["lowerT"] = 280.0
    pushfirst!(rec2["TPMethods"], low)
    path = joinpath(mktempdir(), "hematite-below.json")
    write(path, JSON.json(merge(raw, Dict("substances" => [rec2], "reactions" => Any[]))))
    hem2 = only(build_species(path; verbose = false))
    for T in (250.0, 220.0)
        @test hem2[:ΔₐH⁰](T = T) ≈ H(T) - dH rtol = 1.0e-12
        @test hem2[:S⁰](T = T) ≈ S(T) - dS rtol = 1.0e-12
        @test hem2[:ΔₐG⁰](T = T) ≈ G(T) - dH + T * dS rtol = 1.0e-12
    end
    @test hem2[:ΔₐG⁰](T = 500.0) == G(500.0)
end
