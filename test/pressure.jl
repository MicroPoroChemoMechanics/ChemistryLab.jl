# Pressure: the ideal gas, and the standard state of a condensed species.
#
# A gas is referred to the pure ideal gas at P° = 1 bar, so its activity in an
# ideal mixture is xᵢ P/P° and its chemical potential grows as RT/P, which is
# also its molar volume. A condensed species whose record declares a constant
# molar volume has its standard state at the pressure of the system, its standard
# Gibbs energy and enthalpy moving by V⁰ (P − P°). Every reference below is an
# identity: a Maxwell relation, the Henry scaling of a single gas, the value at
# the standard pressure.

using JSON

@testsection "Pressure" begin

    cem = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-thermofun.json"); verbose = false))

    @testset "the standard pressure is one bar" begin
        @test P_STANDARD === 1.0e5
        @test ustrip(us"Pa", P_STANDARD_Q) === P_STANDARD
    end

    @testset "a gas built without a molar volume has the ideal gas's" begin
        n2 = Species("N2"; aggregate_state = AS_GAS, class = SC_COMPONENT)
        @test !haskey(n2, :V⁰)
        cs = ChemicalSystem([n2], ["N2"])
        st = ChemicalState(cs; T = 310.0u"K", P = 2.0e5u"Pa", n = [0.5u"mol"])
        @test ustrip(us"m^3", volume(st).gas) ≈ 0.5 * R_GAS * 310.0 / 2.0e5 rtol = 1.0e-14
        @test ustrip(us"m^3", volume(st).total) ≈ ustrip(us"m^3", volume(st).gas) rtol = 1.0e-14
        set_pressure!(st, 4.0e5u"Pa")
        @test ustrip(us"m^3", volume(st).gas) ≈ 0.5 * R_GAS * 310.0 / 4.0e5 rtol = 1.0e-14
    end

    @testset "a database gas has V = RT/P at every temperature and pressure" begin
        co2 = cem["CO2"]
        @test aggregate_state(co2) == AS_GAS
        for (T, P) in ((298.15, 1.0e5), (323.15, 3.0e6), (278.15, 5.0e4))
            @test co2[:V⁰](T = T, P = P) ≈ R_GAS * T / P rtol = 1.0e-15
        end
        # And the pressure is in its activity, not in its standard energy.
        @test ForwardDiff.derivative(P -> co2[:ΔₐG⁰](T = 298.15, P = P), 3.0e6) == 0
    end

    @testset "a condensed species: ∂G⁰/∂P = ∂H⁰/∂P = V⁰, S⁰ unmoved" begin
        # `CA` is marked `mv_pvnrt` in CEMDATA18 while carrying a solid's volume:
        # the aggregate state decides, and it is treated as the crystal it is.
        for k in ("Cal", "Portlandite", "ettringite", "CA")
            s = cem[k]
            V = s[:V⁰](T = 298.15, P = 3.0e6)
            @test V == s[:V⁰](T = 298.15, P = 1.0e5)
            @test 0 < V < 1.0e-3
            for P in (5.0e4, 3.0e6, 1.0e8)
                @test ForwardDiff.derivative(Pv -> s[:ΔₐG⁰](T = 310.0, P = Pv), P) ≈ V rtol = 1.0e-14
                @test ForwardDiff.derivative(Pv -> s[:ΔₐH⁰](T = 310.0, P = Pv), P) ≈ V rtol = 1.0e-14
                @test ForwardDiff.derivative(Pv -> s[:S⁰](T = 310.0, P = Pv), P) == 0
            end
            @test s[:ΔₐG⁰](T = 310.0, P = 1.1e6) - s[:ΔₐG⁰](T = 310.0, P = 1.0e5) ≈ V * 1.0e6 rtol = 1.0e-9
        end
    end

    # The solvent follows the equation of state of water in pressure, its volume
    # at P° being the tabulated one. Every relation below is one of a single
    # Gibbs energy, by automatic differentiation.
    # The solvent follows the equation of state of water in temperature and in
    # pressure (test/water_eos_reference.jl for the equation itself); here, what
    # the pressure does to it.
    @testset "the solvent: compressed as the equation of state of water says" begin
        w = cem["H2O@"]
        for T in (283.15, 298.15, 333.15)
            V° = w[:V⁰](T = T, P = P_STANDARD)
            for P in (1.0e6, 3.0e7, 5.0e7)
                V = w[:V⁰](T = T, P = P)
                # ∂G/∂P = V, ∂H/∂P = V − T ∂V/∂T, ∂S/∂P = −∂V/∂T.
                dVdT = ForwardDiff.derivative(Tv -> w[:V⁰](T = Tv, P = P), T)
                @test ForwardDiff.derivative(Pv -> w[:ΔₐG⁰](T = T, P = Pv), P) ≈ V rtol = 1.0e-10
                @test ForwardDiff.derivative(Pv -> w[:ΔₐH⁰](T = T, P = Pv), P) ≈ V - T * dVdT rtol = 1.0e-9
                @test ForwardDiff.derivative(Pv -> w[:S⁰](T = T, P = Pv), P) ≈ -dVdT rtol = 1.0e-8
                # ΔG = ΔH − T ΔS for the changes from P° (the energies are of
                # formation, the entropy absolute), and Cp = ∂H/∂T at the pressure.
                δ(k) = w[k](T = T, P = P) - w[k](T = T, P = P_STANDARD)
                @test δ(:ΔₐG⁰) ≈ δ(:ΔₐH⁰) - T * δ(:S⁰) rtol = 1.0e-10
                @test ForwardDiff.derivative(Tv -> w[:ΔₐH⁰](T = Tv, P = P), T) ≈ w[:Cp⁰](T = T, P = P) rtol = 1.0e-9
                # Compressed, by about 4.5e-10 per pascal near 25 °C.
                @test V < V°
            end
        end
        # Against a constant volume: 10 J/mol at 500 bar.
        V° = w[:V⁰](T = 298.15, P = P_STANDARD)
        Δ(P) = w[:ΔₐG⁰](T = 298.15, P = P) - w[:ΔₐG⁰](T = 298.15)
        @test Δ(5.0e7) < V° * (5.0e7 - P_STANDARD)
        @test 5 < V° * (5.0e7 - P_STANDARD) - Δ(5.0e7) < 15
    end

    @testset "the standard pressure leaves every standard energy where it was" begin
        # At P° the pressure term is an exact zero: the tabulated ΔfG° comes back
        # at each record's own Tst to rounding, as it did without the term.
        raw = JSON.parsefile(datapath("cemdata18-thermofun.json"))
        for r in raw["substances"]
            k = String(r["symbol"])
            haskey(r, "sm_gibbs_energy") && aggregate_state(cem[k]) != AS_AQUEOUS || continue
            Tst = Float64(get(r, "Tst", 298.15))
            @test abs(cem[k][:ΔₐG⁰](T = Tst, P = P_STANDARD) - Float64(r["sm_gibbs_energy"]["values"][1])) < 1.0e-8
        end
    end

    # ∂ln aᵢ/∂P = 1/P for a gas, in every activity model: the Maxwell relation
    # ∂μᵢ/∂P = Vᵢ = RT/P with the ideal gas's volume.
    @testset "every activity model gives a gas ∂ln a/∂P = 1/P" begin
        names = ["H2O@", "Na+", "Cl-", "CO2@", "CO2"]
        cs = ChemicalSystem([cem[k] for k in names], ["H2O@", "Na+", "Cl-", "CO2@"])
        st = ChemicalState(cs)
        for (k, q) in zip(names, (1.0u"kg", 0.1u"mol", 0.1u"mol", 1.0e-3u"mol", 0.2u"mol"))
            set_quantity!(st, k, q)
        end
        n = ustrip.(us"mol", st.n)
        p = ChemistryLab._build_params(st)
        ig = findfirst(==("CO2"), symbol.(cs.species))
        models = (
            DiluteSolutionModel(),
            HKFActivityModel(),
            DaviesActivityModel(),
            TruesdellJonesActivityModel(; parameters = Dict("Na+" => (4.0, 0.075), "Cl-" => (3.5, 0.015))),
            SITActivityModel(; parameters = SITParameters([("Na+", "Cl-") => 0.03])),
            PitzerActivityModel(; parameters = build_pitzer_parameters(datapath("pitzer-reardon1990.toml"))),
        )
        for m in models
            lna = activity_model(cs, m)
            f(P) = lna(n, merge(p, (P = P,)))[ig]
            for P in (5.0e4, 1.0e5, 2.0e6)
                @test ForwardDiff.derivative(f, P) ≈ 1 / P rtol = 1.0e-14
                # A single gas: x = 1, so its activity is P/P°.
                @test f(P) ≈ log(P / P_STANDARD) atol = 1.0e-14
            end
        end
    end

    # A single gas over water: its activity is P/P°, so the dissolved amount
    # follows the pressure. Where the gas activity ignored P, the dissolved CO2
    # moved between 1 and 10 bar only by the pressure dependence of the solute's
    # own standard energy, a hundredth of what Henry's law asks.
    @testset "Henry: the dissolved gas follows its pressure" begin
        cs = ChemicalSystem([cem["H2O@"], cem["CO2@"], cem["CO2"]], ["H2O@", "CO2@"])
        model = DiluteSolutionModel()
        μ = build_potentials(cs, model)
        lnr, g_aq, g_gas = Float64[], Float64[], Float64[]
        for P in (1.0e5, 1.0e6)
            st = ChemicalState(cs; P = P * u"Pa")
            set_quantity!(st, "H2O@", 1.0u"kg")
            set_quantity!(st, "CO2", 0.5u"mol")
            eq, cert = equilibrate_certified(st; model)
            @test cert.optimal
            n = ustrip.(us"mol", eq.n)
            p = ChemistryLab._build_params(eq)
            mu = μ(n, p)
            # The gas carries ln(P/P°) and meets the dissolved CO2 there.
            @test mu[3] ≈ p.ΔₐG⁰overRT[3] + log(P / P_STANDARD) atol = 1.0e-12
            @test mu[2] ≈ mu[3] atol = 1.0e-8
            push!(lnr, log(n[2] / n[1]))
            push!(g_aq, p.ΔₐG⁰overRT[2])
            push!(g_gas, p.ΔₐG⁰overRT[3])
        end
        # Under the ideal model ln a = ln(n/n_w) + const for the solute, so the
        # ratio of dissolved to water follows ln 10 and the two standard energies.
        @test lnr[2] - lnr[1] ≈ log(10) - (g_aq[2] - g_aq[1]) + (g_gas[2] - g_gas[1]) atol = 1.0e-8
        @test g_gas[2] == g_gas[1]
    end

    # Two gases: each meets its solute at its own fugacity xᵢ P. The dual solver
    # files a gas among the pure phases, which holds because a present gas is
    # solved on its own stationarity with the activity of the mixture.
    @testset "a gas mixture: each gas at its fugacity" begin
        cs = ChemicalSystem([cem["H2O@"], cem["CO2@"], cem["CO2"], cem["N2"]], ["H2O@", "CO2@", "N2"])
        model = DiluteSolutionModel()
        μ = build_potentials(cs, model)
        P = 2.0e5
        st = ChemicalState(cs; P = P * u"Pa")
        set_quantity!(st, "H2O@", 1.0u"kg")
        set_quantity!(st, "CO2", 0.5u"mol")
        set_quantity!(st, "N2", 0.5u"mol")
        eq, cert = equilibrate_certified(st; model)
        @test cert.optimal
        n = ustrip.(us"mol", eq.n)
        p = ChemistryLab._build_params(eq)
        mu = μ(n, p)
        x = n[3] / (n[3] + n[4])
        @test mu[3] ≈ p.ΔₐG⁰overRT[3] + log(x) + log(P / P_STANDARD) atol = 1.0e-12
        @test mu[2] ≈ mu[3] atol = 1.0e-8
        @test n[4] ≈ 0.5 rtol = 1.0e-12
    end
end
