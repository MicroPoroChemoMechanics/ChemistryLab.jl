# Real gases: the equation of state of Peng and Robinson (1976).
#
# Every reference but one is an identity the equation obeys whatever its
# constants: the fugacity coefficients are the derivatives of the residual Gibbs
# energy, they satisfy Gibbs-Duhem, their pressure derivative is the partial
# molar volume, they tend to the second virial coefficient as P → 0, and the
# cubic has its triple root at the critical point. The one comparison with
# experiment is the vapor pressure of carbon dioxide at 25 °C.

@testsection "Real gases" begin
    CL = ChemistryLab
    cem = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-thermofun.json"); verbose = false))
    co2 = peng_robinson(cem["CO2"])
    n2 = peng_robinson(cem["N2"]; kij = [:CO2 => 0.1])
    ch4 = peng_robinson(Species("CH4"; aggregate_state = AS_GAS, class = SC_COMPONENT))
    # The mixing description of a list of gases, as a system holding them gets it.
    mixing(sp...) = CL._gas_mixing((idx_gas = collect(eachindex(sp)), species = collect(sp)))

    @testset "the constants, from phreeqc.dat or given" begin
        t = literature_table("ParkhurstAppelo2013", "gas_critical_constants")
        @test properties(co2)[:T_c] == 304.2
        @test properties(co2)[:P_c] ≈ 72.86 * 101325 rtol = 1.0e-15
        @test properties(co2)[:ω] == 0.225
        @test properties(co2)[:kij] == Pair{Symbol, Float64}[]
        @test properties(n2)[:kij] == [:CO2 => 0.1]
        @test length(t.gas) == 8
        # The species of the database is left as it was: the constants are on a copy.
        @test !haskey(properties(cem["CO2"]), :T_c)
        @test cem["CO2"][:ΔₐG⁰](T = 298.15, P = 1.0e5) == co2[:ΔₐG⁰](T = 298.15, P = 1.0e5)
        given = peng_robinson(cem["CO2"]; T_c = 304.128u"K", P_c = 73.773u"bar", ω = 0.2239)
        @test properties(given)[:T_c] == 304.128
        @test properties(given)[:P_c] ≈ 7.3773e6 rtol = 1.0e-15
        @test properties(peng_robinson(cem["CO2"]; T_c = 304.128, P_c = 7.3773e6, ω = 0.2239))[:P_c] == 7.3773e6
        @test_throws ArgumentError peng_robinson(cem["Cal"])
        @test_throws ArgumentError peng_robinson(Species("Ar"; aggregate_state = AS_GAS))
        @test_throws ArgumentError peng_robinson(cem["CO2"]; T_c = 304.2)
        @test_throws ArgumentError peng_robinson(cem["CO2"]; T_c = -1.0, P_c = 7.0e6, ω = 0.2)
        @test_throws DimensionError peng_robinson(cem["CO2"]; T_c = 304.2u"K", P_c = 300.0u"K", ω = 0.2)
    end

    @testset "the cubic has its triple root at the critical point" begin
        # The two constants are the exact roots of the critical conditions, which
        # Peng and Robinson print to five decimals.
        @test round(CL._PR_ΩA; digits = 5) == 0.45724
        @test round(CL._PR_ΩB; digits = 5) == 0.0778
        A, B, Z = CL._PR_ΩA, CL._PR_ΩB, CL._PR_ZC
        @test abs(CL._pr_cubic(Z, A, B)) < 1.0e-16
        @test abs(CL._pr_cubic′(Z, A, B)) < 1.0e-15
        @test abs(6Z - 2 * (1 - B)) < 1.0e-15
        @test round(Z; digits = 4) == 0.3074
        # Every gas at its own critical point, whatever its acentric factor. A
        # triple root moves by the cube root of a perturbation of the cubic, so
        # the rounding of A and B, 1e-16, leaves it good to about 5e-6 only.
        for g in (co2, n2, ch4)
            mix = mixing(g)
            _, Zc = CL._pr_ln_phi(mix, [1.0], mix.T_c[1], mix.P_c[1])
            @test Zc ≈ CL._PR_ZC rtol = 2.0e-5
        end
    end

    @testset "ln φᵢ = ∂(n G_res/RT)/∂nᵢ, Gibbs-Duhem, and the partial molar volume" begin
        mix = mixing(co2, n2, ch4)
        @test mix.kij == [0 0.1 0; 0.1 0 0; 0 0 0]
        R = CL.R_GAS
        # A dilute gas, a dense fluid rich in carbon dioxide at 25 °C, a
        # supercritical one, and a cold one.
        for (T, P, n) in (
                (298.15, 1.0e5, [0.3, 0.5, 0.2]), (298.15, 1.0e7, [0.9, 0.05, 0.05]),
                (350.0, 4.0e7, [0.3, 0.6, 0.1]), (250.0, 5.0e6, [0.5, 0.3, 0.2]),
            )
            lnφ, Z = CL._pr_ln_phi(mix, n ./ sum(n), T, P)
            @test ForwardDiff.gradient(x -> CL._pr_residual_gibbs(mix, x, T, P), n) ≈ lnφ atol = 1.0e-13
            J = ForwardDiff.jacobian(x -> CL._pr_ln_phi(mix, x ./ sum(x), T, P)[1], n)
            @test J ≈ J' atol = 1.0e-12
            @test maximum(abs, n' * J) < 1.0e-12
            # ∂ln φᵢ/∂P = V̄ᵢ/(RT) − 1/P, with V̄ᵢ = ∂V/∂nᵢ of the phase volume Z N R T/P.
            dlnφ = ForwardDiff.derivative(Pv -> CL._pr_ln_phi(mix, n ./ sum(n), T, Pv)[1], P)
            V̄ = ForwardDiff.gradient(x -> CL._pr_phase_volume(mix, x, T, P), n)
            @test R * T .* (dlnφ .+ 1 / P) ≈ V̄ rtol = 1.0e-10
            @test CL._pr_phase_volume(mix, n, T, P) ≈ Z * sum(n) * R * T / P rtol = 1.0e-15
            @test sum(n ./ sum(n) .* dlnφ) ≈ (Z - 1) / P rtol = 1.0e-10
        end
    end

    @testset "the ideal gas as P → 0: the second virial coefficient" begin
        mix = mixing(co2)
        R = CL.R_GAS
        for T in (250.0, 298.15, 400.0)
            a, b = CL._pr_ab(mix, T)
            lnφ, Z = CL._pr_ln_phi(mix, [1.0], T, 1.0e-3)
            @test abs(lnφ[1]) < 1.0e-9
            @test Z ≈ 1 atol = 1.0e-9
            # ln φ = B₂ P/(RT) + O(P²), with B₂ = b − a/(RT) for this equation:
            # at 100 Pa the second-order term is a few 1e-6 of the first, and the
            # rounding of a ln φ of 1e-5 is far below it.
            P = 100.0
            lnφ, _ = CL._pr_ln_phi(mix, [1.0], T, P)
            @test lnφ[1] / P ≈ (b[1] - a[1] / (R * T)) / (R * T) rtol = 1.0e-5
        end
    end

    @testset "the vapor pressure of carbon dioxide at 25 °C" begin
        mix = mixing(co2)
        T = 298.15
        # The two roots of the cubic as a function of the pressure, in the range
        # where there are three: equal energies mark the saturation.
        function roots(P)
            a, b = CL._pr_ab(mix, T)
            A, B = a[1] * P / (CL.R_GAS * T)^2, b[1] * P / (CL.R_GAS * T)
            r = CL._pr_roots(A, B)
            return extrema(r), A, B
        end
        function Δg(P)
            (zl, zv), A, B = roots(P)
            return CL._pr_g_res(zv, A, B) - CL._pr_g_res(zl, A, B)
        end
        lo, hi = 6.1e6, 6.6e6
        @test length(CL._pr_roots(roots(lo)[2], roots(lo)[3])) == 3
        @test Δg(lo) < 0 < Δg(hi)
        for _ in 1:60
            mid = (lo + hi) / 2
            Δg(mid) < 0 ? (lo = mid) : (hi = mid)
        end
        Psat = (lo + hi) / 2
        ref = ustrip(us"Pa", literature_value("Spycher2003", "co2_vapor_pressure_25C"))
        @test abs(Psat / ref - 1) < 2.0e-3
        # The stable root jumps from the vapor to the liquid there, and the
        # fugacity does not: that is what the saturation is.
        below, above = CL._pr_ln_phi(mix, [1.0], T, Psat * (1 - 1.0e-6)), CL._pr_ln_phi(mix, [1.0], T, Psat * (1 + 1.0e-6))
        @test below[2] > 0.4
        @test above[2] < 0.2
        @test below[1][1] ≈ above[1][1] atol = 1.0e-5
    end

    @testset "every activity model gives a real gas its fugacity coefficient" begin
        names = ["H2O@", "Na+", "Cl-", "CO2@", "CO2"]
        sp = [k == "CO2" ? co2 : cem[k] for k in names]
        cs = ChemicalSystem(sp, ["H2O@", "Na+", "Cl-", "CO2@"])
        st = ChemicalState(cs)
        for (k, q) in zip(names, (1.0u"kg", 0.1u"mol", 0.1u"mol", 1.0e-3u"mol", 0.2u"mol"))
            set_quantity!(st, k, q)
        end
        n = ustrip.(us"mol", st.n)
        p = CL._build_params(st)
        ig = findfirst(==("CO2"), symbol.(cs.species))
        mix = mixing(co2)
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
            for P in (1.0e5, 5.0e6, 2.0e7)
                lnφ, Z = CL._pr_ln_phi(mix, [1.0], p.T, P)
                @test f(P) ≈ lnφ[1] + log(P / P_STANDARD) atol = 1.0e-14
                # ∂μ/∂P = V: the molar volume of the real gas, Z R T/P.
                @test ForwardDiff.derivative(f, P) ≈ Z / P rtol = 1.0e-10
            end
        end
    end

    @testset "a mixture of real gases: Gibbs-Duhem on the activities" begin
        cs = ChemicalSystem([cem["H2O@"], cem["CO2@"], co2, n2], ["H2O@", "CO2@", "N2"])
        st = ChemicalState(cs; P = 1.0e7u"Pa")
        for (k, q) in zip(("H2O@", "CO2@", "CO2", "N2"), (1.0u"kg", 0.5u"mol", 0.7u"mol", 0.3u"mol"))
            set_quantity!(st, k, q)
        end
        n = ustrip.(us"mol", st.n)
        p = CL._build_params(st)
        ϵ = hasproperty(p, :ϵ) ? p.ϵ : CL._AMOUNT_FLOOR
        ig = cs.idx_gas
        lna = activity_model(cs, DiluteSolutionModel())
        y = (n[ig] .+ ϵ) ./ sum(n[ig] .+ ϵ)
        lnφ, _ = CL._pr_ln_phi(mixing(co2, n2), y, p.T, p.P)
        @test lna(n, p)[ig] ≈ log.(y) .+ lnφ .+ log(p.P / P_STANDARD) atol = 1.0e-14
        J = ForwardDiff.jacobian(x -> lna(x, p)[ig], n)[:, ig]
        @test J ≈ J' atol = 1.0e-12
        @test maximum(abs, (n[ig] .+ ϵ)' * J) < 1.0e-12
    end

    @testset "a gas phase is real or ideal as a whole" begin
        cs = ChemicalSystem([cem["H2O@"], cem["CO2@"], co2, cem["N2"]], ["H2O@", "CO2@", "N2"])
        @test_throws ArgumentError activity_model(cs, DiluteSolutionModel())
        twice = peng_robinson(cem["CO2"]; kij = [:N2 => 0.05])
        @test_throws ArgumentError mixing(twice, n2)
        # A partner absent from the phase is ignored.
        @test mixing(peng_robinson(cem["CO2"]; kij = [:CH4 => 0.09]), n2).kij == [0 0.1; 0.1 0]
    end

    @testset "the volume of a real gas phase is Z N R T/P" begin
        cs = ChemicalSystem([co2], ["CO2"])
        T, P = 310.0, 1.0e7
        st = ChemicalState(cs; T = T * u"K", P = P * u"Pa", n = [0.5u"mol"])
        _, Z = CL._pr_ln_phi(mixing(co2), [1.0], T, P)
        @test Z < 0.5
        @test ustrip(us"m^3", volume(st).gas) ≈ Z * 0.5 * CL.R_GAS * T / P rtol = 1.0e-14
        @test CL._total_volume(cs, [0.5], T, P) ≈ Z * 0.5 * CL.R_GAS * T / P rtol = 1.0e-14
        # A constraint counting the volume as Σ nᵢ V⁰ᵢ refuses a real gas.
        @test_throws ArgumentError CL._molar_volumes(cs, T * u"K", P * u"Pa")
        # An ideal gas keeps the volume of the ideal gas.
        ideal = ChemicalSystem([cem["CO2"]], ["CO2"])
        @test CL._total_volume(ideal, [0.5], T, P) ≈ 0.5 * CL.R_GAS * T / P rtol = 1.0e-14
        # The accessors read the same equation of state.
        lnφ, _ = CL._pr_ln_phi(mixing(co2), [1.0], T, P)
        @test compressibility_factor(st) == Z
        @test fugacity_coefficients(st) == OrderedDict("CO2" => exp(lnφ[1]))
        st_ideal = ChemicalState(ideal; T = T * u"K", P = P * u"Pa", n = [0.5u"mol"])
        @test compressibility_factor(st_ideal) == 1
        @test fugacity_coefficients(st_ideal) == OrderedDict("CO2" => 1.0)
        @test_throws ArgumentError compressibility_factor(ChemicalState(ChemicalSystem([cem["H2O@"]], ["H2O@"])))
        # In a mixture, at the composition of the phase.
        cs2 = ChemicalSystem([co2, n2], ["CO2", "N2"])
        st2 = ChemicalState(cs2; T = T * u"K", P = P * u"Pa", n = [0.3u"mol", 0.7u"mol"])
        y = [0.3, 0.7] .+ CL._AMOUNT_FLOOR
        lnφ2, Z2 = CL._pr_ln_phi(mixing(co2, n2), y ./ sum(y), T, P)
        @test collect(values(fugacity_coefficients(st2))) ≈ exp.(lnφ2) rtol = 1.0e-14
        @test compressibility_factor(st2) ≈ Z2 rtol = 1.0e-14
        @test ustrip(us"m^3", volume(st2).gas) ≈ Z2 * CL.R_GAS * T / P rtol = 1.0e-12
    end

    @testset "Henry's law at the fugacity" begin
        # Carbon dioxide over water at 50 bar and 25 °C, in excess: under the
        # ideal solution model the dissolved amount is proportional to the
        # fugacity of the gas, so the real gas dissolves φ times what the ideal one
        # does. Water is fixed by the hydrogen balance, so the ratio is exact.
        T, P = 298.15, 5.0e6
        model = DiluteSolutionModel()
        dissolved = Float64[]
        for gas in (cem["CO2"], co2)
            cs = ChemicalSystem([cem["H2O@"], cem["CO2@"], gas], ["H2O@", "CO2@"])
            st = ChemicalState(cs; T = T * u"K", P = P * u"Pa")
            set_quantity!(st, "H2O@", 1.0u"kg")
            set_quantity!(st, "CO2", 5.0u"mol")
            eq, cert = equilibrate_certified(st; model)
            @test cert.optimal
            @test cert.scope === :global_minimum
            n = ustrip.(us"mol", eq.n)
            @test n[3] > 1
            mu = build_potentials(cs, model)(n, CL._build_params(eq))
            @test mu[2] ≈ mu[3] atol = 1.0e-8
            push!(dissolved, n[2])
        end
        lnφ, _ = CL._pr_ln_phi(mixing(co2), [1.0], T, P)
        @test log(dissolved[2] / dissolved[1]) ≈ lnφ[1] atol = 1.0e-8
        @test lnφ[1] < -0.2
    end

    @testset "a mixture of real gases at equilibrium with water" begin
        cs = ChemicalSystem([cem["H2O@"], cem["CO2@"], co2, n2], ["H2O@", "CO2@", "N2"])
        model = DiluteSolutionModel()
        st = ChemicalState(cs; P = 1.0e7u"Pa")
        set_quantity!(st, "H2O@", 1.0u"kg")
        set_quantity!(st, "CO2", 2.0u"mol")
        set_quantity!(st, "N2", 2.0u"mol")
        eq, cert = equilibrate_certified(st; model)
        @test cert.optimal
        # No convexity is proved for a mixture of real gases.
        @test cert.scope !== :global_minimum
        @test any(r -> occursin("Peng and Robinson", r), cert.scope_reasons)
        n = ustrip.(us"mol", eq.n)
        p = CL._build_params(eq)
        mu = build_potentials(cs, model)(n, p)
        y = n[3] / (n[3] + n[4])
        lnφ, _ = CL._pr_ln_phi(mixing(co2, n2), [y, 1 - y], p.T, p.P)
        @test mu[3] ≈ p.ΔₐG⁰overRT[3] + log(y) + lnφ[1] + log(p.P / P_STANDARD) atol = 1.0e-8
        @test mu[2] ≈ mu[3] atol = 1.0e-8
        @test n[4] ≈ 2 rtol = 1.0e-12
    end

    @testset "constants being fitted carry their derivatives" begin
        # The critical constants, the acentric factor and the interaction
        # parameters were stored as `Float64`, and a dual one raised.
        lnφk(k) = CL._pr_ln_phi(mixing(co2, peng_robinson(cem["N2"]; kij = [:CO2 => k])), [0.5, 0.5], 298.15, 5.0e6)[1][1]
        dk = ForwardDiff.derivative(lnφk, 0.1)
        @test isfinite(dk) && dk != 0
        lnφω(w) = CL._pr_ln_phi(mixing(peng_robinson(cem["CO2"]; T_c = 304.2, P_c = 7.38e6, ω = w)), [1.0], 298.15, 5.0e6)[1][1]
        dω = ForwardDiff.derivative(lnφω, 0.225)
        @test isfinite(dω) && dω != 0
    end
end
