# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# The ODE route under a rate law that reads the speciation.
#
# Calcite dissolving under r = k(1 − Ω) over 1e5 s, k = 1e-4 mol/s, the case the
# implicit step is checked on (test_implicit_step.jl). With the partition frozen
# within a step, the rate is constant within it, the stiff method integrates the
# extent explicitly, and a step longer than the relaxation of Ω overshot the
# equilibrium: Rodas5P ended on 487 mol of calcite from 0.05 and reported
# success. The partition is now solved in the right-hand side for such a law.

using LinearAlgebra: I
using OrderedCollections

@testsection "a rate law that reads the speciation" begin
    sp = Dict(symbol(x) => x for x in build_species(datapath("slop98-inorganic-thermofun.json"); verbose = false))
    cs = ChemicalSystem(
        [sp[x] for x in split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Cal")],
        ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"],
    )
    i_Cal = findfirst(x -> symbol(x) == "Cal", cs.species)
    model = DiluteSolutionModel()
    GT = [
        ustrip(us"J/mol", x[:ΔₐG⁰](T = 298.15u"K", P = 1.0e5u"Pa"; unit = true)) / (R_GAS * 298.15)
            for x in cs.species
    ]
    calcite() = Reaction(
        OrderedDict(cs["Cal"] => 1.0), OrderedDict(cs["Ca+2"] => 1.0, cs["CO3-2"] => 1.0);
        symbol = "calcite",
    )
    stoich = Float64.(
        KineticReaction(cs, calcite(), KineticFunc((T, P, t, n, lna, n0) -> 0.0, NamedTuple(), u"mol/s")).stoich
    )
    k = 1.0e-4
    ω_law = KineticFunc(
        (T, P, t, n, lna, n0) -> k * (1 - saturation_ratio(stoich, [lna[symbol(x)] for x in cs.species], GT)),
        NamedTuple(), u"mol/s",
    )
    amount_law = KineticFunc(
        (T, P, t, n, lna, n0) -> 1.0e-6 * max(n["Cal"], zero(eltype(n.data))), NamedTuple(), u"mol/s",
    )
    function state()
        st = ChemicalState(cs)
        set_quantity!(st, "Cal", 5.0e-2u"mol")
        set_quantity!(st, "H2O@", 1.0u"kg")
        set_quantity!(st, "H+", 1.0e-7u"mol")
        set_quantity!(st, "OH-", 1.0e-7u"mol")
        return st
    end
    function problem(law, tend)
        rxn = calcite()
        rxn[:rate] = law
        return KineticsProblem(
            cs, [rxn], state(), (0.0, tend);
            activity_model = model, equilibrium_solver = EquilibriumSolver(cs, model, OptimaOptimizer()),
        )
    end
    Cal_eq = ustrip(us"mol", SciMLBase.solve(DualEquilibriumSolver(cs, model), state()).n[i_Cal])
    ks = KineticsSolver(; ode_solver = Rodas5P(), reltol = 1.0e-8, abstol = 1.0e-12)
    kp = problem(ω_law, 1.0e5)

    @testset "the law is seen to read the partition, and the run reaches the equilibrium" begin
        sol = integrate(kp, ks)
        p = sol.prob.p
        @test p.rates_read_speciation[]
        @test p.rhs_mode[] === :rhs
        @test SciMLBase.successful_retcode(sol)
        cal(u) = u[p.n_be + 1]
        @test cal(sol.u[end]) ≈ Cal_eq rtol = 1.0e-6
        # Feasible along the whole trajectory, not only at its end.
        @test all(u -> -1.0e-12 <= cal(u) <= 5.0e-2 * (1 + 1.0e-8), sol.u)
        @test !any(u -> ChemistryLab._kinetic_state_infeasible(p, u), sol.u)
    end

    @testset "the Jacobian of the right-hand side holds the partition's derivative" begin
        # Frozen, the columns of `bₑ` were identically zero. Solved in the
        # right-hand side, they are ∂r/∂nₑ ⋅ ∂nₑ/∂bₑ, by the chain rule.
        f! = build_kinetics_ode(kp)
        p = build_kinetics_params(kp)
        u = build_u0(kp)
        ChemistryLab.respeciate!(p, u)
        p.rhs_mode[] = :rhs
        u[p.n_be + 1] -= 1.0e-4                 # some calcite dissolved, Ω < 1
        u[1:(p.n_be)] .+= p.Ae * (p.νe' * [1.0e-4])
        J = ForwardDiff.jacobian((du, x) -> f!(du, x, p, 0.0), similar(u), u)
        nb = p.n_be
        @test any(!iszero, J[:, 1:nb])
        be = u[1:nb]
        S = ForwardDiff.jacobian(b -> ChemistryLab._rhs_partition(p, b, 298.15), be)
        ne = ChemistryLab._rhs_partition(p, be, 298.15)
        n0 = StateView(p.n_initial_full, p.species_index)
        function rates(x)
            n = convert(Vector{eltype(x)}, p.n_full)
            n[p.idx_kinetic[1]] = u[nb + 1]
            n[p.idx_equilibrium] .= x
            lna = p.lna_fn(n, p)
            return [kr.rate_fn(p.T, p.P, 0.0, StateView(n, p.species_index), StateView(lna, p.species_index), n0) for kr in p.kin_rxns]
        end
        Jr = ForwardDiff.jacobian(rates, ne)
        M = vcat(p.Ae * p.νe', p.νk', Matrix(1.0I, 1, 1))
        @test J[:, 1:nb] ≈ M * (Jr * S) rtol = 1.0e-8
    end

    @testset "a law on amounts keeps the frozen route, bit for bit" begin
        kp2 = problem(amount_law, 3600.0)
        auto = integrate(kp2, ks)
        @test !auto.prob.p.rates_read_speciation[]
        @test auto.prob.p.rhs_mode[] === :frozen
        frozen = integrate(kp2, ks; speciation = :frozen)
        @test auto.u == frozen.u
    end

    @testset "frozen under the Ω law, the run is not reported a success" begin
        sol = @test_logs (:warn, r"no chemistry can produce") match_mode = :any integrate(kp, ks; speciation = :frozen)
        @test !SciMLBase.successful_retcode(sol)
        @test any(u -> ChemistryLab._kinetic_state_infeasible(sol.prob.p, u), sol.u)
    end

    @testset "a state with no partition makes the right-hand side NaN" begin
        # Negative element totals: no partition meets them, from the last one
        # or from a reconstruction, and the right-hand side is NaN, which the
        # integrator rejects as a step.
        f! = build_kinetics_ode(kp)
        p = build_kinetics_params(kp)
        u = build_u0(kp)
        ChemistryLab.respeciate!(p, u)
        p.rhs_mode[] = :rhs
        u[1:(p.n_be)] .*= -1
        quiet(f) = Base.CoreLogging.with_logger(f, Base.CoreLogging.NullLogger())
        @test quiet(() -> ChemistryLab._rhs_partition(p, u[1:(p.n_be)], 298.15)) === nothing
        du = similar(u)
        quiet(() -> f!(du, u, p, 0.0))
        @test all(isnan, du)
    end

    @testset "the modes are checked" begin
        @test_throws ArgumentError integrate(kp, ks; speciation = :sometimes)
        saved = ChemistryLab._DUAL_AVAILABLE[]
        try
            ChemistryLab._DUAL_AVAILABLE[] = false
            @test_throws ArgumentError integrate(kp, ks; speciation = :rhs)
        finally
            ChemistryLab._DUAL_AVAILABLE[] = saved
        end
    end

    @testset "what the system cannot hold" begin
        p = build_kinetics_params(kp)
        u = build_u0(kp)
        @test !ChemistryLab._kinetic_state_infeasible(p, u)
        u[p.n_be + 1] = 0.05 * (1 + 1.0e-6)     # more calcium than the system has
        @test ChemistryLab._kinetic_state_infeasible(p, u)
        u[p.n_be + 1] = -1.0e-6                  # a negative amount
        @test ChemistryLab._kinetic_state_infeasible(p, u)
        u[p.n_be + 1] = NaN
        @test ChemistryLab._kinetic_state_infeasible(p, u)
    end
end
