# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# Differentiating a kinetic run with respect to its parameters: a rate constant
# captured by the rate law, an initial amount. The references are closed forms,
# or the certified equilibrium at a budget known in closed form; no difference
# quotient anywhere.

@testsection "a kinetic run is differentiated with respect to its parameters" begin
    subs = build_species(datapath("cemdata18-thermofun.json"))
    sp = speciation(subs, ["Cal", "Portlandite"]; aggregate_state = [AS_AQUEOUS])
    cs = ChemicalSystem(sp, CEMDATA_PRIMARIES)
    prim = [p for p in ("Ca+2", "CO3-2", "H2O@", "H+") if p in symbol.(cs.species)]
    names = symbol.(cs.species)
    i_cal = findfirst(==("Cal"), names)
    i_ca = findfirst(==("Ca+2"), names)
    tend = 1.0e4
    solver = KineticsSolver(; ode_solver = Rodas5P(), reltol = 1.0e-12, abstol = 1.0e-14)

    function water_and_calcite(n_cal)
        n = [zero(n_cal) * u"mol" for _ in cs.species]
        st = ChemicalState(cs, n)
        set_quantity!(st, "H2O@", 1.0u"kg")
        set_quantity!(st, "Cal", n_cal * u"mol")
        return st
    end
    function dissolution(rate)
        rxn = Reaction([cs["Cal"]], [cs[p] for p in prim]; symbol = "calcite dissolution")
        rxn[:rate] = KineticFunc(rate, (T = 298.15u"K", P = 1.0e5u"Pa"), u"mol/s")
        return rxn
    end

    @testset "the dual numbers a rate law captures are found" begin
        k = ForwardDiff.Dual{:t}(1.0e-4, 1.0)
        f = (T, P, t, n, lna, n0) -> k * n["Cal"]
        @test ChemistryLab._captured_number_type(f) == typeof(k)
        @test ChemistryLab._captured_number_type((T, P, t, n, lna, n0) -> 1.0) == Float64
        wrapped = (T, P, t, n, lna, n0) -> 2 * f(T, P, t, n, lna, n0)
        @test ChemistryLab._captured_number_type(wrapped) == typeof(k)
    end

    # First order, no equilibrium partition: `n(t) = n₀ e^{−kt}`.
    n0 = 0.5
    function calcite_left(k, n_cal = n0)
        kp = KineticsProblem(
            cs, [dissolution((T, P, t, n, lna, n0) -> k * n["Cal"])],
            water_and_calcite(n_cal), (0.0, tend); equilibrium_solver = nothing,
        )
        return integrate(kp, solver).u[end][1]
    end
    k0 = 1.0e-4
    e = exp(-k0 * tend)

    @testset "a rate constant, to first and second order" begin
        @test calcite_left(k0) ≈ n0 * e rtol = 1.0e-8
        @test ForwardDiff.derivative(calcite_left, k0) ≈ -tend * n0 * e rtol = 1.0e-6
        d2 = ForwardDiff.derivative(k -> ForwardDiff.derivative(calcite_left, k), k0)
        @test d2 ≈ tend^2 * n0 * e rtol = 1.0e-5
    end

    @testset "an initial amount" begin
        @test ForwardDiff.derivative(m -> calcite_left(k0, m), n0) ≈ e rtol = 1.0e-8
    end

    # Under partial equilibrium, at a constant rate: the extent is `k t`, the
    # budget of the equilibrium partition is `b₀ + Aₑ νₑᵀ k t` in closed form,
    # and the partition at the end is the certified equilibrium at that budget.
    # Differentiated through the run, it must be what the equilibrium gives
    # differentiated at that budget.
    @testset "the equilibrium partition of a run" begin
        model = DiluteSolutionModel()
        problem_at(k) = KineticsProblem(
            cs, [dissolution((T, P, t, n, lna, n0) -> k + zero(n["Cal"]))],
            water_and_calcite(n0), (0.0, tend);
            activity_model = model,
            equilibrium_solver = EquilibriumSolver(cs, model, OptimaOptimizer()),
        )
        run_at(k) = (kp = problem_at(k); (kp, integrate(kp, solver)))
        function ca_by_run(k)
            kp, sol = run_at(k)
            return sol.prob.p.n_full[i_ca]
        end
        function ca_by_equilibrium(k)
            kp = problem_at(1.0)
            idx = kp.idx_equilibrium
            p = ChemistryLab.build_kinetics_params(kp)
            be = kp.Ae * Float64[ustrip(us"mol", x) for x in kp.initial_state.n[idx]] .+
                (kp.Ae * transpose(kp.νe)) * [k * tend]
            st = ChemicalState(p.eq_system, Float64[ustrip(us"mol", x) for x in kp.initial_state.n[idx]] .* u"mol")
            eq, cert = solve_certified(p.eq_dual, (st,); b = be)
            @test cert.optimal
            return ustrip(us"mol", eq.n[findfirst(==(i_ca), idx)])
        end
        # A ten-thousandth of the calcite dissolves over the run. On plain numbers
        # the partition of the run is the interior point's, escalated to the
        # certified solve only past 1e-6 mol of imbalance: it agrees with the
        # certified answer to that accuracy (1.4e-6 measured). On dual numbers it
        # is the certified answer, lifted, and so is its derivative.
        kc = 1.0e-8
        @test ca_by_run(kc) ≈ ca_by_equilibrium(kc) rtol = 1.0e-5
        @test ForwardDiff.derivative(ca_by_run, kc) ≈ ForwardDiff.derivative(ca_by_equilibrium, kc) rtol = 1.0e-6
        # The certified replay of the run on dual numbers lifts each instant.
        function ca_by_replay(k)
            kp, sol = run_at(k)
            st = only(speciated_states(sol, kp; times = [tend]))
            return ustrip(us"mol", st.n[i_ca])
        end
        @test ForwardDiff.derivative(ca_by_replay, kc) ≈ ForwardDiff.derivative(ca_by_equilibrium, kc) rtol = 1.0e-6

        # What a run reaches only when a certified solve fails. The continuation
        # of the replay, from the composition certified at one instant, reaches
        # the answer the replay certifies at the next; it says so when it cannot,
        # whether the target budget has no equilibrium or the solver no step.
        kp1, sol1 = run_at(kc)
        p1 = sol1.prob.p
        t1, t2 = tend / 2, tend
        s1, s2 = speciated_states(sol1, kp1; times = [t1, t2])
        idx = kp1.idx_equilibrium
        n1 = [ustrip(us"mol", s1.n[i]) for i in idx]
        n2 = [ustrip(us"mol", s2.n[i]) for i in idx]
        sub = ChemistryLab._equilibrium_subsystem(kp1.system, kp1.idx_equilibrium)
        des = DualEquilibriumSolver(sub, model)
        Tv = ustrip(us"K", temperature(s2)) * u"K"
        Pv = ustrip(us"Pa", pressure(s2)) * u"Pa"
        be2 = collect(sol1(t2)[1:(p1.n_be)])
        cont(d, be) = ChemistryLab._replay_continuation(sol1, kp1, p1, d, sub, copy(n1), n1, s1, t1, t2, Tv, Pv, be)
        proved, n_c, _, _ = cont(des, be2)
        @test proved
        @test n_c ≈ n2 rtol = 1.0e-8
        @test !first(cont(des, -be2))
        stuck = DualEquilibriumSolver(sub, model; maxit = 0)
        quiet = Base.CoreLogging.NullLogger()
        @test !first(Base.CoreLogging.with_logger(() -> cont(stuck, be2), quiet))
        # An intermediate that raises is a failed intermediate, not an error: a
        # guess of the wrong length fails every one, and the step shrinks to none.
        @test !first(ChemistryLab._replay_continuation(sol1, kp1, p1, des, sub, n1[1:(end - 1)], n1, s1, t1, t2, Tv, Pv, be2))
        # The values a differentiated step keeps when its certified solve fails
        # are lifted where they stand: at a certified answer, the derivative of
        # the certified solve.
        e = zeros(length(be2)); e[1] = 1.0e-3
        lifted = ForwardDiff.derivative(x -> ChemistryLab._lifted_partition(p1, n2, Tv, Pv, be2 .+ x .* e), 0.0)
        st2 = ChemicalState(p1.eq_system, n2 .* u"mol"; T = Tv, P = Pv)
        direct = ForwardDiff.derivative(
            x -> [ustrip(us"mol", v) for v in first(solve_certified(p1.eq_dual, (st2,); b = be2 .+ x .* e)).n], 0.0,
        )
        @test lifted ≈ direct rtol = 1.0e-8 atol = 1.0e-14
    end
end
