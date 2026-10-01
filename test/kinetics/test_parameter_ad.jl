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
    end
end
