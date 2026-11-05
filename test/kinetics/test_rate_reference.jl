# The amounts a rate law reads as `n_initial`, when a host code restarts a run at
# each of its steps.
#
# The law `r = k n/n_initial` of calcite dissolving has the closed form
# `n(t) = n(0) exp(−k t/n_ref)` with `n_ref` the reference amount. A run cut in
# two halves, the second restarted from the state the first ended at, follows it
# only if the second half keeps the reference of the first.

using ForwardDiff
using OrdinaryDiffEq

@testsection "The reference amounts of the rate laws" begin
    sp = Dict(symbol(s) => s for s in build_species(datapath("slop98-inorganic-thermofun.json"); verbose = false))
    cs = ChemicalSystem([sp[s] for s in split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Cal")], ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"])
    dissolution() = Reaction(OrderedDict(cs["Cal"] => 1.0), OrderedDict(cs["Ca+2"] => 1.0, cs["CO3-2"] => 1.0); symbol = "calcite")
    k = 1.0e-6                                   # mol/s
    law = KineticFunc((T, P, t, n, lna, n0) -> k * n["Cal"] / n0["Cal"], NamedTuple(), u"mol/s")
    function water_with(cal)
        st = ChemicalState(cs)
        set_quantity!(st, "Cal", cal * u"mol")
        set_quantity!(st, "H2O@", 1.0u"kg")
        set_quantity!(st, "H+", 1.0e-7u"mol")
        set_quantity!(st, "OH-", 1.0e-7u"mol")
        return st
    end
    n0, T = 1.0e-2, 1.0e4                        # k T / n0 = 1
    cal(st) = ustrip(us"mol", moles(st, "Cal"))
    exact(t; nref = n0) = n0 * exp(-k * t / nref)
    ks = KineticsSolver(; ode_solver = Rodas5P(), reltol = 1.0e-10, abstol = 1.0e-14)

    @testset "a run restarted halfway keeps its reference" begin
        rxn = dissolution()
        rxn[:rate] = law
        st0 = water_with(n0)
        kp1 = KineticsProblem(cs, [rxn], st0, (0.0, T / 2))
        mid = state_at(integrate(kp1, ks), kp1, T / 2)
        @test cal(mid) ≈ exact(T / 2) rtol = 1.0e-7
        kept = KineticsProblem(cs, [rxn], mid, (T / 2, T); n_initial = st0)
        @test cal(state_at(integrate(kept, ks), kept, T)) ≈ exact(T) rtol = 1.0e-7
        # without it, the second half counts from the middle, and the law changes
        reset = KineticsProblem(cs, [rxn], mid, (T / 2, T))
        n_reset = cal(state_at(integrate(reset, ks), reset, T))
        @test n_reset ≈ cal(mid) * exp(-k * (T / 2) / cal(mid)) rtol = 1.0e-7
        @test abs(n_reset - exact(T)) > 0.1 * exact(T)
    end

    @testset "derivative with respect to the reference amount" begin
        # n(T) = n0 exp(−k T/n_ref): ∂n(T)/∂n_ref = n(T) k T / n_ref²
        rxn = dissolution()
        rxn[:rate] = law
        function final_amount(nref)
            st0 = water_with(n0)
            # the reference built from its amounts, which carry the derivative
            n = [ustrip(us"mol", x) * one(nref) for x in st0.n]
            n[findfirst(s -> symbol(s) == "Cal", cs.species)] = nref
            kp = KineticsProblem(cs, [rxn], st0, (0.0, T); n_initial = ChemicalState(cs, n .* u"mol"))
            return ustrip(us"mol", state_at(integrate(kp, ks), kp, T).n[findfirst(s -> symbol(s) == "Cal", cs.species)])
        end
        nref = 2.0e-2
        @test final_amount(nref) ≈ exact(T; nref) rtol = 1.0e-7
        @test ForwardDiff.derivative(final_amount, nref) ≈ exact(T; nref) * k * T / nref^2 rtol = 1.0e-6
    end

    @testset "the implicit step reads the reference it is given" begin
        kss = KineticStepSolver(cs, DiluteSolutionModel(), [KineticReaction(cs, dissolution(), law)])
        st0, st_mid = water_with(n0), water_with(0.6 * n0)
        Δt = 10.0
        for (reference, nref) in ((nothing, 0.6 * n0), (st0, n0))
            Δξ = Ref(Float64[])
            st1 = kinetic_step(kss, st_mid, Δt * u"s"; n_initial = reference, parameters = Δξ)
            # backward Euler: the rate at the end of the step
            @test only(Δξ[]) ≈ Δt * k * cal(st1) / nref rtol = 1.0e-8
            # and the adaptive step, whose three steps share the reference
            Δξa = Ref(Float64[])
            sta, used, _ = kinetic_step_adaptive(kss, st_mid, Δt * u"s"; n_initial = reference, parameters = Δξa)
            @test used == Δt
            @test cal(sta) ≈ cal(st_mid) * exp(-k * Δt / nref) rtol = 1.0e-6
        end
    end

    @testset "a replay that starts where the last one ended" begin
        # The same law, the solution now at equilibrium with what dissolves:
        # the speciated end of the second half, replayed from the speciated
        # middle, is the one the full replay finds.
        rxn = dissolution()
        rxn[:rate] = law
        model = DiluteSolutionModel()
        es = EquilibriumSolver(cs, model, OptimaOptimizer())
        st0 = water_with(n0)
        build(st, span) = KineticsProblem(cs, [rxn], st, span; activity_model = model, equilibrium_solver = es, n_initial = st0)
        kp1 = build(st0, (0.0, T / 2))
        mid = only(speciated_states(integrate(kp1, ks), kp1; times = [T / 2]))
        kp2 = build(mid, (T / 2, T))
        sol2 = integrate(kp2, ks)
        fast = only(speciated_states(sol2, kp2; times = [T], start = mid))
        full = only(speciated_states(sol2, kp2; times = [T]))
        amounts(st) = [ustrip(us"mol", x) for x in st.n]
        @test amounts(fast) ≈ amounts(full) rtol = 1.0e-8
        @test cal(fast) ≈ exact(T) rtol = 1.0e-6
        other = ChemicalSystem([sp[s] for s in split("H2O@ H+ OH- Ca+2 CO3-2 Cal")], ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"])
        @test_throws ArgumentError speciated_states(sol2, kp2; times = [T], start = ChemicalState(other))
        # The start a replay falls back on when the others fail: the cast
        # composition, on the budget of the instant and inside the bounds.
        p2 = sol2.prob.p
        be = sol2(T)[1:(p2.n_be)]
        cold = ChemistryLab._cold_replay_start(p2, be)
        @test p2.Ae * cold ≈ be rtol = 1.0e-10
        @test all(>=(ChemistryLab._EQ_GUESS_FLOOR), cold)
    end

    @testset "a reference of an equal system is read" begin
        # A system built again from the same species is another object: its
        # states are read by position, as those of the problem's own.
        twin = ChemicalSystem([sp[s] for s in split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Cal")], ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"])
        @test twin !== cs
        reference = ChemicalState(twin)
        set_quantity!(reference, "Cal", 2 * n0 * u"mol")
        rxn = dissolution()
        rxn[:rate] = law
        kp = KineticsProblem(cs, [rxn], water_with(n0), (0.0, T); n_initial = reference)
        @test cal(state_at(integrate(kp, ks), kp, T)) ≈ exact(T; nref = 2 * n0) rtol = 1.0e-7
    end

    @testset "a reference of another system is refused" begin
        other = ChemicalSystem([sp[s] for s in split("H2O@ H+ OH- Ca+2 CO3-2 Cal")], ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"])
        foreign = ChemicalState(other)
        rxn = dissolution()
        rxn[:rate] = law
        @test_throws ArgumentError KineticsProblem(cs, [rxn], water_with(n0), (0.0, T); n_initial = foreign)
        kss = KineticStepSolver(cs, DiluteSolutionModel(), [KineticReaction(cs, dissolution(), law)])
        @test_throws ArgumentError kinetic_step(kss, water_with(n0), 1.0u"s"; n_initial = foreign)
    end
end
