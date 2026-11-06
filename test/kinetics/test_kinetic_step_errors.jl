# What `kinetic_step` returns and what it refuses: the eliminated step that
# `coupling = :reactions` falls back on for a mineral per reaction, its domain
# `n ≥ 0`, and the `KineticStepError` that replaces any answer that is not
# certified.

using LinearAlgebra

@testsection "Kinetic step: certified answers and explicit refusals" begin

    sp = Dict(symbol(s) => s for s in build_species(datapath("slop98-inorganic-thermofun.json")))
    cs = ChemicalSystem(
        [sp[s] for s in split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Cal")],
        ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"],
    )
    A = Matrix{Float64}(cs.SM.A)
    i_Cal = findfirst(s -> symbol(s) == "Cal", cs.species)
    i_Ca = findfirst(s -> symbol(s) == "Ca+2", cs.species)
    rxn = Reaction(
        OrderedDict(cs["Cal"] => 1.0), OrderedDict(cs["Ca+2"] => 1.0, cs["CO3-2"] => 1.0);
        symbol = "calcite dissolution",
    )
    law(f) = KineticReaction(cs, rxn, KineticFunc(f, NamedTuple(), u"mol/s"))
    amounts(st) = Float64[ustrip(us"mol", x) for x in st.n]

    function water(; calcite = 1.0e-2, calcium = 0.0)
        st = ChemicalState(cs)
        set_quantity!(st, "H2O@", 1.0u"kg")
        set_quantity!(st, "H+", 1.0e-7u"mol")
        set_quantity!(st, "OH-", 1.0e-7u"mol")
        calcite > 0 && set_quantity!(st, "Cal", calcite * u"mol")
        if calcium > 0
            set_quantity!(st, "Ca+2", calcium * u"mol")
            set_quantity!(st, "CO3-2", calcium * u"mol")
        end
        return st
    end

    @testsection "`:reactions` keeps the eliminated form of a mineral per reaction" begin
        kss = KineticStepSolver(cs, DiluteSolutionModel(), [law((T, P, t, n, lna, n0) -> 1.0e-6)])
        @test kss.coupling === :reactions && kss.dual_free !== nothing
        @test_throws "`coupling` must be" KineticStepSolver(cs, DiluteSolutionModel(), [law((T, P, t, n, lna, n0) -> 1.0e-6)]; coupling = :auto)
        # the same step under both, where both are defined: K square, invertible,
        # on pure phases
        @test ChemistryLab._eliminable(reshape([-1.0], 1, 1), [1], [1])
        @test !ChemistryLab._eliminable([-1.0 0.0; 0.0 -1.0; 1.0 1.0], [1, 2, 3], [1, 2, 3])
        @test !ChemistryLab._eliminable(reshape([-1.0], 1, 1), [1], Int[])
        @test !ChemistryLab._eliminable([-1.0 -1.0; 1.0 1.0], [1, 2], [1, 2])
    end

    # The rate `k (c − n_Ca)` vanishes at the calcium of the equilibrium of the
    # same budget, so a step much longer than `1/k` lands on that equilibrium:
    # an independent reference, from the equilibrium solver.
    @testsection "a long step reaches the equilibrium, from either side" begin
        for (start, label) in ((water(), "dissolving"), (water(calcite = 0.0, calcium = 2.0e-3), "precipitating from zero"))
            eq = equilibrate(start; model = DiluteSolutionModel())
            c = amounts(eq)[i_Ca]
            kss = KineticStepSolver(cs, DiluteSolutionModel(), [law((T, P, t, n, lna, n0) -> 1.0e-3 * (c - n["Ca+2"]))])
            cert = Ref{Any}(nothing)
            st = kinetic_step(kss, start, 1.0e9; certificate = cert)
            n = amounts(st)
            @test cert[].optimal
            @test n[i_Cal] ≈ amounts(eq)[i_Cal] rtol = 1.0e-6
            @test maximum(abs, A * n .- A * amounts(start)) <= 1.0e-12
            @test n[i_Cal] >= 0
        end
    end

    @testsection "a rate that does not vanish at zero is refused, and why" begin
        # 1e-3 mol/s for 100 s dissolves ten times the calcite present.
        kss = KineticStepSolver(cs, DiluteSolutionModel(), [law((T, P, t, n, lna, n0) -> 1.0e-3)])
        err = try
            kinetic_step(kss, water(), 100.0)
        catch e
            e
        end
        @test err isa KineticStepError
        @test err.cause === :exhausted
        @test occursin("Cal", err.message) && occursin("must vanish", err.message)
        @test occursin("KineticStepError (exhausted, Δt = 100.0 s)", sprint(showerror, err))
        # the adaptive step halves on it, then reports the last refusal
        err_a = try
            kinetic_step_adaptive(kss, water(), 100.0; max_halvings = 2)
        catch e
            e
        end
        @test err_a isa KineticStepError && err_a.cause === :exhausted
        @test occursin("The last refusal", err_a.message)
    end

    @testsection "a rate with no root is refused on its rate equation" begin
        # A sign: the rate jumps at the equilibrium, and `Δξ = Δt·r` has no
        # solution between the two branches.
        c = amounts(equilibrate(water(); model = DiluteSolutionModel()))[i_Ca]
        kss = KineticStepSolver(cs, DiluteSolutionModel(), [law((T, P, t, n, lna, n0) -> 1.0e-5 * sign(c - n["Ca+2"]))])
        err = try
            kinetic_step(kss, water(), 100.0)
        catch e
            e
        end
        @test err isa KineticStepError && err.cause === :rate
        @test occursin("calcite dissolution", err.message)
    end

    @testsection "an equilibrium that does not certify is refused" begin
        kss = KineticStepSolver(cs, DiluteSolutionModel(), [law((T, P, t, n, lna, n0) -> 1.0e-6)]; maxit = 1)
        err = try
            kinetic_step(kss, water(), 10.0)
        catch e
            e
        end
        @test err isa KineticStepError && err.cause === :equilibrium
    end

    @testsection "a march that cannot meet its tolerance says so" begin
        kss = KineticStepSolver(cs, DiluteSolutionModel(), [law((T, P, t, n, lna, n0) -> 1.0e-6)])
        err = try
            kinetic_step_adaptive(kss, water(), 10.0; reltol = 0.0, abstol = 0.0, max_halvings = 1)
        catch e
            e
        end
        @test err isa KineticStepError && err.cause === :tolerance
    end

    @testsection "the extents found around the equilibrium of the reactions" begin
        # `_kinetic_step_nested`, the route of `:reactions` where the
        # simultaneous solve does not certify and the eliminated form does not
        # apply, against the eliminated form and the equilibrium.
        start = water()
        c = amounts(equilibrate(start; model = DiluteSolutionModel()))[i_Ca]
        relax = law((T, P, t, n, lna, n0) -> 1.0e-3 * (c - n["Ca+2"]))
        kss_r = KineticStepSolver(cs, DiluteSolutionModel(), [relax]; coupling = :reactions)
        kss_s = KineticStepSolver(cs, DiluteSolutionModel(), [relax]; coupling = :species)
        for Δt in (10.0, 1.0e9)
            st_r, cert_r, q_r = ChemistryLab._kinetic_step_nested(kss_r, start, Δt)
            q_s = Ref{Any}(nothing)
            st_s = kinetic_step(kss_s, start, Δt; parameters = q_s)
            @test cert_r.optimal
            @test q_r ≈ kss_r.M * q_s[] rtol = 1.0e-8
            @test amounts(st_r)[i_Cal] ≈ amounts(st_s)[i_Cal] rtol = 1.0e-8
            @test maximum(abs, A * amounts(st_r) .- A * amounts(start)) <= 1.0e-12
        end
        # a rate that does not vanish at zero, and an equilibrium that cannot
        # certify, refused by name
        kss_x = KineticStepSolver(cs, DiluteSolutionModel(), [law((T, P, t, n, lna, n0) -> 1.0e-3)]; coupling = :reactions)
        err = try
            ChemistryLab._kinetic_step_nested(kss_x, water(), 100.0)
        catch e
            e
        end
        @test err isa KineticStepError && err.cause === :exhausted
        kss_e = KineticStepSolver(cs, DiluteSolutionModel(), [relax]; coupling = :reactions, maxit = 1)
        err = try
            ChemistryLab._kinetic_step_nested(kss_e, water(), 10.0)
        catch e
            e
        end
        @test err isa KineticStepError && err.cause === :equilibrium
        # the `:reactions` step whose simultaneous solve does not certify falls
        # back on the eliminated form, the same step: the same answer
        err = try
            kinetic_step(kss_x, water(), 100.0)
        catch e
            e
        end
        @test err isa KineticStepError && err.cause === :exhausted
    end

    @testsection "a stiff rate the simultaneous solve does not certify" begin
        # `k (1 − (n_Ca/c)⁶)`: the simultaneous Newton method of `:reactions`
        # does not certify this step, pinned or free; it falls back on the
        # eliminated form, the same step, and the extents found around the
        # equilibrium give it too.
        start = water()
        c = amounts(equilibrate(start; model = DiluteSolutionModel()))[i_Ca]
        stiff = law((T, P, t, n, lna, n0) -> 1.0 * (1 - (n["Ca+2"] / c)^6))
        kss_r = KineticStepSolver(cs, DiluteSolutionModel(), [stiff]; coupling = :reactions)
        kss_s = KineticStepSolver(cs, DiluteSolutionModel(), [stiff]; coupling = :species)
        @test !ChemistryLab._kinetic_step_augmented(kss_r, start, 10.0; pin_minerals = true)[2].optimal
        @test !ChemistryLab._kinetic_step_augmented(kss_r, start, 10.0; pin_minerals = false)[2].optimal
        q_r, q_s, c_r = Ref{Any}(nothing), Ref{Any}(nothing), Ref{Any}(nothing)
        st_r = kinetic_step(kss_r, start, 10.0; parameters = q_r, certificate = c_r)
        st_s = kinetic_step(kss_s, start, 10.0; parameters = q_s)
        @test c_r[].optimal
        @test q_r[] ≈ kss_r.M * q_s[] rtol = 1.0e-8
        @test amounts(st_r) ≈ amounts(st_s) rtol = 1.0e-8
        st_n, cert_n, q_n = ChemistryLab._kinetic_step_nested(kss_r, start, 10.0)
        @test cert_n.optimal
        @test amounts(st_n)[i_Cal] ≈ amounts(st_s)[i_Cal] rtol = 1.0e-8
        # the adaptive march resolves the pinning of `:reactions` once
        st_a, h, _ = kinetic_step_adaptive(kss_r, start, 10.0)
        @test h > 0 && amounts(st_a)[i_Cal] < amounts(start)[i_Cal]
        @test maximum(abs, A * amounts(st_a) .- A * amounts(start)) <= 1.0e-12
    end

    @testsection "a rate law with corners: the simultaneous solve certifies" begin
        # Parrott–Killoh's rate is the least of three mechanisms, so it has
        # corners. On this clinker the Newton method on the extents of
        # `:species` stops on the first step, and the simultaneous solve of
        # `:reactions` certifies it.
        subs = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
        α_max, blaine = powers_alpha_max(0.4), 380.0u"m^2/kg"
        laws = Dict(
            "C3S" => parrott_killoh_avrami(PK84_PARAMS_C3S, "C3S"; α_max, blaine),
            "C2S" => parrott_killoh_avrami(PK84_PARAMS_C2S, "C2S"; α_max, blaine),
            "C3A" => parrott_killoh_avrami(PK84_PARAMS_C3A, "C3A"; α_max, blaine),
            "C4AF" => parrott_killoh_avrami(PK84_PARAMS_C4AF, "C4AF"; α_max, blaine),
        )
        species = speciation(
            subs, split("C3S C2S C3A C4AF Portlandite Jennite ettringite monosulphate12 C3AH6 C3FH6 H2O@");
            aggregate_state = [AS_AQUEOUS],
        )
        cs_k = ChemicalSystem(species, CEMDATA_PRIMARIES; kinetic_species = laws)
        st0 = ChemicalState(cs_k)
        for (name, kg) in (("C3S", 0.619), ("C2S", 0.165), ("C3A", 0.08), ("C4AF", 0.087), ("H2O@", 0.4))
            set_quantity!(st0, name, kg * u"kg")
        end
        es = EquilibriumSolver(cs_k, DiluteSolutionModel(), OptimaOptimizer())
        krs = KineticsProblem(cs_k, st0, (0.0, 1.0); equilibrium_solver = es).kinetic_reactions
        kss_r = KineticStepSolver(cs_k, DiluteSolutionModel(), krs)
        kss_s = KineticStepSolver(cs_k, DiluteSolutionModel(), krs; coupling = :species)
        Δt = 3456.0
        @test_throws KineticStepError kinetic_step(kss_s, st0, Δt; n_initial = st0)
        c_r, q_r = Ref{Any}(nothing), Ref{Any}(nothing)
        st_r = kinetic_step(kss_r, st0, Δt; n_initial = st0, certificate = c_r, parameters = q_r)
        @test c_r[].optimal
        # the reactivity rows: the extents are those of the amounts
        n_r, n_0 = amounts(st_r), amounts(st0)
        @test transpose(kss_r.K) * (n_r .- n_0) ≈ q_r[] rtol = 1.0e-8
        A_k = Matrix{Float64}(cs_k.SM.A)
        @test maximum(abs, A_k * (n_r .- n_0)) <= 1.0e-10 * maximum(abs, A_k * n_0)
    end

    @testsection "the adaptive step: its error grows as Δt², and the next step follows" begin
        # Calcite relaxing to the equilibrium calcium in about 10³ s. Richardson's
        # estimate measures the local error of an implicit Euler step, `O(Δt²)`:
        # halving the step divides it by four, and the suggested step is
        # `0.9·Δt/√e`.
        start = water()
        c = amounts(equilibrate(start; model = DiluteSolutionModel()))[i_Ca]
        kss = KineticStepSolver(cs, DiluteSolutionModel(), [law((T, P, t, n, lna, n0) -> 1.0e-3 * (c - n["Ca+2"]))])
        e1, e2 = Ref(0.0), Ref(0.0)
        _, h1, next1 = kinetic_step_adaptive(kss, start, 100.0; error_estimate = e1)
        _, h2, _ = kinetic_step_adaptive(kss, start, 50.0; error_estimate = e2)
        @test h1 == 100.0 && h2 == 50.0
        @test e1[] / e2[] ≈ 4 rtol = 0.1
        @test next1 ≈ h1 * clamp(0.9 / sqrt(e1[]), 0.2, 2.0)
    end

    @testsection "a step the domain refuses from its start, and a singular Newton step" begin
        # Calcite absent and a constant dissolution: no extent above zero is in
        # the domain.
        kss = KineticStepSolver(cs, DiluteSolutionModel(), [law((T, P, t, n, lna, n0) -> 1.0e-3)])
        err = try
            kinetic_step(kss, water(calcite = 0.0, calcium = 1.0e-3), 10.0)
        catch e
            e
        end
        @test err isa KineticStepError && err.cause === :exhausted
        # `r = −k n_Cal` with `k Δt = 1`: the Jacobian `1 − k Δt` of the rate
        # equation vanishes, Newton's method has no step, and the residual left
        # is reported.
        kss1 = KineticStepSolver(cs, DiluteSolutionModel(), [law((T, P, t, n, lna, n0) -> -0.01 * n["Cal"])])
        err = try
            kinetic_step(kss1, water(), 100.0)
        catch e
            e
        end
        @test err isa KineticStepError && err.cause === :rate
    end

    @testsection "the active set on a problem with a surface potential" begin
        # A diffuse layer, its potential a parameter of the problem: the active
        # set on the amounts and the potential together against the certified
        # dual solve of the same equilibrium.
        h2o, hp, oh, na, cl = (sp[s] for s in ("H2O@", "H+", "OH-", "Na+", "Cl-"))
        g0(v) = SymbolicFunc(v * u"J/mol")
        RT = ChemistryLab.R_GAS * 298.15
        free = Species("XsOH"; aggregate_state = AS_SURFACE, class = SC_SURFCOMPLEX)
        free[:ΔₐG⁰] = g0(0.0)
        prot = Species("XsOH2+"; aggregate_state = AS_SURFACE, class = SC_SURFCOMPLEX)
        prot[:ΔₐG⁰] = g0(-RT * 7.29 * log(10.0))
        depr = Species("XsO-"; aggregate_state = AS_SURFACE, class = SC_SURFCOMPLEX)
        depr[:ΔₐG⁰] = g0(RT * 8.93 * log(10.0))
        family = SiteFamily(
            "Xs", free, [prot, depr];
            capacity = TotalSiteAmount(2.0e-4), support = SurfaceSupport("hfo", nothing, FixedSurfaceArea(600.0)),
            model = DiffuseLayer(; area = 600.0),
        )
        ds = ChemicalSystem([h2o, hp, oh, na, cl, free, prot, depr], [h2o, hp, na, cl, free]; site_families = [family])
        st = ChemicalState(ds, [55.508u"mol", 1.0e-7u"mol", 1.0e-7u"mol", 1.0e-2u"mol", 1.0e-2u"mol", 2.0e-4u"mol", 1.0e-12u"mol", 1.0e-12u"mol"])
        des = DualEquilibriumSolver(ds, DiluteSolutionModel())
        p = ChemistryLab._build_params(st)
        n0 = amounts(st)
        b = des.A * n0
        surf = ChemistryLab._surface_potential_blocks(des, st, p, n0)
        bl = ChemistryLab._with_surface(ChemistryLab._constraint_blocks(FixedTP(), des, st, p, n0), surf, 0)
        @test bl.nq == 1 && size(bl.Aq) == (size(des.A, 1), 1)
        # beside the extents of one reaction, the rows of the reaction pad it
        kin = (; nq = 1, gq = nothing, hq = nothing, pq = nothing, cq = (x, q, pp) -> q, Aq = zeros(size(des.A, 1) + 1, 1), q0 = [0.0], qscale = [1.0], apply = (T, P, q) -> (T, P))
        both = ChemistryLab._with_surface(kin, surf, 1)
        @test both.nq == 2 && size(both.Aq) == (size(des.A, 1) + 1, 2)
        prob = ChemistryLab._dual_problem(des, p, n0, bl)
        r = ChemistryLab._optima_parametric_solve(prob, b, max.(n0, 1.0e-16), bl.q0, fill(1.0e-16, length(n0)), des.opts)
        @test r.converged
        eq, cert = solve_certified(des, (st,); b = b)
        @test cert.optimal
        @test r.x[6:8] ≈ amounts(eq)[6:8] rtol = 1.0e-6
        # and without OptimaSolver, what it lacks is named
        @test_throws ErrorException ChemistryLab._optima_parametric_solve()
    end

    @testsection "extents projected back into the domain" begin
        # n = n₀ + K q: the extent that takes the amount below zero is moved to
        # the boundary, the others left alone.
        @test ChemistryLab._project_extents([0.5, 0.1], [0.2, 1.0], [-1.0 0.0; 0.0 -1.0]) ≈ [0.2, 0.1]
        @test ChemistryLab._project_extents([0.1], [0.2], reshape([-1.0], 1, 1)) == [0.1]
    end

    @testsection "a potential carried as a dual number is seen by the promotion" begin
        D = ForwardDiff.Dual{Nothing, Float64, 1}
        @test ChemistryLab._number_type_t(Vector{ChemistryLab._PotentialSlot{D}}) === D
        @test ChemistryLab._number_type_t(Vector{ChemistryLab._PotentialSlot{Float64}}) === Float64
        @test ChemistryLab._number_type_t(Vector{D}) === D
        @test ChemistryLab._number_type_t(Matrix{Float64}) === Float64
        @test ChemistryLab._number_type_t(Vector{Any}) === Float64
        # Parameters as an activity model receives them: the answer is a constant
        # of their types, so asking costs nothing.
        params = (; T = 298.15, ψ_site = ChemistryLab._PotentialSlot{D}[nothing], b = [1.0, 2.0])
        @test @inferred(ChemistryLab._number_type_of(params)) === D
        @test (@allocated ChemistryLab._number_type_of(params)) == 0
    end
end
