# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using ChemistryLab
using DynamicQuantities
using ForwardDiff
using OrderedCollections
using OrdinaryDiffEq
using Test

include("../reference_species.jl")

# A slow surface: a site family whose states change only through declared
# reactions, under `sorption_rate`. The references are closed forms and the
# certified equilibrium; nothing is fitted.
#
# The sorbent sits in a solution saturated with portlandite, which holds the
# activities of Ca+2 and H+ whatever the surface takes, as long as the sites are
# few. Under ideal site mixing the occupancy of
#
#     ≡XOH + Ca+2 → ≡XOCa+ + H+,     r = k n(≡XOH) a(Ca+2) (1 − Ω)
#
# then relaxes exponentially: dn_c/dt = k a_Ca (N − n_c) − (k/K) a_H n_c, so
# n_c(t) = n_eq (1 − e^{−λt}) with λ = k (a_Ca + a_H/K) and n_eq = k a_Ca N/λ.

const _SORB_RT = ChemistryLab.R_GAS * 298.15

function _sorb_system(;
        logK = -11.0, n_slow = 1.0e-9, n_fast = 0.0, model = IdealSiteMixing(),
        fast_model = IdealSiteMixing(), shared_support = false, extra_slow_state = false,
        logK_extra = log10(25.0)
    )
    g(v) = SymbolicFunc(v * u"J/mol")
    sf(sym, v) = (s = Species(sym; aggregate_state = AS_SURFACE, class = SC_SURFCOMPLEX); s[:ΔₐG⁰] = g(v); s)
    h2o, hp, oh, ca = reference_species(("H2O@", "H+", "OH-", "Ca+2"))
    pt = reference_species("Portlandite"; db = :cemdata18)
    G_CA = ustrip(us"J/mol", ca[:ΔₐG⁰](T = 298.15u"K", P = 1.0e5u"Pa"; unit = true))
    # ≡XOH + Ca+2 = ≡XOCa+ + H+: G(≡XOCa+) = −RT ln K + G(Ca+2), G(≡XOH) = G(H+) = 0.
    G_bound = -_SORB_RT * log(10.0^logK) + G_CA

    species = Any[h2o, hp, oh, ca, pt]
    primaries = Any[h2o, hp, ca]
    families = SiteFamily[]
    sup_s = SurfaceSupport("slow sorbent", nothing, FixedSurfaceArea(1.0))
    sup_f = shared_support ? sup_s : SurfaceSupport("fast sorbent", nothing, FixedSurfaceArea(1.0))
    if n_slow > 0
        free, bound = sf("XsOH", 0.0), sf("XsOCa+", G_bound)
        states = [bound]
        # A third state, ≡XOH + Ca+2 = ≡XOHCa+2, its constant chosen so that it
        # holds a share of the sites comparable to ≡XOCa+ in this solution.
        extra_slow_state && push!(states, sf("XsOHCa+2", -_SORB_RT * log(10.0^logK_extra) + G_CA))
        push!(species, free, states...)
        push!(primaries, free)
        push!(families, SiteFamily("Xs", free, states; capacity = TotalSiteAmount(n_slow * u"mol"), support = sup_s, model))
    end
    if n_fast > 0
        free, bound = sf("XfOH", 0.0), sf("XfOCa+", G_bound)
        push!(species, free, bound)
        push!(primaries, free)
        push!(families, SiteFamily("Xf", free, [bound]; capacity = TotalSiteAmount(n_fast * u"mol"), support = sup_f, model = fast_model))
    end
    cs = ChemicalSystem([species...], [primaries...]; site_families = families)
    nm = symbol.(cs.species)
    idx(s) = findfirst(==(s), nm)
    n0 = Any[fill(0.0u"mol", length(nm))...]
    n0[idx("H2O@")] = moles_of_water() * u"mol"
    n0[idx("Portlandite")] = 0.1u"mol"
    n_slow > 0 && (n0[idx("XsOH")] = n_slow * u"mol")
    n_fast > 0 && (n0[idx("XfOH")] = n_fast * u"mol")
    return cs, ChemicalState(cs, n0), idx
end

function _sorb_reaction(cs, idx; prefix = "Xs", site_first = true)
    sp(s) = cs.species[idx(s)]
    reactants = site_first ? OrderedDict(sp("$(prefix)OH") => 1, sp("Ca+2") => 1) :
        OrderedDict(sp("Ca+2") => 1, sp("$(prefix)OH") => 1)
    return Reaction(
        reactants, OrderedDict(sp("$(prefix)OCa+") => 1, sp("H+") => 1);
        symbol = "Ca on $prefix", equal_sign = '→'
    )
end

function _sorb_run(cs, state, kr, tend; saveat, reltol = 1.0e-10)
    kp = KineticsProblem(
        cs, [kr], state, (0.0, tend);
        equilibrium_solver = EquilibriumSolver(cs, DiluteSolutionModel(), OptimaOptimizer())
    )
    ks = KineticsSolver(; ode_solver = Rodas5P(), reltol, abstol = 1.0e-20, saveat)
    return kp, integrate(kp, ks)
end

# The activities the surface sees, from the first speciated state of a run.
function _sorb_activities(st, idx)
    n = Float64[ustrip(us"mol", x) for x in st.n]
    des = DualEquilibriumSolver(st.system, DiluteSolutionModel())
    lna = des.lna(n, ChemistryLab._build_params(st))
    return exp(lna[idx("Ca+2")]), exp(lna[idx("H+")])
end

@testsection "A slow surface: sorption_rate and kinetic site families" begin

    logK = -11.0
    K = 10.0^logK

    @testset "the whole family goes to the kinetic side, the free site with it" begin
        cs, st, idx = _sorb_system()
        rxn = _sorb_reaction(cs, idx)
        kr = KineticReaction(cs, rxn, sorption_rate(0.02, cs, rxn))
        @test kr.idx_mineral == idx("XsOH")
        kp = KineticsProblem(
            cs, [kr], st, (0.0, 1.0);
            equilibrium_solver = EquilibriumSolver(cs, DiluteSolutionModel(), OptimaOptimizer())
        )
        @test sort(kp.idx_kinetic) == sort([idx("XsOH"), idx("XsOCa+")])
        # The site row leaves the equilibrium partition with the family.
        sub = ChemistryLab._equilibrium_subsystem(cs, kp.idx_equilibrium)
        @test sub.site_families === nothing
        # And the implicit step classifies it the same way.
        K_mat = ChemistryLab._reactivity_matrix([kr], cs, :auto)
        @test Set(i for i in axes(K_mat, 1) if any(!iszero, K_mat[i, :])) ==
            Set([idx("XsOH"), idx("XsOCa+")])
    end

    @testset "(a) kinetic Langmuir in a buffered solution, against its closed form" begin
        N, k = 1.0e-9, 0.02
        cs, st, idx = _sorb_system(; logK, n_slow = N)
        rxn = _sorb_reaction(cs, idx)
        kr = KineticReaction(cs, rxn, sorption_rate(k, cs, rxn))
        ts = [0.0, 300.0, 1000.0, 3000.0, 6000.0]
        kp, sol = _sorb_run(cs, st, kr, ts[end]; saveat = ts)
        states = speciated_states(sol, kp)
        a_ca, a_h = _sorb_activities(states[1], idx)
        λ = k * (a_ca + a_h / K)
        n_eq = k * a_ca * N / λ
        bound = Float64[ustrip(us"mol", s.n[idx("XsOCa+")]) for s in states]
        free = Float64[ustrip(us"mol", s.n[idx("XsOH")]) for s in states]
        for (t, b, f) in zip(ts, bound, free)
            @test b ≈ n_eq * (1 - exp(-λ * t)) atol = 1.0e-6 * N
            @test b + f ≈ N rtol = 1.0e-10
        end
        # The constants of the case are what make it a test: the occupancy is
        # partial at equilibrium, and the run covers several relaxation times.
        @test 0.2 < n_eq / N < 0.8
        @test λ * ts[end] > 5
    end

    @testset "(b) the long-time limit is the certified equilibrium" begin
        # Many sites this time, two percent of the calcium: the solution moves,
        # so the limit is no longer the closed form, only the equilibrium.
        N = 2.0e-4
        cs, st, idx = _sorb_system(; logK, n_slow = N)
        rxn = _sorb_reaction(cs, idx)
        kr = KineticReaction(cs, rxn, sorption_rate(0.05, cs, rxn))
        kp, sol = _sorb_run(cs, st, kr, 2.0e4; saveat = [0.0, 2.0e4])
        last = speciated_states(sol, kp)[end]
        eq, cert = equilibrate_certified(st; model = DiluteSolutionModel())
        @test cert.optimal
        for s in ("XsOCa+", "XsOH", "Ca+2", "H+", "Portlandite")
            @test ustrip(us"mol", last.n[idx(s)]) ≈ ustrip(us"mol", eq.n[idx(s)]) rtol = 1.0e-8
        end
    end

    @testset "(c) a fast and a slow family reach the equilibrium of both" begin
        N = 1.0e-4
        cs, st, idx = _sorb_system(; logK, n_slow = N, n_fast = N)
        rxn = _sorb_reaction(cs, idx)
        kr = KineticReaction(cs, rxn, sorption_rate(0.05, cs, rxn))
        kp, sol = _sorb_run(cs, st, kr, 2.0e4; saveat = [0.0, 50.0, 2.0e4])
        @test !(idx("XfOH") in kp.idx_kinetic)
        states = speciated_states(sol, kp)
        # The fast family is at equilibrium with the solution at every instant,
        # the slow one only at the end.
        des = DualEquilibriumSolver(cs, DiluteSolutionModel())
        gap(st, pre) = begin
            n = Float64[ustrip(us"mol", x) for x in st.n]
            lna = des.lna(n, ChemistryLab._build_params(st))
            lna[idx("$(pre)OCa+")] + lna[idx("H+")] - lna[idx("$(pre)OH")] - lna[idx("Ca+2")] - logK * log(10)
        end
        for s in states
            @test abs(gap(s, "Xf")) < 1.0e-6
        end
        @test abs(gap(states[2], "Xs")) > 0.1
        @test abs(gap(states[end], "Xs")) < 1.0e-6
        eq, cert = equilibrate_certified(st; model = DiluteSolutionModel())
        @test cert.optimal
        for s in ("XsOCa+", "XfOCa+", "Ca+2")
            @test ustrip(us"mol", states[end].n[idx(s)]) ≈ ustrip(us"mol", eq.n[idx(s)]) rtol = 1.0e-8
        end
    end

    @testset "(d) the implicit step follows its own discrete closed form" begin
        # Backward Euler on a linear relaxation: n_{j+1} = (n_j + Δt λ n_eq)/(1 + λΔt).
        # A thousandth of a millimole of sites: the solution moves by a part in 10⁴.
        N, k, Δt = 1.0e-6, 0.02, 20.0
        cs, st, idx = _sorb_system(; logK, n_slow = N)
        rxn = _sorb_reaction(cs, idx)
        kr = KineticReaction(cs, rxn, sorption_rate(k, cs, rxn))
        kss = KineticStepSolver(cs, DiluteSolutionModel(), [kr])
        # The bare sorbent in the solution the run starts from: the speciated
        # state of the ODE route at t = 0.
        kp0, sol0 = _sorb_run(cs, st, kr, 1.0; saveat = [0.0])
        cur = speciated_states(sol0, kp0)[1]
        a_ca, a_h = _sorb_activities(cur, idx)
        λ = k * (a_ca + a_h / K)
        n_eq = k * a_ca * N / λ
        expected = 0.0
        for j in 1:10
            cur = kinetic_step(kss, cur, Δt; t = (j - 1) * Δt)
            expected = (expected + Δt * λ * n_eq) / (1 + λ * Δt)
            @test ustrip(us"mol", cur.n[idx("XsOCa+")]) ≈ expected rtol = 1.0e-3
        end
    end

    @testset "(e) the derivative with respect to k is the closed form's" begin
        # n_c(t) = n_eq (1 − e^{−k s t}) with s = a_Ca + a_H/K, and n_eq free of k:
        # dn_c/dk = n_eq s t e^{−k s t}.
        N, k0, t1 = 1.0e-9, 0.02, 500.0
        cs, st, idx = _sorb_system(; logK, n_slow = N)
        rxn = _sorb_reaction(cs, idx)
        run_k(k) = _sorb_run(cs, st, KineticReaction(cs, rxn, sorption_rate(k, cs, rxn)), t1; saveat = [0.0, t1])
        function bound_at(k)
            kp, sol = run_k(k)
            j = findfirst(==(idx("XsOCa+")), kp.idx_kinetic)
            return sol.u[end][sol.prob.p.n_be + j]
        end
        kp0, sol0 = run_k(k0)
        a_ca, a_h = _sorb_activities(speciated_states(sol0, kp0)[1], idx)
        s_ = a_ca + a_h / K
        n_eq = a_ca * N / s_
        @test bound_at(k0) ≈ n_eq * (1 - exp(-k0 * s_ * t1)) atol = 1.0e-6 * N
        @test ForwardDiff.derivative(bound_at, k0) ≈ n_eq * s_ * t1 * exp(-k0 * s_ * t1) rtol = 1.0e-5
    end

    @testset "a trace site budget is computed at equilibrium" begin
        # Down to 1e-10 mol of sites beside a kilogram of water; until OptimaSolver
        # 0.8.2 the total of the family was seeded at no less than 1e-6 and the
        # search did not come back down to it (4e-9 mol of sites for 1e-9).
        for N in (1.0e-10, 1.0e-9)
            cs, st, idx = _sorb_system(; logK, n_slow = N)
            eq, cert = equilibrate_certified(st; model = DiluteSolutionModel())
            @test cert.optimal
            n(s) = ustrip(us"mol", eq.n[idx(s)])
            @test n("XsOH") + n("XsOCa+") ≈ N rtol = 1.0e-10
            la = log_activities(eq, DiluteSolutionModel())
            @test la["XsOCa+"] + la["H+"] - la["XsOH"] - la["Ca+2"] ≈ logK * log(10) atol = 1.0e-6
        end
    end

    @testset "what cannot be slow is refused, by name" begin
        cs, st, idx = _sorb_system()
        rxn = _sorb_reaction(cs, idx)
        @test_throws "Write the site species first" sorption_rate(0.02, cs, _sorb_reaction(cs, idx; site_first = false))
        @test_throws "non-negative" sorption_rate(-1.0, cs, rxn)
        ca, oh = cs.species[idx("Ca+2")], cs.species[idx("OH-")]
        pt = cs.species[idx("Portlandite")]
        no_site = Reaction(OrderedDict(pt => 1), OrderedDict(ca => 1, oh => 2); symbol = "pt", equal_sign = '→')
        @test_throws "one surface species on each side" sorption_rate(0.02, cs, no_site)

        # A state of a slow family that no reaction forms would keep its amount.
        cs3, st3, idx3 = _sorb_system(; extra_slow_state = true)
        rxn3 = _sorb_reaction(cs3, idx3)
        kr3 = KineticReaction(cs3, rxn3, sorption_rate(0.02, cs3, rxn3))
        @test_throws "XsOHCa+2" KineticsProblem(cs3, [kr3], st3, (0.0, 1.0))

        # A charged slow family and a fast one on one support would share a
        # potential the equilibrium solve cannot compute.
        ccm = ConstantCapacitance(; C = 1.0, area = 1.0)
        csE, stE, idxE = _sorb_system(; n_fast = 1.0e-9, model = ccm, fast_model = ccm, shared_support = true)
        rxnE = _sorb_reaction(csE, idxE)
        krE = KineticReaction(csE, rxnE, sorption_rate(0.02, csE, rxnE))
        @test_throws "separate supports" KineticsProblem(csE, [krE], stE, (0.0, 1.0))
        # On separate supports, it is accepted.
        csS, stS, idxS = _sorb_system(; n_fast = 1.0e-9, model = ccm, fast_model = ccm)
        rxnS = _sorb_reaction(csS, idxS)
        krS = KineticReaction(csS, rxnS, sorption_rate(0.02, csS, rxnS))
        @test KineticsProblem(csS, [krS], stS, (0.0, 1.0)) isa KineticsProblem
    end

    @testset "the second law: dissipation tells a consistent law from one that is not" begin
        # The same mass-action law twice, once with its reverse rate tied to the
        # constant of the database, once with a reverse constant ten times too
        # small: the second stops where Ω = 10, and from Ω = 1 on it runs against
        # its affinity.
        N, k = 2.0e-4, 0.05
        cs, st, idx = _sorb_system(; logK, n_slow = N)
        rxn = _sorb_reaction(cs, idx)
        good = KineticReaction(cs, rxn, sorption_rate(k, cs, rxn))
        bad_law = (T, P, t, n, lna, n0) ->
        k * n["XsOH"] * exp(lna["Ca+2"]) - (k / (10K)) * n["XsOCa+"] * exp(lna["H+"])
        bad = KineticReaction(cs, rxn, bad_law)
        # Ω passes one near 600 s under the inconsistent law: instants on both
        # sides of it, and at rest.
        ts = [0.0, 20.0, 50.0, 100.0, 300.0, 1000.0, 2000.0, 4000.0, 2.0e4]
        for (kr, consistent) in ((good, true), (bad, false))
            kp, sol = _sorb_run(cs, st, kr, ts[end]; saveat = ts)
            d = dissipation(sol, kp)
            @test size(d.affinity) == (length(ts), 1)
            @test d.affinity == reaction_affinities(sol, kp)
            if consistent
                @test isempty(d.violations)
                @test all(>=(0), d.entropy_production)
                @test abs(d.affinity[end, 1]) < 1.0e-6 * ChemistryLab.R_GAS * 298.15
            else
                @test !isempty(d.violations)
                @test all(v -> v.affinity * v.rate < 0, d.violations)
                # Where the inconsistent law stops: Ω = 10, i.e. 𝒜 = −RT ln 10.
                @test d.affinity[end, 1] ≈ -ChemistryLab.R_GAS * 298.15 * log(10) rtol = 1.0e-4
            end
        end
    end

    @testset "(f) a cycle of three reactions carries no flux at rest" begin
        # Three states on one family, three reactions closing a cycle. Each law
        # in (1 − Ω) stops at the equilibrium of the database, so all three rates
        # vanish together; with one reverse rate off by ten, the cycle keeps
        # turning at rest, and the dissipation sees it.
        N, k = 1.0e-4, 0.05
        cs, st, idx = _sorb_system(; logK, n_slow = N, extra_slow_state = true)
        sp(s) = cs.species[idx(s)]
        r1 = _sorb_reaction(cs, idx)
        r2 = Reaction(
            OrderedDict(sp("XsOH") => 1, sp("Ca+2") => 1), OrderedDict(sp("XsOHCa+2") => 1);
            symbol = "XsOH to XsOHCa+2", equal_sign = '→'
        )
        r3 = Reaction(
            OrderedDict(sp("XsOHCa+2") => 1), OrderedDict(sp("XsOCa+") => 1, sp("H+") => 1);
            symbol = "XsOHCa+2 to XsOCa+", equal_sign = '→'
        )
        consistent = [KineticReaction(cs, r, sorption_rate(k, cs, r)) for r in (r1, r2, r3)]
        # ≡XOHCa+2 → ≡XOCa+ + H+: its constant from the energies built above.
        K3 = exp(
            -(
                ustrip(us"J/mol", sp("XsOCa+")[:ΔₐG⁰](T = 298.15u"K", P = 1.0e5u"Pa"; unit = true)) -
                    ustrip(us"J/mol", sp("XsOHCa+2")[:ΔₐG⁰](T = 298.15u"K", P = 1.0e5u"Pa"; unit = true))
            ) / _SORB_RT
        )
        off = (T, P, t, n, lna, n0) -> k * n["XsOHCa+2"] - (k / (K3 / 10)) * n["XsOCa+"] * exp(lna["H+"])
        inconsistent = [consistent[1], consistent[2], KineticReaction(cs, r3, off)]
        tend = 5.0e4
        for (krs, ok) in ((consistent, true), (inconsistent, false))
            kp = KineticsProblem(
                cs, krs, st, (0.0, tend);
                equilibrium_solver = EquilibriumSolver(cs, DiluteSolutionModel(), OptimaOptimizer())
            )
            sol = integrate(kp, KineticsSolver(; ode_solver = Rodas5P(), reltol = 1.0e-10, abstol = 1.0e-20, saveat = [0.0, tend]))
            d = dissipation(sol, kp)
            flux = maximum(abs, d.rate[end, :])
            if ok
                @test flux < 1.0e-10 * k * N
                @test isempty(d.violations)
            else
                @test flux > 1.0e-3 * k * N
                @test !isempty(d.violations)
            end
        end
    end
end

# A slow surface on a host that dissolves under its own rate law: the sites of
# portlandite, `AreaSiteDensity` over a specific area, their budget following the
# kinetic host, and calcium binding to them slowly,
#
#     ≡Xs + Ca+2 + H₂O → ≡XsOHCa+ + H+,     log K = −9.2,
#
# the free site bare, as a budget following a kinetic host requires.
function _slow_host_system(; n_host0 = 1.0e-3, Γ = 1.0e-5, ssa = 1.0e5, area = BETSurfaceArea(ssa * u"m^2/kg"))
    g(v) = SymbolicFunc(v * u"J/mol")
    sf(sym, v) = (s = Species(sym; aggregate_state = AS_SURFACE, class = SC_SURFCOMPLEX); s[:ΔₐG⁰] = g(v); s)
    h2o, hp, oh, ca = reference_species(("H2O@", "H+", "OH-", "Ca+2"))
    pt = reference_species("Portlandite"; db = :cemdata18)
    G(s) = ustrip(us"J/mol", s[:ΔₐG⁰](T = 298.15u"K", P = 1.0e5u"Pa"; unit = true))
    free, bound = sf("Xs", 0.0), sf("XsOHCa+", -_SORB_RT * log(10.0^-9.2) + G(h2o) + G(ca))
    family = SiteFamily(
        "Xs", free, [bound];
        capacity = AreaSiteDensity(Γ * u"mol/m^2"),
        support = SurfaceSupport("portlandite", "Portlandite", area; coupling = SITES_FOLLOW_HOST),
    )
    cs = ChemicalSystem(
        AbstractSpecies[h2o, hp, oh, ca, pt, free, bound], AbstractSpecies[h2o, hp, ca, free];
        site_families = [family],
    )
    idx(s) = findfirst(==(s), symbol.(cs.species))
    st = ChemicalState(cs)
    set_quantity!(st, "H2O@", moles_of_water() * u"mol")
    set_quantity!(st, "Portlandite", n_host0 * u"mol")
    set_quantity!(st, "Ca+2", 1.0e-4u"mol")
    set_quantity!(st, "OH-", 2.0e-4u"mol")
    return cs, host_consistent_state(st), idx, family
end

# The host dissolving at `k_d n_host`: its amount is n₀ e^{−k_d t}, exactly.
function _host_dissolution(cs, idx; k_d)
    sp(s) = cs.species[idx(s)]
    rxn = Reaction(
        OrderedDict(sp("Portlandite") => 1), OrderedDict(sp("Ca+2") => 1, sp("OH-") => 2);
        symbol = "portlandite", equal_sign = '→',
    )
    st = Float64[get(Dict("Portlandite" => -1.0, "Ca+2" => 1.0, "OH-" => 2.0), symbol(s), 0.0) for s in cs.species]
    return KineticReaction(rxn, (T, P, t, n, lna, n0) -> k_d * n["Portlandite"], idx("Portlandite"), st)
end

function _slow_sorption(cs, idx, k)
    sp(s) = cs.species[idx(s)]
    rxn = Reaction(
        OrderedDict(sp("Xs") => 1, sp("Ca+2") => 1, sp("H2O@") => 1),
        OrderedDict(sp("XsOHCa+") => 1, sp("H+") => 1);
        symbol = "Ca on Xs", equal_sign = '→',
    )
    return KineticReaction(cs, rxn, sorption_rate(k, cs, rxn))
end

function _slow_host_run(cs, st, krs, ts)
    kp = KineticsProblem(
        cs, krs, st, (0.0, ts[end]);
        equilibrium_solver = EquilibriumSolver(cs, DiluteSolutionModel(), OptimaOptimizer()),
    )
    sol = integrate(kp, KineticsSolver(; ode_solver = Rodas5P(), reltol = 1.0e-10, abstol = 1.0e-20, saveat = ts))
    return kp, sol, speciated_states(sol, kp)
end

@testsection "A slow surface on a kinetic host" begin
    cs, st, idx, family = _slow_host_system()
    mol(s, x) = ustrip(us"mol", s.n[idx(x)])
    M = ustrip(us"kg/mol", cs.species[idx("Portlandite")][:M])
    N(n_host) = 1.0e-5 * 1.0e5 * M * n_host
    k_d = 2.0e-4
    ts = [0.0, 1000.0, 3000.0, 6000.0]

    @testset "matter and charge are conserved, and the members follow the budget" begin
        kp, sol, states = _slow_host_run(cs, st, [_host_dissolution(cs, idx; k_d), _slow_sorption(cs, idx, 0.05)], ts)
        @test sort(kp.idx_kinetic) == sort([idx("Portlandite"), idx("Xs"), idx("XsOHCa+")])
        host = [mol(s, "Portlandite") for s in states]
        @test host ≈ 1.0e-3 .* exp.(-k_d .* ts) rtol = 1.0e-7
        @test [mol(s, "Xs") + mol(s, "XsOHCa+") for s in states] ≈ N.(host) rtol = 1.0e-8
        ca(s) = mol(s, "Portlandite") + mol(s, "Ca+2") + mol(s, "XsOHCa+")
        z = [Float64(charge(sp)) for sp in cs.species]
        q(s) = sum(z .* [ustrip(us"mol", x) for x in s.n])
        @test maximum(abs(ca(s) - ca(st)) for s in states) < 1.0e-12
        @test maximum(abs(q(s) - q(st)) for s in states) < 1.0e-12
        @test mol(states[end], "XsOHCa+") > 0
        # The extents no longer hold the members, and the replay is refused.
        @test ChemistryLab.extent_residual(sol, kp, reaction_extents(sol, kp)) < 1.0e-12
        @test_throws "speciated_states" state_at(sol, kp, ts[end])
    end

    @testset "as the rate grows, the family at equilibrium" begin
        # The same host with the family left at equilibrium (a budget following a
        # kinetic host, the fast route), against the slow family at k and 10k:
        # the lag of the occupancy is one over the rate constant.
        _, _, fast = _slow_host_run(cs, st, [_host_dissolution(cs, idx; k_d)], ts)
        gap(k) = begin
            _, _, slow = _slow_host_run(cs, st, [_host_dissolution(cs, idx; k_d), _slow_sorption(cs, idx, k)], ts)
            abs(mol(slow[end], "XsOHCa+") - mol(fast[end], "XsOHCa+"))
        end
        g1, g10 = gap(0.05), gap(0.5)
        @test 0.05 < g10 / g1 < 0.2
    end

    @testset "everything bound returns when the host is gone" begin
        tend = 25 / k_d                         # e^{−25}: 1.4e-11 of the host left
        _, _, states = _slow_host_run(cs, st, [_host_dissolution(cs, idx; k_d), _slow_sorption(cs, idx, 0.05)], [0.0, tend])
        @test mol(states[end], "Xs") + mol(states[end], "XsOHCa+") < 1.0e-9 * N(1.0e-3)
        @test mol(states[end], "Ca+2") ≈ mol(st, "Ca+2") + mol(st, "XsOHCa+") + 1.0e-3 rtol = 1.0e-8
    end

    @testset "the budget term, differentiated" begin
        # dnₖ/dt gains Ṅ nₖ/N when sites vanish: its Jacobian in the members is
        # Ṅ (δₖₘ/N − nₖ/N²), and what returns to bₑ is minus that times the
        # column of what each member holds.
        idx_kin = [idx("Portlandite"), idx("Xs"), idx("XsOHCa+")]
        idx_eq = setdiff(eachindex(cs.species), idx_kin)
        Ae = Float64.(ChemistryLab._constraint_matrix(ChemistryLab._equilibrium_subsystem(cs, idx_eq)))
        sh = only(ChemistryLab._kinetic_site_hosts(cs, idx_kin, idx_eq, Ae))
        @test sh isa ChemistryLab._SlowSiteHost
        n_be = size(Ae, 1)
        n_host, rate = 6.0e-4, -1.0e-7
        p = (; idx_kinetic = idx_kin, n_initial_full = fill(1.0e-3, length(cs.species)), n_be)
        function term(nm)
            n_full = zeros(eltype(nm), length(cs.species))
            n_full[idx("Portlandite")] = n_host
            n_full[idx("Xs")], n_full[idx("XsOHCa+")] = nm
            du = zeros(eltype(nm), n_be + 3)
            ChemistryLab._move_site_budget!(du, sh, p, n_full, [rate, 0.0, 0.0])
            return du
        end
        nm = [2.0e-5, 2.4e-5]
        Ndot = 1.0e-5 * 1.0e5 * M * rate
        Nm = sum(nm)
        J = ForwardDiff.jacobian(term, nm)
        analytic = [Ndot * ((k == m) / Nm - nm[k] / Nm^2) for k in 1:2, m in 1:2]
        @test J[n_be .+ (2:3), :] ≈ analytic rtol = 1.0e-12
        @test J[1:n_be, :] ≈ -sh.content * analytic rtol = 1.0e-12
        # The sites born with a growing host are free, and hold nothing.
        function grow(nm)
            n_full = zeros(eltype(nm), length(cs.species))
            n_full[idx("Portlandite")] = n_host
            n_full[idx("Xs")], n_full[idx("XsOHCa+")] = nm
            du = zeros(eltype(nm), n_be + 3)
            ChemistryLab._move_site_budget!(du, sh, p, n_full, [-rate, 0.0, 0.0])
            return du
        end
        @test grow(nm)[n_be .+ (2:3)] ≈ [-Ndot, 0.0] rtol = 1.0e-14
        @test all(iszero, grow(nm)[1:n_be])
    end

    @testset "what a slow surface on a host still refuses" begin
        # A host the equilibrium moves.
        @test_throws "no rate law controls" KineticsProblem(
            cs, [_slow_sorption(cs, idx, 0.05)], st, (0.0, 1.0);
            equilibrium_solver = EquilibriumSolver(cs, DiluteSolutionModel(), OptimaOptimizer()),
        )
        # A member holding a primary the partition does not keep: what it holds
        # could not be returned there when its site vanishes.
        idx_kin = [idx("Portlandite"), idx("Xs"), idx("XsOHCa+")]
        rows = ChemistryLab._partition_rows(cs, setdiff(eachindex(cs.species), idx_kin))
        no_ca = filter(r -> ChemistryLab.symbol(cs.SM.primaries[r]) != "Ca+2", rows)
        @test_throws "not in the equilibrium partition" ChemistryLab._member_content(cs, family, no_ca)
        # The implicit step, whose reactivity constraints cannot hold a member
        # that moves with its host.
        @test_throws "integrate" KineticStepSolver(
            cs, DiluteSolutionModel(), [_host_dissolution(cs, idx; k_d), _slow_sorption(cs, idx, 0.05)],
        )
    end
end
