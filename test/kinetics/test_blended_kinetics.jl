# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# A blended cement on the coupled kinetic path: the CEM I 52.5 N and slag paste
# of Gruyaert et al. (2010) with 50 % slag, its clinker under Parrott–Killoh,
# its slag glass under the Waller law, the equilibrium solved along the run
# (`scripts/gruyaert2010_kinetics.jl`). The 28-month comparison with the article
# is the page `examples/blended_slag_kinetics.md`; this file holds the run to
# what it must satisfy whatever the data, over one week.

isdefined(@__MODULE__, :gruyaert_run) || include(joinpath(pkgdir(ChemistryLab), "scripts", "gruyaert2010_kinetics.jl"))

@testsection "a blended cement on the coupled kinetic path" begin
    day = 86400.0
    τ = gruyaert_slag_time()
    n = WALLER_PARAMS_FLY_ASH.n
    α(t) = 1 / (1 + (ustrip(us"d", τ) / t)^n)
    # The calibration is the closed form of the law solved for τ: it gives back
    # the one measurement it was set on, and nothing else was fitted.
    @test α(852.0) ≈ only(literature_table("Gruyaert2010", "hydration_degree_slag"; age_days = 852u"d", slag_to_binder = 0.5).alpha_slag) / 100 rtol = 1.0e-12

    # The slag's rate carries the apparent activation energy the article fits
    # at the paste's cement-to-binder ratio (Eq. 3): its Arrhenius ratio between
    # two temperatures is that of E_S, exactly.
    G10q(k) = ustrip(us"J/mol", literature_value("Gruyaert2010", k))
    E_S = G10q("activation_energy_slag_slope") * 0.5 + G10q("activation_energy_slag_intercept")
    w = gruyaert_rates(; slag = 0.5, τ_slag = τ)["BFS"]
    ix = Dict("BFS" => 1)
    at_T(T) = w(T, 1.0e5, 86400.0, StateView([0.6], ix), StateView([0.0], ix), StateView([1.0], ix))
    @test at_T(303.15) / at_T(293.15) ≈ exp(-E_S / R_GAS * (1 / 303.15 - 1 / 293.15)) rtol = 1.0e-12

    slag_sp = gruyaert_slag_species()
    cs = gruyaert_kinetic_system(slag_sp)
    r = gruyaert_run(cs; slag = 0.5, τ_slag = τ, slag_species = slag_sp, days = 7)
    @test SciMLBase.successful_retcode(r.sol)
    # The slag follows the law it was given: at the reference temperature and
    # fineness of the law, the integrated degree is its closed form.
    for d in (1.0, 2.0, 7.0)
        @test gruyaert_degree(r, ["BFS"], d * day) ≈ α(d) rtol = 1.0e-3
    end
    # The clinker reacts, the faster phases first.
    @test gruyaert_degree(r, ["C3S"], 7day) > gruyaert_degree(r, ["C2S"], 7day) > 0

    # The replayed states hold the elements the paste was mixed with, kinetic
    # species and equilibrium together, and carry hydrates.
    A = Float64.(cs.SM.A)
    b0 = A * ustrip.(us"mol", r.kp.initial_state.n)
    states = speciated_states(r.sol, r.kp; times = [2day, 7day])
    for st in states
        @test A * ustrip.(us"mol", st.n) ≈ b0 rtol = 1.0e-9 atol = 1.0e-12
    end
    # Each replayed instant is certified, the first included, which has no
    # certified neighbor to walk from: at 2 days neither start certifies it,
    # and the replay ends on the full search, where it once fell back to an
    # interior-point composition at a pH of 15.3.
    @test all(st -> isapprox(pH(st), 12.82; atol = 0.01), states)
    p = r.sol.prob.p
    sub = ChemistryLab._equilibrium_subsystem(r.kp.system, r.kp.idx_equilibrium)
    des = DualEquilibriumSolver(sub, activity_model(p.eq_solver))
    for (t, st) in zip((2day, 7day), states)
        be = Float64.(r.sol(t)[1:(p.n_be)])
        at = ChemicalState(sub, st.n[r.kp.idx_equilibrium]; T = 293.15u"K")
        @test optimality_certificate(des, at; b = be).optimal
    end
    i_ch = findfirst(s -> symbol(s) == "Portlandite", cs.species)
    @test ustrip(us"mol", states[2].n[i_ch]) > ustrip(us"mol", states[1].n[i_ch]) > 0
end
