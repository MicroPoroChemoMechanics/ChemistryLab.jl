# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

module OptimaSolverExt

using ChemistryLab
import ChemistryLab:
    _AMOUNT_FLOOR,
    _LOG_UNDERFLOW,
    EquilibriumProblem,
    EquilibriumSolver,
    ChemicalState,
    _build_params,
    _build_n0,
    _solution_transform,
    _update_derived!
using OptimaSolver: OptimaOptimizer, DualNewtonProblem, DualNewtonOptions,
    SolutionPhase, dual_newton_solve, dual_newton_tangent, kkt_certificate, lp_start
# Read at a dual answer exactly as the certificate reads it.
using OptimaSolver: current_g, current_h, _degenerate_conservation_rows
using SciMLBase
using LinearAlgebra: dot, mul!
using DynamicQuantities
using ForwardDiff

# ── OptimizationProblem helpers (NoAD — OptimaOptimizer handles gradients) ────

# The element-conservation matrix and vector are passed through the parameters.
# Without them `OptimaSolver` falls back to rebuilding `A` by finite differences
# on the constraint function, which caps the achievable feasibility at ~1e-6
# whatever tolerance is requested — `A` is known exactly, so it is handed over.
# Only in the linear parameterization: in log space the constraint is
# A·exp(x) = b, which is not linear in the optimization variables, so handing
# over `A` would be wrong there.
"""
    _hessian_diagonal(μ, q) -> (hf, n) -> hf

Exact diagonal of ∇²G, handed to the back-end alongside `A` and `b`.

`G(n) = nᵀ μ(n)` gives `∂G/∂nᵢ = μᵢ + Σⱼ nⱼ ∂μⱼ/∂nᵢ`, and the sum vanishes by
Gibbs–Duhem at fixed `T`, `P`. So `∇²G = ∂μ/∂n` — the Hessian *is* the Jacobian
of the chemical potentials, which `ForwardDiff` returns exactly.

It is worth the derivative evaluation. The entries are `1/nᵢ` for a solute and
exactly zero for a pure phase, and both fallbacks in the back-end get one of the
two wrong. Approximating this does not merely slow Newton down: understating the
curvature on a trace ion makes the line search reject its step, the barrier is
exhausted, and the solve stops on a point that is not the minimum — water then
comes out with `[H⁺]/[OH⁻] = 3.8` instead of 1.
"""
function _hessian_diagonal(μ, q)
    return function (hf, n)
        J = ForwardDiff.jacobian(nn -> μ(nn, q), n)
        @inbounds for i in eachindex(hf)
            hf[i] = max(J[i, i], zero(eltype(hf)))
        end
        return hf
    end
end

function _with_constraints(q::NamedTuple, A, b, μ)
    base = merge(q, (A = A, b = b))
    return ChemistryLab.EXACT_HESSIAN[] ?
        merge(base, (hdiag = _hessian_diagonal(μ, q),)) : base
end

function _build_optima_opt_prob(ep::EquilibriumProblem, μ, ::Val{:linear})
    f_gibbs(x, q) = dot(x, μ(x, q))
    # ∂G/∂nᵢ = μᵢ exactly: the remaining term Σⱼ nⱼ ∂μⱼ/∂nᵢ vanishes by
    # Gibbs–Duhem at fixed T, P. Handing the gradient over rather than leaving
    # the back-end to differentiate `dot(n, μ(n))` matters, because that
    # cancellation is only exact in theory. Evaluated, the solvent term alone is
    # 55 mol × 0.018 ≈ 1 against a μ of order 10, so the back-end was steering
    # on a gradient wrong by some ten percent — and a Gibbs energy is stationary
    # precisely where its gradient is, so the answer it settled on was not the
    # equilibrium. Water came out with [H⁺]/[OH⁻] = 4.
    g_gibbs!(g, x, q) = (g .= μ(x, q))
    cons!(res, x, _) = mul!(res, ep.A, x) .-= ep.b
    optf = SciMLBase.OptimizationFunction{true}(f_gibbs; grad = g_gibbs!, cons = cons!)
    return SciMLBase.OptimizationProblem(
        optf, ep.u0, _with_constraints(ep.p, ep.A, ep.b, μ);
        lb = ep.lb, ub = ep.ub,
        lcons = zeros(size(ep.A, 1)),
        ucons = zeros(size(ep.A, 1)),
    )
end

function _build_optima_opt_prob(ep::EquilibriumProblem, μ, ::Val{:log})
    f_gibbs(x, q) = (n = exp.(x); dot(n, μ(n, q)))
    # In `x = ln n` the gradient of G is `n ∘ μ`, by the chain rule on the one of
    # the linear route above, and it is handed over for the same reason: the
    # derivative of `dot(n, μ(n))` carries the term `Jᵀn`, which is zero only
    # where the model satisfies the Gibbs–Duhem relation, and steering on it
    # settled on another composition than the equilibrium.
    g_gibbs!(g, x, q) = (n = exp.(x); g .= n .* μ(n, q))
    cons!(res, x, _) = (n = exp.(x); mul!(res, ep.A, n); res .-= ep.b)
    optf = SciMLBase.OptimizationFunction{true}(f_gibbs; grad = g_gibbs!, cons = cons!)
    return SciMLBase.OptimizationProblem(
        optf, log.(ep.u0), ep.p;
        lb = log.(ep.lb), ub = log.(ep.ub),
        lcons = zeros(size(ep.A, 1)),
        ucons = zeros(size(ep.A, 1)),
    )
end

# ── solve(EquilibriumSolver{OptimaOptimizer}, ChemicalState) ──────────────────

"""
    SciMLBase.solve(esolver::EquilibriumSolver{<:Function, <:OptimaOptimizer},
                   state::ChemicalState; ϵ=1e-16) -> ChemicalState

Solve a chemical equilibrium problem using an `OptimaOptimizer` solver.
Loaded automatically when `using OptimaSolver` is active.
"""
# The first parameter is bounded as the struct bounds it, and this is what makes
# the method reachable at all once `OptimizationIpoptExt` is loaded. That
# extension defines `solve(::EquilibriumSolver, ::ChemicalState)` for every
# back end, and this signature used to read `EquilibriumSolver{F,
# <:OptimaOptimizer, V} where {F, V}`: with `F` and `V` free of the bounds the
# struct declares, Julia does not rank it as more specific than the bare
# `EquilibriumSolver`, and the generic method won. Every `OptimaOptimizer` solve
# of a session that had also loaded Ipopt then went through the generic
# `OptimizationProblem` -- no exact conservation matrix, the gradient of
# `dot(n, μ(n))` by automatic differentiation instead of `μ` itself -- which is
# the path the comments above measure as wrong. On a 109-species cement it
# returned a start 2e-3 mol away from this method's, which the dual solve could
# not certify, and a certified solve that costs 0.7 s here cost 34 to 78 s
# there. The documentation loads Ipopt, so every page paid it.
function SciMLBase.solve(
        esolver::EquilibriumSolver{<:Function, <:OptimaOptimizer},
        state::ChemicalState;
        ϵ::Float64 = _AMOUNT_FLOOR,
        b = nothing,
        certificate = nothing,
    )
    # A problem carrying dual numbers, in its state, its budget, its data or its
    # activity model, takes the implicit-function route: primal solve, then the
    # derivatives at the answer. No solver is asked to iterate on dual numbers.
    n0 = max.(_build_n0(state), ϵ)
    p = _build_params(state; ϵ = ϵ)
    ChemistryLab._has_dual_inputs(n0, b, p, esolver.model) &&
        return ChemistryLab._solve_dual(esolver, state, ϵ; b = b)

    # `b` given explicitly is Leal's φ(b): minimize G subject to A n = b, with
    # `state` supplying only the starting guess and the T, P conditions. The
    # element totals then come from the caller — the ODE state of a kinetics
    # run — instead of being derived from a composition that may not carry them.
    A = ChemistryLab._constraint_matrix(state.system)
    prob = isnothing(b) ?
        EquilibriumProblem(A, esolver.μ, n0; p = p) :
        EquilibriumProblem(A, esolver.μ, n0; b = collect(b), p = p)
    opt_prob = _build_optima_opt_prob(prob, esolver.μ, esolver.variable_space)

    # The polish decides on the answer, so the interior point's own return code
    # is not a reason to raise under `STRICT_CONVERGENCE` when there is one.
    polish = ChemistryLab._POLISH[] && ChemistryLab._DUAL_AVAILABLE[] &&
        ChemistryLab._dual_applicable(state.system)
    run() = ChemistryLab._check_converged(
        SciMLBase.solve(opt_prob, esolver.solver; esolver.kwargs...),
        "equilibrium solve",
    )
    sol = polish ? ChemistryLab._relaxed_convergence(run) : run()
    transform = _solution_transform(esolver.variable_space)

    state_eq = copy(state)
    for (i, nᵢ) in enumerate(sol.u)
        state_eq.n[i] = max(transform(nᵢ), ϵ) * u"mol"
    end
    _update_derived!(state_eq)

    return ChemistryLab._finish_backend_solve(esolver, state, state_eq; ϵ = ϵ, b = b, certificate = certificate)
end

# ── __init__: register default solver (high priority — always overrides) ──────
#
# Measured on the cement equilibria this package targets, once the exact
# conservation matrix is handed over (see `_with_constraints` above):
# OptimaSolver reaches a 4e-14 element balance and is 3 to 26 times faster than
# Ipopt. Pass a solver explicitly to override the choice.

_default_optima_solver() = OptimaOptimizer(;
    use_fd_hessian = !ChemistryLab.EXACT_HESSIAN[],
    nullspace_step = ChemistryLab.NULLSPACE_STEP[],
)

# ── the certified KKT solver ──────────────────────────────────────────────────
#
# `ChemistryLab.DualEquilibriumSolver` supplies the chemistry — the conservation
# matrix, the reference potentials, the activity model, and which species are
# strictly positive, which may vanish, and which is the solvent. The algorithm
# is `OptimaSolver`'s, and these four methods are the join.

# A phase whose composition the substitution cannot recover (a sublattice
# model) is inverted by Newton's method, and the members it may lack entirely
# are declared: see `SolutionPhase` in OptimaSolver.
function _solution_phase(ph)
    return SolutionPhase(
        ph.members, ph.j_ref;
        always_present = ph.always_present, mole_fraction = ph.mole_fraction,
        split_starts = get(ph, :split_starts, ()), newton = get(ph, :newton, false),
        bounded_members = get(ph, :bounded_members, Int[]),
        local_h = get(ph, :local_h, nothing), invert = get(ph, :invert, nothing),
    )
end

function ChemistryLab._optima_dual_problem(
        A, g, lna, phases, idx_bounded, params,
        gq = nothing, hq = nothing, cq = nothing,
        q0 = Float64[], qscale = Float64[], Aq = nothing,
        always_active = Int[], conservation_rows = nothing,
    )
    return DualNewtonProblem(
        A, g, lna;
        phases = [_solution_phase(ph) for ph in phases],
        idx_bounded = idx_bounded, params = params,
        gq = gq, hq = hq, cq = cq, q0 = q0, qscale = qscale,
        Aq = Aq === nothing ? zeros(Float64, size(A, 1), length(q0)) : Aq,
        always_active = always_active,
        conservation_rows = conservation_rows === nothing ?
            (1:size(A, 1)) : conservation_rows,
    )
end

# The options of a `DualEquilibriumSolver`, as OptimaSolver takes them.
_dual_newton_options(o) = DualNewtonOptions(;
    tol = o.tol, maxit = o.maxit,
    max_active_updates = o.max_active_updates,
    si_tol = o.si_tol, inner_tol = o.inner_tol,
    inner_maxit = o.inner_maxit, inner_fall_bound = o.inner_fall_bound,
    lenient_line_search = o.lenient_line_search, verbose = o.verbose,
)

function ChemistryLab._optima_dual_solve(prob, b, x0, o)
    return dual_newton_solve(prob, b, x0; opts = _dual_newton_options(o))
end

ChemistryLab._optima_lp(prob, b) = lp_start(prob, b)

# The solutes of `solutes`, each below the activity floor `ϵ`, at the amount the
# potentials of `res` give them: `ϵ·exp(r)` where `r = uᵢ − ∇fᵢ > 0`, the rest
# unchanged. `u = −Aᵀy` with the multipliers of the solve, and `∇f` read at the
# answer, as the certificate reads both. A species a vanished component forces to
# zero is left alone. See `ChemistryLab._complete_floored_solutes`.
function ChemistryLab._optima_complete_floored(prob, res, solutes, ϵ)
    x = copy(res.x)
    u = -(transpose(prob.A) * res.y)
    ∇f = current_g(prob, res.q) .+ current_h(prob, x, res.q)
    rows = _degenerate_conservation_rows(prob, prob.A * x)
    for i in solutes
        any(abs(prob.A[k, i]) > 0 for k in rows) && continue
        r = u[i] - ∇f[i]
        r > 0 && (x[i] = ϵ * exp(min(r, -_LOG_UNDERFLOW)))
    end
    return x
end

# The answer `x` of `prob` on the values of its data, lifted to the dual numbers
# `b` or the data carry: the implicit-function theorem at the answer, with the
# active set frozen. See `ChemistryLab._lift_equilibrium`.
ChemistryLab._optima_tangent(prob, b, x; q = nothing, floor = 1.0e-25, primal = nothing) =
    dual_newton_tangent(prob, b, x; q = q, floor = floor, primal = primal)

function ChemistryLab._optima_kkt_certificate(
        prob, x, b, floor, tol, si_tol, q = nothing,
    )
    return q === nothing ?
        kkt_certificate(prob, x, b; floor = floor, tol = tol, si_tol = si_tol) :
        kkt_certificate(prob, x, b; floor = floor, tol = tol, si_tol = si_tol, q = q)
end

function __init__()
    ChemistryLab.register_solver_factory!(_default_optima_solver)
    ChemistryLab._DUAL_AVAILABLE[] = true
    return ChemistryLab._DEFAULT_SOLVER_FACTORY[] = _default_optima_solver
end

end # module OptimaSolverExt
