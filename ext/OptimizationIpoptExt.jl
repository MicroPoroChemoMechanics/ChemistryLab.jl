# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

module OptimizationIpoptExt

using ChemistryLab
import ChemistryLab:
    _AMOUNT_FLOOR,
    EquilibriumProblem,
    EquilibriumSolver,
    ChemicalState,
    AbstractActivityModel,
    DiluteSolutionModel,
    _build_params,
    _build_n0,
    _solution_transform,
    _update_derived!
using Optimization
using OptimizationIpopt
using SciMLBase
using LinearAlgebra: dot, mul!
using DynamicQuantities
using ForwardDiff

# ── OptimizationProblem conversions ──────────────────────────────────────────

# Ipopt works with the Hessian of the Lagrangian. Handed a first-order backend,
# OptimizationBase builds `SecondOrder(backend, backend)` itself and warns that it
# did, once per solve; the certified cascade runs this back end many times per
# equilibrium, and the documentation pages printed the warning by the hundred.
# Asking for that same pair explicitly changes the computation in nothing.
const _SECOND_ORDER_AD = Optimization.OptimizationBase.DifferentiationInterface.SecondOrder(
    Optimization.AutoForwardDiff(), Optimization.AutoForwardDiff()
)

"""
    SciMLBase.OptimizationProblem(ep::EquilibriumProblem, ::Val{:linear}; kwargs...)

Convert an `EquilibriumProblem` to an `OptimizationProblem` in linear variable space.
The objective minimizes G = n⋅μ(n) subject to the mass conservation constraint A*n = b.
"""
function SciMLBase.OptimizationProblem(ep::EquilibriumProblem, ::Val{:linear}; kwargs...)
    Gibbs_energy(x, p) = dot(x, ep.μ(x, p))
    conservation_constraints(res, x, _) = mul!(res, ep.A, x) .-= ep.b

    optf = OptimizationFunction(
        Gibbs_energy,
        _SECOND_ORDER_AD;
        cons = conservation_constraints,
    )

    return OptimizationProblem(
        optf,
        ep.u0,
        ep.p;
        lb = ep.lb,
        ub = ep.ub,
        lcons = zeros(size(ep.A, 1)),
        ucons = zeros(size(ep.A, 1)),
        kwargs...,
    )
end

"""
    SciMLBase.OptimizationProblem(ep::EquilibriumProblem, ::Val{:log}; kwargs...)

Convert an `EquilibriumProblem` to an `OptimizationProblem` in log variable space.
Solves for x = log(n), which automatically enforces positivity and is more robust
for systems spanning many orders of magnitude.
"""
function SciMLBase.OptimizationProblem(ep::EquilibriumProblem, ::Val{:log}; kwargs...)
    Gibbs_energy_log(x, p) = (n = exp.(x); dot(n, ep.μ(n, p)))
    conservation_constraints_log(res, x, _) = (n = exp.(x); mul!(res, ep.A, n); res .-= ep.b)

    optf = OptimizationFunction(
        Gibbs_energy_log,
        _SECOND_ORDER_AD;
        cons = conservation_constraints_log,
    )

    return OptimizationProblem(
        optf,
        log.(ep.u0),
        ep.p;
        lb = log.(ep.lb),
        ub = log.(ep.ub),
        lcons = zeros(size(ep.A, 1)),
        ucons = zeros(size(ep.A, 1)),
        kwargs...,
    )
end

# ── solve(EquilibriumProblem, solver) ─────────────────────────────────────────

"""
    SciMLBase.solve(ep::EquilibriumProblem, solver; variable_space=Val(:linear), kwargs...)

Solve an `EquilibriumProblem`. The solution is transformed back to physical (mole) space.
"""
function SciMLBase.solve(
        ep::EquilibriumProblem,
        solver;
        variable_space = Val(:linear),
        kwargs...,
    )
    opt_prob = SciMLBase.OptimizationProblem(ep, variable_space)
    sol = ChemistryLab._check_converged(
        SciMLBase.solve(opt_prob, solver; kwargs...),
        "equilibrium solve",
    )
    transform = _solution_transform(variable_space)
    sol.u .= transform.(sol.u)
    return sol
end

# ── solve(EquilibriumSolver, ChemicalState) ───────────────────────────────────

"""
    SciMLBase.solve(esolver::EquilibriumSolver, state::ChemicalState; ϵ=1e-16) -> ChemicalState

Solve a chemical equilibrium problem from an initial `ChemicalState`.

The conservation matrix `A` is `conservation_matrix(state.system)`, which is `SM.A` unless a site family follows its host.
The initial mole vector `n0` and thermodynamic parameters `ΔₐG⁰/RT`
are extracted from `state` at its current `T` and `P`.

Returns a new `ChemicalState` with equilibrium mole amounts,
sharing the same `ChemicalSystem` as the input.

# Arguments

  - `esolver`: the `EquilibriumSolver` to use.
  - `state`: initial state — defines `T`, `P`, and initial composition.
  - `ϵ`: regularization floor (default: `1e-16`).

# Examples
```julia
solver = EquilibriumSolver(cs, DiluteSolutionModel(), IpoptOptimizer(); abstol=1e-10)
state0 = ChemicalState(cs, n0; T=298.15u"K", P=1u"bar")
state_eq = solve(solver, state0)
```

The logarithmic variable space (`variable_space = Val(:log)`) refines a solved
state; started from amounts held at the floor `ϵ`, it returns the start, since
the gradient in `log n` of such a species is of order `ϵ`.

This back end minimizes `n⋅μ(n)`, whose minimum is the equilibrium only where the
activities satisfy the Gibbs–Duhem relation. With OptimaSolver loaded and an
aqueous system, its answer is polished by the dual Newton, which solves the
conditions of equilibrium themselves, and `certificate` (a `Ref`) receives the
certificate of the answer returned. Without it, a model that breaks the relation
is refused: see `ChemistryLab._require_gibbs_duhem`.
"""
function SciMLBase.solve(
        esolver::EquilibriumSolver,
        state::ChemicalState;
        ϵ::Float64 = _AMOUNT_FLOOR,
        b = nothing,
        certificate = nothing,
        polish::Bool = ChemistryLab._POLISH[],
    )
    # A problem carrying dual numbers, in its state, its budget, its data or its
    # activity model, takes the implicit-function route: primal solve, then the
    # derivatives at the answer. No solver is asked to iterate on dual numbers.
    n0 = max.(_build_n0(state), ϵ)
    p = _build_params(state; ϵ = ϵ)
    ChemistryLab._has_dual_inputs(n0, b, p, esolver.model) &&
        return ChemistryLab._solve_dual(esolver, state, ϵ; b = b, polish = polish)

    # Nothing to polish the answer with: the scalar minimized must then be the
    # Gibbs energy, which needs the Gibbs–Duhem relation. Checked whatever
    # `_POLISH` says, since a start asked of this back end is no better for it.
    ChemistryLab._DUAL_AVAILABLE[] && ChemistryLab._dual_applicable(state.system) ||
        ChemistryLab._require_gibbs_duhem(esolver, state.system, p)

    # `b` given explicitly is Leal's φ(b): minimize G subject to A n = b, with
    # `state` supplying only the starting guess and the T, P conditions. The
    # element totals then come from the caller — the ODE state of a kinetics
    # run — instead of being derived from a composition that may not carry them.
    A = ChemistryLab._constraint_matrix(state.system)
    prob = isnothing(b) ?
        EquilibriumProblem(A, esolver.μ, n0; p = p) :
        EquilibriumProblem(A, esolver.μ, n0; b = collect(b), p = p)
    # The polish decides on the answer, so the back end's own return code is not
    # checked when there is one: it is neither a warning nor, under
    # `STRICT_CONVERGENCE`, an error, and a polish that fails says so itself. The
    # problem is solved here rather than by `solve(::EquilibriumProblem, …)`,
    # which checks it, and not under `_relaxed_convergence`, whose scoped value
    # around the solve costs seconds of compilation (see `OptimaSolverExt`).
    polish = polish && ChemistryLab._DUAL_AVAILABLE[] &&
        ChemistryLab._dual_applicable(state.system)
    raw = SciMLBase.solve(SciMLBase.OptimizationProblem(prob, esolver.variable_space), esolver.solver; esolver.kwargs...)
    sol = polish ? raw : ChemistryLab._check_converged(raw, "equilibrium solve")
    sol.u .= _solution_transform(esolver.variable_space).(sol.u)

    state_eq = copy(state)
    for (i, nᵢ) in enumerate(sol.u)
        state_eq.n[i] = max(nᵢ, ϵ) * u"mol"
    end
    _update_derived!(state_eq)

    return ChemistryLab._finish_backend_solve(
        esolver, state, state_eq; ϵ = ϵ, b = b, certificate = certificate, polish = polish,
    )
end

# ── Default Ipopt solver factory ──────────────────────────────────────────────

function _default_ipopt_solver()
    return IpoptOptimizer(
        acceptable_tol = 1.0e-12,
        dual_inf_tol = 1.0e-12,
        acceptable_iter = 1000,
        constr_viol_tol = 1.0e-12,
        warm_start_init_point = "no",
        # SILENT, because this solver is one back end of a multi-start cascade
        # that may run it hundreds of times in one call. At Ipopt's default
        # `print_level = 5` each of those prints a banner and a full iteration
        # table: measured on the `quickstart` page, a single equilibrium produced
        # 110 lines of output of which about 80 were Ipopt's, burying the three
        # lines the example was written to show.
        #
        # `sb = "yes"` suppresses the banner as well; `print_level` alone leaves
        # it. A caller who wants the trace builds their own optimizer —
        # `IpoptOptimizer(print_level = 5)` — and passes it explicitly, which is
        # the case where seeing it is a choice rather than an accident.
        # `additional_options`, because `IpoptOptimizer` has no `print_level`
        # field: it forwards this dictionary to Ipopt verbatim.
        additional_options = Dict{String, Any}("print_level" => 0, "sb" => "yes"),
    )
end

# ── __init__: register default solver (low priority) ─────────────────────────

function __init__()
    ChemistryLab.register_solver_factory!(_default_ipopt_solver)
    return if isnothing(ChemistryLab._DEFAULT_SOLVER_FACTORY[])
        ChemistryLab._DEFAULT_SOLVER_FACTORY[] = _default_ipopt_solver
    end
end

end # module OptimizationIpoptExt
