# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using Base.ScopedValues: ScopedValue, with
using DynamicQuantities
using ForwardDiff
using LinearAlgebra: pinv
using OrderedCollections
using SciMLBase

# ── EquilibriumSolver ─────────────────────────────────────────────────────────

"""
    struct EquilibriumSolver{F<:Function, S, V<:Val}

Encapsulates all fixed ingredients of a chemical equilibrium calculation:
the potential function, the SciML solver, and the variable space.

Construct once, call repeatedly with different `ChemicalState` inputs.

# Fields

  - `μ`: chemical potential closure `μ(n, p) -> Vector`, in the number type of `n` and `p`.
  - `solver`: any Optimization.jl-compatible solver (e.g. `IpoptOptimizer()`).
  - `variable_space`: variable space — `Val(:linear)` or `Val(:log)`.
  - `kwargs`: solver keyword arguments forwarded to `solve`.

# Examples
```julia
julia> cs = ChemicalSystem([
           Species("H2O"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLVENT),
           Species("H+";  aggregate_state=AS_AQUEOUS, class=SC_AQSOLUTE),
           Species("OH-"; aggregate_state=AS_AQUEOUS, class=SC_AQSOLUTE),
       ]);

julia> solver = EquilibriumSolver(cs, DiluteSolutionModel(), IpoptOptimizer());

julia> solver isa EquilibriumSolver
true
```
"""
struct EquilibriumSolver{F <: Function, S, V <: Val, M <: AbstractActivityModel}
    μ::F                    # potential closure — built once from cs and model
    solver::S               # SciML-compatible solver
    variable_space::V       # Val(:linear) or Val(:log)
    kwargs::Base.Pairs      # forwarded to solve
    model::M                # the activity model μ was built from
end

"""
    EquilibriumSolver(cs, model, solver; variable_space=Val(:linear), kwargs...)

Construct an `EquilibriumSolver` from a `ChemicalSystem`, an activity model,
and a SciML solver.

The potential function is built once at construction time from `cs` and `model`.
Repeated calls to `solve` with different `ChemicalState` inputs reuse it.

# Arguments

  - `cs`: the `ChemicalSystem` defining species and conservation matrix.
  - `model`: an `AbstractActivityModel` (e.g. `DiluteSolutionModel()`).
  - `solver`: any Optimization.jl solver.
  - `variable_space`: `Val(:linear)` (default) or `Val(:log)`.
  - `kwargs...`: forwarded to the underlying `solve` call (tolerances, verbosity...).

A solver is not meant to be used by two threads at once: its optimizer may keep
the last answer as the next start (`OptimaOptimizer`). Give each thread its own.
"""
function EquilibriumSolver(
        cs::ChemicalSystem,
        model::AbstractActivityModel,
        solver::S;
        variable_space::V = Val(:linear),
        kwargs...,
    ) where {S, V <: Val}
    _refuse_state_keywords(kwargs, "EquilibriumSolver")
    μ = build_potentials(cs, model)     # built once — captures indices and constants
    # `model` is kept alongside `μ`, which is opaque once compiled. A coupled
    # kinetics run has to rebuild the solver for the equilibrium SUB-system, and
    # without this it had no way to know which activity model to rebuild with —
    # it silently fell back to the problem's own default and ran a dilute solve
    # on a pore solution the user had asked to treat with HKF.
    return EquilibriumSolver{typeof(μ), S, V, typeof(model)}(
        μ, solver, variable_space, kwargs, model
    )
end

"""
    activity_model(solver::EquilibriumSolver) -> AbstractActivityModel

The activity model the solver's potential function was built from.
"""
activity_model(solver::EquilibriumSolver) = solver.model

# ── Temperature and pressure belong to the state ──────────────────────────────

const _STATE_KEYWORDS = (:T, :P, :temperature, :pressure)

"""
    _refuse_state_keywords(kwargs, caller)

Throw if `kwargs` names a temperature or a pressure.

The solves take the temperature and the pressure from the `ChemicalState`, and
forward the keywords they do not know to the optimizer's options. So a
`T = 293.15u"K"` handed to one of them was not an error: it reached the
optimizer, was ignored there, and the solve ran at the state's own temperature
without a word. On a C-S-H pore solution measured at 20 °C that is a 0.17
difference in pH — the change in `pKw` between 20 and 25 °C — reported as if it
were the model's.
"""
function _refuse_state_keywords(kwargs, caller::AbstractString)
    bad = [k for k in keys(kwargs) if k in _STATE_KEYWORDS]
    isempty(bad) || throw(
        ArgumentError(
            "$caller has no keyword " * join(string.(bad), ", ") * ": the " *
                "temperature and the pressure belong to the state. Set them with " *
                "`ChemicalState(cs; T = ..., P = ...)`, `set_temperature!` or " *
                "`set_pressure!`; given here they would reach the optimizer's " *
                "options and be ignored.",
        ),
    )
    return nothing
end

# ── Internal helper: build p from ChemicalState ───────────────────────────────

"""
    _build_params(state::ChemicalState; ϵ=1e-16) -> NamedTuple

Extract dimensionless parameters from a `ChemicalState`.
`ΔₐG⁰overRT` is evaluated at the current `T` and `P` of the state.
Units are stripped — compatible with ForwardDiff dual numbers.

# Returned fields

  - `ΔₐG⁰overRT`: vector of standard Gibbs energies divided by RT (dimensionless).
  - `T`: temperature in K (plain number, Dual-safe).
  - `P`: pressure in Pa (plain number, Dual-safe).
  - `ϵ`: regularization floor (default `1e-16`).
  - `ϵa`: the floor of the activities, `min(ϵ, _ACTIVITY_FLOOR)`; see
    [`_ACTIVITY_FLOOR`](@ref).

`T` and `P` are included so that temperature-dependent activity models
(e.g. [`HKFActivityModel`](@ref) with `temperature_dependent=true`) can
recompute their parameters inside the potential closure.
"""
function _build_params(state::ChemicalState; ϵ::Float64 = _AMOUNT_FLOOR)
    T = temperature(state)
    P = pressure(state)
    R = Constants.R
    RT = R * T                  # keeps units — division below strips them

    # ustrip without forced Float64 conversion — preserves Dual if T is Dual
    # In the common number type of the potentials: one species whose data are
    # being differentiated makes the whole vector dual, not a `Vector{Real}`.
    ΔₐG⁰overRT = _promoted(
        [
            ustrip(s[:ΔₐG⁰](T = T, P = P; unit = true) / RT)
                for s in state.system.species
        ]
    )

    T_K = ustrip(us"K", T)   # Quantity{Dual} → Dual, Float64 → Float64
    P_Pa = ustrip(us"Pa", P)

    # Inside a search run on the values (`_STRIP_TAGS`), the data are values too.
    return (
        ΔₐG⁰overRT = _unscoped(ΔₐG⁰overRT), T = _unscoped(T_K), P = _unscoped(P_Pa),
        ϵ = ϵ, ϵa = min(ϵ, _ACTIVITY_FLOOR),
    )
end

"""
    _build_n0(state::ChemicalState) -> Vector

Extract the dimensionless mole vector from a `ChemicalState`.
Type is inferred from the state — compatible with ForwardDiff dual numbers.
"""
function _build_n0(state::ChemicalState)
    return ustrip.(us"mol", state.n)    # no Float64 cast — type follows eltype(state.n)
end

# ── Default solver factory ────────────────────────────────────────────────────

"""
    _DEFAULT_SOLVER_FACTORY

Internal `Ref{Union{Nothing, Function}}` — populated by extension `__init__`
functions to register a default solver factory.

  - `OptimizationIpoptExt.__init__`: registers only if nothing is set (low priority).
  - `OptimaSolverExt.__init__`: always overrides (high priority).

Result: OptimaSolver wins whenever loaded, regardless of load order.
"""
const _DEFAULT_SOLVER_FACTORY = Ref{Union{Nothing, Function}}(nothing)

"""
    _SOLVER_FACTORIES

Every back end an extension has registered, in load order. The default back end is
[`_DEFAULT_SOLVER_FACTORY`](@ref), not the first entry; [`equilibrate_certified`](@ref)
uses all of them as starting points, because neither back end dominates the other —
see [`solve_certified`](@ref).

The search tries them in this order, and the order is not neutral: a start that
certifies ends the search, so which one comes first decides both the time and, on a
hard case, the verdict. Putting the default back end first was measured and not
kept. On the documentation's cements it saved about an eighth of the time overall,
but it slowed two pages down and it lost the certificate of a composite binder at
seven days (element balance 1.0 against 7.4e-13), which the search starting from
Ipopt finds.

"""
const _SOLVER_FACTORIES = Function[]

"""
    _DUAL_AVAILABLE

Whether the KKT solver and its certificate are loaded. Set by
`OptimaSolverExt.__init__`, and false otherwise: the dual Newton lives in
`OptimaSolver`, so with only the Ipopt extension loaded there is no certifying
route and `equilibrate` must take its plain path rather than raise.

Registering a back end is not the same question. `OptimizationIpopt` registers a
factory, which makes it a usable STARTING POINT for the certified route, but it
cannot certify anything by itself.
"""
const _DUAL_AVAILABLE = Ref(false)

"""
    register_solver_factory!(f)

Called from an extension's `__init__` to make its back end available to
[`equilibrate_certified`](@ref). Idempotent.
"""
function register_solver_factory!(f::Function)
    f in _SOLVER_FACTORIES || push!(_SOLVER_FACTORIES, f)
    return _SOLVER_FACTORIES
end

# ── Differentiating an equilibrium ────────────────────────────────────────────

"""
    _primal(state::ChemicalState) -> ChemicalState

The same state with every dual number replaced by its value.
"""
function _primal(state::ChemicalState)
    n = [_plain(ustrip(us"mol", nᵢ)) * u"mol" for nᵢ in state.n]
    T = _plain(ustrip(us"K", temperature(state))) * u"K"
    P = _plain(ustrip(us"Pa", pressure(state))) * u"Pa"
    return ChemicalState(state.system, n; T = T, P = P)
end

"""
    _equilibrium_sensitivity(A, H, gθ, bdot, nstar;
                             pinned = falses(n), pinnable = falses(n), maxpin = 8) -> Vector

Sensitivity of an equilibrium composition, from the optimality conditions.

Equilibrium is the Gibbs minimization of [Leal2017](@citet), the problem Reaktoro
solves: `min G(n)` subject to `A n = b`, `n >= 0`, whose first-order conditions
are `grad G(n) - A' y - z = 0`, `A n = b`, `n_i z_i = 0`, with `y` the element
potentials and `z >= 0` the stability multipliers.

Differentiating them gives, on the FREE set (species present, `z_i = 0`),

```
    | H   A' | | ndot |   | -dgradG |
    | A   0  | | ydot | = |   bdot  |
```

with `ndot = 0` on the complement and `H = grad^2 G` at the solution. One
factorization serves every partial derivative, and the answer is exact — no
finite difference, no step size.

# Which species are held at zero

`pinned` names the species the caller knows to be absent, the pure phases at
their bound. A member of a mixing phase is never absent while the phase exists,
and an aqueous trace in particular keeps its conservation law however small its
amount: pinning it would erase the perturbation of the component it carries.

No back-end returns `z`, so a pure phase near its bound may still be left free.
Among the species marked `pinnable` — the pure phases — one that is negligible
on the scale of the system yet takes a leading share of the response is pinned,
and the system solved again; on calcite and CO₂ in water with a gas phase
declared, the absent gas left free took the whole perturbation, satisfying
`A ndot = bdot` to 4e-16 and meaning nothing. Each pass pins at least one
species, so the loop terminates.

# How the system is solved

A trace species has a curvature near `1/n`, up to `1e300`, beside a
conservation block of order one, so the matrix is equilibrated by rows and
columns before a rank-revealing solve; a pseudo-inverse of the unscaled matrix
discards the conservation equations silently. `H` is singular by construction, a
pure phase having unit activity and hence a zero row, which the saddle-point form
handles and a method inverting `H` does not. Stationarity and the imposed
component perturbation are checked on the answer, and a failure of either
raises rather than returning a sensitivity that does not solve the problem.
"""
function _equilibrium_sensitivity(
        A, H, gθ, bdot, nstar;
        pinned::AbstractVector{Bool} = falses(length(nstar)),
        pinnable::AbstractVector{Bool} = falses(length(nstar)),
        maxpin::Int = 8,
    )
    ns = length(nstar)
    scale = maximum(abs, nstar)
    pinned = BitVector(pinned)
    ndot = zeros(ns)
    for _ in 0:maxpin
        free = findall(!, pinned)
        fill!(ndot, 0.0)
        ndot[free] .= _scaled_kkt_solve(A, H, gθ, bdot, free)
        big = maximum(abs, ndot)
        offenders = [
            i for i in free
                if pinnable[i] && nstar[i] < 1.0e-6 * scale && abs(ndot[i]) > 0.1 * big
        ]
        isempty(offenders) && break
        pinned[offenders] .= true
    end
    norm(A * ndot - bdot, Inf) <= 1.0e-8 * max(1.0, norm(bdot, Inf)) || error(
        "equilibrium sensitivity does not conserve the imposed component perturbation"
    )
    return ndot
end

"""
    _scaled_kkt_solve(A, H, gθ, bdot, free) -> Vector

The species part of the solution of the sensitivity system restricted to the
species `free`, equilibrated by rows and columns and solved by a truncated
singular value decomposition; raises when the answer does not satisfy the
equations to `1e-8` of their scale.
"""
function _scaled_kkt_solve(A, H, gθ, bdot, free)
    m = size(A, 1)
    K = [H[free, free] A[:, free]'; A[:, free] zeros(m, m)]
    rhs = vcat(-gθ[free], bdot)
    all(isfinite, K) && all(isfinite, rhs) ||
        error("nonfinite equilibrium sensitivity system")
    scaled = copy(K)
    load = copy(rhs)
    columns = ones(size(K, 2))
    for _ in 1:32
        rs = map(x -> iszero(x) ? 1.0 : inv(x), vec(maximum(abs, scaled; dims = 2)))
        scaled .*= rs
        load .*= rs
        cs = map(x -> iszero(x) ? 1.0 : inv(x), vec(maximum(abs, scaled; dims = 1)))
        scaled .*= cs'
        columns .*= cs
    end
    F = svd(scaled)
    cutoff = maximum(F.S) * maximum(size(scaled)) * eps(Float64)
    y = F.V * [σ > cutoff ? v / σ : 0.0 for (σ, v) in zip(F.S, F.U' * load)]
    norm(scaled * y - load, Inf) <= 1.0e-8 * max(1.0, norm(load, Inf)) ||
        error("equilibrium sensitivity does not satisfy stationarity")
    return (columns .* y)[1:length(free)]
end

"""
    _pure_phase_indices(cs) -> Vector{Int}

The species of `cs` the dual solver treats as pure phases: every species outside
the aqueous phase that belongs to no solid solution and no site family. These are
the ones a sensitivity may hold at zero when absent.
"""
function _pure_phase_indices(cs::ChemicalSystem)
    mixing = Set(vcat(cs.ss_groups..., cs.site_groups...))
    return [i for (i, s) in enumerate(cs.species) if aggregate_state(s) != AS_AQUEOUS && !(i in mixing)]
end

"""
    STRICT_CONVERGENCE

Whether a non-converged equilibrium solve raises (`true`) or warns (`false`,
the default).

The default is *not* strict, and deliberately so: the back end's convergence
flag does not by itself indicate whether a point is the minimum on these
problems. It reports `MaxIters` on points that are numerically excellent — pure
water comes back flagged while giving `[H⁺]/[OH⁻] = 1.000003` — and reports
success on points that are not the minimum. Raising on the flag alone would
reject good answers and would still miss the bad ones, so it is offered as an
opt-in for callers who want the strictest possible reading.

This is the **session default**, set by the caller and read by every solve.
Internally the package sometimes has to suspend it — while it is computing a
starting point rather than an answer — and it does so through
[`_relaxed_convergence`](@ref), which is scoped to one task and leaves this
`Ref` untouched. Read the effective value with [`_strict_convergence`](@ref),
never this `Ref` directly.
"""
const STRICT_CONVERGENCE = Ref(false)

"""
    _STRICT_OVERRIDE

A per-task override of [`STRICT_CONVERGENCE`](@ref), `nothing` when none is in
force.

It exists because the two readers of that setting want different things. A
caller sets it once for a session; the package suspends it for the duration of a
start search. Saving the `Ref`, writing it, and restoring it in a `finally` does
the second correctly only while no two of them overlap — and two solves on two
tasks do overlap. Measured, with A entering first and leaving first: B lost the
relaxation inside its own region, and the flag was left **set** after both had
finished, so every later solve in that session inherited it.

A `ScopedValue` has exactly the semantics the `finally` was imitating — dynamic
extent, inherited by child tasks, invisible to siblings — and cannot be left
behind, because there is nothing to restore.
"""
const _STRICT_OVERRIDE = ScopedValue{Union{Nothing, Bool}}(nothing)

"""
    _AUTO_SPLIT

Whether a phase declared `instances = :auto` may be given its second instance
during the solves of this task: `true`, except inside a solve that builds a
starting point. Such a solve (the ideal pre-solve, the ideal-mixing pre-solve)
hands its answer to a search in the system it was given, so its answer must stay
in that system. Given the second instance, it returned a state with two more
species than the search it was meant for, and the dual solve built from both
indexed past the end of the conservation matrix.
"""
const _AUTO_SPLIT = ScopedValue(true)

"""
    _STRIP_TAGS

The tags of the dual numbers a certified search is run without, outermost first.

A derivative is taken by solving on the values and lifting the answer by the
implicit-function theorem (`_lift_equilibrium`). The values have to be those of
every input the duals of the differentiation reach: the state and the budget,
which are stripped before the search, but also the thermodynamic data, the
parameters of the activity and mixing models, and the targets of a constraint,
which live in the species and in closures the search builds as it goes. While a
search runs in this scope, `_build_params`, the activity closures and the blocks
of a constraint return their values without these tags. A tuple and not one tag,
because a differentiation nested in another strips its own level on top of the
outer one's.
"""
const _STRIP_TAGS = ScopedValue{Tuple}(())

# A tag is whatever ForwardDiff was given: a `Tag` type, or any other value.
_strip_tag(x, tg) = x
# At every level: a dual of another tag may carry duals of `tg` inside, as the
# Jacobian a solver takes of an activity closure whose model is differentiated.
function _strip_tag(x::ForwardDiff.Dual{T}, tg) where {T}
    T === tg && return ForwardDiff.value(x)
    v = _strip_tag(ForwardDiff.value(x), tg)
    ps = map(d -> _strip_tag(d, tg), ForwardDiff.partials(x).values)
    return ForwardDiff.Dual{T}(v, ps...)
end
_strip_tag(x::AbstractArray, tg) =
    eltype(x) <: ForwardDiff.Dual && ForwardDiff.tagtype(eltype(x)) === tg ? ForwardDiff.value.(x) :
    (eltype(x) <: Number && isconcretetype(eltype(x)) && !(eltype(x) <: ForwardDiff.Dual)) ? x :
    _promoted(map(v -> _strip_tag(v, tg), x))
_strip_tags(x, tags::Tuple) = foldl((v, tg) -> _strip_tag(v, tg), tags; init = x)
_unscoped(x) = _strip_tags(x, _STRIP_TAGS[])

# The state without the duals of tag `Tg`.
function _strip_state(state::ChemicalState, Tg)
    n = [_strip_tag(ustrip(us"mol", nᵢ), Tg) * u"mol" for nᵢ in state.n]
    T = _strip_tag(ustrip(us"K", temperature(state)), Tg) * u"K"
    P = _strip_tag(ustrip(us"Pa", pressure(state)), Tg) * u"Pa"
    return ChemicalState(state.system, n; T = T, P = P)
end

# An activity closure that, built inside the scope, returns its values.
_scoped_lna(lna) = (tags = _STRIP_TAGS[]; isempty(tags) ? lna : (n, p) -> _strip_tags(lna(n, p), tags))

# The blocks of a constraint, inside the scope, returning their values.
function _scoped_blocks(bl)
    tags = _STRIP_TAGS[]
    isempty(tags) && return bl
    s(f) = f === nothing ? nothing : (args...) -> _strip_tags(f(args...), tags)
    return merge(bl, (gq = s(bl.gq), hq = s(bl.hq), cq = s(bl.cq), q0 = _strip_tags(bl.q0, tags)))
end

"""
    _strict_convergence() -> Bool

The effective strict-convergence setting: the innermost
[`_STRICT_OVERRIDE`](@ref) in force, or the session's
[`STRICT_CONVERGENCE`](@ref) when there is none.
"""
_strict_convergence() = something(_STRICT_OVERRIDE[], STRICT_CONVERGENCE[])

"""
    _relaxed_convergence(f)

Run `f` with strict convergence suspended, for this task only.

Used wherever the package is computing a **starting point** rather than a
result: a back end that reports `MaxIters` on a candidate must not raise, or the
`catch` around it silently loses that candidate and a caller asking for strict
results gets a worse search than one who did not.
"""
_relaxed_convergence(f) = with(f, _STRICT_OVERRIDE => false)

"""
    _EXPLORING_STARTS

Set while a multi-start route is computing or trying **starting points**, so the
diagnostics of a candidate are not reported as diagnostics of the answer.

A start that does not converge is ordinary and expected: `equilibrate_certified`
runs every back end from several compositions precisely because none of them
works on every problem, and it keeps whichever answer the certificate proves.
Left unguarded, a call that ends `optimal = true` still printed "returned
`MaxIters`" and "did not certify optimality" from candidates along the way, which
reads as a failed solve and is not one.

The verdict on the *answer* is untouched: `equilibrate_certified` warns or raises
on its own certificate after the search, and `verbose = true` still reports every
rejected start.

Scoped to the task that set it, for the reason given under
[`_STRICT_OVERRIDE`](@ref) — and here the leak was the quieter of the two: a
flag left set silences the non-convergence warnings of every later solve in the
session, so a genuine failure stops announcing itself.
"""
const _EXPLORING_STARTS = ScopedValue(false)

"""
    _exploring_starts(f)

Run `f` with [`_EXPLORING_STARTS`](@ref) set, for this task only. Nesting is
safe, and so is overlapping with another task's search.
"""
_exploring_starts(f) = with(f, _EXPLORING_STARTS => true)

"""
    _POLISH

Whether the answer of a back end (Ipopt, the interior point of OptimaSolver, any
optimizer behind an [`EquilibriumSolver`](@ref)) is polished by the dual Newton
of the system before it is returned: `true` unless suspended.

A back end minimizes a scalar, `n⋅μ(n)`, or steers on a gradient, and the dual
Newton solves the conditions of equilibrium themselves, `μ(n) = −Aᵀy` on the
species present, with the activities of the model whatever they derive from.
The two coincide when the activities are the gradient of one Gibbs energy,
homogeneous of degree one: then `n⋅μ(n)` is that energy and its gradient is `μ`.
For the extended Debye–Hückel models in general use they are not, and the
minimum of `n⋅μ(n)` is another composition than the equilibrium. Polished, every
route returns the composition its certificate describes, and the derivatives
lifted at it are those of the map it returns.

It is the default of the `polish` keyword of a back end's `solve`. The searches
of this package, which ask a back end only for a starting point they polish
themselves, pass `polish = false` rather than run the solve under
[`_unpolished`](@ref): a scoped value around a solve is inferred through, and on
the first cement equilibrium of a session that cost seconds of compilation.
"""
const _POLISH = ScopedValue(true)

"""
    _unpolished(f)

Run `f` with [`_POLISH`](@ref) off, for this task only: every back-end solve
inside it returns its own answer. For one solve, pass `polish = false` instead.
"""
_unpolished(f) = with(f, _POLISH => false)

"""
    NONCONVERGED :: Threads.Atomic{Int}

Running count of equilibrium solves that returned a non-success retcode.

With `STRICT_CONVERGENCE[] = false` — the default — a non-converged solve is a
`@warn` at `maxlog = 1` and its result is used anyway. Over the thousands of
steps of a coupled kinetics run that is one warning for an arbitrary number of
bad speciations, and it never reached the failure count reported by
[`integrate`](@ref), which only saw solves that actually *threw*.

Reset it with `ChemistryLab.NONCONVERGED[] = 0` before a run and read it after.
**The caller does this, not the package**: nothing in `src/` resets or reports
the counter, so a figure read without resetting first covers the whole session
and not the run. (This docstring claimed [`integrate`](@ref) did it; it never
has.)

It is an `Atomic` rather than a `Ref` only so that increments cannot be lost when
solves run on several tasks, which `+= 1` on a `Ref` does not guarantee. Reading
and writing it are unchanged.
"""
const NONCONVERGED = Threads.Atomic{Int}(0)

"""
    EXACT_HESSIAN

Whether to hand the back-end the exact Gibbs Hessian diagonal `∂μ/∂n` computed
by `ForwardDiff`, instead of letting it approximate one. Default `false`.

It is off by default only because it currently trades one defect for another,
and both are measured:

| | pure water | calcite + CO₂ |
|:--|--:|--:|
| `false` (back-end approximation) | `[H⁺]/[OH⁻] = 3.78` ✗ | worst ×19.8 |
| `true` (exact `∂μ/∂n`) | `[H⁺]/[OH⁻] = 1.000003` ✓ | worst ×3751 ✗ |

Turn it on for an **aqueous-only** system, where it makes the water
autoprotolysis come out right. Leave it off when a pure phase is present: the
exact curvature of such a phase is zero, the interior-point iteration then
stalls essentially at its starting point, and the whole speciation is wrong.

That stall is the open problem. It is not a matter of iterating longer (the
answer is identical at 300 and at 200 000 iterations) nor of the near-singular
Hessian entry (capping the inverse curvature over five orders of magnitude
changes nothing, because the solve stops before the barrier has decayed). The
C++ Optima this back-end is ported from offers a `Nullspace` linear solver in
addition to the `Rangespace` one implemented here — the latter carries the
warning that it suits diagonal Hessians only, and it is the one that inverts
`H`. Porting the nullspace path, which never inverts `H`, is the identified
next step.
"""
const EXACT_HESSIAN = Ref(false)

"""
    NULLSPACE_STEP

Whether the back-end computes its Newton step by the nullspace method rather
than the Schur complement. Default `true`.

The Schur complement forms `S = A H⁻¹ Aᵀ`, so it needs `H` invertible — and a
pure phase has unit activity, hence exactly zero curvature. The nullspace method
writes `dn = dnₚ + Z dz` with `Z` a basis of `null(A)` and solves
`(Zᵀ H Z) dz = −Zᵀ(ex + H dnₚ)`, in which `H` appears only as a product. This is
the route the C++ Optima takes by default, and its `Rangespace` counterpart —
the Schur complement — is documented there as suitable for invertible diagonal
Hessians only.

It is what makes the water autoprotolysis come out right: pure water gives
`[H⁺]/[OH⁻] = 1.0` and `pKw = 13.9994`, against 3.78 and 13.9897 through the
Schur complement, with no change to mixed solid/aqueous systems.
"""
const NULLSPACE_STEP = Ref(true)

"""
    _check_converged(sol, what) -> sol

Return `sol`, raising or warning according to [`STRICT_CONVERGENCE`](@ref) when
the optimizer did not converge. Before this existed, neither extension looked at
the return code and a non-converged iterate was written into the state as though
it were the equilibrium.
"""
function _check_converged(sol, what::AbstractString)
    SciMLBase.successful_retcode(sol) && return sol
    Threads.atomic_add!(NONCONVERGED, 1)
    if !_strict_convergence()
        _EXPLORING_STARTS[] || @warn "$what returned `$(sol.retcode)`; the \
               composition may not be an equilibrium. Set \
               `ChemistryLab.STRICT_CONVERGENCE[] = true` to raise instead." maxlog = 1
        return sol
    end
    throw(
        ErrorException(
            """
            $what did not converge: the optimizer returned `$(sol.retcode)`.

            Raise `max_iter`, loosen `tol`, or start from a better-conditioned \
            composition.
            """
        ),
    )
end

"""
    _solve_dual(esolver, state, ϵ; b = nothing) -> ChemicalState

Equilibrium of a problem carrying dual numbers: in the state (its amounts, its
temperature, its pressure), in the budget `b`, in the standard potentials of its
species or in the parameters of its activity model.

No optimization solver is asked to iterate on dual numbers — most cannot, and
Ipopt, a C library, does not operate on them. The equilibrium is solved by
`esolver` on the values of the outermost level of duals, polished by the dual
Newton (see [`_POLISH`](@ref)), and the answer is lifted by the
implicit-function theorem at it, as the certified route lifts its own
(`_lift_equilibrium`): exact at every level of a nested differentiation, with
the active set read off the answer. The conditions lifted are those of the dual
Newton, `μ(n) = −Aᵀy` on the species present, and the polish is what makes them
the conditions the returned answer satisfies; lifted at an unpolished answer of
a back end that minimized `n⋅μ(n)`, they described another map than the one it
returned. A pure phase holding no more than the certificate's floor is absent
(`10ϵ` at an answer the polish was suspended for).

Without the certified solver of the system (OptimaSolver not loaded, or no
aqueous phase), only the amounts of the state may carry duals, one level, and the
derivative is that of the optimality conditions of the unconstrained problem
(`_attach_sensitivity`).

Called from the back-end `solve` methods, which dispatch on the solver type;
making this a method of `solve` dispatching on the *state* would be ambiguous
with them.
"""
function _solve_dual(
        esolver::EquilibriumSolver, state::ChemicalState, ϵ::Float64; b = nothing,
        polish::Bool = _POLISH[],
    )
    D = _input_number_type(state, b; model = esolver.model)
    D <: ForwardDiff.Dual || throw(ArgumentError("_solve_dual: nothing to differentiate."))
    if _DUAL_AVAILABLE[] && _dual_applicable(state.system)
        Tg = ForwardDiff.tagtype(D)
        eq_v = with(_STRIP_TAGS => (_STRIP_TAGS[]..., Tg)) do
            # Rebuilt in the scope, so that its potentials return values.
            es = EquilibriumSolver(
                state.system, esolver.model, esolver.solver;
                variable_space = esolver.variable_space, esolver.kwargs...,
            )
            SciMLBase.solve(
                es, _strip_state(state, Tg); ϵ = ϵ, b = b === nothing ? nothing : _strip_tag(collect(b), Tg),
                polish = polish,
            )
        end
        des = DualEquilibriumSolver(state.system, esolver.model)
        bd = b === nothing ? des.A * _build_n0(state) : collect(b)
        floor = polish ? _CERTIFICATE_FLOOR : max(_CERTIFICATE_FLOOR, 10ϵ)
        eq_d, _ = _lift_equilibrium(des, state, eq_v, bd; ϵ = ϵ, strip_tag = Tg, floor = floor)
        return eq_d
    end
    R = _amount_number_type(state)
    (R <: ForwardDiff.Dual && !(ForwardDiff.valtype(R) <: ForwardDiff.Dual) && !_carries_duals(b)) || throw(
        ArgumentError(
            "differentiating `equilibrate(state, solver)` with respect to anything but " *
                "the amounts of the state, or more than once, needs the certified solver " *
                "of the system: load OptimaSolver, and give the system an aqueous phase " *
                "with `H2O@`. Without it, only one level of duals in the state is lifted.",
        ),
    )
    state_v = _primal(state)
    eq_v = SciMLBase.solve(esolver, state_v; ϵ = ϵ, b = b === nothing ? nothing : _plain.(b))
    nstar = Float64[ustrip(us"mol", nᵢ) for nᵢ in eq_v.n]
    # A budget given in values does not move: its derivative is zero.
    return _attach_sensitivity(state, nstar, esolver.μ, ϵ; b = b === nothing ? nothing : R.(collect(b)))
end

# Whether a back end's solve carries dual numbers it cannot iterate on: in the
# amounts or the temperature (`n0`), the budget, the data (`p`) or the activity
# model, less the levels a solve on values has stripped (`_STRIP_TAGS`).
_has_dual_inputs(n0, b, p, model) =
    eltype(n0) <: ForwardDiff.Dual || _carries_duals(b) || eltype(p.ΔₐG⁰overRT) <: ForwardDiff.Dual ||
    _type_strip(_captured_number_type(model), _STRIP_TAGS[]) <: ForwardDiff.Dual


"""
    _attach_sensitivity(state, nstar, μ, ϵ; b = nothing) -> ChemicalState

Differentiate the equilibrium map implicitly at a composition already found, and
return it carrying the dual parts.

Split out of [`_solve_dual`](@ref) so the same derivative can be attached to an
answer obtained by any route — in particular to a certified one, which is solved
in real arithmetic by construction. Pushing duals through the iteration itself
would be wrong anyway: an active set has a discrete component, so the map
`b ↦ n*(b)` is smooth only piecewise, and the derivative belongs at the solution
with the active set frozen.
"""
function _attach_sensitivity(
        state::ChemicalState{C, S, Q, R}, nstar, μ, ϵ::Float64; b = nothing,
    ) where {C, S, Q, R <: ForwardDiff.Dual}
    esolver = (; μ = μ)
    state_v = _primal(state)

    A = Float64.(_constraint_matrix(state.system))
    # The physical log activities are differentiated, not the floor the primal
    # optimizer regularizes with: in a certified state an aqueous trace may sit
    # below `ϵ`, and its curvature at the floor is not its curvature.
    positive = filter(>(0), nstar)
    floor_s = isempty(positive) ? ϵ : min(ϵ, max(floatmin(Float64), minimum(positive) / 10))
    p_v = _build_params(state_v; ϵ = floor_s)
    H = ForwardDiff.jacobian(n -> esolver.μ(n, p_v), nstar)

    # The parameter enters through the potentials and through the element
    # amounts; both partial derivatives are read off the dual parts.
    p_d = _build_params(state; ϵ = floor_s)
    μ_d = esolver.μ(nstar, p_d)
    n0_d = _build_n0(state)

    # The pure phases at their bound are absent; everything in a mixing phase,
    # an aqueous trace included, stays free.
    pure = falses(length(nstar))
    pure[_pure_phase_indices(state.system)] .= true
    pinned = pure .& (nstar .<= 10ϵ)

    npart = ForwardDiff.npartials(R)
    ns = length(nstar)
    ndot = Matrix{Float64}(undef, ns, npart)
    for k in 1:npart
        gθ = Float64[ForwardDiff.partials(μᵢ, k) for μᵢ in μ_d]
        bdot = isnothing(b) ?
            A * Float64[ForwardDiff.partials(nᵢ, k) for nᵢ in n0_d] :
            Float64[ForwardDiff.partials(bᵢ, k) for bᵢ in b]
        @views ndot[:, k] .= _equilibrium_sensitivity(
            A, H, gθ, bdot, nstar; pinned = pinned, pinnable = pure,
        )
    end

    Tag = ForwardDiff.tagtype(R)
    n_dual = [
        ForwardDiff.Dual{Tag}(nstar[i], ForwardDiff.Partials(ntuple(k -> ndot[i, k], npart)))
            for i in 1:ns
    ]
    return ChemicalState(
        state.system, n_dual .* u"mol";
        T = temperature(state), P = pressure(state),
    )
end

"""
    _finish_backend_solve(esolver, state, eq; ϵ, b = nothing, certificate = nothing,
                          polish = _POLISH[]) -> ChemicalState

The answer `eq` a back end returned from `state`, polished by the dual Newton
when `polish` holds (its default is [`_POLISH`](@ref)), OptimaSolver is loaded
and the system has an aqueous phase with `H2O@`, and returned as it is otherwise. `certificate`, a `Ref`,
receives the certificate of the answer returned, or `nothing` when none was
computed.

A polish that does not certify keeps the better of the two answers, ranked as
the certified search ranks them, and says so as a non-converged solve does: a
warning, or an error under [`STRICT_CONVERGENCE`](@ref).
"""
function _finish_backend_solve(
        esolver::EquilibriumSolver, state::ChemicalState, eq::ChemicalState;
        ϵ::Float64 = _AMOUNT_FLOOR, b = nothing, certificate = nothing,
        polish::Bool = _POLISH[],
    )
    certificate === nothing || (certificate[] = nothing)
    (polish && _DUAL_AVAILABLE[] && _dual_applicable(state.system)) || return eq
    des = DualEquilibriumSolver(state.system, esolver.model)
    bv = b === nothing ? des.A * _build_n0(state) : collect(b)
    eqp, cert = _exploring_starts(() -> solve_certified(des, (eq,); b = bv, ϵ = ϵ))
    if eqp === nothing || !cert.optimal
        raw = optimality_certificate(des, eq; b = bv, ϵ = ϵ)
        eqp, cert = eqp === nothing ? (eq, raw) : _keep_better(eqp, cert, eq, raw)
        Threads.atomic_add!(NONCONVERGED, 1)
        _strict_convergence() && throw(
            ErrorException(
                "the answer of the back end could not be polished into a " *
                    "certified equilibrium: stationarity $(cert.stationarity), " *
                    "balance $(cert.balance) mol."
            ),
        )
        _EXPLORING_STARTS[] || @warn """the answer of the back end could not be \
        polished into a certified equilibrium; the better of the two is \
        returned. Audit it with `optimality_certificate`.""" maxlog = 1
    end
    certificate === nothing || (certificate[] = cert)
    return eqp
end

"""
    _require_gibbs_duhem(esolver, system, p)

Refuse a back end that minimizes `n⋅μ(n)` on a model whose activities do not
satisfy the Gibbs–Duhem relation, when no dual Newton is there to polish its
answer.

The gradient of `n⋅μ(n)` is `μ + Jᵀn`, with `J = ∂μ/∂n`, and `Jᵀn = 0` is the
Gibbs–Duhem relation, `Σᵢ nᵢ ∂μᵢ/∂nⱼ = 0`. Where it holds, the minimum of
`n⋅μ(n)` is the equilibrium; where it fails, it is another composition, which
nothing would then correct. It is checked at a composition holding every
species, a kilogram of solvent with a tenth of a mole of each solute and one
mole of everything else, where the terms that break it are all present.
"""
function _require_gibbs_duhem(esolver::EquilibriumSolver, system::ChemicalSystem, p)
    n = ones(length(system))
    if !isempty(system.idx_solvent)
        jw = only(system.idx_solvent)
        n[jw] = 1 / ustrip(us"kg/mol", system.species[jw][:M])
        n[system.idx_solutes] .= 0.1
    end
    J = ForwardDiff.jacobian(x -> esolver.μ(x, p), n)
    defect, column = _gibbs_duhem_defect(J, n)
    defect <= _SCOPE_ASYMMETRY && return nothing
    throw(
        ArgumentError(
            "$(typeof(esolver.model).name.name) does not satisfy the Gibbs–Duhem " *
                "relation on this system (a relative defect of " *
                "$(round(defect, sigdigits = 2)) on the column of " *
                "$(symbol(system.species[column]))), so the minimum of n⋅μ(n) that this " *
                "back end computes is not the equilibrium. Load OptimaSolver, whose dual " *
                "Newton solves the conditions of equilibrium themselves and polishes the " *
                "answer, or choose a model derived from one Gibbs energy: " *
                "DiluteSolutionModel, PitzerActivityModel, HKFActivityModel with a common " *
                "ion size and Ḃ = Kₙ = 0, DaviesActivityModel without neutral solutes or " *
                "with bₙ = 0."
        ),
    )
end

# ── equilibrate ──────────────────────────────────────────────────────────────

"""
    equilibrate(state::ChemicalState, solver; model=..., variable_space=..., ϵ=...) -> ChemicalState
    equilibrate(state::ChemicalState; kwargs...) -> ChemicalState

Compute the chemical equilibrium state by minimizing the Gibbs free energy.

**Two-argument form** (solver explicit, always available once an extension is loaded):

```julia
using Optimization, OptimizationIpopt
state_eq = equilibrate(state, IpoptOptimizer())

using OptimaSolver
state_eq = equilibrate(state, OptimaOptimizer())
```

**One-argument form** — solves by every available route and returns the answer
[`optimality_certificate`](@ref) accepts; its `scope` says what that proves (a
global minimum, a KKT point, or a self-consistent speciation):

```julia
state_eq = equilibrate(state)                  # certified where it can be
state_eq = equilibrate(state; certify = false) # single back end, polished
cert = Ref{Any}()
state_eq = equilibrate(state; certificate = cert)   # cert[] === nothing if none
```

The certified route is the default because a single back end does not certify
every case here: measured on calcite dissolving in water, the interior point
returns a composition whose **charge balance is met to 3 % only**, because
the fraction-to-boundary rule caps its step and the residual stops moving. The
dual Newton gets that case to 1e-12 but fails to admit a supersaturated phase on
a low-water cement. Offering both and keeping a proved answer certifies all ten
cases of the reference battery; either alone certifies at most nine.

Use [`equilibrate_certified`](@ref) when the certificate itself is wanted, or
pass a `Ref` as `certificate`, which receives it, or `nothing` when no
certificate was computed. `certify` has three values:

  - `nothing` (the default): the certified search where it applies — OptimaSolver
    loaded, an aqueous phase with `H2O@` — and the single back end otherwise;
  - `true`: the certified search, and an `ArgumentError` where it does not apply,
    rather than an answer that was never certified;
  - `false`: the single back end, its answer polished by the dual Newton where
    the certified search would apply.

Every back end's answer is polished in that way, the two-argument form's
included, so that all routes return a composition satisfying the same conditions
of equilibrium (see `ChemistryLab._POLISH`). Without OptimaSolver, a back end
that minimizes `n⋅μ(n)`, such as Ipopt, refuses an activity model that breaks the
Gibbs–Duhem relation, since its minimum is then not the equilibrium.

When both extensions are loaded, `OptimaSolverExt` provides the default single
back end.

# Arguments

  - `state`: initial `ChemicalState` — defines the system, T, P, and composition.
  - `solver`: any SciML-compatible solver (e.g. `IpoptOptimizer()`, `OptimaOptimizer()`).
  - `model`: activity model (default: `DiluteSolutionModel()`).
  - `variable_space`: `Val(:linear)` (default) or `Val(:log)`.
  - `ϵ`: regularization floor for mole amounts (default: `1e-16`).
  - `certificate`: a `Ref` that receives the certificate of the answer returned,
    or `nothing` when none was computed.
  - `kwargs...`: forwarded to the underlying solver. The temperature and the
    pressure are not among them: they are the state's, and a `T` or a `P` given
    here is refused rather than passed on and ignored.
"""
function equilibrate(
        state::ChemicalState,
        solver;
        model::AbstractActivityModel = DiluteSolutionModel(),
        variable_space::Val = Val(:linear),
        ϵ::Float64 = _AMOUNT_FLOOR,
        certificate::Union{Nothing, Base.RefValue} = nothing,
        kwargs...,
    )
    _refuse_state_keywords(kwargs, "equilibrate")
    esolver = EquilibriumSolver(
        state.system, model, solver;
        variable_space = variable_space,
        kwargs...,
    )
    return SciMLBase.solve(esolver, state; ϵ = ϵ, certificate = certificate)
end

function equilibrate(
        state::ChemicalState;
        certify::Union{Nothing, Bool} = nothing,
        model::AbstractActivityModel = DiluteSolutionModel(),
        constraint::EquilibriumConstraint = FixedTP(),
        certificate::Union{Nothing, Base.RefValue} = nothing,
        kwargs...,
    )
    _refuse_state_keywords(kwargs, "equilibrate")
    f = _DEFAULT_SOLVER_FACTORY[]
    if isnothing(f)
        error(
            "equilibrate without explicit solver requires an extension.\n" *
                "Add `using Optimization, OptimizationIpopt` or `using OptimaSolver`, " *
                "or call `equilibrate(state, solver; ...)` explicitly.",
        )
    end
    applicable = _DUAL_AVAILABLE[] && _dual_applicable(state.system)
    certify === true && !applicable && throw(
        ArgumentError(
            "certify = true asks for a certified equilibrium, and the certified " *
                "search needs OptimaSolver loaded" * (_DUAL_AVAILABLE[] ? "" : " — it is not") *
                " and a system with an aqueous phase and `H2O@`. Leave `certify` " *
                "at its default to take the single back end where it does not apply.",
        )
    )
    if (certify !== false || !(constraint isa FixedTP)) && applicable
        eq, cert = equilibrate_certified(
            state; model = model, constraint = constraint, kwargs...,
        )
        certificate === nothing || (certificate[] = cert)
        return eq
    end
    constraint isa FixedTP || throw(
        ArgumentError(
            "a constraint other than `FixedTP` is imposed inside the dual " *
                "solver's own system, so it needs `OptimaSolver` loaded" *
                (_DUAL_AVAILABLE[] ? "" : " — it is not") *
                " and a system with an aqueous phase and `H2O@`.",
        )
    )
    return equilibrate(state, f(); model = model, certificate = certificate, kwargs...)
end
