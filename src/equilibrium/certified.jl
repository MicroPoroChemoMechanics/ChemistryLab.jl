# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# ── certified.jl ──────────────────────────────────────────────────────────────
#
# An equilibrium that comes with a proof, or says it has none.
#
# The back ends do not agree, and the disagreement is not small. Measured on the
# same problems, with the element balance judged row by row against each row's
# own budget:
#
#   | case                        | interior point | dual Newton |
#   |-----------------------------|----------------|-------------|
#   | calcite in water, 10-40 °C  | 3.0e-2         | 3.0e-12     |
#   | calcite + 1 mmol CO2        | 1.0e-3         | 6.3e-13     |
#   | calcite + 50 mmol CO2       | 1.3e-16        | 8.5e-14     |
#   | pure water                  | 1.3e-16        | 1.3e-16     |
#   | CEM I paste, w/c 0.45, 0.60 | 1.5e-14        | 2.0e-14     |
#   | CEM I paste, w/c 0.30       | 7.2e-16        | 3.3e-10, not certified |
#
# The interior point is wrong in the second digit of the charge balance on a
# calcite solution, because the fraction-to-boundary rule caps its step and the
# residual stops moving; the dual Newton fails to admit a supersaturated phase on
# a low-water cement. Neither is the right default on its own.
#
# What settles it is that the problem is convex whenever the mixing terms are, so
# the KKT conditions are sufficient and `optimality_certificate` DECIDES rather
# than ranks. Given a decision procedure, solving by every available route and
# keeping a proved answer is exact, not heuristic.

"""
    _MAX_RESTARTS

How many times [`equilibrate_certified`](@ref) may restart from its own answer
before giving up. One round is what the measured cases need; the bound exists so
a case that improves by a hair every round cannot loop.
"""
const _MAX_RESTARTS = 3

"""
    _keep_better(eq, cert, eq2, cert2) -> (eq, cert)

Keep the better of two answers: a certificate of optimality beats none, and
otherwise the smaller KKT error wins. A round that buys nothing changes nothing,
which is what lets the restart loop run without ever making the answer worse.

The optimality flag is compared **first**, in both directions. Ranking on the
residuals alone would let an uncertified point displace a certified one, and no
residual is worth trading a proof for.

Among uncertified answers the comparison is [`_kkt_error`](@ref) — the worst of
the three residuals — and not the stationarity alone, which is the same ranking
[`solve_certified`](@ref) uses internally. The distinction is not academic: an
answer stationary to 2.4e-3 whose element balance was off by 6.7 mol beat every
candidate the continuation produced, because those were stationary to only 1e-2
while conserving mass. The search computed a usable answer and discarded it.
"""
function _keep_better(eq, cert, eq2, cert2)
    cert2.optimal == cert.optimal || return cert2.optimal ? (eq2, cert2) : (eq, cert)
    a, a2 = _within_domain(eq), _within_domain(eq2)
    a == a2 || return a2 ? (eq2, cert2) : (eq, cert)
    return _kkt_error(cert2) < _kkt_error(cert) ? (eq2, cert2) : (eq, cert)
end

"""
    _within_domain(eq) -> Bool

Whether the answer is inside the domain the model is written for: an aqueous
system whose solvent has been eaten by the solids is not.

This is a ranking question, not a diagnosis. [`_check_solvent`](@ref) already
reports such a state on the answer that is returned; what it could not do is stop
one from being *chosen*. Measured on a CEM I with eight solid solutions: a wandering
iterate came back with the solvent at `x_w = 0.033` and an element balance of
71 mol, the continuation's answer stood at 160, and the smaller number won — so
the search settled on a composition in which the water had gone into the solids,
and every route that conserved mass was discarded behind it.

Ranking on residuals alone cannot separate the two: both are large, and one is
merely larger. The distinction that matters is not how big the residual is but
whether the point is a composition this model can describe at all.
"""
function _within_domain(eq::ChemicalState)
    isempty(eq.system.idx_solvent) && return true
    return solvent_fraction(eq) >= SOLVENT_FRACTION_FLOOR
end

# The ranking is exercised on placeholders in `test/certified_equilibrium.jl`,
# where the answers are stand-ins rather than compositions. Nothing that is not a
# state carries a solvent to have lost, so the question does not arise.
_within_domain(::Any) = true

"""
    _repair_start(eq, model, bfix, ϵ) -> Union{ChemicalState, Nothing}

Build a starting point that makes the phases the certificate says are missing
actually present. `nothing` when there are none.

This is what turns a diagnosis into a repair. `optimality_certificate` reports a
positive worst supersaturation when a phase sits at the lower bound while the
solution is supersaturated with respect to it — a genuine KKT failure on a convex
problem, so the active set is wrong and the answer is not the answer.
[`saturation_indices`](@ref) says *which* phase, and the obvious move is then to
put it in and solve again.

The failure it exists for is a **phase swap**, which an active-set loop that
admits one phase at a time cannot perform. Measured on a CEM I paste, reported
from one machine while another certified the same source: all 0.02515 mol of
magnesium sat in brucite with `hydrotalcite` absent and supersaturated by 5.58
log units, and admitting the hydrotalcite requires dissolving the brucite
entirely and taking aluminum back from the hydrogarnet in the same step. The
solve was otherwise impeccable — stationarity 1.5e-16, element balance 1.8e-14 —
which is exactly what a converged-onto-the-wrong-active-set answer looks like.

The amount each missing phase is given is what the recipe could make of it,
`min_c b_c / A_cs` over the components it consumes, scaled by `_REPAIR_FRACTION`.
That is a chemical bound, not a guess at the answer: a carbonate in a system with
1e-9 mol of carbon is offered 1e-9 mol and no more. What matters is only that the
phase starts well away from the boundary, since being *at* the boundary is what
the active-set loop cannot recover from. The element balance of the result is not
this function's business — `b` is fixed once by the caller and every start is
projected onto it, so a start that over-spends the budget costs nothing.
"""
function _repair_start(eq::ChemicalState, model, bfix, ϵ::Float64)
    cs = eq.system
    si = saturation_indices(eq, model; ϵ = ϵ)
    n = Float64[ustrip(us"mol", x) for x in eq.n]
    A = _constraint_matrix(cs)
    # A solid-solution end-member is not a pure phase and its own index does not
    # decide anything: the activity of a member is `ln x`, so an index taken at
    # the bound reports that its mole fraction is small, not that the phase should
    # form. The phase has its own criterion, applied below.
    ss = Set(cs.idx_ssendmembers)
    missing_phases = [
        i for i in cs.idx_crystal
            if !(i in ss) && n[i] <= 10ϵ && get(si, symbol(cs.species[i]), -Inf) > 1.0e-4
    ]

    n2 = copy(n)
    cap_of(i) = minimum(
        (
            bfix[c] / A[c, i]
                for c in eachindex(bfix) if A[c, i] > 0 && bfix[c] > 0
        );
        init = Inf,
    )
    for i in missing_phases
        cap = cap_of(i)
        isfinite(cap) && cap > 0 || continue
        n2[i] = max(n2[i], _REPAIR_FRACTION * cap)
    end

    # ── the phases that are not pure ──
    #
    # An ideal solid solution can only exist at `xᵢ = 10^SIᵢ`, with `SIᵢ` the index
    # of the PURE end-member: that is what `μᵢ⁰ + RT ln xᵢ = νᵢ` says, one equation
    # per member. A composition needs `Σᵢ xᵢ = 1`, so the phase is exactly
    # saturated when `Ω = Σᵢ 10^SIᵢ = 1` and should form above it. Those numbers
    # are available here: `saturation_indices` reports `SIᵢ − log₁₀ xᵢ`, and adding
    # `ln aᵢ / ln 10` back removes the mole fraction — and with it the `0/0` a
    # phase sitting entirely at the floor would otherwise produce.
    #
    # Without this a search that has lost a solid solution cannot get it back:
    # every member reports a small index because its mole fraction is small, no
    # member is offered, and the round does nothing. On a CEM I with eight solid
    # solutions that is the difference between an answer whose element balance is
    # 13.6 mol and one certified to 1.1e-14.
    lna = log_activities(eq, model; ϵ = ϵ)
    inv_ln10 = inv(log(10))
    for (grp, phase) in zip(cs.ss_groups, cs.solid_solutions)
        all(n[i] <= 10ϵ for i in grp) || continue        # already present
        sip = [
            si[symbol(cs.species[i])] + lna[symbol(cs.species[i])] * inv_ln10
                for i in grp
        ]
        m = maximum(sip)
        isfinite(m) || continue
        Ω = sum(10^(v - m) for v in sip) * 10^m
        Ω > 1 || continue                                 # cannot form, whatever mixing
        x = [10^(v - m) for v in sip]
        x ./= sum(x)                                      # the composition it would take
        cap = minimum((cap_of(i) for i in grp); init = Inf)
        isfinite(cap) && cap > 0 || continue
        for (j, i) in enumerate(grp)
            n2[i] = max(n2[i], _REPAIR_FRACTION * cap * x[j])
        end
    end

    n2 == n && return nothing
    return ChemicalState(
        cs, n2 .* u"mol"; T = temperature(eq), P = pressure(eq),
    )
end

"""
    _LazyStarts

The starting points offered to [`solve_certified`](@ref), each back end solved
only when the search actually asks for it, and cached once it has been.

`solve_certified` returns at the **first** start that certifies, so solving
every registered back end before the search begins pays for answers the search
may never look at.

Measured honestly, the gain depends entirely on whether the first start
certifies. On the CEM I of `docs/src/examples/cem1_solid_solutions.md` it does
**not**, the later starts are genuinely consumed, and this buys nothing: 107.9 s
against 105.1 s eager, within noise. The saving appears only where the first
back end settles the problem, which is the common case for the smaller systems.
So this is a correctness-of-effort change — never compute a start no one reads —
and not a cement optimization; the expensive cement case needs a different
answer.

Cached, and that is not optional: the route search offers the same starts again
after the ideal pre-solve and again after the continuation, so recomputing them
each time would more than undo the gain.

Iterating yields each back end's answer in registration order, skipping any that
threw, and finally `tail` — the caller's own state, which is the only start left
if every back end failed — unless `offer_tail` is false.
"""
struct _LazyStarts{F, S}
    solve_one::F                     # factory -> a state, or `nothing` if it threw
    factories::Vector{Function}
    tail::S                          # the caller's own state, offered last
    cache::Vector{S}
    tried::Base.RefValue{Int}
    offer_tail::Bool
end

Base.IteratorSize(::Type{<:_LazyStarts}) = Base.SizeUnknown()
Base.eltype(::Type{_LazyStarts{F, S}}) where {F, S} = S

function Base.iterate(s::_LazyStarts, i::Int = 1)
    # Solve just far enough to answer for element `i`, and no further.
    while i > length(s.cache) && s.tried[] < length(s.factories)
        s.tried[] += 1
        r = s.solve_one(s.factories[s.tried[]])
        r === nothing || push!(s.cache, r)
    end
    i <= length(s.cache) && return (s.cache[i], i + 1)
    i == length(s.cache) + 1 && s.offer_tail && return (s.tail, i + 1)
    return nothing
end

"""
    _DeferredStarts(make)

Starts computed only when the search reaches them, and once: `make()` returns an
iterable of starting states, or `nothing`, and is called the first time the
iterator is iterated; the search offers its starts again after a continuation or
a repair, and the second pass reads what the first made. For a start that is a
solve in its own right, such as [`_ideal_mixing_start`](@ref
ChemistryLab._ideal_mixing_start), which a search that certifies from an earlier
start never needs.
"""
mutable struct _DeferredStarts{F}
    make::F
    made::Any
end
_DeferredStarts(make) = _DeferredStarts(make, nothing)

Base.IteratorSize(::Type{<:_DeferredStarts}) = Base.SizeUnknown()

function Base.iterate(d::_DeferredStarts, st = nothing)
    if st === nothing
        d.made === nothing && (d.made = something(d.make(), ()))
        r = iterate(d.made)
    else
        r = iterate(d.made, only(st))
    end
    r === nothing && return nothing
    return (r[1], (r[2],))
end

"""
    _ideal_start(state, model, bfix, ϵ, constraint, verbose; kwargs...)
        -> Union{ChemicalState, Nothing}

A certified answer to the same problem under **ideal** activities, to be used as a
starting point. `nothing` when that solve does not certify either, or when
`model` is already the ideal one.

The easier question is the useful one here. Without activity coefficients the
residual does not depend on the composition through a second, non-linear path, so
the solve is far better conditioned and certifies where the non-ideal model does
not; and the phases it finds are the same ones — they differ in amount, not in
identity — so the non-ideal solve that starts from it begins with the correct
active set instead of discovering it.

That discovery is what was fragile. On a CEM I at `w/c = 0.5` with the eight
distinct CEMDATA18 solid solutions, eighty phases sit at the bound in the cold
state and the active-set search decides its route on comparisons of nearly equal
quantities: `100/sum(oxides)` summed over a `Dict` and over an `OrderedDict`
differ by one ulp, and that was enough to choose between an equilibrium certified
to 1.1e-14 and a failure with an element balance of 71 mol. Started from the ideal
answer, both reach the same certified composition.

Any failure of the inner solve is swallowed: this builds a starting point, and a
caller who asked for a result is entitled to the outer verdict rather than to an
error raised inside a heuristic.
"""
function _ideal_start(
        state::ChemicalState, model, bfix, ϵ::Float64, constraint,
        verbose::Bool; kwargs...,
    )
    model isa DiluteSolutionModel && return nothing
    return try
        # A start for a search in this system, so never in an enlarged one.
        eq0, cert0 = with(_AUTO_SPLIT => false) do
            equilibrate_certified(
                state; model = DiluteSolutionModel(), b = bfix, ϵ = ϵ,
                constraint = constraint, verbose = false, autostart = true, kwargs...,
            )
        end
        cert0.optimal ? eq0 : nothing
    catch err
        verbose && @info "the ideal pre-solve did not run" err
        nothing
    end
end

"""
    _linear_program(des, state, bfix) -> Union{NamedTuple, Nothing}

The linear program over the pure phases of the problem `equilibrate_certified`
is about to solve: minimize `gᵀn` subject to `A n = b`, `n ≥ 0`, with the
standard potentials, the conservation matrix and the budget of the dual solve
itself (OptimaSolver's `lp_start`). `nothing` when it could not be set up.

The program is what the equilibrium becomes when every activity is one, and it
answers two questions at the cost of a few hundred pivots.

  - **Can the budget be met at all?** The minimization has no solution when no
    non-negative amounts of the declared species reproduce `b`, and the program
    says so with a Farkas vector `z`, `Aᵀz ≥ 0` and `bᵀz < 0`, verified: a
    combination of the balances that every species raises and the budget lowers.
    The cascade would otherwise spend its restarts, repairs and continuation on a
    problem without an answer.
  - **Where to start.** Its vertex is the assemblage of pure phases of lowest
    standard energy, and its multipliers give every other species an amount; see
    `_lp_lifted_state`.
"""
function _linear_program(des::DualEquilibriumSolver, state::ChemicalState, bfix)
    return try
        p = _build_params(state)
        n0 = Float64[ustrip(us"mol", x) for x in state.n]
        prob = _dual_problem(des, p, n0)
        (; start = _optima_lp(prob, bfix), prob, params = p)
    catch err
        nothing
    end
end

"""
    _lp_lifted_state(des, state, lp, ϵ) -> ChemicalState

The start the linear program gives: its vertex for the species it holds, and for
the others the amount the multipliers `y` of the program give them in their own
phase. With `u = −Aᵀy`, a species absent from the vertex has the activity
`aⱼ = exp(uⱼ − gⱼ)`, which is at most one, and its amount follows from what
the activity of its phase is:

  - a solute: a molality, so `aⱼ` mol/kg of the water at the vertex;
  - a member of a solid solution present at the vertex: a mole fraction, so `aⱼ`
    times the amount of the phase there;
  - a pure phase, or a member of a solid solution absent from the vertex: zero.
    The phase is not part of the assemblage, and the solve decides whether it
    enters.

A dead species (one whose component the budget lacks) stays at zero.

No member of a phase starts below `ϵ`, the floor of the search, as none does in
a state as given. The potentials of the program know nothing of a trace
species: at the pH of a cement they give H⁺ an activity near 1e-100, far below
anything the search resolves, and a logarithmic iteration climbs back from such
an amount slowly if at all.

The amounts are those of an ideal phase at the potentials of the program, which
is what the dual solve needs to start from; they do not meet the budget exactly,
which it does not need either, `b` being fixed for the whole call.

The first form gave every absent species `aⱼ` mol whatever its phase, capped at
one mole. On a cement that is dozens of species near a mole, and pure phases the
program found undersaturated among them. Measured from the cast state of four
cements, the start was then 11 to 66 mol off a 4-mol budget, and on two of them
nothing certified from it, so the search lost 9 and 15 s before the state as
given certified. Placed in its phase, the start is 1.3 to 1.8 mol off, and all
four certify from it at the first attempt, in 0.1 to 0.2 s.
"""
function _lp_lifted_state(des::DualEquilibriumSolver, state::ChemicalState, lp, ϵ::Real)
    s, prob = lp.start, lp.prob
    u = -(transpose(prob.A) * s.y)
    dead(j) = isnan(s.reduced_costs[j])
    activity(j) = dead(j) ? 0.0 : exp(clamp(u[j] - prob.g[j], _LOG_UNDERFLOW, 0.0))
    n = [s.x[j] > 0 ? s.x[j] : 0.0 for j in eachindex(prob.g)]
    # The solutes, per kilogram of the water the vertex holds.
    if des.j_solvent > 0
        jw = des.idx_aq[des.j_solvent]
        kg = n[jw] * ustrip(us"kg/mol", state.system.species[jw][:M])
        for j in des.idx_aq
            j == jw || n[j] > 0 || dead(j) || (n[j] = max(activity(j) * kg, ϵ))
        end
    end
    # The members of each solid solution, as fractions of the phase; none if the
    # vertex holds none of it.
    for grp in des.ss_groups
        N = sum(n[j] for j in grp)
        N > 0 || continue
        for j in grp
            n[j] > 0 || dead(j) || (n[j] = max(activity(j) * N, ϵ))
        end
    end
    return ChemicalState(state.system; T = state.T[1], P = state.P[1], n = n .* u"mol")
end

"""
    _refuse_infeasible_budget(des, state, bfix, lp, ϵ) -> (state, certificate)

The answer to a budget no amounts of the declared species can meet: a copy of
the state as given, and its certificate with `optimal = false`, `budget_feasible = false`,
`route = :infeasible` and `unplaceable`, the reason in words, read from the
Farkas vector of the linear program.
"""
function _refuse_infeasible_budget(des::DualEquilibriumSolver, state, bfix, lp, ϵ)
    z = lp.start.farkas
    A = des.A
    prim = [symbol(sp) for sp in des.system.SM.primaries]
    rowname(k) = k <= length(prim) ? prim[k] : "site-coupling row $(k - length(prim))"
    zmax = maximum(abs, z)
    rows = [k for k in eachindex(z) if abs(z[k]) > 1.0e-12 * zmax]
    amount(k) = round(bfix[k]; sigdigits = 4)
    why = if length(rows) == 1
        k = only(rows)
        z[k] > 0 ?
            "the budget asks for $(amount(k)) mol of $(rowname(k)), and every declared species holds it with a non-negative coefficient" :
            "the budget asks for $(amount(k)) mol of $(rowname(k)), and no declared species holds a positive amount of it"
    else
        Az = transpose(A) * z
        tie = [symbol(des.system.species[j]) for j in eachindex(Az) if abs(Az[j]) <= 1.0e-12 * zmax]
        "the combination " *
            join(["$(round(z[k] / zmax; sigdigits = 3)) × $(rowname(k)) (budget $(amount(k)) mol)" for k in rows], " + ") *
            " is negative for the budget and non-negative for every declared species" *
            (isempty(tie) ? "" : ", zero only for $(join(tie, ", "))")
    end
    cert = optimality_certificate(des, state; b = bfix, ϵ = ϵ)
    cert = merge(
        cert,
        (; optimal = false, budget_feasible = false, unplaceable = why, route = :infeasible, n_dual_solves = 0),
    )
    msg = "equilibrate_certified: no non-negative amounts of the declared species meet the " *
        "element budget: $why. There is no equilibrium to search for; the state is returned as given."
    _strict_convergence() && error(msg)
    @warn msg maxlog = 1
    return (copy(state), cert)
end

"""
    _ideal_mixing_system(cs) -> Union{ChemicalSystem, Nothing}

`cs` with every sublattice phase replaced by ideal mixing of the same
end-members, the species, their order and every other phase unchanged; `nothing`
when `cs` has no sublattice phase.
"""
_mixes_on_sites(ph) = model(ph) isa Union{SublatticeModel, CompoundEnergyModel}
_has_site_mixing(cs::ChemicalSystem) =
    cs.solid_solutions !== nothing && any(_mixes_on_sites, cs.solid_solutions)

function _ideal_mixing_system(cs::ChemicalSystem)
    ss = cs.solid_solutions
    site_model = _mixes_on_sites
    _has_site_mixing(cs) || return nothing
    ideal = [
        site_model(ph) ?
            SolidSolutionPhase{eltype(ph.end_members), IdealSolidSolutionModel}(
                ph.name, ph.end_members, IdealSolidSolutionModel(), ph.instances,
                ph.max_instances, ph.declared,
            ) : ph
            for ph in ss
    ]
    return ChemicalSystem(
        (f === :solid_solutions ? ideal : getfield(cs, f) for f in fieldnames(typeof(cs)))...,
    )
end

"""
    _ideal_mixing_start(state, model, bfix, ϵ, constraint, verbose; kwargs...)
        -> Union{ChemicalState, Nothing}

A certified answer to the same problem with every sublattice phase mixing its
end-members ideally, as a state of `state.system`, to start the search from.
`nothing` when the system has no sublattice phase, or when that solve does not
certify.

The two problems differ by the site-mixing term alone, and the one without it is
the solve the package has always done. Measured on a CEM I paste with the CNASH
gel of [Myers2014](@citet) (CEMDATA18, its eq. C.1 activity model): certified in 20.5 s
from the state as given, in 0.67 s from this start, which costs 0.25 s itself.
"""
function _ideal_mixing_start(
        state::ChemicalState, model, bfix, ϵ::Float64, constraint,
        verbose::Bool; kwargs...,
    )
    cs = _ideal_mixing_system(state.system)
    cs === nothing && return nothing
    # A starting point, not a result: the strict flag is for the answer the
    # caller receives, and a refusal here only means no start.
    st = ChemicalState(cs; T = state.T[1], P = state.P[1], n = state.n)
    eq0, cert0 = with(_STRICT_OVERRIDE => false, _AUTO_SPLIT => false) do
        equilibrate_certified(
            st; model = model, b = bfix, ϵ = ϵ, constraint = constraint,
            verbose = false, autostart = true, kwargs...,
        )
    end
    return cert0.optimal ? ChemicalState(state.system; T = state.T[1], P = state.P[1], n = eq0.n) : nothing
end

"""
    _REPAIR_FRACTION

What fraction of the amount the recipe could make of a missing phase
[`_repair_start`](@ref) puts in. Far enough off the boundary for the active-set
loop to work with, small enough not to pretend it knows the answer.
"""
const _REPAIR_FRACTION = 0.1

"""
    _repair_round(eq, cert, model, bfix, ϵ, solve_from, verbose)
        -> (eq, cert, improved)

One round of "the certificate names a missing phase, so put it in and solve
again". Returns the better of the two answers and whether the round bought
anything, so the caller can stop as soon as it does not.

`solve_from` is the search, injected: given a starting composition it returns
`(state, certificate)`. Passing it in rather than closing over the caller's
solver is what makes this round testable on its own — the situation it exists
for, a back end that converges onto the wrong active set, needs a system of some
135 species to arise, while the round's logic needs eight.
"""
function _repair_round(eq, cert, model, bfix, ϵ::Float64, solve_from, verbose::Bool)
    fixed = _repair_start(eq, model, bfix, ϵ)
    fixed === nothing && return (eq, cert, false)
    verbose && @info "repairing a missing phase" worst_si = cert.worst_supersaturation
    eq2, cert2 = solve_from(fixed)
    improved = cert2.optimal || _kkt_error(cert2) < _kkt_error(cert)
    kept, kept_cert = _keep_better(eq, cert, eq2, cert2)
    return (kept, kept_cert, improved)
end

"""
    equilibrate_split(state; model, b, maxpasses = 3, share = 0.5, kwargs...)
        -> (state, certificate)

The certified equilibrium of a system whose mixing phases may **unmix**, found by
giving a phase that wants to split a second composition to split into.

# The problem this solves

Inside a miscibility gap the Gibbs minimum of a mixing phase is two coexisting
compositions, not one. A formulation carrying one amount per species describes
that by declaring the phase twice — `SolidSolutionPhase(...; instances = 2)` —
and when the element balance **pins** the phase's overall composition inside the
gap, declaring it is enough: the minimization separates the two instances onto
the common-tangent pair by itself, and the certificate proves it. Measured on a
calcite/magnesite binary, where 0.025 mol of each fixes x̄ = 1/2 whatever the
energetics say, the instances land on the binodal to within 1e-3 and in the
proportions the lever rule asks for.

This function is for the case where nothing pins it. In a cement paste the AFm
composition is free — the sulfate has ettringite to go to and the hydroxide is
abundant — so **two instances started at the same composition stay there**. The
symmetric state satisfies every first-order condition jointly, so it is a
stationary point of the minimization, and no descent direction leads away from it
however unstable it is. The composition is typically *metastable* rather than
unstable — outside the spinodal, inside the binodal — where reaching the pair
needs a finite jump and not a gradient step.

# What it does

Michelsen's stability analysis [Michelsen1982](@cite), which is what the
certificate already runs on every present mixing phase, does not only answer
*whether* a phase splits: the trial composition that minimizes the tangent-plane
distance is an estimate of the **incipient phase**, and that is what seeds the
second instance here. The pass is then repeated until the certificate accepts or
stops improving.

The seed comes from the analysis of the **full system** — the trial composition
is computed with the chemical potentials the pore solution actually has — and not
from the mixing model alone. [`common_tangent`](@ref) gives the binodal of the
*isolated* binary, which is a different pair and is the wrong place to start
from: seeding there was measured to collapse straight back to one composition.

# What it returns, and what it promises

The best certified answer found, or — if none certifies — the last one, exactly
as [`equilibrate_certified`](@ref) would. A pass is kept only when the
certificate's **KKT error** improves — the worst of stationarity, element
balance, supersaturation and the constraint residual, not one of them — so the
result is never worse than the answer without splitting.

`share` is how much of the phase's amount is moved into the incipient instance on
each pass; `maxpasses` bounds the work.

Under [`STRICT_CONVERGENCE`](@ref) the passes are searches, not results: they run
with strictness suspended, and the answer they end on is judged strictly — an
error if it does not certify. Until 0.25.1 the flag was honored by the first
pass, which is by construction the one expected not to certify, so the function
raised before it had seeded anything.

!!! note "This is where convexity has already been given up"
    A phase that unmixes has a concave mixing energy, so `G` is not convex and
    `cert.optimal` no longer proves a *global* minimum — it proves a KKT point
    whose present phases are additionally stable against splitting, which is
    strictly more than stationarity gives. `SolidSolutionPhase` refuses such a
    model unless `instances > 1` is asked for, so a system reaching this function
    was built deliberately.

See also: [`equilibrate_certified`](@ref), [`common_tangent`](@ref),
[`miscibility_split`](@ref).
"""
function equilibrate_split(
        state::ChemicalState;
        model::AbstractActivityModel = DiluteSolutionModel(),
        b = nothing,
        maxpasses::Int = 3,
        share::Float64 = 0.5,
        kwargs...,
    )
    0 < share < 1 || throw(
        ArgumentError("`share` must lie strictly between 0 and 1, got $share."),
    )
    # Every pass is a search. Under the strict flag the first one -- the pass
    # this function exists to improve on -- would raise before any seeding, so
    # strictness is suspended while they run and applied to the final answer,
    # as `equilibrate_certified` does with its own starting routes.
    best_eq, best_cert = _relaxed_convergence() do
        _split_passes(state, model, b, maxpasses, share; kwargs...)
    end
    if best_cert !== nothing && !best_cert.optimal && _strict_convergence()
        error(
            "equilibrate_split: no pass produced a certifiable equilibrium: " *
                "stationarity $(best_cert.stationarity), " *
                "$(_balance_text(best_cert)), worst supersaturation " *
                "$(best_cert.worst_supersaturation). " *
                "`ChemistryLab.STRICT_CONVERGENCE[]` is set, so this raises rather " *
                "than returning an answer that is not an equilibrium.",
        )
    end
    _check_solvent(best_eq)
    return best_eq, best_cert
end

function _split_passes(state, model, b, maxpasses, share; kwargs...)
    eq, cert = equilibrate_certified(state; model = model, b = b, kwargs...)
    (cert === nothing || cert.optimal) && return eq, cert

    cs = state.system
    twin, untwin = _instance_pairs(cs)
    isempty(twin) && return eq, cert

    best_eq, best_cert = eq, cert
    for _ in 1:maxpasses
        trials = get(best_cert, :split_trials, nothing)
        (trials === nothing || isempty(trials)) && break

        n = Float64[ustrip(us"mol", x) for x in best_eq.n]
        _seed_split!(n, trials, twin, untwin, share) || break

        seeded = ChemicalState(cs, n .* u"mol"; T = temperature(best_eq), P = pressure(best_eq))
        # `autostart = false`: the seed IS the information, and the route search
        # would discard it for a start of its own choosing.
        eq2, cert2 = equilibrate_certified(
            seeded; model = model, b = b, autostart = false, kwargs...,
        )
        # Ranked by the package's own KKT error and not by one residual of it.
        # `worst_supersaturation` alone would accept a pass that improved the
        # saturation indices while losing moles of an element — the very failure
        # `_kkt_error` exists to prevent.
        _kkt_error(cert2) < _kkt_error(best_cert) || break
        best_eq, best_cert = eq2, cert2
        best_cert.optimal && break
    end
    return best_eq, best_cert
end

"""
    _instance_pairs(cs) -> (twin, untwin)

The species index of each `#2` twin, by the index of the species it copies, and
the inverse map.

Both directions, because Michelsen's trial can be reported on **either** instance
of a pair and the pair has to be recovered from whichever it names. Empty when
the system was not declared with `instances = 2`, which is how
[`equilibrate_split`](@ref) knows there is nowhere to split into.
"""
function _instance_pairs(cs)
    names = String[String(symbol(s)) for s in cs.species]
    twin = Dict{Int, Int}()
    untwin = Dict{Int, Int}()
    for (i, nm) in enumerate(names)
        j = findfirst(==(nm * "#2"), names)
        j === nothing && continue
        twin[i] = j
        untwin[j] = i
    end
    return twin, untwin
end

"""
    _seed_split!(n, trials, twin, untwin, share) -> Bool

Move material from the fuller instance of each flagged phase into the emptier
one, **at the composition Michelsen's analysis asks for**. Mutates `n` and
returns whether anything moved.

# The element budget is not touched, and that is the whole design

The transfer takes the incipient composition `t.x` out of one instance and puts
the same `t.x` into the other, so the total of every end-member across the pair
is unchanged and `A n` is exactly what it was. The obvious alternative — remove
at the DONOR's composition and add at the trial's — moves the same number of
moles but a different mixture, so it silently rewrites the element budget the
solver is about to be measured against. Two end-members of one binary are
different substances; `C4AH13` and `monosulphate12` do not have the same sulfur.

`move` is then bounded so that no end-member of the donor goes negative, which
is what makes `share` a request rather than a command.

# From the fuller into the emptier

Not the other way round, and not from whichever instance the trial happened to
name. The material is typically all in one of the two, and moving a share of an
instance holding 1.5e-4 mol while its twin holds 5.1e-2 is a perturbation of
three parts in a thousand — a seed that cannot move the answer is
indistinguishable from no seed at all.

# The receiver is emptied first

The two instances are copies of the same end-members, so moving the receiver's
whole content back into the donor leaves the budget as it was, and the receiver
then holds the trial composition alone. Without it the trial material mixes with
what the receiver already holds. Measured on two instances of a carbonate
binary that both sat at x = 0.117: the receiver came out at 0.209 instead of the
0.973 asked for, no pass reached the pair, and which instance was the donor was
decided by the rounding of two equal amounts. With the receiver emptied, the
seed is the same whichever instance gives, and the first pass certifies.
"""
function _seed_split!(
        n::AbstractVector{Float64}, trials, twin::Dict{Int, Int},
        untwin::Dict{Int, Int}, share::Float64,
    )
    moved = false
    seen = Set{Vector{Int}}()
    for (_, t) in trials
        # Recover the pair of instances from the one the trial names. The two
        # member lists are parallel — same end-member order — so `t.x` indexes
        # either of them.
        base = if all(haskey(twin, i) for i in t.members)
            collect(t.members)
        elseif all(haskey(untwin, i) for i in t.members)
            [untwin[i] for i in t.members]
        else
            continue                        # not an instanced phase
        end
        base in seen && continue            # both instances flag the same pair
        push!(seen, base)
        other = [twin[i] for i in base]

        n_base = sum(n[i] for i in base)
        n_other = sum(n[i] for i in other)
        total = max(n_base, n_other)
        total > 0 || continue
        from, to = n_base >= n_other ? (base, other) : (other, base)
        # The receiver is emptied into the donor first, so that it holds the
        # trial composition alone once the transfer is made.
        for (j, i) in enumerate(to)
            n[from[j]] += n[i]
            n[i] = 0.0
        end
        total = sum(n[i] for i in from)

        move = share * total
        for (j, i) in enumerate(from)
            t.x[j] > 0 || continue
            move = min(move, n[i] / t.x[j])
        end
        move > 0 || continue

        for (j, i) in enumerate(from)
            n[i] -= move * t.x[j]
        end
        for (j, i) in enumerate(to)
            n[i] += move * t.x[j]
        end
        moved = true
    end
    return moved
end

"""
    equilibrate_path(state, budgets; model, kwargs...) -> (states, certificates)

A **sequence** of certified equilibria, each one started from the last that
certified.

`budgets` is any iterable of element budgets — the vectors `equilibrate_certified`
takes as `b`. The first is solved from `state`; every later one is solved from the
previous answer, and a previous answer is reused only once the certificate has
accepted it. Where none has yet, `state` is used again.

# Why this exists

A cement equilibrium is hard to start cold and easy to start warm, and the gap is
not marginal. Measured on the 135-species paste of `scripts/ionic_hydration.jl`:

| | |
|:--|--:|
| cold start, the full multi-start cascade | 15.2 s |
| warm start from a neighboring answer | 0.19 s |

Eighty to one. So a sweep that rebuilds its state at every point pays the cold
price at every point, and — worse — can fail at one while both of its neighbors
certify, which is a *starting point* and not an infeasibility. Both blended-binder
sweeps in this package's documentation had such a point before they were written
this way.

# What it does and does not change

For a **convex** problem the minimum is unique, so walking to it cannot change
*what* is found — only whether the search finds it. That premise is checked
rather than assumed: [`SolidSolutionPhase`](@ref) refuses a mixing model whose
energy has a spinodal, so a system that was constructed at all is convex unless
the refusal was explicitly waived or the phase given two instances. Then this
becomes a genuine choice of branch, because inside a gap the starting point
decides which lobe the answer lands in — see [`common_tangent`](@ref).

The certificate still decides every point. A refused point is returned like any
other, with its certificate, and does **not** become the next start.

# Examples

```julia
budgets = [budget_at(f) for f in 0.0:0.05:0.30]
states, certs = equilibrate_path(fresh, budgets; model = HKFActivityModel())
all(c.optimal for c in certs)      # every point proved, not merely converged
```

See also: [`equilibrate_certified`](@ref), [`common_tangent`](@ref).
"""
function equilibrate_path(
        state::ChemicalState, budgets;
        model::AbstractActivityModel = DiluteSolutionModel(), kwargs...,
    )
    states = ChemicalState[]
    certs = Any[]
    warm = nothing
    for b in budgets
        eq, cert = equilibrate_certified(
            something(warm, state); model = model, b = b, kwargs...,
        )
        cert.optimal && (warm = eq)
        push!(states, eq)
        push!(certs, cert)
    end
    return states, certs
end

"""
    with_instances(state, cs) -> ChemicalState

`state` carried over to `cs`, a system [`with_instances`](@ref) built from its
own: every species keeps its amount, and the copies it did not have start at
zero, so the element budget is unchanged. Temperature and pressure are kept.
"""
function with_instances(state::ChemicalState, cs::ChemicalSystem)
    have = Dict(symbol(sp) => x for (sp, x) in zip(state.system.species, state.n))
    unit = first(state.n) * 0
    n = [get(have, symbol(sp), unit) for sp in cs.species]
    for (sp, x) in have
        haskey(cs.dict_species, sp) || throw(
            ArgumentError("with_instances: the species \"$sp\" of the state is not in the new system."),
        )
    end
    return ChemicalState(cs, n; T = temperature(state), P = pressure(state))
end

"""
    equilibrate_certified(state; model, ϵ, b, verbose, autostart, lp_start, dual,
                          fallback_model, fallback_on) -> (state, certificate)

Equilibrium composition together with a proof of its global optimality, obtained
by solving from every registered back end and keeping the answer
[`optimality_certificate`](@ref) proves optimal.

`state` is not modified: the answer is a new state, and the same `state` can be
solved again, or under another model, without a copy.

# The linear program over the pure phases comes first

With every activity set to one, the equilibrium is a linear program:
`minimize gᵀn subject to A n = b, n ≥ 0`. It is solved before anything else
(OptimaSolver's `lp_start`), and it does two things.

  - **It refuses an impossible budget.** When no non-negative amounts of the
    declared species meet `b`, the program proves it with a combination of the
    balances that every species raises and the budget lowers. The state is then
    returned as given, with `optimal = false`, `budget_feasible = false`,
    `route = :infeasible` and `unplaceable`, the reason in words; nothing is
    searched, and `STRICT_CONVERGENCE` raises.
  - **It gives the first start.** Its vertex, with every species it leaves out
    given the amount its multipliers give it in its own phase, is where the
    search begins: on two cold cements, 0.05 s and 0.26 s against 11.3 s and
    10.1 s without it, to the same composition.

`lp_start = false` skips it, and so does `autostart = false`. The certificate
reports `budget_feasible`, the start the answer came from as `route`
(`:lp_start`, `:state`, `:ideal_mixing`, `:ideal`, `:continuation`, `:restart`
or `:repair`) and the number of dual solves the search ran, `n_dual_solves`.

# The starting point is found, not asked for

When no back end certifies from the state as given, an initial approximation is
computed by continuation — [`homotopy_initial_state`](@ref) — and every back end
is run again from it. This is what makes a realistic cement solvable without the
caller knowing anything about the answer: from the cold state of a CEM I paste
(all the mass in the reactants, every product at the `ϵ` floor) no route reaches
the optimum, and with the continuation the same call certifies.

It costs nothing in the ordinary case, because it only runs when nothing else
certified. `autostart = false` declines it, which is what the coupled kinetic
step does: there the caller already supplies the previous instant as a warm
start, and a handful of extra solves inside an implicit ODE step would be paid
at every step.

# Options of the solvers

`dual` is a `NamedTuple` of keywords for the [`DualEquilibriumSolver`](@ref) the
route certifies with — `maxit`, `max_active_updates`, `inner_tol`, `inner_maxit`,
`tol`, `si_tol` — for instance `dual = (; maxit = 1000)`. Until 0.25.2 nothing
reached it. `tol` and `si_tol` are also the thresholds of the certificate, so
loosening them loosens the proof. The other keywords go to the interior-point
starts: `variable_space` sets their formulation, Ipopt takes the common
arguments of Optimization.jl (`maxiters`, `reltol`, `maxtime`, `verbose`), which
OptimizationIpopt maps to its options, and `OptimaOptimizer` ignores the rest.

`certificate.optimal == true` is a **proof** of a global minimum when the log
activities are the gradient of one Gibbs energy and that energy is proved convex
over the whole domain: ideal mixing, convex solid solutions, and the aqueous
models of a convexity bound (the ideal dilute model, Debye–Hückel with one ion
size and Davies on ions of charge two at most). `certificate.scope` says which
case holds: `:global_minimum`; `:local_minimum` when convexity is not proved
(a concave mixing phase, Pitzer, a convexity bound above one) and the Hessian
over the directions that conserve matter is positive definite at the answer, its
smallest eigenvalue being `reduced_curvature`; `:kkt_point` when that curvature
is not positive, or a constraint leaves the sufficiency unestablished;
`:self_consistent` when the activities are not the gradient of one energy, as
with the B-dot and Davies models at their default settings, and the answer is a
speciation consistent with its own activities.

When no route yields a proof, the answer with the smallest KKT error is returned,
its certificate says so, and a warning names the residual. That is the honest
outcome, and it is not the same thing as a failure: on a low-water cement the
returned composition satisfies the element balance to 1e-15 and has a
supersaturated phase left out, which the certificate reports as
`worst_supersaturation > 0`.

Requires `OptimaSolver` (the dual Newton lives there). Systems without an aqueous
phase, or without `H2O@`, cannot use the dual route; for those, this falls back to
the plain [`equilibrate`](@ref) and returns `nothing` as the certificate.

# A second activity model, when the first does not hold

`fallback_model` names an activity model to use when `model` does not give a
certified answer (`fallback_on = :refusal`, the default), or also when its
certified answer lies past the ionic strength `model` is stated valid for
(`fallback_on = :out_of_range`, see [`activity_model_range`](@ref)). The fallback
is solved from the first answer when it certified, from `state` otherwise. When
it certifies, its answer is returned, a warning names both models and the
reason, and the certificate carries `fallback_used = true`, `fallback_reason`
and the first certificate as `primary_certificate`. When it does not, the first
answer is returned with `fallback_used = false`. Without a fallback, nothing
changes.

```julia
eq, cert = equilibrate_certified(state; model = cemdata18_activity_model(),
                                 fallback_model = HKFActivityModel(),
                                 fallback_on = :out_of_range)
```

# A phase that may unmix: `instances = :auto`

A solid solution declared `SolidSolutionPhase(...; instances = :auto)` is solved
with one composition. When the certificate finds it wanting to split, the system
is rebuilt with a second instance ([`with_instances`](@ref)), the answer is
carried over, and the passes of [`equilibrate_split`](@ref) look for the pair.
Their answer is kept when its KKT error is smaller: it is then a state of the
enlarged system, and the certificate names the phases in `instances_added`.
Under [`STRICT_CONVERGENCE`](@ref) the first solve is a search and the final
answer is judged.

# Example

```julia
using ChemistryLab, OptimaSolver
eq, cert = equilibrate_certified(state)
cert.optimal          # true — proved globally optimal
cert.balance          # element balance residual, in mol
cert.balance_relative # the same, relative to what each row holds
cert.worst_supersaturation   # negative: every absent phase undersaturated
```
"""
function equilibrate_certified(state::ChemicalState; kwargs...)
    (_AUTO_SPLIT[] && _has_auto_instances(state.system)) || return _certified_with_fallback(state; kwargs...)
    # A phase declared `instances = :auto` may be solved first with one
    # composition inside its gap, which cannot certify. The search stops as soon
    # as a phase asks to split (see `solve_certified`), without the restarts that
    # would only confirm the refusal, and the second instance is tried here; the
    # answer is judged once that is done.
    eq, cert = _relaxed_convergence() do
        _after_first_search(state, _certified_with_fallback(state; kwargs...)...; kwargs...)
    end
    if cert !== nothing && !cert.optimal && _strict_convergence()
        error(
            "equilibrate_certified: no certifiable equilibrium, a second instance of the " *
                "`instances = :auto` phases included: stationarity $(cert.stationarity), " *
                "$(_balance_text(cert)). `ChemistryLab.STRICT_CONVERGENCE[]` is set, " *
                "so this raises rather than returning an answer that is not an equilibrium.",
        )
    end
    return eq, cert
end

# What follows the first search of a system with an `instances = :auto` phase:
# its answer when it certified or no such phase asks to split; otherwise the
# second instance, and when that does not help either, the ordinary search that
# the early stop cut short, so that `:auto` is never worse than one instance.
function _after_first_search(state, first_eq, first_cert; kwargs...)
    (first_cert === nothing || first_cert.optimal) && return first_eq, first_cert
    _wants_auto_split(state.system, first_cert) || return first_eq, first_cert
    quiet(f) = Base.CoreLogging.with_logger(f, Base.CoreLogging.NullLogger())
    split_eq, split_cert = quiet(() -> _auto_split(first_eq, first_cert; kwargs...))
    hasproperty(split_cert, :instances_added) && return split_eq, split_cert
    return _certified_with_fallback(state; kwargs..., split_early = false)
end

_has_auto_instances(cs) = cs.solid_solutions !== nothing &&
    any(ss -> _max_instances(ss) > _instances(ss), cs.solid_solutions)

_max_instances(ss) = hasproperty(ss, :max_instances) ? Int(getproperty(ss, :max_instances)) : _instances(ss)

# The declarations of `cs` that `cert` finds wanting to split and that were
# declared `instances = :auto`, by name.
function _auto_phases_to_split(cs, cert)
    grow = String[]
    phases = cs.solid_solutions
    trials = get(cert, :split_trials, nothing)
    (phases === nothing || trials === nothing) && return grow
    for (_, t) in trials
        k = findfirst(g -> first(t.members) in g, cs.ss_groups)
        k === nothing && continue
        d = only(ss for ss in phases if name(ss) == _declared(phases[k]))
        _max_instances(d) > _instances(d) && !(name(d) in grow) && push!(grow, name(d))
    end
    return grow
end

# The split test says something only about an answer that has otherwise reached
# the solution: on one that has not, a phase can look unstable for no other
# reason. So a phase is given its second instance only when the answer is
# stationary and balanced, and the split is its worst violation. Read on every
# certificate that failed, as it first was, the rule split the AFt of a slag
# cement paste on a first answer with an element balance of 0.17, and the search
# on the enlarged system then failed where one composition certifies.
function _converged_but_unstable(cert)
    all(k -> hasproperty(cert, k), (:stationarity, :balance, :worst_violation_split, :worst_supersaturation)) ||
        return false
    return cert.stationarity <= 1.0e-8 && _judged_balance(cert) <= 1.0e-8 &&
        cert.worst_violation_split >= cert.worst_supersaturation
end

_wants_auto_split(cs, cert) = _converged_but_unstable(cert) && !isempty(_auto_phases_to_split(cs, cert))

"""
    _auto_split(eq, cert; model, b, kwargs...) -> (state, certificate)

What `instances = :auto` asks for. When the certificate finds a phase declared
`:auto` wanting to split, the system is rebuilt with a second instance of it
([`with_instances`](@ref)), the answer is carried over with the copies at zero,
and the passes of [`equilibrate_split`](@ref) run on the enlarged system, seeded
from the composition the stability test gives. Their answer is kept only when its
KKT error is smaller, so the result is never worse than the first answer; it then
carries `instances_added`, the phases that were given a second instance.
"""
function _auto_split(
        eq, cert; model::AbstractActivityModel = DiluteSolutionModel(), b = nothing,
        maxpasses::Int = 3, share::Float64 = 0.5, kwargs...,
    )
    (cert === nothing || cert.optimal) && return eq, cert
    cs = eq.system
    grow = _auto_phases_to_split(cs, cert)
    isempty(grow) && return eq, cert
    declaration(n) = only(ss for ss in cs.solid_solutions if name(ss) == n)
    T = ustrip(us"K", temperature(eq))
    enlarged = with_instances(cs, (n => _max_instances(declaration(n)) for n in grow)...; T = T)
    eq2, cert2 = _split_passes(with_instances(eq, enlarged), model, b, maxpasses, share; kwargs...)
    (cert2 !== nothing && _kkt_error(cert2) < _kkt_error(cert)) || return eq, cert
    return eq2, merge(cert2, (; instances_added = grow))
end

function _certified_with_fallback(
        state::ChemicalState;
        model::AbstractActivityModel = DiluteSolutionModel(),
        fallback_model::Union{Nothing, AbstractActivityModel} = nothing,
        fallback_on::Symbol = :refusal,
        kwargs...,
    )
    fallback_on in (:refusal, :out_of_range) || throw(
        ArgumentError("fallback_on must be :refusal or :out_of_range; got :$fallback_on"),
    )
    fallback_model === nothing && return _equilibrate_certified(state; model = model, kwargs...)

    # The first pass is silent and never raises: whether it failed is decided
    # here, and said once, with the fallback's outcome.
    quiet(f) = Base.CoreLogging.with_logger(f, Base.CoreLogging.NullLogger())
    strict = _strict_convergence()
    eq, cert = with(_STRICT_OVERRIDE => false) do
        quiet(() -> _equilibrate_certified(state; model = model, kwargs...))
    end
    cert === nothing && return (eq, cert)
    reason = !cert.optimal ? :refusal :
        (fallback_on === :out_of_range && cert.within_activity_range === false) ? :out_of_range :
        nothing
    reason === nothing && return (eq, merge(cert, (; fallback_used = false)))

    # From the certified answer when there is one, which is then only outside
    # the range of its model, and from the caller's state otherwise.
    start = cert.optimal ? eq : state
    eq2, cert2 = with(_STRICT_OVERRIDE => false) do
        quiet(() -> _equilibrate_certified(start; model = fallback_model, kwargs...))
    end
    why = reason === :refusal ?
        "it did not certify (stationarity $(cert.stationarity), $(_balance_text(cert)))" :
        "its answer lies at an ionic strength of $(round(cert.ionic_strength; sigdigits = 3)) mol/kg, " *
        "past the $(cert.activity_range) mol/kg the model is stated valid for"
    if cert2 !== nothing && cert2.optimal
        @warn "equilibrate_certified: the answer is that of the fallback activity model " *
            "$(nameof(typeof(fallback_model))), because with $(nameof(typeof(model))) $why." maxlog = 1
        return (eq2, merge(cert2, (; fallback_used = true, fallback_reason = reason, primary_certificate = cert)))
    end
    msg = "equilibrate_certified: $(nameof(typeof(model))) was not kept because $why, and the " *
        "fallback $(nameof(typeof(fallback_model))) did not certify either (stationarity " *
        "$(cert2 === nothing ? NaN : cert2.stationarity)); returning the first answer."
    strict && !cert.optimal && error(msg)
    @warn msg maxlog = 1
    return (eq, merge(cert, (; fallback_used = false, fallback_reason = reason, fallback_certificate = cert2)))
end

function _equilibrate_certified(
        state::ChemicalState;
        model::AbstractActivityModel = DiluteSolutionModel(),
        b = nothing,
        ϵ::Float64 = _AMOUNT_FLOOR,
        verbose::Bool = false,
        constraint::EquilibriumConstraint = FixedTP(),
        parameters::Union{Nothing, Base.RefValue} = nothing,
        autostart::Bool = true,
        lp_start::Bool = true,
        dual::NamedTuple = NamedTuple(),
        split_early::Bool = true,
        kwargs...,
    )
    _refuse_state_keywords(kwargs, "equilibrate_certified")
    if !_DUAL_AVAILABLE[]
        error(
            "equilibrate_certified needs `OptimaSolver`: the KKT solver and the " *
                "certificate live there. Add `using OptimaSolver` — and " *
                "optionally `using Optimization, OptimizationIpopt` as a second " *
                "starting route, since neither back end certifies every case.",
        )
    end

    if !_dual_applicable(state.system)
        constraint isa FixedTP || throw(
            ArgumentError(
                "a constraint other than `FixedTP` needs the dual route, which " *
                    "requires an aqueous phase and `H2O@` among the species: the " *
                    "prescribed property is an unknown of that solver's system, " *
                    "and the interior-point back ends have nowhere to put it.",
            )
        )
        # No dual route: return the plain answer and say there is no proof.
        eq = equilibrate(state; model = model, ϵ = ϵ, certify = false, kwargs...)
        return (eq, nothing)
    end

    # Duals take the implicit-function route, dispatched on the state's element
    # type rather than tested for. See `_certified_primal_then_derivative`.
    dual_route = _certified_dual_route(
        _input_number_type(state, b; model, constraint), state, model, b, ϵ, verbose, constraint,
        parameters, (; dual = dual, kwargs...),
    )
    dual_route === nothing || return dual_route

    des = DualEquilibriumSolver(state.system, model; verbose = verbose, dual...)

    # `b` is fixed ONCE, from the state as given. Letting each start define its
    # own would pose a different problem for each: a start that violates the
    # balance — the interior point does, by 3e-6 mol on this class of problem —
    # shifts the component totals by exactly its own infeasibility, and the dual
    # solve then certifies the answer to the shifted problem. Measured, that gave
    # two "certified" compositions 0.2 % apart on dissolved calcium, which on a
    # convex problem with one minimum can only mean two different problems.
    bfix = b === nothing ?
        des.A * Float64[ustrip(us"mol", x) for x in state.n] :
        Float64.(collect(b))

    # The linear program over the pure phases, solved before anything else: it
    # proves a budget infeasible when it is, and otherwise gives a start. See
    # `_linear_program`.
    lp = (autostart && lp_start && constraint isa FixedTP && isempty(state.system.site_groups)) ?
        _linear_program(des, state, bfix) : nothing
    if lp !== nothing && lp.start.status === :infeasible
        return _refuse_infeasible_budget(des, state, bfix, lp, ϵ)
    end

    s = _CertifiedSearch(des, state, model, bfix, ϵ, constraint, parameters, verbose, dual, kwargs, split_early)
    starts = _first_starts(s, lp, autostart)
    eq, cert = _search(s, starts)
    split_next = _stops(s, cert)
    if autostart && !cert.optimal && !split_next && !(model isa DiluteSolutionModel)
        eq, cert = _from_ideal_answer(s, eq, cert, starts)
    end
    # What the automatic start did, in words, for the diagnostic of a refusal.
    note = autostart ? "not reached (the first route certified)" : "declined (autostart = false)"
    if autostart && !cert.optimal && !split_next
        eq, cert, note = _from_continuation(s, eq, cert, starts)
        eq, cert = _restarts(s, eq, cert)
        eq, cert = _repairs(s, eq, cert, starts)
    end
    !cert.optimal && !split_next && _uncertified(eq, cert, model, note)
    _check_solvent(eq)
    eq, cert = _reported(des, eq, cert, true; ϵ, floor = _CERTIFICATE_FLOOR, constraint)
    return (eq, merge(cert, (; route = _route(s, eq), n_dual_solves = length(s.memo), budget_feasible = true)))
end

# ── The stages of the search ──────────────────────────────────────────────────

# One certified search: the problem posed once (the solver `des` and the budget
# `bfix`), the dual solves already made (`memo`), and the route each start came
# by (`route_of`). The same cached starts are offered again after the
# continuation, after each restart and in each repair round, and a start already
# solved is not solved twice: `memo` returns what it gave the first time, which is
# the same answer bit for bit, since `des`, `bfix`, `ϵ` and `constraint` are fixed
# for the whole call. See `solve_certified`.
#
# The fields are of abstract types on purpose: each stage is then compiled once,
# and what it calls is compiled when it is first called. Typed by parameters,
# the stages let inference reach the fallbacks a certified solve never takes,
# and the first equilibrium of a cement compiled 20 s longer (172.8 against
# 152.4 s); with these fields, 158.6 s. A dynamic call costs nothing beside a
# solve.
struct _CertifiedSearch
    des::DualEquilibriumSolver
    state::ChemicalState
    model::AbstractActivityModel
    bfix::Vector{Float64}
    ϵ::Float64
    constraint::EquilibriumConstraint
    parameters::Union{Nothing, Base.RefValue}
    verbose::Bool
    dual::NamedTuple
    kwargs::Any
    stop::Any
    memo::IdDict{Any, Any}
    route_of::IdDict{Any, Symbol}
end

function _CertifiedSearch(des, state, model, bfix, ϵ, constraint, parameters, verbose, dual, kwargs, split_early)
    system = state.system
    # A phase declared `instances = :auto` that asks to split cannot certify with
    # one composition from any start, so the search stops there and the answer
    # goes to `equilibrate_certified`, which gives the phase its second instance.
    stop = split_early && _AUTO_SPLIT[] && _has_auto_instances(system) ?
        (c -> !c.optimal && _wants_auto_split(system, c)) : nothing
    return _CertifiedSearch(
        des, state, model, bfix, ϵ, constraint, parameters, verbose, dual, kwargs, stop,
        IdDict{Any, Any}(), IdDict{Any, Symbol}(),
    )
end

_stops(s::_CertifiedSearch, cert) = s.stop !== nothing && s.stop(cert)

# Every back end's answer from `from`, and `from` itself — the only start
# available if they all threw.
#
# Strict convergence is suspended for the duration, through
# `_relaxed_convergence` and so for this task alone: what this computes is a
# STARTING POINT, not a result. Left set, a back end that
# reports `MaxIters` raises, the `catch` below swallows it, and the search
# silently loses that candidate — so a caller asking for strict results gets
# a *worse* search than a caller who did not. Measured on a CEM I paste where
# the interior point ends on `MaxIters`: with the flag set the route returned
# an element balance of 27.6 mol, and with it clear the very same call
# returned 1.8e-14. The result is still judged strictly, at the end of
# `_equilibrate_certified`, which is where the flag belongs.
function _starts_from(
        s::_CertifiedSearch, from::ChemicalState, what::AbstractString, route::Symbol;
        factories::Vector{Function} = copy(_SOLVER_FACTORIES), offer_tail::Bool = true,
    )
    s.route_of[from] = route
    solve_one = function (f)
        r = _relaxed_convergence() do
            try
                # A candidate is a start: polishing it here would run the
                # dual Newton the search runs on it anyway, twice.
                _exploring_starts() do
                    esolver = EquilibriumSolver(s.state.system, s.model, f(); s.kwargs...)
                    SciMLBase.solve(esolver, from; ϵ = s.ϵ, b = s.bfix, polish = false)
                end
            catch err
                s.verbose && @info "$what rejected" backend = f err
                nothing
            end
        end
        r === nothing || (s.route_of[r] = route)
        return r
    end
    return _LazyStarts(solve_one, factories, from, typeof(from)[], Ref(0), offer_tail)
end

# The certified solve over `starts`. Every candidate the search tries is a
# candidate, and a candidate that does not converge is what the search exists
# for: its diagnostics stay quiet, and the verdict on the answer is pronounced
# once, on the certificate. The search decides on the verdicts; what the proof
# covers is reported once, for the answer returned.
_search(s::_CertifiedSearch, starts) = _exploring_starts() do
    solve_certified(
        s.des, starts; b = s.bfix, ϵ = s.ϵ, constraint = s.constraint, parameters = s.parameters,
        memo = s.memo, stop = s.stop, report = false,
    )
end

# The starts of the first search: the vertex of the linear program when it has
# one, ideal mixing for a phase mixing on sites, and the state as given.
function _first_starts(s::_CertifiedSearch, lp, autostart::Bool)
    state_starts = _starts_from(s, s.state, "start", :state)
    starts = state_starts
    lp_first = false
    first_start = nothing

    # The vertex of the linear program, with the species it leaves out given the
    # amounts its multipliers give them in their phases, comes first: the answer
    # of the default back end from it, and nothing else. See `_lp_lifted_state`,
    # whose docstring records the slowdown the first form of the lifting caused.
    #
    # Only that one candidate, because it is the one that pays. Measured on seven
    # solves of five cements, it certified at once on the cold ones, where
    # nothing from the state as given does. The other two candidates the point
    # would give never certified before a start from the state did. The
    # interior point from it never certified at all, and cost 1.5 to 10 s. The
    # point itself certified twice, both times after a start from the state had.
    # On a warm start from a neighboring answer, those two cost 10 s the state
    # did not need, and on a paste nothing certifies from, 10 s before the
    # continuation.
    #
    # Against the same calls with `lp_start = false`, on two cements
    # (cement107, and a CEM I with the CNASH gel): cold 0.005 and 0.03 times the
    # time, warm (from the answer) 0.33 and 0.18, a neighbor (1 % more water,
    # from the answer) 0.32 and 0.02, every composition the same to 3e-10.
    default = _DEFAULT_SOLVER_FACTORY[]
    if lp !== nothing && lp.start.status === :optimal && default !== nothing
        lifted = _lp_lifted_state(s.des, s.state, lp, s.ϵ)
        first_start = _starts_from(
            s, lifted, "start from the linear program", :lp_start;
            factories = Function[default], offer_tail = false,
        )
        starts = Iterators.flatten((first_start, starts))
        lp_first = true
    end

    # Ideal mixing as a stepping stone for a sublattice phase: see
    # `_ideal_mixing_start`. After the start of the linear program and before the
    # state as given, and solved only if the search reaches it: it is a whole
    # certified solve of its own, and from a neighboring answer, or wherever the
    # start of the linear program certifies, it is not needed. Measured on a CEM I
    # paste with the CNASH gel on its sites, ten warm restarts on neighboring
    # budgets: computed first, it was 47 % of their time.
    if autostart && _has_site_mixing(s.state.system)
        mixed = _DeferredStarts() do
            m = _ideal_mixing_start(s.state, s.model, s.bfix, s.ϵ, s.constraint, s.verbose; dual = s.dual, s.kwargs...)
            m === nothing ? nothing : _starts_from(s, m, "start from ideal mixing", :ideal_mixing)
        end
        starts = lp_first ? Iterators.flatten((first_start, mixed, state_starts)) :
            Iterators.flatten((mixed, state_starts))
    end
    return starts
end

# The ideal model as a stepping stone.
#
# A start near the answer is what this problem needs, and the cheapest good
# one is the answer to an easier question: the same minimization under ideal
# activities, which has no activity coefficients to make the residual depend
# on the composition and certifies where the non-ideal model does not. Its
# assemblage is the right one -- the phases present differ from the non-ideal
# answer by their amounts, not by their identity -- so the non-ideal solve
# starts with the correct active set instead of discovering it. See
# `_ideal_start`, whose docstring records the ulp sensitivity it fixed.
#
# Without it, a cold 109-species CEM IV paste without ash does not certify
# from either back end (dual balances 4.5 and 0.12); from the ideal answer it
# certifies at 3.6e-15, pH 13.444.
#
# Only when nothing else certified, so the ordinary case pays nothing, and
# guarded against recursion: the inner call is already ideal.
function _from_ideal_answer(s::_CertifiedSearch, eq, cert, starts)
    ideal = _ideal_start(s.state, s.model, s.bfix, s.ϵ, s.constraint, s.verbose; dual = s.dual, s.kwargs...)
    ideal === nothing && return eq, cert
    return _keep_better(
        eq, cert,
        _search(s, Iterators.flatten((_starts_from(s, ideal, "start from the ideal answer", :ideal), starts)))...,
    )
end

# An automatic initial approximation, computed rather than asked for.
#
# Only when nothing above certified, so the common case pays nothing for it.
# A realistic cement does not converge from the state as given — all the mass
# in the reactants, every product at the `ϵ` floor — and the caller should
# not have to know that, nor supply a chemically informed guess.
# `homotopy_initial_state` walks the solute amount up from a dilute system,
# which costs a handful of extra solves and needs nothing from the caller.
# It says what it did, in words, for the diagnostic of a refusal: when a
# solve fails on one machine and not another, the first thing anyone needs to
# know is whether the continuation ran at all and whether it helped — and
# asking for that should not require a second run with `verbose = true`.
function _from_continuation(s::_CertifiedSearch, eq, cert, starts)
    # Walked under the IDEAL model, deliberately, whatever `model` is: the
    # non-ideal ones do not walk (the a = 0 Debye-Huckel runs away to
    # I = 18 mol/kg, its coefficients falling with I raising solubility
    # raising I). The ideal endpoint is then a good start for `model`,
    # which is what the back-end loop below does with it.
    guess = homotopy_initial_state(s.state; ϵ = s.ϵ, verbose = s.verbose)
    guess === nothing && return eq, cert,
        "the continuation produced no usable start: every rung was " *
        "refused, or the system has no aqueous solvent to walk"
    before = _kkt_error(cert)
    eq, cert = _keep_better(
        eq, cert,
        _search(s, Iterators.flatten((_starts_from(s, guess, "start from the continuation", :continuation), starts)))...,
    )
    note = cert.optimal ?
        "the continuation certified it" :
        (
            _kkt_error(cert) < before ?
            "the continuation improved the KKT error from $before to " *
            "$(_kkt_error(cert)) without certifying" :
            "the continuation ran and its answer was no better than " *
            "$before, so it was discarded"
        )
    return eq, cert, note
end

# Restart from the answer. The continuation ends on a composition that is
# nearly the equilibrium but not certifiably so, and one more solve from
# there closes the gap — measured on a CEM I paste under the per-species
# Debye-Huckel model, stationarity 9.9e-7 (uncertified) becomes 1.5e-16
# with the worst absent phase 1.4e-5 below saturation. It is the same
# observation that motivates the continuation, applied once more: a start
# near the answer is what this problem needs, and the best one available
# is the answer already in hand.
#
# Bounded, and it stops as soon as a round buys nothing, so a genuinely
# hard case costs a fixed handful of solves rather than looping.
function _restarts(s::_CertifiedSearch, eq, cert)
    for _ in 1:_MAX_RESTARTS
        cert.optimal && break
        eq2, cert2 = _search(s, _starts_from(s, eq, "restart from the answer", :restart))
        improved = cert2.optimal || _kkt_error(cert2) < _kkt_error(cert)
        eq, cert = _keep_better(eq, cert, eq2, cert2)
        improved || break
    end
    return eq, cert
end

# Act on what the certificate says. A positive worst supersaturation
# names a phase that should be present and is not, which no amount of
# restarting from the same active set will fix: see `_repair_start`.
# The accumulated starts go back in with the repaired composition.
# Measured on the paste, a solve from the repair start alone reaches the
# right assemblage -- 74.19 cm3, all the magnesium back in the
# hydrotalcite -- and still fails its certificate on an unrelated trace
# component: the recipe's 1e-9 mol of carbon is lost, leaving an element
# balance of exactly 1e-9 against a tolerance of 1e-10. Ranked on the
# worst residual, that answer loses to the very point it was meant to
# replace. Handing the search the repaired composition *and* the
# candidates it already had keeps the chemistry of the one and the trace
# components of the others.
function _repairs(s::_CertifiedSearch, eq, cert, starts)
    repair_search(f) = _search(s, Iterators.flatten((_starts_from(s, f, "repair start", :repair), starts)))
    for _ in 1:_MAX_RESTARTS
        cert.optimal && break
        eq, cert, improved = _repair_round(eq, cert, s.model, s.bfix, s.ϵ, repair_search, s.verbose)
        improved || break
    end
    return eq, cert
end

# `STRICT_CONVERGENCE[]` is honored here, not only on the interior-point
# retcode. A caller who sets it is asking that a non-converged solve
# never pass as a result, and an uncertified answer from this route is
# exactly that: it can violate the element balance by moles and still
# come back looking like an ordinary `ChemicalState` — measured, a paste
# returned with a balance off by 6.7 mol, every hydrate at zero and a
# table of amounts that reads as a result. A warning is the right default
# (the answer is still the best one found, and `optimality_certificate`
# audits it), but under the strict flag it must raise.
function _uncertified(eq, cert, model, note)
    msg = "no route produced a certifiable equilibrium: stationarity " *
        "$(cert.stationarity), $(_balance_text(cert)), worst " *
        "supersaturation $(cert.worst_supersaturation). Automatic initial " *
        "approximation: $note" * _activity_range_hint(eq, model)
    _strict_convergence() && error(
        msg * ". `ChemistryLab.STRICT_CONVERGENCE[]` is set, so this raises " *
            "rather than returning an answer that is not an equilibrium. " *
            "Audit it with `optimality_certificate`; " *
            "`homotopy_initial_state(state; verbose = true)` reports each rung."
    )
    @warn msg * "; returning the answer with the smallest KKT error — audit it with `optimality_certificate`" maxlog = 1
    return nothing
end

# Which start the answer came from, found by identity among the dual solves.
function _route(s::_CertifiedSearch, eq)
    for (s0, hit) in s.memo
        first(hit) === eq && return get(s.route_of, s0, :other)
    end
    return :other
end

"""
    _activity_range_hint(eq, model) -> String

A sentence for the refusal of `equilibrate_certified` when the answer it gives up
on has an ionic strength past the range the manual states for `model`
([`activity_model_range`](@ref)), and an empty string otherwise.

Measured on blended cement pastes taken to full reaction, the refusal can be the
activity model's rather than the solver's. The Debye–Hückel limiting law,
`HKFActivityModel(å = 0)`, has its `log γ` still falling
at several mol/kg, and the search does not conclude; with an ion size per ion the
same budgets certify. Without this sentence the message names a stationarity and
a balance, and nothing in it points at the activity model.
"""
function _activity_range_hint(eq::ChemicalState, model::AbstractActivityModel)
    r = activity_model_range(model)
    r === nothing && return ""
    # No aqueous phase, or an answer whose solvent is gone: nothing to say here,
    # and `_check_solvent` speaks for the second.
    I = try
        _plain(ionic_strength(eq))
    catch
        return ""
    end
    (isfinite(I) && I > r) || return ""
    tail = if model isa HKFActivityModel && model.å !== nothing
        "; with an ion size per ion, `HKFActivityModel()`, the same budget " *
            "may certify, and the activity model is then the thing to report"
    else
        "; a model stated for higher ionic strengths, such as " *
            "`PitzerActivityModel` with parameters for these ions, is the next " *
            "one to try, or a budget that releases less salt"
    end
    return ". The ionic strength it stopped at, $(round(I; sigdigits = 3)) mol/kg, " *
        "is past the $(r) mol/kg the manual states for " *
        "$(nameof(typeof(model))), so the failure may be the activity model's " *
        "rather than the solver's" * tail
end

"""
    _check_solvent(eq)

Say so when the answer has no solution left to be an answer about.

A certificate proves that a composition minimizes the Gibbs energy of the
problem as posed. It says nothing about whether the problem was posed inside the
model's domain, and there is one way to leave it that produces an ordinary-looking
`ChemicalState`: let the solids take all the water. Every aqueous quantity is
then computed per kilogram of a solvent that is not there.

Measured on a sealed cement paste below its stoichiometric water demand, at
w/c = 0.28: the free water goes to 6e-9 mol, the solvent falls to a **fifth** of
its own aqueous phase, and the ionic strength is reported as 409 mol/kg by a
Debye-Huckel model valid to about one. Nothing in the certificate objects,
because nothing is wrong with the minimization — the mix simply does not contain
enough water to be a solution chemistry problem.

Warned, not raised, under the default flag: the composition of the *solids* in
such a solve still carries the mass-balance information a caller may legitimately
want, and it is the caller who knows whether that is what they asked for. Under
`STRICT_CONVERGENCE[]` it raises, like any other answer that is not one.
"""
function _check_solvent(eq::ChemicalState)
    x_w = solvent_fraction(eq)
    x_w >= SOLVENT_FRACTION_FLOOR && return nothing
    msg = "the aqueous phase has effectively vanished: the solvent is only " *
        "$(round(x_w; sigdigits = 3)) of it in mole fraction, against a floor of " *
        "$(SOLVENT_FRACTION_FLOOR). Molality, ionic strength, activity and pH are " *
        "defined per kilogram of solvent and are meaningless here, and the " *
        "activity model is being evaluated far outside its range. The solids " *
        "have taken the water: this is a mix below its stoichiometric water " *
        "demand, not an equilibrium the model can describe"
    _strict_convergence() && error(
        msg * ". `ChemistryLab.STRICT_CONVERGENCE[]` is set, so this raises."
    )
    @warn msg maxlog = 1
    return nothing
end

"""
    _dual_applicable(system) -> Bool

Whether [`DualEquilibriumSolver`](@ref) can be built for `system`: it needs an
aqueous phase, and `H2O@` among the species, because it parameterizes the interior
variables by the solvent's chemical potential.
"""
function _dual_applicable(system::ChemicalSystem)
    isempty(system.idx_aqueous) && return false
    return haskey(system.dict_species, "H2O@")
end


"""
    _certified_dual_route(state, model, b, ϵ, verbose, constraint, parameters, kwargs)

`nothing` for real-valued data; the certified answer with its derivative attached
when the amounts, the temperature, the pressure or the budget `b` carry
`ForwardDiff.Dual` numbers. The derivative is that of the answer the search
certified, lifted one level of duals at a time (`_lift_equilibrium`), so nested
differentiations are each exact and kept apart by their tags.

Dispatched on the element type, positionally, so the choice is the type system's
and neither path pays for the other. Two paths are needed for a mathematical
reason, not for want of a generic element type: making the component totals
generic would let duals flow into the solver, and what came back would be the
derivative of the **algorithm** — an active set decided by sign tests, a line
search with branches, an iteration count that varies with the data — rather than
the derivative of the **solution**. The map `b ↦ n*(b)` is smooth only piecewise,
and on each piece the implicit function theorem gives its derivative at the
solution with the active set frozen. That is how `OptimaSolver` computes its own
`Sensitivity`, and how Optima does upstream.
"""
_amount_number_type(state::ChemicalState) = _number_type(eltype(state.n))
# The number type of a solve's inputs: the state, the budget, the activity model,
# the constraint, the thermodynamic data and the mixing models of the system.
# Inside a search run on the values (`_STRIP_TAGS`) the tags of the scope are
# not counted: the data and the models still carry them, but every output the
# search reads has them removed.
function _input_number_type(
        state::ChemicalState, b; model = nothing, constraint = nothing, captured = nothing,
    )
    D = promote_type(Float64, _amount_number_type(state))
    _carries_duals(b) && (D = promote_type(D, mapreduce(typeof, promote_type, (x for x in b if x isa ForwardDiff.Dual))))
    D = promote_type(
        D, captured !== nothing ? captured :
            promote_type(
                model === nothing ? Float64 : _captured_number_type(model),
                _captured_number_type(_MixingTerms(state.system)),
                # The capacities of the site families, which the conservation
                # matrix carries when sites follow their host.
                _captured_number_type(state.system.site_families),
            ),
    )
    constraint === nothing || (D = promote_type(D, _captured_number_type(constraint)))
    D = promote_type(D, eltype(_build_params(state).ΔₐG⁰overRT))
    return _type_strip(D, _STRIP_TAGS[])
end

_type_strip(::Type{T}, ::Tuple) where {T} = T
_type_strip(::Type{ForwardDiff.Dual{Tg, V, N}}, tags::Tuple) where {Tg, V, N} =
    Tg in tags ? _type_strip(V, tags) : ForwardDiff.Dual{Tg, _type_strip(V, tags), N}
_number_type(::Type{<:DynamicQuantities.AbstractQuantity{T}}) where {T} = T
_number_type(::Type{T}) where {T <: Real} = T

_certified_dual_route(::Type{<:Real}, state, model, b, ϵ, verbose, constraint, parameters, kwargs) =
    nothing

function _certified_dual_route(
        ::Type{D}, state, model, b, ϵ, verbose, constraint,
        parameters, kwargs,
    ) where {D <: ForwardDiff.Dual}
    # The search runs on the values of the outermost level of duals (`Tg`), the
    # state and budget stripped here and every other input inside the scope, so
    # that a differentiation nested in another is solved on the outer one's
    # duals and lifted from there: `_plain` would strip every level and lose the
    # inner derivatives.
    Tg = ForwardDiff.tagtype(D)
    qv = Ref{Any}(Float64[])
    eq_v, cert = with(_STRIP_TAGS => (_STRIP_TAGS[]..., Tg)) do
        _equilibrate_certified(
            _strip_state(state, Tg); model = model, ϵ = ϵ, verbose = verbose,
            constraint = constraint, parameters = qv,
            b = b === nothing ? nothing : _strip_tag(collect(b), Tg), kwargs...,
        )
    end
    cert === nothing && return (eq_v, cert)
    des = DualEquilibriumSolver(state.system, model; verbose = verbose, get(kwargs, :dual, NamedTuple())...)
    bd = b === nothing ? des.A * _build_n0(state) : collect(b)
    eq_d, q_d = _lift_equilibrium(des, state, eq_v, bd; ϵ = ϵ, constraint = constraint, q = qv[], strip_tag = Tg)
    parameters === nothing || (parameters[] = q_d)
    return (eq_d, cert)
end
