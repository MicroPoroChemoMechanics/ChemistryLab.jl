# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# ── dual_solver.jl ────────────────────────────────────────────────────────────
#
# Chemical equilibrium by exact Newton on the KKT system, certified.
#
# THE ALGORITHM IS NOT HERE. It lives in `OptimaSolver.dual_newton_solve`, where
# it belongs: the active set on the bound-constrained variables, the degeneracy
# test on the rows of `A`, the two-level Newton and the certificate are all
# statements about a convex program, not about chemistry. What this file supplies
# is the four things that ARE chemistry:
#
#   1. the conservation matrix `A` and the reference potentials `g = Δ_a G⁰/RT`;
#   2. the activity model, as the callback `h` with `∇f = g + h(n)`;
#   3. which species are strictly positive at any solution — the aqueous ones,
#      parameterized by `ln n` — and which may vanish — the pure phases;
#   4. which species plays the reference role: the solvent, whose activity is a
#      mole fraction and therefore bounded above, so that its stationarity cannot
#      be inverted and must be carried by the outer system.
#
# The mathematics is documented in the manual under
# "Proving that an answer is the answer", and the general form in OptimaSolver's
# theory page.

using LinearAlgebra

"""
    DualEquilibriumSolver(system, model; tol, maxit, max_active_updates, si_tol, verbose)

Equilibrium by Newton on the KKT system, in element-potential space.

Where [`EquilibriumSolver`](@ref) minimizes `G` by an interior-point method over
all species and stops on `MaxIters` for a cement, this solves the stationarity
conditions directly and returns a result that
[`optimality_certificate`](@ref) can **prove** optimal — the Gibbs problem being
convex, the KKT conditions are sufficient.

Writing `u = −Aᵀy` for the element potentials, an aqueous species obeys the
mass-action law `aᵢ = exp(uᵢ − gᵢ)` and a pure phase is present exactly when
`uᵢ = gᵢ`, absent when undersaturated: the classical phase-stability criterion.

Starting from the answer of an [`EquilibriumSolver`](@ref) is the intended use —
the interior-point method reaches a neighborhood, this reaches the conditions.

See also: [`optimality_certificate`](@ref), [`speciated_states`](@ref).
"""
struct DualEquilibriumSolver{L, M <: AbstractActivityModel}
    system::ChemicalSystem
    lna::L
    model::M
    idx_aq::Vector{Int}
    idx_pure::Vector{Int}
    j_solvent::Int
    ss_groups::Vector{Vector{Int}}   # end-members of each declared solid solution
    site_groups::Vector{Vector{Int}} # members of each surface site family
    A::Matrix{Float64}
    opts::NamedTuple
end

function DualEquilibriumSolver(
        system::ChemicalSystem,
        model::AbstractActivityModel = DiluteSolutionModel();
        tol::Float64 = 1.0e-10,
        maxit::Int = 200,
        max_active_updates::Int = 200,
        si_tol::Float64 = 1.0e-8,
        verbose::Bool = false,
    )
    idx_aq = [i for (i, s) in enumerate(system.species) if aggregate_state(s) == AS_AQUEOUS]
    idx_pure = [i for (i, s) in enumerate(system.species) if aggregate_state(s) != AS_AQUEOUS]
    isempty(idx_aq) && throw(
        ArgumentError(
            "DualEquilibriumSolver needs an aqueous phase: the element potentials are " *
                "determined by the mass-action laws of the aqueous species."
        )
    )
    jw = something(findfirst(i -> symbol(system.species[i]) == "H2O@", idx_aq), 0)
    jw == 0 && throw(
        ArgumentError(
            "DualEquilibriumSolver needs `H2O@` among the species: the solvent's " *
                "activity is a mole fraction, hence bounded above, so its stationarity " *
                "cannot be inverted and is carried by the outer system instead."
        )
    )
    # A solid solution is a MIXING phase, exactly like the aqueous one: its
    # end-members have composition-dependent activities, so none of them is ever
    # exactly absent while the phase exists, and the active set belongs at the
    # level of the phase. They are therefore removed from the bound-constrained
    # set, where they would have been treated as pure phases.
    ss = system.solid_solutions
    ss_groups = Vector{Int}[]
    if ss !== nothing
        byname = Dict(symbol(sp) => i for (i, sp) in enumerate(system.species))
        for phase in ss
            idx = [byname[symbol(em)] for em in phase.end_members if haskey(byname, symbol(em))]
            length(idx) == length(phase.end_members) && push!(ss_groups, idx)
        end
    end
    # A surface site family is a mixing phase too, and for a sharper reason: its
    # members' activities are site fractions, so none is ever exactly absent
    # while the family has a budget, and the budget is pinned by a conservation
    # row. Left in the bound-constrained set they would be pure phases of unit
    # activity — the site mixing never computed, the balance unsatisfiable.
    site_groups = copy(system.site_groups)
    in_mixing = Set(vcat(ss_groups..., site_groups...))
    idx_pure = [i for i in idx_pure if !(i in in_mixing)]

    return DualEquilibriumSolver(
        system, activity_model(system, model), model,
        idx_aq, idx_pure, jw, ss_groups, site_groups, Float64.(system.SM.A),
        (; tol, maxit, max_active_updates, si_tol, verbose),
    )
end

activity_model(des::DualEquilibriumSolver) = des.model

"""
    _dual_problem(des, p, n0) -> DualNewtonProblem

Package the chemistry as the convex program `OptimaSolver` solves. Built per
solve because the reference potentials `Δ_a G⁰/RT` depend on temperature and
pressure.
"""
function _dual_phases(des::DualEquilibriumSolver, n0)
    # One mixing phase for the aqueous solution — always present, the solvent as
    # its reference — and one more per declared solid solution, whose presence
    # the tangent-plane test decides.
    #
    # The reference is the member whose stationarity the OUTER system carries
    # instead of inverting, so it has to be the one the phase is mostly made of.
    # For the aqueous phase that is the solvent, and not by convention: water is
    # the only species there whose activity is a mole fraction, bounded above, so
    # its law cannot be inverted at all. Inside a solid solution every member is
    # a mole fraction and any of them could serve — but the choice is not
    # indifferent. The outer unknown is `ln x_ref`, and picking a member that
    # happens to be minor sends that unknown towards −∞ and the Jacobian with it.
    # So the reference is chosen by magnitude, from the composition the caller
    # supplied, which is what the solvent already is for the aqueous phase.
    # `split_starts` is empty here and carried all the same: the vector's element
    # type is fixed by its FIRST entry, so an aqueous phase without the field
    # would make every solid solution that has one unpushable. And the aqueous
    # phase cannot unmix anyway -- there is one solvent.
    phases = [
        (
            members = des.idx_aq, j_ref = des.j_solvent,
            always_present = true, mole_fraction = false,
            split_starts = Vector{Vector{Float64}}(),
        ),
    ]
    models = _ss_models(des)
    for (k, grp) in pairs(des.ss_groups)
        j_ref = argmax(@view n0[grp])
        push!(
            phases,
            (
                members = grp, j_ref = j_ref,
                always_present = false, mole_fraction = true,
                split_starts = _split_starts(get(models, k, nothing), length(grp)),
            ),
        )
    end

    # One more per surface site family, and `always_present = true` is not a
    # convenience here. With `false`, the seeding test asks whether the members
    # already hold more than 1e-9 mol; from a cold start they all sit at the
    # floor, the phase is never activated, and the conservation row
    # `Σ n_member = N_t` with `N_t > 0` is infeasible from the first Newton
    # iteration. `true` also exempts it from the drop rule, which is right: a
    # phase whose total is pinned must not be dropped because its members are
    # small.
    #
    # The reference is the **free site**, the first member by construction. It is
    # the state a fresh surface is mostly in, which is the same criterion the
    # solvent meets for the aqueous phase, and it is stable: a family cannot run
    # out of free sites without the occupied ones taking their place.
    #
    # A site family does not unmix, so no split start is offered.
    for grp in des.site_groups
        push!(
            phases,
            (
                members = grp, j_ref = 1,
                always_present = true, mole_fraction = true,
                split_starts = Vector{Vector{Float64}}(),
            ),
        )
    end
    return phases
end

"""
    _ss_models(des) -> Dict{Int, Any}

The mixing model of each entry of `des.ss_groups`, by position.

`ss_groups` is built by walking `system.solid_solutions` and **skipping** any
declaration whose end-members are not all in the species list, so the two lists
are the same length only when nothing was skipped. Rebuilding the correspondence
the same way is the only way to be sure a model is matched to its own group; a
positional zip would silently pair a phase with someone else's model the first
time a declaration is dropped.
"""
function _ss_models(des::DualEquilibriumSolver)
    out = Dict{Int, Any}()
    ss = des.system.solid_solutions
    ss === nothing && return out
    byname = Dict(symbol(sp) => i for (i, sp) in enumerate(des.system.species))
    k = 0
    for phase in ss
        idx = [byname[symbol(em)] for em in phase.end_members if haskey(byname, symbol(em))]
        length(idx) == length(phase.end_members) || continue
        k += 1
        out[k] = model(phase)
    end
    return out
end

"""
    _split_starts(model, nmembers) -> Vector{Vector{Float64}}

Where to look for the other lobe of a phase that may unmix, as mole fractions
over its members.

The tangent-plane search in `OptimaSolver` probes the **corners** of the
composition simplex and refines by successive substitution, which converges to
the stationary point nearest its start. A corner is usually on the right side of
the barrier and sometimes is not; where it is not, the iteration walks back to
the phase's own composition and reports nothing, though splitting would lower the
energy.

Measured on the AFm sulfate/hydroxide binary with the published Redlich-Kister
parameters (spinodal [0.631, 0.914], binodal [0.4999, 0.9700]):

| phase sits at | corners alone | with the binodal handed over |
|:--|--:|--:|
| x = 0.5268, where a CEM I settles | +3.99e-02, found | the same verdict |
| x = 0.95 | +2.4e-16, **missed** | +1.23e-01, trial x = 0.444 |
| x = 0.98, outside the binodal | stable | stable |

Both flagged compositions are metastable, so being metastable is not by itself
what defeats the corners — sitting in the lobe they lead back into is, and that
is not knowable in advance. Hence a start supplied unconditionally rather than a
rule for when to supply one.

[`common_tangent`](@ref) computes the binodal from the mixing model alone, in
microseconds and with no reference to the rest of the system, and the search then
refines it with the chemical potentials the system actually has — which is the
pair that matters. That division is the whole point: the model's binodal is a
good *place to look*, not the answer. Extra starts can only raise the maximum the
search returns, so they never take a verdict away and, as the last row shows, do
not invent one.

Empty for anything but a binary with a gap, which is the only case
`common_tangent` is defined for.
"""
function _split_starts(model, nmembers::Int)
    out = Vector{Vector{Float64}}()
    (model === nothing || nmembers != 2) && return out
    pair = try
        common_tangent(model)
    catch
        nothing
    end
    pair === nothing && return out
    for x in pair
        (isfinite(x) && 0 < x < 1) || continue
        push!(out, Float64[1 - x, x])
    end
    return out
end

function _dual_problem(des::DualEquilibriumSolver, p, n0, blocks = nothing)
    phases = _dual_phases(des, n0)
    bl = blocks === nothing ?
        _constraint_blocks(FixedTP(), des, nothing, p, n0) : blocks
    return _optima_dual_problem(
        des.A, Float64.(p.ΔₐG⁰overRT), des.lna, phases, des.idx_pure, p,
        bl.gq, bl.hq, bl.cq, bl.q0, bl.qscale, bl.Aq, Int[], nothing,
    )
end

"""
    SciMLBase.solve(des::DualEquilibriumSolver, state; b = nothing, ϵ = 1e-16)
        -> ChemicalState

Equilibrium composition. `state` supplies the temperature, the pressure and the
starting guess; `b` the component totals, defaulting to those of `state`.

The amounts are returned as computed, without a floor: an amount of `e⁻³⁰⁰` IS
the mass-action answer for a species that is not there, and raising it to `ϵ`
falsifies its activity by hundreds of `RT` units — which is enough to make
[`optimality_certificate`](@ref) report a residual of 74 on a composition solved
to `5e-12`.
"""
function SciMLBase.solve(
        des::DualEquilibriumSolver,
        state::ChemicalState;
        b = nothing,
        ϵ::Float64 = 1.0e-16,
        constraint::EquilibriumConstraint = FixedTP(),
        parameters::Union{Nothing, Base.RefValue} = nothing,
        surface_potential::Symbol = :auto,
    )
    surface_potential in (:auto, :unknown, :eliminated) || throw(
        ArgumentError(
            "surface_potential must be :auto, :unknown or :eliminated; " *
                "got :$surface_potential"
        ),
    )
    p = _build_params(state; ϵ = ϵ)
    n0 = Float64[ustrip(us"mol", x) for x in state.n]
    bv = b === nothing ? des.A * n0 : Float64.(collect(b))

    blocks = _constraint_blocks(constraint, des, state, p, n0)
    if surface_potential !== :eliminated
        surf = _surface_potential_blocks(des, state, p, n0)
        surface_potential === :unknown && surf === nothing && throw(
            ArgumentError(
                "surface_potential = :unknown was asked for, and no site family of " *
                    "this system needs one. Only a `DiffuseLayer` does."
            ),
        )
        blocks = _compose_blocks(blocks, surf)
    end
    res = _optima_dual_solve(_dual_problem(des, p, n0, blocks), bv, n0, des.opts)

    res.converged || begin
        Threads.atomic_add!(NONCONVERGED, 1)
        # Silent while a multi-start route is trying candidates: one of them not
        # converging is what the search is for, and the verdict belongs to the
        # certificate of the answer, not to a candidate.
        _EXPLORING_STARTS[] || @warn """the dual equilibrium solve did not certify \
        optimality; audit it with `optimality_certificate`.""" maxlog = 1
    end

    # The parameters the constrained solve FOUND — the temperature an adiabatic
    # solve reached, the titrant amount a prescribed pH required. They are part of
    # the answer, not a diagnostic, so the caller can ask for them.
    parameters === nothing || (parameters[] = copy(res.q))

    T_out, P_out = blocks.apply(temperature(state), pressure(state), res.q)
    return ChemicalState(
        des.system, [nᵢ * u"mol" for nᵢ in res.x];
        T = T_out, P = P_out,
    )
end

"""
    optimality_certificate(des, state; b = nothing, ϵ = 1e-16, floor = 1e-25)
        -> (; stationarity, balance, worst_supersaturation, n_interior,
             n_absent_component, param_residual, worst_violation_split,
             split_phases, split_trials, optimal)

Check the KKT conditions at a composition, independently of how it was obtained.

For a convex problem these conditions are sufficient, so `optimal = true` is a
proof of **global** optimality. Use it to audit any solver — including
[`EquilibriumSolver`](@ref), whose interior-point iteration reports `MaxIters` on
a cement equilibrium and cannot say whether the point it returns is the answer.

The three quantities are the stationarity of the interior species, the component
balance, and the worst saturation index among absent phases (negative when every
one of them is undersaturated, as optimality requires).

`worst_violation_split` extends that last test to the phases that are **present**:
Michelsen's tangent-plane distance, which asks whether a mixing phase would lower
the Gibbs energy by separating into two compositions. It is `-Inf` when no phase
could be tested, negative when every one of them is stable, and positive when one
wants to unmix — `split_phases` then names them and `split_trials` carries, per
phase, the composition it wants to split into. That composition is what
[`equilibrate_split`](@ref) seeds a second instance with; it exists nowhere else,
being a property of the full system and not of the mixing model alone.
"""
function optimality_certificate(
        des::DualEquilibriumSolver, state::ChemicalState;
        b = nothing, ϵ::Float64 = 1.0e-16, floor::Float64 = 1.0e-25,
        constraint::EquilibriumConstraint = FixedTP(),
        q = nothing,
    )
    p = _build_params(state; ϵ = ϵ)
    n = Float64[ustrip(us"mol", x) for x in state.n]
    bv = b === nothing ? des.A * n : Float64.(collect(b))

    # The certificate has to audit the problem that was SOLVED, and a constraint
    # is part of that problem. Rebuilt with `FixedTP` — which is what this did
    # until 0.16 — the audit misses two things at once: the conservation rows
    # lose their `Aq q` term, so a prescribed activity or pH is measured against
    # a budget short by exactly the titrant amount; and `hq` is skipped, so a
    # constraint that shifts a chemical potential is measured against the
    # unshifted one and can never certify however right it is.
    blocks = _constraint_blocks(constraint, des, state, p, n)
    qv = blocks.nq == 0 ? nothing :
        (q === nothing ? Float64.(collect(blocks.q0)) : Float64.(collect(q)))

    c = _optima_kkt_certificate(
        _dual_problem(des, p, n, blocks), n, bv, floor,
        des.opts.tol, des.opts.si_tol, qv,
    )
    return (;
        stationarity = c.stationarity, balance = c.feasibility,
        # The unscaled stationarity, in RT units. `stationarity` is divided by the
        # size of the potentials it is built from — they are of order 10²-10³, so
        # an absolute threshold on their residual would ask for thirteen digits of
        # cancellation — and this is the raw figure for reporting.
        stationarity_abs = c.stationarity_abs,
        worst_supersaturation = c.worst_violation, n_interior = c.n_interior,
        n_absent_component = c.n_forced_zero,
        # Zero on the unconstrained route, so it costs nothing there and is the
        # constraint's own residual when there is one.
        param_residual = hasproperty(c, :param_residual) ? c.param_residual : 0.0,
        # Michelsen's split verdict on the PRESENT mixing phases, forwarded
        # rather than dropped. `equilibrate_split` seeds a second instance from
        # `split_trials`, and it is the only place that composition exists: the
        # trial is a property of the full system, fixed jointly with the solution
        # the phase sits in, so a caller cannot recompute it from the mixing
        # model alone. Read through `hasproperty` because a certificate also
        # arrives from a back end that does not run the test.
        worst_violation_split = hasproperty(c, :worst_violation_split) ?
            c.worst_violation_split : -Inf,
        split_phases = hasproperty(c, :split_phases) ? c.split_phases : Int[],
        split_trials = hasproperty(c, :split_trials) ? c.split_trials :
            Dict{Int, NamedTuple{(:members, :x), Tuple{Vector{Int}, Vector{Float64}}}}(),
        optimal = c.optimal,
    )
end

"""
    _kkt_error(cert) -> Float64

How far a composition is from satisfying the KKT conditions: the worst of the
three residuals the certificate reports, in one number.

All of them, and not the stationarity alone. A composition can be stationary to
1e-3 while violating mass conservation by **moles** — measured, an answer with
stationarity 2.4e-3, an element balance off by 6.7 mol and a phase supersaturated
by 45 — and that is not a near-answer, it is not an answer to this problem at
all. Ranking on stationarity alone lets such a point beat a candidate that
conserves mass, which is how a multi-start search can discard the good answer it
just computed.

The constraint's own residual counts too, when there is one. Leaving it out would
rank a candidate that minimizes the Gibbs energy while violating the very
equation that makes it a *constrained* answer above one that satisfies both —
the same hole OptimaSolver records for the kinetic step, where "a march that
should have stopped at saturation dissolved everything and was proved optimal".
Read with `hasproperty` because a certificate also arrives from the kinetic
route, which builds its own.

`worst_supersaturation` is clamped at zero because a negative value is not an
error: it means every absent phase is undersaturated, as optimality requires.
"""
_kkt_error(cert) = max(
    cert.stationarity, cert.balance, max(cert.worst_supersaturation, 0.0),
    hasproperty(cert, :param_residual) ? cert.param_residual : 0.0,
)

"""
    solve_certified(des, starts; b = nothing, ϵ = 1e-16, floor = 1e-25)
        -> (state, certificate)

Solve from each starting composition in `starts` and return the first answer
[`optimality_certificate`](@ref) **proves** optimal. If none is proved, return the
one with the smallest KKT error, together with its certificate, so the caller sees
what it is getting.

# Why trying more than one start is the rigorous thing to do, not a fudge

The certificate is an **oracle**: for a convex problem it does not rank answers, it
decides them. Given an oracle, running several routes and keeping a proved answer
is not guesswork — the proof is the same proof whichever route produced the point,
and nothing about it depends on having predicted the winner. What would be a fudge
is picking a route by taste and reporting its output unproved.

It is also necessary, because no single route dominates. Measured on an LC³
equilibrium at three degrees of reaction, with the dual Newton started from each
interior-point back end in turn:

| degree of reaction | from `OptimaOptimizer` | from `IpoptOptimizer` |
|---|---|---|
| 0.05 | certified, 1.8e-12 | certified, 9.1e-13 |
| 0.25 | **not certified**, 7.2 | certified, 4.2e-12 |
| 1.00 | certified, 1.2e-11 | **not certified**, 3.5e-3 |

Each back end solves what the other misses. With both offered, all three are
proved.

# Example

```julia
using Optimization, OptimizationIpopt      # for IpoptOptimizer

starts = [
    SciMLBase.solve(EquilibriumSolver(cs, model, OptimaOptimizer()), st; b = b),
    SciMLBase.solve(EquilibriumSolver(cs, model, IpoptOptimizer()), st; b = b),
]
eq, cert = solve_certified(des, starts; b = b)
cert.optimal || @warn "no start produced a certifiable answer" cert
```

The starts are supplied by the caller rather than built here, so this adds no
dependency: whichever back ends are loaded are the ones available.
"""
function solve_certified(
        des::DualEquilibriumSolver, starts;
        b = nothing, ϵ::Float64 = 1.0e-16, floor::Float64 = 1.0e-25,
        constraint::EquilibriumConstraint = FixedTP(),
        parameters::Union{Nothing, Base.RefValue} = nothing,
    )
    best = nothing
    best_cert = nothing
    best_err = Inf
    best_q = Float64[]
    best_ok = false
    for s0 in starts
        # Each candidate's own parameters, captured here rather than written
        # straight into the caller's `Ref`. Written straight through, the `Ref`
        # would end up holding the LAST candidate's parameters while the state
        # returned is the best-by-error one — so a prescribed-pH scan would
        # report a titrant amount belonging to a different composition.
        qref = Ref(Float64[])
        eq = SciMLBase.solve(
            des, s0; b = b, ϵ = ϵ, constraint = constraint, parameters = qref,
        )
        # The certificate is evaluated at the T and P the constrained solve
        # FOUND, which `eq` carries — not at the ones the start had. `∇f` depends
        # on both, so certifying against the start's conditions would measure the
        # stationarity of a different problem. `q` is part of what the solve
        # found in exactly the same way, and is passed for the same reason.
        cert = optimality_certificate(
            des, eq; b = b, ϵ = ϵ, floor = floor,
            constraint = constraint, q = qref[],
        )
        if cert.optimal
            parameters === nothing || (parameters[] = qref[])
            return (eq, cert)
        end
        # Ranked on the KKT error, but only among answers the model can describe:
        # a composition whose solvent has been taken by the solids is outside the
        # formulation, and letting it win on a smaller residual hides every
        # candidate that conserved mass behind it. See `_within_domain`.
        err = _kkt_error(cert)
        ok = _within_domain(eq)
        if (ok && !best_ok) || ((ok == best_ok) && err < best_err)
            best_err = err
            best = eq
            best_cert = cert
            best_q = qref[]
            best_ok = ok
        end
    end
    parameters === nothing || (parameters[] = best_q)
    return (best, best_cert)
end

# ── hooks filled in by the OptimaSolver extension ────────────────────────────
#
# The algorithm is `OptimaSolver`'s, and that package is a weak dependency here,
# so the three entry points are indirected through functions the extension
# overrides. Called without it loaded, they say what is missing.

_optima_dual_problem(args...) = _need_optima()
_optima_dual_solve(args...) = _need_optima()
_optima_kkt_certificate(args...) = _need_optima()

_need_optima() = error(
    "DualEquilibriumSolver needs OptimaSolver ≥ 0.3: the KKT solver and its " *
        "certificate live there. Add `using OptimaSolver`."
)
