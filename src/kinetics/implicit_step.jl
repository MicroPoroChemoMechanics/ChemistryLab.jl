# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using LinearAlgebra

# ── implicit_step.jl ──────────────────────────────────────────────────────────
#
# A kinetic step as Leal et al. (2017) pose it, and as Reaktoro implements it: ONE
# fully implicit problem per step, not an operator split.
#
#     min G(n)   s.t.   [A ; Kᵀ]·n + [0 ; −I]·Δξ = [b₀ ; ξ₀]      (linear)
#                       Δξ − Δt·M·r(n) = 0,   M = KᵀK             (nonlinear)
#                       n ≥ 0
#
# The reaction extents are unknowns of the SAME system as the amounts and the
# element potentials, and the rate is evaluated at the unknown end-of-step state —
# backward Euler. Two consequences follow, and both are the point:
#
#   * the equilibrium species and the kinetic ones are solved TOGETHER, so there
#     is no frozen speciation in a right-hand side and no lag between the two;
#   * the reactivity constraint `Kᵀn − Δξ = ξ₀` is LINEAR, so it joins the
#     conservation block beside the elements and the charge. The algebraic cost of
#     kinetics is the number of REACTIONS, not the number of species.
#
# `M = KᵀK` and not the identity: the constraint defines `Δξ = Kᵀ(n − n₀)`, and
# with `n − n₀ = Kξ` that is `KᵀK ξ`. Writing `Δξ = Δt·r` instead is wrong by a
# factor `M`.
#
# This coexists with `integrate`, which puts the same chemistry in an `ODEProblem`
# and hands it to SciML. That route has adaptive time-stepping and a choice of
# stiff solvers; this one has an exact end-of-step equilibrium and no splitting
# error. Neither dominates, and they are checked against each other.

"""
    KineticStepSolver(system, model, reactions; opts...)

A solver for one fully implicit kinetic step, following [Leal2017](@citet).

`reactions` are [`KineticReaction`](@ref) objects. Every species not carried by a
declared reaction is at equilibrium; every declared reaction is kinetic and must
have a rate law, which `KineticReaction` already enforces.

`coupling` says how the kinetic species are held: `:reactions`, one reactivity
row per reaction in the minimization (Leal's form); `:species`, the kinetic
species eliminated and their amounts set by the extents, the extents found by a
Newton method around an equilibrium of the others. `:reactions`, the default,
solves the composition and the extents together; where that does not certify,
it finds the same step by Newton's method on the extents around the equilibrium
they leave, in the eliminated form when the two pose the same step (every
kinetic species a pure phase, `K` square and invertible: a mineral per
reaction). The simultaneous solve is the faster and crosses the corners of a
rate law that switches mechanism; the Newton method on the extents holds where
rates are far from saturation.

Requires `OptimaSolver`.

# Example

```julia
kss = KineticStepSolver(cs, DiluteSolutionModel(), [kr_calcite])
st1 = kinetic_step(kss, st0, 60.0u"s")
```

See also [`kinetic_step`](@ref), [`integrate`](@ref) for the `ODEProblem` route.
"""
struct KineticStepSolver{D, R <: AbstractVector}
    dual::D
    reactions::R
    K::Matrix{Float64}          # species × reactions, restricted to the partition
    M::Matrix{Float64}          # KᵀK (`:reactions` coupling only)
    idx_kin::Vector{Int}        # the kinetic partition, in system order
    coupling::Symbol            # :reactions or :species
    free::Vector{Int}           # `:species` only: the species left to the minimum
    dual_free::Any              # `:species` only: the solver over those species
end

function KineticStepSolver(
        system::ChemicalSystem, model::AbstractActivityModel, reactions;
        kinetic_species = :auto,
        coupling::Symbol = :reactions,
        kwargs...,
    )
    coupling in (:reactions, :species) || throw(
        ArgumentError(
            "`coupling` must be `:reactions` (the default: one constraint " *
                "per reaction, Leal's form) or `:species` (one per kinetic species, " *
                "the assemblage imposed by the stoichiometry); got :$coupling."
        )
    )
    isempty(reactions) && throw(
        ArgumentError(
            "no kinetic reaction given: with none, the step is a plain " *
                "equilibrium — call `equilibrate` instead."
        )
    )
    ns = length(system.species)
    for kr in reactions
        length(kr.stoich) == ns || throw(
            DimensionMismatch(
                "a reaction's `stoich` has $(length(kr.stoich)) " *
                    "entries for a system of $ns species."
            )
        )
    end
    K = _reactivity_matrix(reactions, system, kinetic_species)
    idx_kin = [i for i in 1:ns if any(!iszero, @view K[i, :])]
    des = DualEquilibriumSolver(system, model; kwargs...)
    eliminable = _eliminable(K, idx_kin, des.idx_pure)

    # `:species` ELIMINATES the pinned species instead of constraining them.
    #
    # Holding them by a linear row leaves their stationarity in the system, and
    # its multiplier must reach the mineral's own chemical potential — of order
    # 10²-10³ in RT units. The Newton then sits at a fixed point its line search
    # cannot leave: measured, an element balance of 6.1e-7 mol whatever `maxit`,
    # `tol` or the number of active-set updates, against 2e-11 for `:reactions`.
    #
    # Their amounts are an explicit affine function of the extents, so they need
    # not be unknowns at all. Removing them from the system and subtracting their
    # element content from the budget leaves a plain, well-conditioned equilibrium
    # over what remains, and the extents are found by a Newton on `nr` unknowns
    # outside it. That is an outer loop, deliberately: it is the price of imposing
    # an assemblage, and `nr` is a handful where the composition is dozens.
    # Built for `:species`, and for `:reactions` wherever the eliminated form
    # poses the same step, which the `:reactions` step falls back on.
    free, dual_free = if coupling === :species || eliminable
        fr = setdiff(1:ns, idx_kin)
        # `_equilibrium_subsystem`, not a bare `ChemicalSystem`. Rebuilding the
        # free side from its species and the parent's primary NAMES dropped
        # every piece of metadata the parent carried: the solid solutions —
        # which is the defect `kinetics_problems.jl` documents at length,
        # end-members silently becoming separate pure phases with no mixing
        # entropy — and the site families, whose loss
        # leaves their members in the sub-system as `AS_SURFACE` species
        # belonging to no family, which `ChemicalSystem` refuses outright. So
        # `coupling = :species` could not be combined with a surface at all.
        #
        # It also filters the primaries to those actually on this side, where
        # `_primary_symbols` handed over the parent's whole list, and it refuses
        # a family split across the partition by name instead of by its
        # downstream orphan symptom.
        (fr, DualEquilibriumSolver(_equilibrium_subsystem(system, fr), model; kwargs...))
    else
        (Int[], nothing)
    end

    return KineticStepSolver(
        des, collect(reactions), K, transpose(K) * K, idx_kin, coupling,
        free, dual_free,
    )
end


# Whether the two couplings pose the same step, so that `:reactions` may fall
# back on the eliminated one. They do when every kinetic species is a pure phase and `K`,
# restricted to them, is square and invertible: `Kᵀ(n − n₀) = Δt·KᵀK·r` then
# gives `n − n₀ = Δt·K·r`, which is what `:species` imposes. The eliminated step
# solves a plain equilibrium inside a Newton on the extents, which certifies
# where the augmented problem does not: measured on a cement paste with four
# AFt/AFm minerals and a C-S-H solution, with or without a diffuse layer.
function _eliminable(K, idx_kin, idx_pure)
    length(idx_kin) == size(K, 2) || return false
    all(in(idx_pure), idx_kin) || return false
    return rank(K[idx_kin, :]) == size(K, 2)
end

"""
    _warm_x0(model, state, n0, ϵ) -> Vector{Float64}

The equilibrated starting guess for an implicit step, **or `n0` when it cannot be
computed**.

# Why the guess is equilibrated at all

The dual solver admits a mixing phase by a tangent-plane test inside its
active-set loop, and from a composition where the phase is entirely absent that
admission fails. Measured on C₃S dissolving into a C-S-H solid solution: a cold
start left one end-member at 2.7e-9 with a stationarity residual of 6.5, while
the same problem started from an equilibrated guess certified at 3.6e-13 with all
four end-members present. The failure belongs to the cold start, not to the
kinetics — a plain equilibrium fails the same way.

`autostart = false`, because the continuation fallback of
[`equilibrate_certified`](@ref) is for a caller who has no starting point. Here
`state` **is** one — the previous instant of the integration — and a handful of
extra solves would be paid at every implicit step. Same principle as Reaktoro's
coupled use, where the solver is reused across instants.

# Why a failure here is not a failure of the step

What this computes is a STARTING POINT, not a result, so strict convergence is
suspended for it — the same distinction the continuation makes. Left set, the
strict flag would turn an uncertified guess into a raise, the `catch` would
swallow it, and the step would fall back on `n0` anyway: a caller asking for
strict results would silently get a *worse* start than a caller who did not.

And whatever else goes wrong, the fallback **is** the cold start. A guess that
cannot be computed must not take the step down with it: the step is judged on its
own result, not on the quality of where it began. Both suspensions are scoped to
this task, so a trajectory integrated beside another cannot take its flags away —
nor leave them behind.
"""
function _warm_x0(model, state, n0, ϵ)
    return _relaxed_convergence() do
        _exploring_starts() do
            try
                eq = first(
                    equilibrate_certified(state; model = model, ϵ = ϵ, autostart = false)
                )
                Float64[ustrip(us"mol", x) for x in eq.n]
            catch
                n0
            end
        end
    end
end

"""
    kinetic_step(kss, state, Δt; t = 0.0, ϵ = 1e-16, parameters = nothing)
        -> ChemicalState

One implicit kinetic step of duration `Δt` from `state`.

`t` is the absolute time handed to rate laws that depend on it. Pass a `Ref` as
`parameters` to receive the reaction extents `Δξ` the step found.

`n_initial`, a state of the same system, gives the amounts the rate laws read as
their `n_initial` argument; by default they read `state`, the start of the step.
A march of steps passes the state its history started from, so that a degree of
reaction `1 − n/n_initial` is counted from there: from `state`, a Parrott–Killoh
law would start every step at `α = 0`.

`warm_start` equilibrates the starting GUESS when the system carries solid
solutions, leaving the component totals untouched. Keep it on: a mixing phase is
admitted by a tangent-plane test, and from a composition where the phase is
absent the step does not reach its answer. Measured on C₃S dissolving into a
C-S-H solution, a cold start did not certify: with OptimaSolver 0.7.4 it left one
end-member at 2.7e-9 with a stationarity residual of 6.5, and with the exact outer
Jacobian since it admits every end-member and stops at a KKT error of 1.75.

The rate laws are evaluated at the **end-of-step** composition, which is what
makes the step implicit and stable on a stiff system. `Δt` is the caller's choice;
[`kinetic_step_adaptive`](@ref) chooses it from a local error estimate instead.

A step whose answer is not certified is not returned: a [`KineticStepError`](@ref)
is thrown instead, naming the condition that failed (the balance, the rate
equations, a reactant exhausted by a rate that does not vanish) and what to change.

!!! warning "One step far beyond the relaxation time can land on the other root"
    The step is a nonlinear problem, and for a rate law that vanishes at
    equilibrium it has a second solution: the composition with the mineral wholly
    dissolved. That one satisfies the element balance and the reactivity row while
    violating `Δξ − Δt·M·r(n) = 0` by the entire extent, so the certificate
    refuses it and the step throws.

    Which root the Newton finds depends on its path, not on the chemistry.
    Measured on calcite under `r = k(1 − Ω)` with `k = 10⁻⁵ mol/s`: with the
    outer Jacobian of OptimaSolver 0.7.4, formed by differences, steps of `10⁴ s`
    and `10⁶ s` converged to the right root on Julia 1.12 and to the other one on
    1.13, and so did `10⁵ s` on 1.13.1. With the exact outer Jacobian all three
    reach the equilibrium root and certify, on 1.13.1. Steps well inside the
    relaxation time certify on either, and the certificate stays the arbiter.

    [`kinetic_step_adaptive`](@ref) is the guarantee whatever the path: it
    refuses an uncertified step and halves until one certifies. With 0.7.4 it
    came back in seven steps on the same case; now it takes the interval in one.
"""
function kinetic_step(
        kss::KineticStepSolver, state::ChemicalState, Δt;
        t = 0.0, ϵ::Float64 = _AMOUNT_FLOOR,
        warm_start::Bool = true,
        pin_minerals = :auto,
        n_initial::Union{Nothing, ChemicalState} = nothing,
        parameters::Union{Nothing, Base.RefValue} = nothing,
        certificate::Union{Nothing, Base.RefValue} = nothing,
    )
    # On plain numbers only, and said at once. On dual numbers the step would
    # differentiate the iterations of its Newton method, not the step, and it
    # failed further down on a plain buffer; the derivatives of a kinetic run
    # are those `integrate` carries along the trajectory.
    D = promote_type(
        _amount_number_type(state), _reference_number_type(n_initial), _captured_number_type(Δt),
        (_captured_number_type(kr.rate_fn) for kr in kss.reactions)...,
    )
    D <: ForwardDiff.Dual && throw(
        ArgumentError(
            "kinetic_step works on plain numbers: the state, the step or a rate law carries " *
                "dual numbers. Differentiate a kinetic run through `integrate`, whose " *
                "trajectory carries the derivatives.",
        ),
    )
    Δt_s = Δt isa Number && !(Δt isa DynamicQuantities.AbstractQuantity) ?
        Float64(Δt) : ustrip(us"s", Δt)
    out, cert, q = kss.coupling === :species ?
        _kinetic_step_eliminated(kss, state, Δt_s; t = t, ϵ = ϵ, n_initial = n_initial) :
        _kinetic_step_reactions(kss, state, Δt_s, pin_minerals; t, ϵ, warm_start, n_initial)
    certificate === nothing || (certificate[] = cert)
    parameters === nothing || (parameters[] = q)
    return out
end

# The `:reactions` step, pinned and free resolved by margin when asked to
# (`pin_minerals = :auto`), and refused when neither certifies.
function _kinetic_step_reactions(kss, state, Δt_s, pin_minerals; t, ϵ, warm_start, n_initial)
    # Whether to hold the kinetic minerals in the active set is decided by the
    # certificate, not by taste, because neither answer works on both cases.
    # Ranked by margin (`_certificate_error`), not by whether one certificate
    # happens to pass: two answers straddling the tolerance would otherwise be
    # decided by rounding, hence by the machine. On the case that exposed it, the
    # free branch has run to equilibrium and violates its own rate equation by
    # the whole extent while the pinned one satisfies it to 1e-13.
    pins = pin_minerals === :auto ? (false, true) : (pin_minerals,)
    best = nothing
    best_err = Inf
    for pin in pins
        r = _kinetic_step_augmented(
            kss, state, Δt_s; t = t, ϵ = ϵ, warm_start = warm_start, pin_minerals = pin,
            n_initial = n_initial,
        )
        err = _certificate_error(r[2])
        (best === nothing || err < best_err) && (best = r; best_err = err)
    end
    best[2].optimal && return best
    # The simultaneous solve did not certify: the same step, its extents found by
    # Newton's method around the equilibrium they leave. On rates far from
    # saturation, Ω of order 10⁶ in a pore solution mixed with an exposure
    # solution, the simultaneous Newton method does not converge and this one
    # does. Where `K` holds a pure phase per reaction the two couplings pose the
    # same step, and the eliminated form, whose kinetic amounts are explicit, is
    # the better conditioned of the two (`_eliminable`).
    if kss.dual_free !== nothing
        out, cert, q = _kinetic_step_eliminated(kss, state, Δt_s; t = t, ϵ = ϵ, n_initial = n_initial)
        # Its extents are those of the kinetic species, `n − n₀ = K q`; the
        # extents of the reactions are `Δξ = Kᵀ(n − n₀) = M q`.
        return out, cert, kss.M * q
    end
    return _kinetic_step_nested(kss, state, Δt_s; t = t, ϵ = ϵ, n_initial = n_initial)
end

function _kinetic_step_augmented(
        kss::KineticStepSolver, state::ChemicalState, Δt_s::Float64;
        t = 0.0, ϵ::Float64 = _AMOUNT_FLOOR, warm_start::Bool = true,
        pin_minerals::Bool = false, n_initial::Union{Nothing, ChemicalState} = nothing,
    )
    des = kss.dual
    p = _build_params(state; ϵ = ϵ)
    n0 = Float64[ustrip(us"mol", x) for x in state.n]

    K, M = kss.K, kss.M
    nr = size(K, 2)
    m = size(des.A, 1)

    # Leal's form: `[A ; Kᵀ] n + [0 ; −I] Δξ = [A n₀ ; Kᵀ n₀]`, one row per
    # reaction. The extents are functionals of the composition, and every
    # species the partition excludes is left to the minimization. (The
    # `:species` coupling, one row per kinetic species, is the elimination
    # returned above.)
    A_aug = vcat(des.A, transpose(K))
    Aq = vcat(zeros(m, nr), -Matrix(1.0I, nr, nr))
    # The element rows hold the budget of the state; the reactivity rows count
    # the extents from the amounts the solve starts from, each at least `ϵ`, the
    # lower bound of every species in the active set. A kinetic mineral absent
    # from the state is then held at that bound by its row instead of being asked
    # for zero, below it; its extent differs from the state's by `ϵ` at most.
    b_aug = vcat(des.A * n0, transpose(K) * max.(n0, ϵ))

    T_K = p.T
    P_Pa = p.P
    idx = des.system.dict_species
    rates = _step_rate_closure(kss, des, T_K, P_Pa, t, _rate_reference(n_initial, n0, des.system))

    # The scale of each extent is the one the explicit rate predicts, the only
    # one available before the step is taken, floored, for a reaction starting at
    # zero rate, at the rounding of the budget: a fraction of what the system
    # holds, not an amount in moles.
    r0 = rates(n0, p)
    S = maximum(abs, b_aug; init = 0.0)
    qscale = max.(abs.(Δt_s .* (M * r0)), eps() * S)

    # `Δξ − Δt·M·r(n) = 0`. The `M` is not decoration: `Kᵀn − Δξ = ξ₀` defines
    # `Δξ = Kᵀ(n − n₀)`, and with `n − n₀ = Kξ` that is `KᵀK ξ`.
    #
    # Each residual is divided by the amount its reaction acts on, `Σᵢ |Kᵢⱼ| nᵢ`,
    # as an element balance is judged relative to what its row holds. The solver
    # judges these residuals, in the iteration and in the certificate, against an
    # absolute tolerance: unscaled, a step whose extent is below that tolerance
    # satisfies the equation where it starts. Measured on a slow site family (10⁻⁶ mol of sites,
    # 2.3 × 10⁻¹⁰ mol/s), the step returned `Δξ = 10⁻¹⁵` instead of `9 × 10⁻¹⁰`,
    # certified with a residual of `9.4 × 10⁻¹⁰`, the whole extent.
    #
    # Not by the extent itself, `|Δt·M·r₀|`, which is the obvious scale and the
    # wrong one, for the reason `kinetic_step_adaptive` gives for its own
    # tolerance: it vanishes with the step and with the rate while the noise of
    # the solve does not. Measured on calcite dissolving into water (the coupling
    # tutorial), the certificate then failed near equilibrium on every step the
    # adaptive march tried above one second: 706 steps reached 842 s, where the
    # march reaches 10⁵ s in 19.
    rscale = _extent_scale(K, n0, qscale)
    cq = (x, q, params) -> (q .- Δt_s .* (M * rates(x, params))) ./ rscale

    # A kinetic mineral is determined by its reactivity row, so it must not be
    # subject to an active-set sign test.
    #
    # This is not defensive: on a mineral that dissolves ENTIRELY at equilibrium
    # the drop rule removes it — its saturation index is far positive, it is on
    # its way out — and the row pinning it then cannot be satisfied. Measured on
    # C₃S dissolving into a C-S-H solution, the step ran to full equilibrium with
    # `C3S = 0` exactly and `Δξ = 0.51` against `Δt·r = 0.02`. On calcite the same
    # bug is invisible, because only a thousandth of it dissolves and the active
    # set keeps it.
    pinned = pin_minerals ? [i for i in kss.idx_kin if i in des.idx_pure] : Int[]

    # Only the ORIGINAL element and charge rows may be judged degenerate. The
    # rows appended here mean something else — a reactivity constraint, or an
    # amount pinned by the extents — and a product that starts absent gives such
    # a row a single positive entry with a zero right-hand side, exactly the
    # shape the degeneracy criterion reads as "this component is absent". It then
    # pins that row's multiplier and declares the species dead, so the product
    # can never form: measured, the stationarity residual sat at 458 and the
    # extents came out 4.3 times short.
    #
    # A surface whose model needs them carries its potentials as parameters
    # beside the extents, with the closure of each as its equation, as the
    # equilibrium solve poses them (`_surface_potential_blocks`): left inside the
    # activity model, the potential couples every member of the surface and the
    # solve did not converge.
    kin = (;
        nq = nr, gq = nothing, hq = nothing, pq = nothing, cq = cq, Aq = Aq,
        q0 = zeros(nr), qscale = qscale, apply = (T, P, q) -> (T, P),
    )
    bl = _with_surface(kin, _surface_potential_blocks(des, state, p, n0), nr)
    prob = _optima_dual_problem(
        A_aug, Float64.(p.ΔₐG⁰overRT), des.lna, _dual_phases(des, n0, p), des.idx_pure, p,
        bl.gq, bl.hq, bl.cq, bl.q0, bl.qscale, bl.Aq, pinned, 1:m, _holdings_rows(des),
    )
    # The GUESS is equilibrated first; `b_aug` stays the true initial budget.
    #
    # This is not a speed-up, it is what makes a solid solution work. The dual
    # solver admits a mixing phase by a tangent-plane test inside its active-set
    # loop, and from a composition where the phase is entirely absent that
    # admission fails: measured on C₃S dissolving into a C-S-H solid solution,
    # a cold start left one end-member at 2.7e-9 with a stationarity residual of
    # 6.5, while the same problem started from an equilibrated guess certified at
    # 3.6e-13 with all four end-members present. The failure belongs to the cold
    # start, not to the kinetics — a plain equilibrium fails the same way.
    x0 = if warm_start && !isempty(des.ss_groups) && _DUAL_AVAILABLE[]
        _warm_x0(des.model, state, n0, ϵ)
    else
        n0
    end
    # The cold start, which `_warm_x0` also falls back to, as the state holds it,
    # as in the equilibrium solve: an exact zero is a pure phase on its bound, and
    # a solute at zero starts at a fixed fraction of the budget (OptimaSolver
    # reads the start against the multipliers, with no amount in moles). Lifted
    # to `ϵ`, an absent product of a reaction read as held where no multiplier of
    # the solution reads it, and started free.
    x0 === n0 && (x0 = max.(n0, 0.0))

    # Seed the pinned species at the amount the EXPLICIT step predicts, not at
    # the floor. A solid product that starts absent would otherwise enter the
    # Newton at 1e-12 while its answer is 1e-3, nine orders away, and the solve
    # stalls near 1e-8 instead of reaching 1e-12 — measured on two C3A pathways.
    # The prediction costs nothing: the rates at the initial composition are
    # already evaluated for `qscale`.
    #
    # The states of a slow site family are kinetic without being pure phases, so
    # they are never pinned; a complex starting bare is seeded the same way.
    seeded = union(pinned, (i for i in kss.idx_kin if i in des.system.idx_surface))
    if !isempty(seeded)
        x0 = copy(x0)
        for i in seeded
            pred = n0[i] + sum(K[i, j] * Δt_s * r0[j] for j in 1:nr)
            pred > x0[i] && (x0[i] = pred)
        end
    end

    res = _optima_dual_solve(prob, b_aug, x0, des.opts)
    res.converged || Threads.atomic_add!(NONCONVERGED, 1)

    # A solute left below the activity floor while its potentials ask for more
    # is given the amount they give it, as the equilibrium solve does
    # (`_complete_floored_solutes`). Without it, measured on two C3A pathways with
    # gypsum, eight such solutes failed the one-sided test of the certificate by
    # 1.1e-4 once the outer Jacobian was exact and the solve took fewer
    # iterations to converge than the solutes took to climb.
    x = res.converged ? _complete_floored_solutes(des, prob, res, _activity_floor(p)) : res.x

    # The certificate must be taken on the AUGMENTED problem, not on the
    # unconstrained equilibrium. A kinetically held mineral is supersaturated by
    # construction — that is what "held back" means — so testing it against the
    # unconstrained KKT conditions reports a stationarity residual of several RT
    # and calls a correct answer wrong. The reactivity rows carry their own
    # multipliers, and those are exactly what makes the mineral's stationarity
    # satisfiable at the amount the constraint demands.
    cert = _normalize_certificate(
        _optima_kkt_certificate(
            prob, x, b_aug, _CERTIFICATE_FLOOR, des.opts.tol, des.opts.si_tol, res.q,
        ),
    )
    out = ChemicalState(des.system, x .* u"mol"; T = temperature(state), P = pressure(state))
    return out, cert, res.q[1:nr]
end

"""
    _kinetic_step_nested(kss, state, Δt; t, ϵ, n_initial, maxit = 50)

The `:reactions` step with its extents found by Newton's method outside the
equilibrium: at given extents `Δξ` the composition is the equilibrium of the
whole system under the element rows and the reactivity rows
`Kᵀn = Kᵀn₀ + Δξ`, and `Δξ` solves `Δξ − Δt·M·r(n(Δξ)) = 0`. The same step as the
simultaneous solve, posed so that a rate far from saturation cannot take the
Newton method of the composition with it. The Jacobian `dn/dΔξ` is that of the
equilibrium, by the implicit-function theorem (`_optima_tangent`).
"""
function _kinetic_step_nested(
        kss::KineticStepSolver, state::ChemicalState, Δt_s::Float64;
        t = 0.0, ϵ::Float64 = _AMOUNT_FLOOR,
        n_initial::Union{Nothing, ChemicalState} = nothing, maxit::Int = 50,
    )
    des = kss.dual
    p = _build_params(state; ϵ = ϵ)
    n0 = Float64[ustrip(us"mol", x) for x in state.n]
    K, M = kss.K, kss.M
    nr = size(K, 2)
    m = size(des.A, 1)
    tol = des.opts.tol
    A_aug = vcat(des.A, transpose(K))
    bA = des.A * n0
    bK = transpose(K) * n0
    # The start of the active set holds every amount at its lower bound `ϵ` or
    # above; its reactivity rows count the extents from there, so that a kinetic
    # mineral absent from the state starts on its bound with its row met. The
    # dual Newton method and the certificate, whose pure phases are bounded at
    # zero, read the budget of the state.
    lift = vcat(zeros(m), transpose(K) * (max.(n0, ϵ) .- n0))
    rates = _step_rate_closure(kss, des, p.T, p.P, t, _rate_reference(n_initial, n0, des.system))
    # The kinetic pure phases are held by their reactivity rows, not by a sign
    # test, in the solve and in its certificate.
    pinned = [i for i in kss.idx_kin if i in des.idx_pure]
    # The potentials of the surfaces that need them are the only parameters.
    bl = _with_surface(_constraint_blocks(FixedTP(), des, state, p, n0), _surface_potential_blocks(des, state, p, n0), 0)
    Aq = vcat(bl.Aq, zeros(nr, bl.nq))
    problem(q0) = _optima_dual_problem(
        A_aug, Float64.(p.ΔₐG⁰overRT), des.lna, _dual_phases(des, n0, p), des.idx_pure, p,
        bl.gq, bl.hq, bl.cq, q0, bl.qscale, Aq, pinned, 1:m, _holdings_rows(des),
    )
    lb = fill(ϵ, length(n0))
    # The equilibrium inside is solved a thousand times tighter than the rate
    # equations are judged: resumed from the last answer, a solve stops as soon
    # as it meets its tolerance, and the noise it leaves, amplified by a rate far
    # from saturation, would otherwise hold the Newton method on the extents at
    # that tolerance. Newton's method converges quadratically, so the margin
    # costs an iteration at most.
    inner_opts = _with_tolerance(des.opts, tol / 1000)
    # The last answer, the start of the next solve; and its certificate.
    last_x = Ref(max.(n0, ϵ))
    last_y = Ref{Union{Nothing, Vector{Float64}}}(nothing)
    last_q = Ref(collect(Float64, bl.q0))
    last_cert = Ref{Any}(nothing)
    # The equilibrium at the budget `b`, on the values: the dual Newton method
    # from the last answer, then from the answer of the active set.
    function equilibrium(b)
        prob = problem(last_q[])
        res = _optima_dual_solve(prob, b, last_x[], inner_opts; y0 = last_y[])
        if !res.converged
            st = _optima_parametric_solve(prob, b .+ lift, max.(n0, ϵ), last_q[], lb, inner_opts)
            prob = problem(st.q)
            res = _optima_dual_solve(prob, b, max.(st.x, 0.0), inner_opts; y0 = st.y)
        end
        x = res.converged ? _complete_floored_solutes(des, prob, res, _activity_floor(p)) : res.x
        last_cert[] = _normalize_certificate(
            _optima_kkt_certificate(prob, x, b, _CERTIFICATE_FLOOR, tol, des.opts.si_tol, res.q),
        )
        if last_cert[].optimal
            last_x[], last_y[], last_q[] = x, collect(Float64, res.y), collect(Float64, res.q)
        end
        return prob, x, res.q
    end
    function compose(q)
        prob, x, qs = equilibrium(vcat(bA, bK .+ ForwardDiff.value.(q)))
        n = eltype(q) <: ForwardDiff.Dual ?
            _optima_tangent(prob, vcat(bA, bK .+ q), x; q = isempty(qs) ? nothing : qs, floor = _CERTIFICATE_FLOOR).x : x
        return q .- Δt_s .* (M * rates(n, p)), n
    end
    admissible(q) = last_cert[] !== nothing && last_cert[].optimal

    r0 = rates(n0, p)
    q = Δt_s .* (M * r0)
    rscale = _rate_scale(K, n0, q, maximum(abs, bA; init = 0.0), tol)
    S = maximum(abs, bA; init = 0.0)
    nt = _extent_newton(q -> (_residual_and_jacobian(compose, q)..., last_cert[]), admissible, pinned, eps() * S, q, rscale, tol, maxit)
    nt.started || throw(_equilibrium_error(Δt_s, last_cert[]))
    q, F = nt.q, nt.F
    x = nt.n
    cert = _with_rate_residual(nt.cert, _rate_error_value(F, rscale), tol)
    cert.optimal || throw(_newton_error(kss, Δt_s, nt, cert, rscale))
    out = ChemicalState(des.system, x .* u"mol"; T = temperature(state), P = pressure(state))
    return out, cert, q
end

"""
    _reactivity_matrix(reactions, system, kinetic_species) -> Matrix

The matrix `K` of the reactivity constraint `Kᵀn − Δξ = ξ₀`: one column per
reaction, carrying that reaction's stoichiometry restricted to the **kinetic
partition** — the species whose amounts the reactions alone decide.

The restriction is the modeling decision of the whole scheme. A rate law says how
fast a mineral reacts; it says nothing about where matter released into solution
ends up, which is the equilibrium problem's business. So an aqueous species must
be left OUT of `K`: it re-speciates, and including it makes `Kᵀn` a functional of
the products, which stops pinning the mineral. Measured on calcite dissolving at
1 µmol/s for 1000 s with the full stoichiometry in `K`, the amount came out
5.7e-5 mol below `n₀ − kΔt`, the constraint having absorbed the carbonate's
speciation into the extent. The answer solved its own equations exactly and meant
something else.

`kinetic_species = :auto` takes the partition to be **every non-aqueous
participant of a declared reaction**. That covers both styles of model:

  - **dissolution into solution** — `Cal → Ca²⁺ + CO₃²⁻`. Only `Cal` is kinetic,
    so `K = [−1]`, `M = KᵀK = 1`, `Δξ = Δt·r`, and the released matter speciates
    freely under element conservation.
  - **solid to solid** — `C₃S + H₂O → C-S-H + CH`, the form
    [Lavergne2018](@citet) use. Every solid product is kinetic too, so the
    assemblage is **imposed by the stoichiometry** rather than found by
    minimizing `G`. That is a different model, deliberately, and it is the one to
    write when the hydrate assemblage is part of what you are prescribing.

Species occupying a **surface site** are excluded from the automatic partition
when they only take part in a reaction another species controls: their family's
site budget is a conservation row of the equilibrium problem, and moving one
member to the kinetic side would move that row with it. When a reaction
controls a member, as [`sorption_rate`](@ref) requires, the adsorption is the
slow step and the whole family is kinetic, as in [`integrate`](@ref).

Pass an explicit list of species (symbols or indices) to override the partition.

Several reactions may share a mineral — the whole point of a reaction-centric
scheme, and what [Lavergne2018](@cite) needs for C₃A reacting through both an
ettringite and a monosulfoaluminate pathway. Their columns must be **linearly
independent**, though: two reactions differing only in products the partition
excludes would give the same row twice, `Kᵀn` could not tell their extents apart,
and the linear system would be singular. That is checked here, with the offending
reactions named.
"""
# A slow family whose budget follows a kinetic host: its members move with the
# host as well as with their reactions, which the reactivity constraints of one
# implicit step, `Kᵀn − Δξ = ξ₀`, cannot express.
function _refuse_slow_host_families(system::ChemicalSystem, kinetic)
    system.site_families === nothing && return nothing
    for (f, grp) in zip(system.site_families, system.site_groups)
        surface_support(f).coupling === SITES_FOLLOW_HOST && any(in(kinetic), grp) && throw(
            ArgumentError(
                "SiteFamily \"$(name(f))\" is slow and its budget follows its host: its " *
                    "states move with the host as well as with their reactions, which the " *
                    "implicit step does not pose. Integrate it with `integrate`.",
            )
        )
    end
    return nothing
end

function _reactivity_matrix(reactions, system::ChemicalSystem, kinetic_species)
    ns = length(system.species)
    syms = symbol.(system.species)
    aqueous = Set(system.idx_aqueous)
    # A species occupying a surface site is **not** kinetic just because it is
    # not aqueous. Its family's budget is fixed and its states redistribute as
    # fast as the aqueous speciation does, so it belongs on the equilibrium side
    # — and putting it on the kinetic one would take its conservation row with
    # it, leaving the equilibrium a surface with no sites.
    #
    # Unless a reaction CONTROLS a member: the adsorption is then the slow step,
    # and the whole family goes to the kinetic side, as `integrate` puts it
    # (`_kinetic_partition`, which also refuses what cannot be slow).
    surface = Set(system.idx_surface)

    kin = if kinetic_species === :auto
        s = Set{Int}()
        for kr in reactions, i in 1:ns
            kr.stoich[i] == 0 && continue
            (i in aqueous || i in surface) || push!(s, i)
        end
        for i in _kinetic_partition(system, reactions)
            i in surface && push!(s, i)
        end
        _refuse_slow_host_families(system, s)
        s
    else
        Set(
            i isa Integer ? Int(i) : begin
                    j = findfirst(==(String(i)), syms)
                    j === nothing && throw(
                        ArgumentError("`$i` is not a species of this system.")
                    )
                    j
                end for i in kinetic_species
        )
    end

    isempty(kin) && throw(
        ArgumentError(
            "no kinetic species: every participant of every declared " *
                "reaction is aqueous, so there is nothing whose amount the " *
                "rate laws decide. Name the controlled species with " *
                "`kinetic_species = [...]`."
        )
    )

    K = zeros(Float64, ns, length(reactions))
    for (j, kr) in enumerate(reactions), i in kin
        K[i, j] = kr.stoich[i]
    end

    for j in axes(K, 2)
        all(iszero, @view K[:, j]) && throw(
            ArgumentError(
                "reaction $j touches no kinetic species, so its extent " *
                    "is not determined by the composition. Its solid " *
                    "participants are missing from the partition."
            )
        )
    end

    nr = size(K, 2)
    if nr > 1 && rank(K) < nr
        # Name the pair that collides, so the message is actionable.
        culprits = String[]
        for a in 1:nr, bb in (a + 1):nr
            rank(K[:, [a, bb]]) < 2 && push!(culprits, "$a and $bb")
        end
        throw(
            ArgumentError(
                "the reactivity constraints are linearly dependent" *
                    (
                    isempty(culprits) ? "" : " (reactions " *
                        join(culprits, ", ") * ")"
                ) *
                    ": `Kᵀn` cannot tell those extents apart and the linear " *
                    "system would be singular. Two reactions may share a mineral " *
                    "— that is the point of a reaction-centric scheme — but they " *
                    "must differ in some species of the kinetic partition. If " *
                    "they differ only in aqueous products, add those products to " *
                    "`kinetic_species`, or model the pathway as one reaction " *
                    "with a combined rate.",
            )
        )
    end
    return K
end

"""
    _step_rate_closure(kss, des, T, P, t, n0) -> (x, params) -> Vector

The reaction rates at composition `x`, in mol/s, one per declared reaction.

Built once per step so the `StateView`s and the initial amounts are not rebuilt
at every residual evaluation — the residual and its Jacobian are evaluated once
per Newton iteration, on dual numbers.
"""
# The amounts the rate laws read as `n_initial`: those of `reference` when the
# caller gives one, those the step starts from otherwise.
_rate_reference(::Nothing, n0, system) = n0
function _rate_reference(reference::ChemicalState, n0, system)
    _check_same_species(reference, system, "n_initial")
    return Float64[ustrip(us"mol", x) for x in reference.n]
end

_reference_number_type(::Nothing) = Float64
_reference_number_type(reference::ChemicalState) = _amount_number_type(reference)

function _step_rate_closure(kss::KineticStepSolver, des, T, P, t, n0)
    index = Dict(symbol(s) => i for (i, s) in enumerate(des.system.species))
    n0_sv = StateView(n0, index)
    return function (x, params)
        xv = collect(x)
        n_sv = StateView(xv, index)
        lna_sv = StateView(collect(des.lna(xv, params)), index)
        # In the number type of the composition: the Jacobian of the step is
        # taken by forward mode through these rates.
        return _promoted([kr.rate_fn(T, P, t, n_sv, lna_sv, n0_sv) for kr in kss.reactions])
    end
end


# ── adaptive stepping ────────────────────────────────────────────────────────
#
# The step of `kinetic_step` is the caller's. The scheme is backward Euler, so
# the error is first order in `Δt`, and a step ten times too large returns an
# answer ten times less accurate with nothing to say so; this estimates it.
#
# The estimate is Richardson's, on the extents: one step of `Δt` against two of
# `Δt/2`. For a first-order method the difference IS the error of the coarse step
# to leading order, and `2·(ξ_fine − ξ_coarse)` estimates the error of the fine
# one. Three implicit solves per accepted step, which is the price.

"""
    kinetic_step_adaptive(kss, state, Δt; reltol, abstol, ...) -> (state, Δt_used, Δt_next)

One kinetic step with its length chosen from a local error estimate.

Takes one step of `Δt` and two of `Δt/2`, compares the extents, and accepts the
finer pair when the estimated error is within tolerance — halving and retrying
otherwise. Returns the accepted state, the step actually taken, and a suggestion
for the next one.

The estimate is Richardson's for a first-order method: `ξ_fine − ξ_coarse` is the
error of the coarse step to leading order. The error controlled is therefore on
the **reaction extents**, which is what a kinetic step advances; the composition
follows from them through an exact equilibrium solve, so it carries no separate
error of its own.

`Δt` is the step to attempt, `max_halvings` bounds the retries, and the tolerance
is `abstol + reltol·sⱼ` per reaction, where `sⱼ` is the **amount the extent acts
on** — `Σᵢ |Kᵢⱼ| nᵢ`, the kinetic minerals of reaction `j` at the start of the
step. `n_initial` is that of [`kinetic_step`](@ref), the same for the three steps:
`state` by default, not the middle of the step for the second half.

Scaling by the extent instead, `reltol·|Δξⱼ|`, does not work, and the reason is
worth stating because it is the obvious thing to write. `Δξ ∝ Δt`, so that
tolerance vanishes with the step while the equilibrium solve's own noise does
not: the measured error then behaves as `noise / (reltol·Δt)` and **grows** as the
step shrinks. Measured, the controller halved to the floor without advancing and
the final answer got worse as the tolerance was tightened — 2.5e-5 at `1e-3`
against 9.2e-5 at `1e-5`. An ODE integrator's `reltol` multiplies the solution,
not the increment, for exactly this reason.

# Example

```julia
st, dt_used, dt_next = kinetic_step_adaptive(kss, st0, 3600.0u"s"; reltol = 1e-4)
```

See also [`kinetic_step`](@ref) for a step of a length you choose, and
[`integrate`](@ref) for the `ODEProblem` route, which brings SciML's own
step-size control and its choice of stiff solvers.
"""
function kinetic_step_adaptive(
        kss::KineticStepSolver, state::ChemicalState, Δt;
        t = 0.0, reltol::Float64 = 1.0e-4, abstol::Float64 = 1.0e-14,
        max_halvings::Int = 12, ϵ::Float64 = _AMOUNT_FLOOR,
        warm_start::Bool = true, pin_minerals = :auto,
        n_initial::Union{Nothing, ChemicalState} = nothing,
        parameters::Union{Nothing, Base.RefValue} = nothing,
        error_estimate::Union{Nothing, Base.RefValue} = nothing,
    )
    Δt_s = Δt isa DynamicQuantities.AbstractQuantity ? ustrip(us"s", Δt) : Float64(Δt)
    Δt_s > 0 || throw(ArgumentError("`Δt` must be positive, got $Δt_s s."))
    # One reference for the three steps: the second half-step starts from the
    # middle, and taking its `n_initial` from there would compare, through
    # Richardson's difference, two steps of different rate laws.
    reference = something(n_initial, state)

    # `:auto` is resolved ONCE, here, and the same choice is used for the coarse
    # step and for both half-steps. Letting each of the three resolve it on its
    # own compares two different formulations through Richardson's difference,
    # and the error estimate is then meaningless: measured, the march accepted a
    # step that dissolved the whole mineral.
    # Resolved by MARGIN, not by whether one branch's certificate happens to
    # pass — see `_certificate_error`. A threshold decision between two answers
    # that straddle the tolerance is settled by rounding, and the whole march
    # then depends on the machine.
    pin_fixed = pin_minerals
    if pin_minerals === :auto && kss.coupling === :reactions
        errs = [
            _certificate_error(
                _kinetic_step_augmented(
                    kss, state, Δt_s; t = t, ϵ = ϵ, warm_start = warm_start,
                    pin_minerals = pin, n_initial = reference,
                )[2]
            ) for pin in (false, true)
        ]
        # Ties, and the case where neither produced a certificate, go to the
        # pinned formulation: it is the one that enforces the reactivity row.
        pin_fixed = !(errs[1] < errs[2])
    end

    common = (;
        ϵ = ϵ, warm_start = warm_start, pin_minerals = pin_fixed, n_initial = reference,
    )

    # The scale the tolerance is relative to: the amount each reaction acts on.
    n_now = Float64[ustrip(us"mol", x) for x in state.n]
    scales = [
        max(sum(abs(kss.K[i, j]) * n_now[i] for i in axes(kss.K, 1)), 1.0e-30)
            for j in axes(kss.K, 2)
    ]

    h = Δt_s
    last_error = nothing
    for _ in 0:max_halvings
        q_c = Ref(Float64[])
        q_1 = Ref(Float64[])
        c_1 = Ref{Any}(nothing)
        q_2 = Ref(Float64[])
        c_2 = Ref{Any}(nothing)
        # A step that does not certify is refused like one out of tolerance.
        fine = try
            kinetic_step(kss, state, h; t = t, parameters = q_c, common...)
            mid = kinetic_step(
                kss, state, h / 2; t = t, parameters = q_1, certificate = c_1, common...,
            )
            kinetic_step(
                kss, mid, h / 2; t = t + h / 2, parameters = q_2, certificate = c_2,
                common...,
            )
        catch err
            err isa KineticStepError || rethrow()
            last_error = err
            h /= 2
            continue
        end
        last_error = nothing
        ξ_fine = q_1[] .+ q_2[]

        err = maximum(
            abs(ξ_fine[j] - q_c[][j]) / (abstol + reltol * scales[j])
                for j in eachindex(ξ_fine)
        )
        error_estimate === nothing || (error_estimate[] = err)

        # A step is accepted on the estimate AND on the proof, never on the
        # estimate alone.
        #
        # Richardson's difference measures the disagreement between two
        # resolutions, so it is blind to an error the two SHARE. Measured on
        # calcite over `10⁵ s`: the coarse step and both half-steps each dissolve
        # the entire mineral, their extents agree to `5×10⁻¹¹`, and the estimator
        # reports `9×10⁻⁶` — a step it should have refused, accepted as excellent.
        # The certificate is not fooled, because a composition that dissolved
        # everything violates `Δξ − Δt·M·r(n) = 0` by the whole extent.
        proved = (c_1[] === nothing || c_1[].optimal) &&
            (c_2[] === nothing || c_2[].optimal)

        if err <= 1 && proved
            parameters === nothing || (parameters[] = ξ_fine)
            # First order, so the step scales as the tolerance ratio itself,
            # capped at doubling — the usual safety factor of 0.9.
            Δt_next = h * min(2.0, max(0.2, 0.9 / max(err, 1.0e-10)))
            return (fine, h, Δt_next)
        end
        h /= 2
    end

    last_error === nothing || throw(
        KineticStepError(
            last_error.cause, Δt_s / 2^max_halvings, last_error.certificate,
            "no step was certified after $max_halvings halvings, down to " *
                "Δt = $(Δt_s / 2^max_halvings) s. The last refusal: " * last_error.message,
        )
    )
    throw(
        KineticStepError(
            :tolerance, Δt_s / 2^max_halvings, nothing,
            "no step of this march was within tolerance after $max_halvings halvings, " *
                "down to Δt = $(Δt_s / 2^max_halvings) s: the tolerance is below what the " *
                "equilibrium solve itself delivers (its own tolerance is $(kss.dual.opts.tol)).",
        )
    )
end


"""
    _kinetic_step_eliminated(kss, state, Δt; ...)

The `:species` step: the pinned species are computed from the extents rather than
solved for, and the extents come from a Newton on `nr` unknowns around a plain
equilibrium over the species that remain.

`F(Δξ) = Δξ − Δt·r(n(Δξ))`, where `n(Δξ)` puts the imposed amounts where the
stoichiometry says and lets the minimization place the rest under the budget those
amounts leave. Backward Euler, as `kinetic_step` is, since `r` is read at the
composition the step produces.
"""
function _kinetic_step_eliminated(
        kss::KineticStepSolver, state::ChemicalState, Δt_s::Float64;
        t = 0.0, ϵ::Float64 = _AMOUNT_FLOOR,
        n_initial::Union{Nothing, ChemicalState} = nothing,
        maxit::Int = 50,
    )
    des = kss.dual
    system = des.system
    p = _build_params(state; ϵ = ϵ)
    n0 = Float64[ustrip(us"mol", x) for x in state.n]
    A = des.A
    b0 = A * n0
    K, free, pinned = kss.K, kss.free, kss.idx_kin
    nr = size(K, 2)
    tol = des.opts.tol

    index = Dict(symbol(sp) => i for (i, sp) in enumerate(system.species))
    n0_sv = StateView(n0, index)
    nref_sv = StateView(_rate_reference(n_initial, n0, system), index)
    # The certificate of the last composition, written by `compose` below: a
    # reference rather than a binding the closure reassigns, which Julia boxes.
    last_cert = Ref{Any}(nothing)

    # The composition at a given set of extents, and the residual it leaves, in
    # the number type of the extents: evaluated on dual numbers, the certified
    # equilibrium of the free species carries its derivative with respect to its
    # budget (`solve_certified`), and the residual its Jacobian. An extent that
    # takes a kinetic species below zero is outside the domain of the step, which
    # refuses it (`exhausted`); the amount is not clamped, so that the Jacobian
    # stays exact on the boundary, where an absent mineral starts to precipitate.
    kinetic_amounts(q) = n0[pinned] .+ K[pinned, :] * q
    function compose(q)
        nf = Vector{eltype(q)}(n0)
        nf[pinned] .= kinetic_amounts(q)
        b_free = b0 .- A[:, pinned] * nf[pinned]
        st_free = ChemicalState(
            kss.dual_free.system, n0[free] .* u"mol";
            T = temperature(state), P = pressure(state),
        )
        # Silent on a solve that does not certify: its certificate is the verdict,
        # read by the admissibility test and reported by `KineticStepError`.
        eq, cert = _exploring_starts() do
            solve_certified(kss.dual_free, (st_free,); b = b_free, ϵ = ϵ, report = false)
        end
        last_cert[] = cert
        nf[free] .= ustrip.(us"mol", eq.n)
        lna_sv = StateView(collect(des.lna(nf, p)), index)
        n_sv = StateView(nf, index)
        rates = _promoted([kr.rate_fn(p.T, p.P, t, n_sv, lna_sv, nref_sv) for kr in kss.reactions])
        return (q .- Δt_s .* rates), nf
    end
    # The kinetic species the extents take below zero, beyond the rounding of
    # the budget.
    S_b = maximum(abs, b0; init = 0.0)
    exhausted(q) = pinned[kinetic_amounts(q) .< -eps() * S_b]
    rj = q -> (_residual_and_jacobian(compose, q)..., last_cert[])
    # A trial is admissible when the equilibrium of the free species certifies
    # and no kinetic species is taken below zero.
    admissible(q) = last_cert[] !== nothing && last_cert[].optimal && isempty(exhausted(q))

    # Newton on `nr` unknowns, from the explicit prediction. Every rate law reads
    # the same activities at n0: they are evaluated once, not once per reaction.
    lna0_sv = StateView(collect(des.lna(n0, p)), index)
    q = Float64[Δt_s * kss.reactions[j].rate_fn(p.T, p.P, t, n0_sv, lna0_sv, nref_sv) for j in 1:nr]
    # The rate equation of each reaction is judged relative to the amount it acts
    # on, or to its predicted extent where that is larger (`_extent_scale`).
    rscale = _rate_scale(K, n0, q, maximum(abs, b0; init = 0.0), tol)
    # Newton's method on the extents, brought back into the domain `n ≥ 0` of
    # the kinetic species as the bounds of the active set are: each step
    # shortened to the first species inside the domain it would take out, the
    # species already on their bound that it pushes out projected back onto it
    # (the least-norm correction of the extents).
    nt = _extent_newton(rj, admissible, pinned, eps() * S_b, q, rscale, tol, maxit)
    q, F = nt.q, nt.F
    if !nt.started
        exh = exhausted(q)
        throw(
            isempty(exh) ? _equilibrium_error(Δt_s, last_cert[]) :
                _exhausted_error(kss, Δt_s, first(exh), q, F, last_cert[])
        )
    end
    n_full = nt.n
    cert = _with_rate_residual(nt.cert, _rate_error_value(F, rscale), tol)
    cert.optimal || throw(_newton_error(kss, Δt_s, nt, cert, rscale))

    out = ChemicalState(system, n_full .* u"mol"; T = temperature(state), P = pressure(state))
    return out, cert, q
end


"""
    _normalize_certificate(c) -> NamedTuple

One shape for the certificate whichever route produced it.

`kinetic_step` fills its `certificate` reference from `OptimaSolver`'s
`kkt_certificate` on the `:reactions` route and from
[`optimality_certificate`](@ref) on the `:species` one, and those name the
element-balance field differently — `feasibility` against `balance`. Returning
two shapes from one keyword is a trap for the caller, so both are mapped onto the
names `optimality_certificate` uses.
"""
function _normalize_certificate(c)
    hasproperty(c, :balance) && return c
    return (;
        stationarity = c.stationarity, balance = c.feasibility_abs,
        balance_relative = c.feasibility_rel,
        stationarity_abs = hasproperty(c, :stationarity_abs) ?
            c.stationarity_abs : c.stationarity,
        worst_supersaturation = c.worst_violation,
        param_residual = hasproperty(c, :param_residual) ? c.param_residual : 0.0,
        n_interior = c.n_interior, n_absent_component = c.n_forced_zero,
        optimal = c.optimal,
    )
end


"""
    _certificate_error(c) -> Float64

One number ranking how far a step is from satisfying its own conditions:
stationarity, element balance, the worst supersaturation among absent phases, and
the residual of the rate equation. `Inf` when no certificate was produced.

Used to choose between the pinned and free formulations of a step. A ranking is
needed rather than the boolean `optimal`, because two answers may straddle the
tolerance and then the choice is decided by rounding — the same code choosing
differently on another machine. Ranked, the gap is fifteen orders of magnitude.
"""
function _certificate_error(c)
    c === nothing && return Inf
    e = max(
        Float64(c.stationarity),
        Float64(_judged_balance(c)),
        max(Float64(c.worst_supersaturation), 0.0),
    )
    hasproperty(c, :param_residual) && (e = max(e, Float64(c.param_residual)))
    return isfinite(e) ? e : Inf
end

# The amount each reaction acts on, `Σᵢ |Kᵢⱼ| n₀ᵢ`, or its predicted extent where
# that is larger: the scale its rate equation is judged against, as an element
# balance is judged against what its row holds. `pred` is already floored at the
# rounding of the budget, so that a reaction starting at zero rate on absent
# minerals keeps a positive scale.
_extent_scale(K, n0, pred) = max.(transpose(abs.(K)) * n0, abs.(pred))

"""
    KineticStepError(cause, Δt, certificate, message)

Thrown by [`kinetic_step`](@ref) when no step of `Δt` seconds is certified, in
place of a composition that would not satisfy its own equations. `cause` says
which condition failed:

  - `:equilibrium`: the equilibrium of the species the reactions do not hold does
    not certify at the budget of the step;
  - `:exhausted`: a reaction consumes more of a kinetic species than the state
    holds within `Δt`, because its rate does not vanish as that species runs out;
  - `:rate`: the rate equations `Δξ = Δt·r(n)` are not met, which happens when a
    rate law is not differentiable at the composition reached;
  - `:tolerance` ([`kinetic_step_adaptive`](@ref) only): every step was certified
    and none within the tolerance of the march.

`message` names the reaction or the species involved and what to change.
[`kinetic_step_adaptive`](@ref) catches it and halves the step.
"""
struct KineticStepError <: Exception
    cause::Symbol
    Δt::Float64
    certificate::Union{Nothing, NamedTuple}
    message::String
end

Base.showerror(io::IO, e::KineticStepError) =
    print(io, "KineticStepError (", e.cause, ", Δt = ", e.Δt, " s): ", e.message)

# The halvings of an extent, or of a Newton step on the extents, before the
# eliminated step gives up: 2⁻³⁰ of the prediction is below the rounding of any
# extent that matters.
const _EXTENT_HALVINGS = 30

_reaction_name(kr) = String(kr.reaction.symbol)

# A kinetic species `i` that a reaction consumes while the state holds none of it.
function _exhausted_error(kss, Δt_s, i, q, F, cert)
    K = kss.K
    js = [j for j in axes(K, 2) if !iszero(K[i, j])]
    names = join((_reaction_name(kss.reactions[j]) for j in js), ", ")
    r = maximum(abs((q[j] - F[j]) / Δt_s) for j in js)
    return KineticStepError(
        :exhausted, Δt_s, cert,
        "$(names) consumes $(symbol(kss.dual.system.species[i])) faster than the state " *
            "holds it: its rate stays at $(r) mol/s while the amount goes to zero, so no " *
            "step keeps it non-negative. A rate law must vanish when its reactant is " *
            "absent; a floor on the reactive area (nucleation) belongs to precipitation only.",
    )
end

# The rate equations not met, named by the reaction with the largest residual.
function _rate_error(kss, Δt_s, cert, F, rscale)
    j = argmax(abs.(F) ./ rscale)
    return KineticStepError(
        :rate, Δt_s, cert,
        "the rate equation of $(_reaction_name(kss.reactions[j])) is not met: its residual is " *
            "$(abs(F[j]) / rscale[j]) of the amount it acts on. A rate law that is not " *
            "differentiable at the composition reached (a `max`, `min` or threshold) stops " *
            "Newton's method: smooth it, or take a shorter Δt.",
    )
end

# The parameters of the step with the potentials of its surfaces beside the
# extents, where a surface needs them.
_with_surface(kin, ::Nothing, nr) = kin
_with_surface(kin, surf, nr) = _compose_blocks(kin, merge(surf, (Aq = vcat(surf.Aq, zeros(nr, surf.nq)),)))

# The extents `q` moved so that the kinetic amounts `n₀ + K q` taken below zero
# return to it: the least-norm correction `δ` with `(K δ)ᵢ = −(n₀ + K q)ᵢ` on those
# rows, the others untouched where `K` lets them be.
function _project_extents(q, n0k, Kk)
    nk = n0k .+ Kk * q
    out = findall(<(0), nk)
    isempty(out) && return q
    return q .- pinv(Kk[out, :]) * nk[out]
end

# The residual of the rate equations and its Jacobian with respect to the
# extents `q`, from one evaluation of `compose` on dual numbers under a tag of its
# own. `compose(q)` returns the residual and the composition, in the number type
# of `q`.
function _residual_and_jacobian(compose, q)
    nr = length(q)
    tag = ForwardDiff.Tag(compose, Float64)
    qd = [
        ForwardDiff.Dual{typeof(tag)}(q[j], ForwardDiff.Partials(ntuple(k -> Float64(k == j), nr)))
            for j in 1:nr
    ]
    Fd, nd = compose(qd)
    return ForwardDiff.value.(Fd), [ForwardDiff.partials(Fd[r], j) for r in 1:nr, j in 1:nr],
        ForwardDiff.value.(nd), [ForwardDiff.partials(nd[i], j) for i in eachindex(nd), j in 1:nr]
end

# The scale each rate equation of the Newton method on the extents is judged
# against: the amount its reaction acts on, `Σᵢ |Kᵢⱼ| n₀ᵢ`, or its predicted extent
# where that is larger, and never below `eps·max|b| / tol`. A residual below the
# rounding of the budget, `eps·max|b|`, cannot be resolved, and a reaction on a
# mineral the state does not hold has no amount of its own to be judged against:
# its residual is met at the rounding of the budget.
_rate_scale(K, n0, pred, S, tol) = max.(transpose(abs.(K)) * n0, abs.(pred), eps() * S / tol)

# The largest residual of the rate equations, each relative to its scale.
_rate_error_value(F, rscale) = maximum(abs.(F) ./ rscale; init = 0.0)

# A certificate of the equilibrium at the extents found, with the residual of
# the rate equations: optimal when both are.
_with_rate_residual(cert, err, tol) =
    merge(cert, (param_residual = max(err, get(cert, :param_residual, 0.0)), optimal = cert.optimal && err <= tol))

"""
    _extent_newton(rj, admissible, pinned, floor, q, rscale, tol, maxit)

Newton's method on the extents `q` of a step, shared by its two couplings. `rj(q)`
returns the residual of the rate equations, its Jacobian, the composition, its
Jacobian with respect to `q` and the certificate of the equilibrium at `q`;
`admissible(q)`, read after `rj(q)`, whether that equilibrium certified. The amounts of the kinetic species `pinned` must
stay non-negative, as the bounds of an active set are kept: each Newton step is
shortened to the first of them inside the domain it would take out, and those
already on their bound (at most `floor`) that it pushes out are projected back
onto it, by the least-norm correction of the extents on the linear model of
their amounts. The step is then halved until its trial is admissible and lowers
the residual, each rate equation judged relative to `rscale`. The start, the
explicit prediction, is halved toward zero until it is admissible
(`started = false` when it never is). Returns the extents, the residual, the
composition and the certificate of the iterate it ends on, which need no solve
again, the last Newton direction, the kinetic species that direction pushes out
of the domain from their bound, and `started`.
"""
function _extent_newton(rj, admissible, pinned, floor, q, rscale, tol, maxit)
    F, J, n, G, cert = rj(q)
    halvings = 0
    while !admissible(q) && halvings < _EXTENT_HALVINGS
        q = q ./ 2
        F, J, n, G, cert = rj(q)
        halvings += 1
    end
    δ = zero(q)
    admissible(q) || return (; q, F, n, cert, δ, held = Int[], started = false, blocked = true)
    # Whether the last search found no admissible trial: the equilibrium does not
    # certify at the extents the rate equations ask for.
    blocked = false
    for _ in 1:maxit
        _rate_error_value(F, rscale) <= tol && break
        δ = try
            J \ (-F)
        catch err
            err isa LinearAlgebra.SingularException || rethrow()
            break
        end
        all(isfinite, δ) || break
        nk, Gk = n[pinned], G[pinned, :]
        dk = Gk * δ
        β = minimum((-nk[k] / dk[k] for k in eachindex(nk) if nk[k] > floor && nk[k] + dk[k] < 0); init = 1.0)
        λ = 1.0
        accepted = false
        blocked = true
        for _ in 1:_EXTENT_HALVINGS
            qn = _project_extents(q .+ (λ * β) .* δ, nk .- Gk * q, Gk)
            Fn, Jn, nn, Gn, cn = rj(qn)
            ok = admissible(qn)
            ok && (blocked = false)
            if ok && _rate_error_value(Fn, rscale) < _rate_error_value(F, rscale)
                q, F, J, n, G, cert = qn, Fn, Jn, nn, Gn, cn
                accepted = true
                break
            end
            λ /= 2
        end
        accepted || break
    end
    # The kinetic species on their bound that the last direction pushes out,
    # where the rate equation of a reaction acting on them is not met: their rate
    # does not vanish there.
    unmet = abs.(F) ./ rscale .> tol
    Gk = G[pinned, :]
    held = [
        i for (k, i) in enumerate(pinned)
            if n[i] <= floor && (Gk * δ)[k] < 0 && any(unmet[j] && !iszero(Gk[k, j]) for j in eachindex(unmet))
    ]
    return (; q, F, n, cert, δ, held, started = true, blocked)
end

# Why the Newton method on the extents stopped short: a kinetic species its rate
# keeps consuming at zero, an equilibrium that does not certify at the extents
# asked for, or a rate equation it cannot solve.
function _newton_error(kss, Δt_s, nt, cert, rscale)
    isempty(nt.held) || return _exhausted_error(kss, Δt_s, first(nt.held), nt.q, nt.F, cert)
    nt.blocked && return KineticStepError(
        :equilibrium, Δt_s, cert,
        "the equilibrium does not certify at the extents the rate equations ask for, " *
            "so Newton's method on the extents cannot reach them. Take a shorter Δt, or " *
            "check the starting state with `equilibrate_certified`.",
    )
    return _rate_error(kss, Δt_s, cert, nt.F, rscale)
end

# The equilibrium at the start of the step does not certify.
_equilibrium_error(Δt_s, cert) = KineticStepError(
    :equilibrium, Δt_s, cert,
    "the equilibrium of the species the reactions do not hold does not certify at the " *
        "budget of the step, even at a vanishing extent. Check the starting state with " *
        "`equilibrate_certified`.",
)

# The options of a solve with the tolerance `tol`, the others unchanged.
_with_tolerance(o::NamedTuple, tol) = merge(o, (tol = tol,))
