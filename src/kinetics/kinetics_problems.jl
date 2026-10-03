# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using DynamicQuantities
using LinearAlgebra
using SciMLBase

# ── KineticsProblem ───────────────────────────────────────────────────────────

"""
    struct KineticsProblem{CS, CAL, ES, AM}

Encapsulates a kinetics simulation following Leal et al. (2017).

The ODE state vector `u` is structured as:
  - Without re-speciation: `u = [nₖ₁, …, nₖ_K, ξ₁, …, ξ_M, [T | Q]]`
  - With re-speciation:    `u = [bₑ₁, …, bₑ_C, nₖ₁, …, nₖ_K, ξ₁, …, ξ_M, [T | Q]]`

where `bₑ` are the element amounts in the equilibrium partition, `nₖ` the moles
of kinetic species and `ξ` the extents of the kinetic reactions. The trailing
slot is present only with a calorimeter: the temperature for
[`SemiAdiabaticCalorimeter`](@ref), the accumulated heat for
[`IsothermalCalorimeter`](@ref).

# Fields

  - `system`: [`ChemicalSystem`](@ref).
  - `kinetic_reactions`: vector of [`KineticReaction`](@ref) objects.
  - `initial_state`: [`ChemicalState`](@ref) providing initial moles, T, P.
  - `tspan`: `(t_start, t_end)` time interval [s].
  - `calorimeter`: `nothing`, [`IsothermalCalorimeter`](@ref),
    or [`SemiAdiabaticCalorimeter`](@ref).
  - `activity_model`: [`AbstractActivityModel`](@ref) for log-activities.
  - `equilibrium_solver`: solver for re-speciation, or `nothing`.
  - `idx_kinetic`: indices of kinetic species in `system.species`.
  - `idx_equilibrium`: indices of equilibrium species.
  - `ν`: stoichiometric matrix (M × N) = `SM.N'` restricted to kinetic reactions.
  - `νe`, `νk`: partitions of `ν` for equilibrium / kinetic species.
  - `Ae`: formula matrix restricted to equilibrium species (C × Nₑ).

See also: [`integrate`](@ref), [`KineticsSolver`](@ref).
"""
struct KineticsProblem{
        CS <: ChemicalSystem,
        KR <: AbstractVector,
        CAL,
        ES,
        AM <: AbstractActivityModel,
    }
    system::CS
    kinetic_reactions::KR
    initial_state::ChemicalState
    tspan::Tuple{Float64, Float64}
    calorimeter::CAL
    activity_model::AM
    equilibrium_solver::ES
    # ── Pre-computed partitions (Leal 2017, Eq. 53) ──
    idx_kinetic::Vector{Int}
    idx_equilibrium::Vector{Int}
    ν::Matrix{Float64}          # (M × N) stoichiometric matrix of kinetic reactions
    νe::Matrix{Float64}         # (M × Nₑ) equilibrium columns
    νk::Matrix{Float64}         # (M × K)  kinetic columns
    Ae::Matrix{Float64}         # (C × Nₑ) formula matrix, equilibrium partition
end

"""
    KineticsProblem(cs, kinetic_reactions, initial_state, tspan; ...) -> KineticsProblem

Construct a [`KineticsProblem`](@ref) from an explicit list of reactions.

Each element of `kinetic_reactions` must be either a [`KineticReaction`](@ref) or
a [`Reaction`](@ref) with a `:rate` entry in its properties. Reaction objects
are automatically wrapped via
`KineticReaction(cs, rxn)`.

    KineticsProblem(cs, initial_state, tspan; ...) -> KineticsProblem

Construct from a [`ChemicalSystem`](@ref) that has `kinetic_species` declared
(reactions and rates auto-generated via the `kinetic_species` keyword).

# Arguments

  - `cs`: [`ChemicalSystem`](@ref).
  - `kinetic_reactions`: `AbstractVector` of [`KineticReaction`](@ref) or
    [`Reaction`](@ref) objects carrying a `:rate` property.
  - `initial_state`: [`ChemicalState`](@ref) providing initial moles, T, P.
  - `tspan`: `(t0, tf)` time interval. Plain `Real` → [s]; `Quantity` → converted.
  - `calorimeter`: `nothing` (no thermal coupling),
    [`IsothermalCalorimeter`](@ref), or [`SemiAdiabaticCalorimeter`](@ref).
  - `activity_model`: activity model for log-activity computation (default: dilute).
  - `equilibrium_solver`: `nothing` (no re-speciation) or an [`EquilibriumSolver`](@ref).

# Examples

```julia
# From explicit reactions
rxn = Reaction(OrderedDict(sp("C3S") => 1.0, sp("H2O@") => 3.33),
               OrderedDict(sp("Jennite") => 0.167, sp("Portlandite") => 1.5))
rxn[:rate] = parrott_killoh(PK_PARAMS_C3S, "C3S"; α_max)
kp = KineticsProblem(cs, [rxn], state0, (0.0, 7 * 86400.0))

# From kinetic_species in ChemicalSystem
cs = ChemicalSystem(species, primaries;
    kinetic_species = Dict("C3S" => pk_C3S, "C2S" => pk_C2S))
kp = KineticsProblem(cs, state0, (0.0, 7 * 86400.0))
```
"""
function _build_kinetics_problem(
        system::ChemicalSystem,
        kin_rxns::AbstractVector{<:KineticReaction},
        initial_state::ChemicalState,
        tspan::Tuple;
        calorimeter = nothing,
        activity_model::AbstractActivityModel = DiluteSolutionModel(),
        equilibrium_solver = nothing,
    )
    n_sp = length(system.species)
    idx_kin = unique!(Int[kr.idx_mineral for kr in kin_rxns])
    idx_eq = setdiff(1:n_sp, idx_kin)
    n_rxn = length(kin_rxns)

    # Stoichiometric matrix ν (M × N) — Leal Eq. 44
    ν = zeros(Float64, n_rxn, n_sp)
    for (i, kr) in enumerate(kin_rxns)
        ν[i, :] .= kr.stoich
    end

    # Partition (Leal Eq. 53): ν = [νₑ  νₖ]
    νe = ν[:, idx_eq]
    νk = ν[:, idx_kin]

    # Formula matrix for equilibrium partition: Aₑ = CSM.A[:, idx_eq]
    # Conservation matrix of the equilibrium partition, taken from the partition
    # SUB-SYSTEM and in the basis the equilibrium solve uses (`SM.A`, with
    # respect to the primaries — not the canonical element matrix `CSM.A`).
    # Conserving the primaries is equivalent to conserving the elements, but the
    # row counts differ, and so do the parent's and the sub-system's: `bₑ` must
    # be built on exactly the matrix the solve is posed on, or every step fails
    # on a dimension mismatch.
    Ae = Float64.(_constraint_matrix(_equilibrium_subsystem(system, idx_eq)))

    return KineticsProblem{
        typeof(system), typeof(kin_rxns), typeof(calorimeter),
        typeof(equilibrium_solver), typeof(activity_model),
    }(
        system,
        kin_rxns,
        initial_state,
        (Float64(safe_ustrip(us"s", tspan[1])), Float64(safe_ustrip(us"s", tspan[2]))),
        calorimeter,
        activity_model,
        equilibrium_solver,
        collect(Int, idx_kin),
        collect(Int, idx_eq),
        ν, νe, νk, Ae,
    )
end

# 4-argument form: explicit reaction list (Reaction or KineticReaction)
function KineticsProblem(
        system::ChemicalSystem,
        kinetic_reactions::AbstractVector,
        initial_state::ChemicalState,
        tspan::Tuple;
        calorimeter = nothing,
        activity_model::AbstractActivityModel = DiluteSolutionModel(),
        equilibrium_solver = nothing,
    )
    kin_rxns = [
        r isa KineticReaction ? r : KineticReaction(system, r)
            for r in kinetic_reactions
    ]
    return _build_kinetics_problem(
        system, kin_rxns, initial_state, tspan;
        calorimeter, activity_model, equilibrium_solver,
    )
end

# 3-argument form: reactions from ChemicalSystem.reactions (kinetic_species API)
function KineticsProblem(
        system::ChemicalSystem,
        initial_state::ChemicalState,
        tspan::Tuple;
        calorimeter = nothing,
        activity_model::AbstractActivityModel = DiluteSolutionModel(),
        equilibrium_solver = nothing,
    )
    isempty(system.idx_kinetic) &&
        throw(
        ArgumentError(
            "ChemicalSystem has no kinetic species. " *
                "Pass kinetic_species to the ChemicalSystem constructor, " *
                "or use the 4-argument form KineticsProblem(cs, reactions, state, tspan)."
        )
    )
    kin_rxns = [KineticReaction(system, rxn) for rxn in system.reactions]
    return _build_kinetics_problem(
        system, kin_rxns, initial_state, tspan;
        calorimeter, activity_model, equilibrium_solver,
    )
end

"""
    _with_equilibrium_solver(kp::KineticsProblem, es) -> KineticsProblem

Return `kp` with its equilibrium solver replaced by `es`, or `kp` itself when
`es` is `nothing`. Used by `integrate` so that a solver passed on the
[`KineticsSolver`](@ref) reaches the ODE; a solver already set on the problem
wins when both are given, and the mismatch is reported.
"""
function _with_equilibrium_solver(kp::KineticsProblem, es)
    isnothing(es) && return kp
    if !isnothing(kp.equilibrium_solver)
        kp.equilibrium_solver === es || @warn(
            "an equilibrium solver is set on both the KineticsProblem and the " *
                "KineticsSolver; the one on the problem is used."
        )
        return kp
    end
    return KineticsProblem{
        typeof(kp.system), typeof(kp.kinetic_reactions), typeof(kp.calorimeter),
        typeof(es), typeof(kp.activity_model),
    }(
        kp.system, kp.kinetic_reactions, kp.initial_state, kp.tspan,
        kp.calorimeter, kp.activity_model, es,
        kp.idx_kinetic, kp.idx_equilibrium, kp.ν, kp.νe, kp.νk, kp.Ae,
    )
end

# ── build_u0 ─────────────────────────────────────────────────────────────────

"""
    build_u0(kp::KineticsProblem) -> Vector{Float64}

Build the initial ODE state vector.

Structure of `u`:
  - Without re-speciation: `u = [nₖ₁, …, nₖ_K, ξ₁, …, ξ_M]`
  - With re-speciation:    `u = [bₑ₁, …, bₑ_C, nₖ₁, …, nₖ_K, ξ₁, …, ξ_M]`
  - Semi-adiabatic adds `T` at the end: `u = [..., T₀]`

The extents of reaction `ξ` are carried alongside the kinetic moles, integrating
`dξ/dt = r`. They are redundant with `nₖ` — the two satisfy
`nₖ = nₖ(0) + νₖᵀ ξ` — but they are what makes the **non-kinetic** amounts
available inside the residual, through `n = n(0) + νᵀ ξ`. Without them, a rate
law gating on a species that is not itself kinetic would read a value frozen at
`t = 0`. They also make [`reaction_extents`](@ref) and [`state_at`](@ref) exact
rather than quadrature-limited.

The calorimeter's slot stays last and is addressed from the end of the vector,
so it is unaffected by the presence of `ξ`.
"""
# ── the number type of a run ─────────────────────────────────────────────────
#
# A run is differentiated with respect to whatever carries dual numbers: the
# amounts, temperature or pressure of the initial state, the constants of a
# calorimeter, the time span, and the parameters a rate law captures (a rate
# constant handed in as a dual by the function being differentiated). The state
# of the integrator and every buffer the run writes into its result are of the
# number type that covers them all. Anything narrower either raises or, worse,
# drops the derivative: the ODE interface promotes the state only when it finds
# the duals in the parameter object, and those of a rate law live in closures it
# does not look into.

"""
    _kinetics_number_type(kp) -> Type

The number type of a run of `kp`: `Float64`, or the dual type covering the
initial state, the time span, the calorimeter and every rate law.
"""
function _kinetics_number_type(kp::KineticsProblem)
    R = promote_type(Float64, _amount_number_type(kp.initial_state), typeof(float(kp.tspan[1])))
    for kr in kp.kinetic_reactions
        R = promote_type(R, _captured_number_type(kr.rate_fn), _captured_number_type(kr.heat_per_mol))
    end
    return promote_type(R, _captured_number_type(kp.calorimeter))
end

function build_u0(kp::KineticsProblem; R::Type = _kinetics_number_type(kp))
    n_mol = R[
        ustrip(us"mol", kp.initial_state.n[i])
            for i in eachindex(kp.system.species)
    ]
    # Kinetic species moles
    nk0 = n_mol[kp.idx_kinetic]
    # Extents of reaction, zero by definition at t = tspan[1]
    ξ0 = zeros(R, length(kp.kinetic_reactions))

    u0 = if isnothing(kp.equilibrium_solver)
        vcat(nk0, ξ0)
    else
        # Element amounts in equilibrium partition: bₑ = Aₑ nₑ
        ne0 = n_mol[kp.idx_equilibrium]
        be0 = kp.Ae * ne0
        vcat(be0, nk0, ξ0)
    end

    # One trailing slot per calorimeter, and only one: the temperature for the
    # semi-adiabatic device, the accumulated heat for the isothermal one. Both are
    # `u[end]`, and `p.has_T` / `p.has_Q` say which.
    if kp.calorimeter isa SemiAdiabaticCalorimeter
        push!(u0, R(safe_ustrip(us"K", kp.calorimeter.T0)))
    elseif kp.calorimeter isa IsothermalCalorimeter
        push!(u0, zero(R))
    end

    return u0
end

# ── build_kinetics_params ────────────────────────────────────────────────────

# A vector whose elements are each of their own type -- one compiled
# thermodynamic function per species, one rate law per kinetic reaction -- has
# no concrete element type and cannot have one. SciMLBase reads such a field of
# the ODE parameters as a mistake and warns on every run, advising a tuple, which
# would compile a method per system. The wrapper is deliberately not an
# `AbstractArray`, which is what SciMLBase inspects; it changes nothing to how the
# elements are reached -- iteration and indexing are the vector's -- and only
# keeps the warning out of the output.
struct _Heterogeneous{V <: AbstractVector}
    fns::V
end
Base.iterate(s::_Heterogeneous, st...) = iterate(s.fns, st...)
Base.length(s::_Heterogeneous) = length(s.fns)
Base.getindex(s::_Heterogeneous, i) = s.fns[i]
Base.eltype(::Type{_Heterogeneous{V}}) where {V} = eltype(V)

"""
    build_kinetics_params(kp::KineticsProblem; ϵ=1e-30) -> NamedTuple

Build the immutable parameter tuple `p` passed to the ODE function.

Key fields: `T`, `P`, `ϵ`, `lna_fn`, `kin_rxns`, `species_index`,
`n_initial_full`, `n_full`, `cp_fns`, `rates_buf`, index ranges
`n_be`, `n_nk`, `idx_kinetic`, `idx_equilibrium`, `νe`, `νk`, `Ae`.
"""
function build_kinetics_params(kp::KineticsProblem; ϵ::Float64 = 1.0e-30, R::Type = _kinetics_number_type(kp))
    state = kp.initial_state
    # Every value the run writes into its result is of the number type `R` of
    # the run (`_kinetics_number_type`). The starting guesses of the
    # re-speciation are the exception, and stay `Float64` on purpose: a starting
    # point carries no derivative, the derivative of the partition coming from
    # the implicit-function theorem at the certified answer.
    T_K = R(ustrip(us"K", temperature(state)))
    P_Pa = R(ustrip(us"Pa", pressure(state)))

    lna_fn = activity_model(kp.system, kp.activity_model)

    # Species name → index dict (built once, shared by all StateViews)
    species_index = Dict{String, Int}()
    for (i, sp) in enumerate(kp.system.species)
        species_index[phreeqc(formula(sp))] = i
        sym = ChemistryLab.symbol(sp)
        !isempty(sym) && (species_index[sym] = i)
    end

    n_sp = length(kp.system.species)
    n_initial_full = R[ustrip(us"mol", state.n[i]) for i in 1:n_sp]
    n_full = copy(n_initial_full)

    cp_fns = _Heterogeneous([haskey(sp, :Cp⁰) ? sp[:Cp⁰] : nothing for sp in kp.system.species])
    h_fns = _Heterogeneous([haskey(sp, :ΔₐH⁰) ? sp[:ΔₐH⁰] : nothing for sp in kp.system.species])

    kin_rxns = _Heterogeneous(kp.kinetic_reactions)
    rates_buf = zeros(R, length(kin_rxns))

    eq_sys, eq_sub, n_eq_init = if isnothing(kp.equilibrium_solver)
        nothing, nothing, Float64[]
    else
        sys_e = _equilibrium_subsystem(kp.system, kp.idx_equilibrium)
        es = kp.equilibrium_solver
        # The sub-solver must be rebuilt because `sys_e` is a different system
        # from `kp.system`, so the compiled `μ` does not carry over. It is
        # rebuilt with the model the USER's solver was built from — not with
        # `kp.activity_model`, which defaults to `DiluteSolutionModel()` and
        # would silently downgrade an HKF run on a cement pore solution.
        es_model = activity_model(es)
        if typeof(es_model) !== typeof(kp.activity_model)
            @warn """
            The equilibrium solver and the kinetics problem carry different \
            activity models. The equilibrium sub-solve uses the solver's \
            ($(nameof(typeof(es_model)))), while the log-activities handed to the \
            rate laws use the problem's ($(nameof(typeof(kp.activity_model)))). \
            Pass `activity_model = $(nameof(typeof(es_model)))()` to \
            `KineticsProblem` to make them agree.""" maxlog = 1
        end
        (
            sys_e,
            EquilibriumSolver(
                sys_e, es_model, es.solver;
                variable_space = es.variable_space, es.kwargs...
            ),
            Float64[_plain(n_initial_full[i]) for i in kp.idx_equilibrium],
        )
    end

    # State layout sizes
    n_be = isnothing(kp.equilibrium_solver) ? 0 : size(kp.Ae, 1)
    n_nk = length(kp.idx_kinetic)
    n_rxn_state = length(kp.kinetic_reactions)
    has_T = kp.calorimeter isa SemiAdiabaticCalorimeter
    has_Q = kp.calorimeter isa IsothermalCalorimeter

    # The heat under partial equilibrium: `−dH/dt` over the whole composition,
    # the equilibrium partition followed through its sensitivity to `bₑ` (see
    # `_equilibrium_heat_rate`). It needs the enthalpy of every species, since a
    # species without one would drop out of the sum and take its heat with it.
    heat_eq = (has_T || has_Q) && n_be > 0
    heat_eq && _refuse_missing_enthalpy(kp.system, h_fns)
    n_eq = length(kp.idx_equilibrium)

    # The element content of every species, and the totals the system was given:
    # no amount of a kinetic species can hold more of an element than there is.
    # See `_kinetic_state_infeasible`.
    elements = sort!(collect(setdiff(union((keys(atoms_charge(sp)) for sp in kp.system.species)...), (:Zz,))))
    E_all = [Float64(get(atoms_charge(sp), el, 0)) for el in elements, sp in kp.system.species]
    B_el = E_all * Float64[_plain(x) for x in n_initial_full]

    # Calorimeter parameters (semi-adiabatic)
    cal = kp.calorimeter
    Cp_calo = cal isa SemiAdiabaticCalorimeter ? R(safe_ustrip(us"J/K", cal.Cp)) : zero(R)
    T_env = cal isa SemiAdiabaticCalorimeter ? R(safe_ustrip(us"K", cal.T_env)) : T_K
    heat_loss_fn = cal isa SemiAdiabaticCalorimeter ? cal.heat_loss : identity

    return (
        T = T_K,
        P = P_Pa,
        ϵ = ϵ,
        has_Q = has_Q,
        lna_fn = lna_fn,
        kin_rxns = kin_rxns,
        species_index = species_index,
        n_initial_full = n_initial_full,
        n_full = n_full,
        cp_fns = cp_fns,
        h_fns = h_fns,
        rates_buf = rates_buf,
        # Index layout
        n_be = n_be,
        n_nk = n_nk,
        n_rxn_state = n_rxn_state,
        has_T = has_T,
        idx_kinetic = kp.idx_kinetic,
        idx_equilibrium = kp.idx_equilibrium,
        # Leal partitions
        νe = Float64.(kp.νe),
        νk = Float64.(kp.νk),
        Ae = Float64.(kp.Ae),
        # Calorimeter
        Cp_calo = Cp_calo,
        T_env = T_env,
        heat_loss_fn = heat_loss_fn,
        # The sensitivity `∂nₑ/∂bₑ` of the equilibrium partition, refreshed by
        # `respeciate!`.
        heat_eq = heat_eq,
        heat_S = Ref(zeros(R, n_eq, n_be)),
        heat_ready = Ref(false),
        # The partition and the element amounts the sensitivity was taken at, and
        # the heat of the part of the last re-speciation it did not predict.
        heat_n = Ref(zeros(R, n_eq)),
        heat_b = Ref(zeros(R, n_be)),
        heat_jump = Ref(zero(R)),
        # In a semi-adiabatic cell, the shift of the partition with temperature,
        # `∂nₑ/∂T`, and the temperature it was taken at.
        heat_dndT = Ref(zeros(R, n_eq)),
        heat_T = Ref(T_K),
        # Equilibrium — Leal et al. (2017) §5. The re-speciation φ(bₑ) is a
        # minimization over the EQUILIBRIUM PARTITION ONLY, at frozen kinetic
        # amounts. Running it over the whole system would let the kinetic
        # minerals equilibrate instantaneously, which is exactly what a kinetic
        # description exists to prevent.
        eq_system = eq_sys,
        eq_solver = eq_sub,
        # The certifying solver over the SAME partition, when one can be built.
        #
        # Re-speciation used to go through the interior point alone, and that is
        # the unreliable path: on calcite under `r = k(1 − Ω)` it returned a
        # partition violating the element balance by 467 mol. (The trajectory of
        # that run, a reaction extent of −457 mol, had another cause, the
        # speciation frozen within a step under a rate law that reads it: see
        # `build_kinetics_ode`.)
        #
        # With the certified route the certificate DECIDES, so a partition that
        # does not conserve matter is not accepted in the first place. The
        # partitions it cannot treat, without an aqueous phase or without
        # `H2O@`, are screened out beforehand; any other failure to build it is
        # an error, reported rather than replaced by the interior point. It is
        # built with the model of the user's solver, as the interior point is
        # above: the problem's own model defaults to the dilute one.
        #
        # Its search is the one for warm-started sequences: a solute may fall to
        # its potential in one sweep, and the line search does not ask a
        # candidate for a converged inner solve the current point lacks. Each
        # step starts from the previous one, where the defaults spent most of the
        # run: 176 s against 10 s for three hours of a CEM I paste, on the same
        # trajectory to the last bit. The certificate is unchanged.
        eq_dual = (
                isnothing(kp.equilibrium_solver) || !_DUAL_AVAILABLE[] ||
                !_dual_applicable(eq_sys)
            ) ? nothing : DualEquilibriumSolver(
                eq_sys, activity_model(kp.equilibrium_solver);
                inner_fall_bound = Inf, lenient_line_search = true,
            ),
        n_eq_init = n_eq_init,
        n_eq_buf = similar(n_eq_init),
        n_eq_buf2 = similar(n_eq_init),
        T_q = Ref(T_K * u"K"),
        P_q = Ref(P_Pa * u"Pa"),
        state_ref = Ref{ChemicalState}(state),
        eq_failures = Ref(0),
        # Worst |Aₑ n − bₑ|∞ over the run. This, not the optimizer's retcode, is
        # the criterion with physical meaning: a solve can stall short of its
        # tolerance and still satisfy the element balance to machine precision.
        eq_worst_residual = Ref(0.0),
        eq_worst_abs = Ref(0.0),
        eq_worst_abs_acc = Ref(0.0),
        on_accepted = Ref(false),
        # Set once a speciation exists, so `respeciate!` can warm-start from it.
        eq_warm = Ref(false),
        # Where the right-hand side reads the equilibrium partition: `:frozen`,
        # as the last accepted step left it, or `:rhs`, at the state it is
        # evaluated at (see `build_kinetics_ode`). `integrate` sets it, from
        # whether a rate law reads the partition (`_rates_read_speciation`).
        rhs_mode = Ref(:frozen),
        rates_read_speciation = Ref(false),
        # The last partition solved in `:rhs` mode, keyed by the values of `bₑ`
        # and of the temperature it was solved at.
        rhs_cache = Ref{Any}(nothing),
        n_rhs = similar(n_full),
        # Feasibility of a kinetic state: the element content of the kinetic
        # species and the element totals of the system.
        E_kin = E_all[:, kp.idx_kinetic],
        B_el = B_el,
        # Whether the LAST respeciation had to fall back on the reconstruction
        # because the warm start was in the wrong basin. An assemblage switch is
        # not a single-step event -- a phase takes several steps to exhaust --
        # so this says which guess to try FIRST on the next call. See
        # `_respeciate_solve!`.
        eq_switching = Ref(false),
    )
end

# ── Equilibrium partition sub-system ─────────────────────────────────────────

"""
    _equilibrium_subsystem(system, idx_equilibrium) -> ChemicalSystem

The chemical system restricted to the equilibrium partition, as the partitioned
formulation of [Leal2017](@citet) requires.

Its formula matrix is exactly `system.CSM.A[:, idx_equilibrium]`, in the same
species order, so the element amounts `bₑ` carried by the ODE state are handed
to it unchanged. The primaries are those of the parent system that survive the
restriction — the kinetic minerals never do, not being in the partition.
"""
function _equilibrium_subsystem(system::ChemicalSystem, idx_equilibrium)
    sub_species = system.species[idx_equilibrium]
    # The parent's primaries are the row labels of its stoichiometric matrix,
    # `system.SM.primaries` — not `idx_components`, which is a class-based index
    # and returns a different, much smaller set. Getting this wrong leaves the
    # sub-system with only a couple of conservation constraints for a dozen
    # species: the minimization is then wildly under-determined and returns an
    # assemblage with essentially no water left.
    comp_names = Set(symbol(sp) for sp in system.SM.primaries)
    prim = [sp for sp in sub_species if symbol(sp) in comp_names]

    # Carry the solid solutions over. Dropping them — as this did until 0.8.2 —
    # is silent and total: the parent may declare CSHQ, AFm or Hydrogarnet, and
    # the partition the equilibrium is actually solved on knows nothing of them,
    # so their end-members are treated as separate pure phases and the mixing
    # entropy never enters the Gibbs energy. Measured on alite and belite with
    # the four CSHQ end-members, the run was bit-identical with and without the
    # declaration, and produced no C-S-H at all: the silicon stayed in solution
    # and the portlandite came out at 4.52 mol against 2.93 with a Jennite
    # end-member.
    #
    # A solution survives only if ALL its end-members are in the partition. One
    # whose members are split between the kinetic and equilibrium sides is not a
    # phase the equilibrium can mix, and passing it truncated would be worse than
    # dropping it.
    sub_names = Set(symbol(sp) for sp in sub_species)
    ss = system.solid_solutions
    sub_ss = if ss === nothing
        nothing
    else
        kept = [
            phase for phase in ss
                if all(em -> symbol(em) in sub_names, phase.end_members)
        ]
        isempty(kept) ? nothing : kept
    end

    # And the surface site families, on the same all-or-nothing rule and for the
    # same reason: a family the partition does not know about is a set of
    # species with no shared budget and no mixing, which is not a surface.
    #
    # A family split between the two partitions is refused rather than dropped.
    # Dropping it would leave its members in `sub_species` as `AS_SURFACE`
    # species belonging to no family, which `ChemicalSystem` refuses anyway —
    # but with a message about orphans rather than about the split that caused
    # them. Sites equilibrate fast by construction in this release, so a split is
    # a declaration error, and it is worth saying which family and why.
    families = system.site_families
    sub_families = if families === nothing
        nothing
    else
        for f in families
            present = count(sp -> symbol(sp) in sub_names, site_members(f))
            0 < present < length(site_members(f)) && throw(
                ArgumentError(
                    "SiteFamily \"$(name(f))\" is split between the kinetic and " *
                        "equilibrium partitions ($present of " *
                        "$(length(site_members(f))) members on the equilibrium " *
                        "side). Its members share one site budget, so they have to " *
                        "be solved together; declare the whole family as fast, or " *
                        "none of it.",
                )
            )
        end
        kept = [
            f for f in families
                if all(sp -> symbol(sp) in sub_names, site_members(f))
        ]
        # A budget that follows its host needs the host in the same partition:
        # a host whose amount a rate law moves would change the site budget
        # between two re-speciations, which nothing here accounts for.
        for f in kept
            sup = surface_support(f)
            sup.coupling === SITES_FOLLOW_HOST || continue
            sup.host in sub_names || throw(
                ArgumentError(
                    "SiteFamily \"$(name(f))\" follows host \"$(sup.host)\", which " *
                        "is a kinetic species of this problem. A site budget that " *
                        "follows a host whose amount a rate law controls is not " *
                        "supported: put the host in the equilibrium partition, or " *
                        "keep the support at SITES_FIXED.",
                )
            )
        end
        isempty(kept) ? nothing : kept
    end

    # A coupled family may carry its site row on a bare site component, which is
    # a primary of the parent and not a species; it goes with its family.
    if sub_families !== nothing
        for pr in system.SM.primaries
            symbol(pr) in sub_names && continue
            any(f -> _is_bare_site(pr, f.site), sub_families) && push!(prim, pr)
        end
    end

    return ChemicalSystem(
        sub_species, isempty(prim) ? sub_species : prim;
        solid_solutions = sub_ss,
        site_families = sub_families,
    )
end

"""
    _EQ_GUESS_FLOOR

Lower floor applied to the starting guess of the equilibrium sub-solve, chosen
strictly above the `1e-16` lower bound that `EquilibriumProblem` imposes.

An interior-point method started on its own bound stalls short of its tolerance
and returns `MaxIters`. Loosening that tolerance is *not* the fix: measured
against Reaktoro on the calcite reference case, the worst species error grows
from 4.3 % at `tol = 1e-10` to 38 % at `1e-8` and 252 % at `1e-7`. Moving the
guess inside keeps the tight tolerance and removes the stalls.
"""
const _EQ_GUESS_FLOOR = 1.0e-10

"""
    EQ_RESIDUAL_TOL

Largest per-element balance violation, relative to that element's own total,
that still counts as a converged speciation. Above it the point is not handed
on as a warm start and the solve is retried from a budget-clipped guess.
"""
const EQ_RESIDUAL_TOL = 1.0e-8

"""
    _RETRY_ABS_TOL

Element-balance violation in moles above which a replayed speciation is solved a
second time from a guess carrying no active set. Chosen well above machine
precision and well below anything chemically meaningful.
"""
const _RETRY_ABS_TOL = 1.0e-6

"""
    _CONTINUATION_STEPS

How many bisection steps a replay may take between the last certified instant and
one it cannot certify directly. Each step halves the jump in the component
totals, so eight of them reduce it by a factor 256.
"""
const _CONTINUATION_STEPS = 8

"""
    RESTORE_MAXIT

Alternating-projection sweeps allowed when restoring the feasibility of an
in-run guess. Exposed because the right value is a trade: the projection
converges linearly and a cement can need tens of thousands of sweeps, while this
runs at every re-speciation of a run.
"""
const RESTORE_MAXIT = Ref(200)

# ── enthalpy tracking ────────────────────────────────────────────────────────

"""
    system_enthalpy(p, u, T) -> Float64

`Σᵢ nᵢ ΔₐH⁰ᵢ(T)` [J] over the composition of the ODE state `u`: the kinetic
minerals read from `u` itself, the equilibrium partition from `p.n_full`.

# Why the heat cannot come from the kinetic reactions here

`heat_rate` sums `rᵢ (−ΔᵣH⁰ᵢ)` over the KINETIC reactions, which is right when
those reactions produce the hydrates. Under partial equilibrium they do not: they
dissolve the anhydrous phases into ions, and the hydrates are precipitated by the
Gibbs minimization, whose heat that sum cannot see. Measured on an ordinary
Portland cement, counting only the dissolution put the semi-adiabatic temperature
rise at 207 K where the test gives some tens of kelvin.

The enthalpy of the whole system has no such blind spot. It is a state function,
so the heat released between two states at the same temperature is their
difference — reactants, ions and hydrates all counted once, with no reaction
stoichiometry to write down. This is Eq. (17)–(21) of Lavergne et al. (2018).
"""
function system_enthalpy(p, u, T)
    # The kinetic amounts are read from the ODE state, NOT from `p.n_full`.
    #
    # `respeciate!` writes only the equilibrium partition into that buffer; the
    # anhydrous amounts in it were last written by the right-hand side, which the
    # stiff solver also evaluates at Jacobian probes and rejected steps. Summing
    # the buffer as it stands therefore pairs an accepted equilibrium with some
    # neighboring point's clinker, and the resulting enthalpy is not a state of
    # the trajectory at all: the recorded heat came out NON-MONOTONE, 936 J/g at
    # one day and 631 J/g at two, which no calorimeter has ever measured.
    kin = p.idx_kinetic
    H = zero(promote_type(eltype(u), typeof(T), eltype(p.n_full)))
    @inbounds for (i, h_fn) in enumerate(p.h_fns)
        isnothing(h_fn) && continue
        j = findfirst(==(i), kin)
        nᵢ = j === nothing ? p.n_full[i] : max(u[p.n_be + j], p.ϵ)
        H += nᵢ * h_fn(; T = T, unit = false)
    end
    return H
end

# ── the heat under partial equilibrium ───────────────────────────────────────
#
# Under partial equilibrium the kinetic reactions only dissolve the anhydrous
# phases into ions, and the hydrates are precipitated by the minimization: the
# heat of the kinetic reactions leaves the precipitation out (on an ordinary
# Portland cement, a semi-adiabatic rise of 207 K). The heat is the rate at which
# the enthalpy of the WHOLE composition falls at fixed temperature,
#
#     q̇ = −dH/dt = −Σᵢ ΔₐH⁰ᵢ(T) dnᵢ/dt ,
#
# the kinetic amounts moving as the ODE moves them, and the equilibrium partition
# as the minimization moves it when `bₑ` changes: `dnₑ/dt = S dbₑ/dt`, with
# `S = ∂nₑ/∂bₑ` from the optimality conditions at the partition in hand
# (`_equilibrium_sensitivity`). Integrated at fixed temperature this is the
# enthalpy difference `heat_release` computes from certified states.

"""
    _refuse_missing_enthalpy(system, h_fns)

Refuse a calorimeter under partial equilibrium when a species carries no
enthalpy of formation: its term would drop out of `−dH/dt`, and the heat with it.
"""
function _refuse_missing_enthalpy(system, h_fns)
    missing_h = [symbol(sp) for (sp, h) in zip(system.species, h_fns) if isnothing(h)]
    isempty(missing_h) && return nothing
    shown = join(missing_h[1:min(5, length(missing_h))], ", ")
    throw(
        ArgumentError(
            "a calorimeter under partial equilibrium takes its heat from the enthalpy " *
                "of the whole composition, and $(length(missing_h)) species carry no " *
                "ΔₐH⁰ ($shown$(length(missing_h) > 5 ? ", …" : "")). Give them one, or " *
                "leave them out of the system."
        )
    )
end

"""
    _heat_sensitivity!(p, n_e, be)

Refresh `∂nₑ/∂bₑ` at the partition `n_e` just computed by `respeciate!` for the
element amounts `be`, from the optimality conditions of that minimization. One
Hessian of the potentials serves every column.

Between two re-speciations the partition is followed linearly, `nₑ + S Δbₑ`.
The re-speciation lands elsewhere wherever that is not the whole story: where
the assemblage changes within the step, a phase appearing for instance, and in a
semi-adiabatic cell where the temperature moved the equilibrium. The enthalpy of
that difference is heat the integrated source did not see.

In a semi-adiabatic cell the partition is also followed in temperature,
`nₑ + S Δbₑ + (∂nₑ/∂T) ΔT`, and the heat it takes up as it shifts is part of the
heat capacity of the cell ([`_equilibrium_shift_capacity`](@ref)). The shift is
taken from the same optimality conditions with the right-hand side of the
Gibbs–Helmholtz relation, `∂(μᵢ/RT)/∂T = −ΔₐH⁰ᵢ/RT²`, over the enthalpies the heat
is counted with. The capacity is then a quadratic form in them, positive where
the optimality conditions are well posed: the stability of an equilibrium, which
heating shifts towards where it absorbs heat. A difference quotient of the
potentials instead carries the temperature dependence of the activity
coefficients, which the enthalpies do not, and it gave the C100 mortar of
Lavergne et al. (2018) a capacity small or negative, and temperatures no water
table covers. Where the form is not positive all the same, the shift is left
to the jump below, in the prediction and in the capacity alike.

The difference between the prediction and the re-speciation is recorded in
`p.heat_jump` (J, positive when released) for the step callback to add to the
calorimeter's state, so that the heat follows the enthalpy of the partition at
every accepted step. The first re-speciation, of the state the run starts from,
is the reference and releases nothing.
"""
function _heat_sensitivity!(p, n_e, be)
    # The reference is a PROVED equilibrium, or nothing moves. The in-run
    # partition is warm-started and not certified, and an assemblage one phase
    # off is worth more than a step's heat: taken as it came, successive
    # partitions of the C100 mortar of Lavergne et al. (2018) differed by
    # 120 kJ, the jumps became tens of kelvin and the Arrhenius terms ran away.
    # An unproved partition leaves the last proved one, and its sensitivity, in
    # place; the next proved one then settles the whole difference.
    n_e = _proved_partition(p, n_e, be)
    n_e === nothing && return nothing
    T = ustrip(us"K", p.T_q[])
    if p.heat_ready[]
        S, n0, b0 = p.heat_S[], p.heat_n[], p.heat_b[]
        ΔT = p.has_T ? T - p.heat_T[] : zero(T)
        jump = zero(eltype(p.heat_n[]))
        for (j, idx) in enumerate(p.idx_equilibrium)
            predicted = n0[j] + p.heat_dndT[][j] * ΔT
            for k in eachindex(b0)
                predicted += S[j, k] * (be[k] - b0[k])
            end
            jump -= p.h_fns[idx](; T = T, unit = false) * (n_e[j] - predicted)
        end
        # Accumulated, not assigned: `respeciate!` may solve twice in one call
        # (the warm guess, then the reconstruction), and the second is measured
        # from the first, so the two parts add up to the whole.
        p.heat_jump[] += jump
    end
    p.heat_n[] .= n_e
    p.heat_b[] .= be
    p.heat_T[] = T
    # `S = ∂nₑ/∂bₑ` and the shift with temperature, from OptimaSolver's tangent
    # at the proved partition (`_partition_sensitivity`), in the number type of
    # the run: a run differentiated with respect to its parameters carries the
    # derivative of `S` as well.
    S, dndT = _partition_sensitivity(p, n_e, be, T)
    p.heat_S[] .= S
    if p.has_T
        # The capacity it implies is a quadratic form, positive when the
        # optimality conditions are well posed. When it is not, the shift is not
        # followed at all, in the prediction as in the heat capacity, and the
        # re-speciation's jump carries it: counting it in one and not the other
        # creates or destroys heat.
        C = sum(p.h_fns[idx](; T = T, unit = false) * dndT[j] for (j, idx) in enumerate(p.idx_equilibrium))
        p.heat_dndT[] .= isfinite(C) && C > 0 ? dndT : zero(dndT)
    end
    p.heat_ready[] = true
    return nothing
end

"""
    _partition_sensitivity(p, n_e, be, T) -> (S, dndT)

`S = ∂nₑ/∂bₑ` at the proved partition `n_e` for the budget `be`, and, in a
semi-adiabatic cell, `∂nₑ/∂T` along the Gibbs–Helmholtz direction
`∂(μᵢ⁰/RT)/∂T = −ΔₐH⁰ᵢ/RT²` (`nothing` otherwise), both from the
implicit-function theorem at the answer with its active set frozen
(OptimaSolver's `dual_newton_tangent`). The direction of the temperature is that
of the standard potentials alone, over the enthalpies the heat is counted with:
the activity coefficients are held, as the heat is.

The active set is read at the certificate's floor: an interior-point partition
leaves an absent pure phase at a negligible amount rather than at zero, and a
pure phase left free there keeps its row of the optimality conditions, as if it
coexisted with the rest. Until 0.28.2 this was a saddle-point system solved by a
truncated singular value decomposition, in `Float64`, with absent phases pinned
by a threshold of its own.
"""
function _partition_sensitivity(p, n_e, be, T)
    des = p.eq_dual
    des === nothing && throw(
        ArgumentError(
            "the heat of a partial equilibrium needs the certified solver of the " *
                "partition: OptimaSolver loaded, an aqueous phase and `H2O@` in it.",
        ),
    )
    st = ChemicalState(p.eq_system, n_e .* u"mol"; T = p.T_q[], P = p.P_q[])
    pv = _build_params(st; ϵ = p.ϵ)
    nb = length(be)
    np = p.has_T ? nb + 1 : nb
    V = promote_type(eltype(n_e), eltype(be), typeof(T), eltype(pv.ΔₐG⁰overRT))
    Tg = typeof(ForwardDiff.Tag(_partition_sensitivity, V))
    D = ForwardDiff.Dual{Tg, V, np}
    unit(k) = ForwardDiff.Partials{np, V}(ntuple(i -> V(i == k), np))
    bd = D[D(V(be[k]), unit(k)) for k in 1:nb]
    g = if p.has_T
        gT = [-p.h_fns[idx](; T = T, unit = false) / (R_GAS * T^2) for idx in p.idx_equilibrium]
        D[
            D(V(pv.ΔₐG⁰overRT[j]), ForwardDiff.Partials{np, V}(ntuple(i -> i == np ? V(gT[j]) : zero(V), np)))
                for j in eachindex(gT)
        ]
    else
        pv.ΔₐG⁰overRT
    end
    pd = merge(pv, (ΔₐG⁰overRT = g,))
    blocks = _constraint_blocks(FixedTP(), des, st, pd, n_e)
    prob = _dual_problem(des, pd, n_e, blocks)
    t = _optima_tangent(prob, bd, n_e; floor = _CERTIFICATE_FLOOR)
    ne = length(n_e)
    S = [ForwardDiff.partials(t.x[j], k) for j in 1:ne, k in 1:nb]
    dndT = p.has_T ? [ForwardDiff.partials(t.x[j], np) for j in 1:ne] : nothing
    return S, dndT
end

"""
    _equilibrium_shift_capacity(p, T) -> Real

`Σᵢ ΔₐH⁰ᵢ(T) ∂nᵢ/∂T` [J/K] over the equilibrium partition: the heat it takes up
per kelvin as it shifts, in the heat capacity of a semi-adiabatic cell. The
shift is kept only where this is positive at the reference temperature
(`_heat_sensitivity!`), and the clamp covers what the change of the enthalpies
with temperature can do to it away from there.
"""
function _equilibrium_shift_capacity(p, T)
    dndT = p.heat_dndT[]
    c = zero(T)
    for (j, idx) in enumerate(p.idx_equilibrium)
        c += p.h_fns[idx](; T = T, unit = false) * dndT[j]
    end
    return max(c, zero(c))
end

"""
    _proved_partition(p, n_e, be) -> Union{Vector, Nothing}

`n_e` if the certificate proves it the equilibrium for `be`, else the certified
solve started from it if that proves one, else `nothing`. Without a certifying
solver for the partition, `n_e` as it is.
"""
function _proved_partition(p, n_e, be)
    p.eq_dual === nothing && return n_e
    st = ChemicalState(p.eq_system, n_e .* u"mol"; T = p.T_q[], P = p.P_q[])
    try
        optimality_certificate(p.eq_dual, st; b = be).optimal && return n_e
        eq_c, cert = solve_certified(p.eq_dual, (st,); b = be, ϵ = p.ϵ)
        cert.optimal && return [ustrip(us"mol", x) for x in eq_c.n]
    catch
        # An audit or a solve that raises proves nothing, and nothing moves.
    end
    return nothing
end

"""
    _equilibrium_heat_rate(p, du, T) -> Real

`q̇ = −Σᵢ ΔₐH⁰ᵢ(T) dnᵢ/dt` [W], the kinetic amounts from `du`, the equilibrium
partition from `S dbₑ/dt`. Generic in the number type of `du` and `T`.
"""
function _equilibrium_heat_rate(p, du, T)
    S = p.heat_S[]
    q = zero(promote_type(eltype(du), typeof(T)))
    for (j, idx) in enumerate(p.idx_equilibrium)
        dn = zero(eltype(du))
        for k in 1:(p.n_be)
            dn += S[j, k] * du[k]
        end
        q -= p.h_fns[idx](; T = T, unit = false) * dn
    end
    for (j, idx) in enumerate(p.idx_kinetic)
        q -= p.h_fns[idx](; T = T, unit = false) * du[p.n_be + j]
    end
    return q
end

"""
    _cell_heat_capacity(p, T) -> Float64

The heat capacity of the semi-adiabatic cell at the composition `p.n_full`: the
vessel, `Σᵢ nᵢ Cp°ᵢ(T)` and the shift of the equilibrium partition, the
denominator of `dT/dt`, for converting a heat into a temperature step.
"""
function _cell_heat_capacity(p, T)
    c = p.Cp_calo
    for (i, cp_fn) in enumerate(p.cp_fns)
        isnothing(cp_fn) && continue
        c += p.n_full[i] * cp_fn(; T = T, unit = false)
    end
    p.heat_eq && p.heat_ready[] && (c += _equilibrium_shift_capacity(p, T))
    return c
end

"""
    _apply_heat_jump!(p, u) -> Bool

Add the heat of the last re-speciation's unpredicted part (`p.heat_jump`) to the
calorimeter's state `u[end]`: to the heat of an isothermal cell, or as a
temperature step to a semi-adiabatic one. Returns whether `u` changed.
"""
function _apply_heat_jump!(p, u)
    p.heat_eq || return false
    q = p.heat_jump[]
    p.heat_jump[] = zero(q)
    iszero(q) && return false
    if p.has_Q
        u[end] += q
    elseif p.has_T
        u[end] += q / _cell_heat_capacity(p, u[end])
    end
    return true
end

# ── respeciate! ──────────────────────────────────────────────────────────────

"""
    respeciate!(p, u) -> Bool

Solve the equilibrium sub-problem once, and write the result into the running
composition `p.n_full`. Returns `true` when a solve actually happened.

This is the second half of the operator-splitting step: the ODE advances the
kinetic minerals with the speciation held frozen, then this function
re-equilibrates the equilibrium partition under the element amounts the ODE has
just produced.

The element amounts `bₑ` carried by the state vector are the constraint of that
sub-problem (Leal et al. 2017, Eq. 54). `solve` conserves `A·n`, so what has to
be handed to it is a composition whose element totals are exactly `bₑ` — here
the previous speciation, projected onto `bₑ` through the pseudo-inverse of
`Aₑ`. Handing over `p.n_full` unchanged, as an earlier version did, discards
`bₑ` entirely and leaves the element balance to drift.
"""
function respeciate!(p, u)
    p.n_be > 0 || return false

    # A semi-adiabatic cell carries the temperature in the state, and the
    # partition is an equilibrium at THAT temperature, not at the initial one.
    p.has_T && (p.T_q[] = u[end] * us"K")

    # The right-hand side solves the partition itself: the accepted one is the
    # same solve, kept as the warm start of the next and for the heat.
    p.rhs_mode[] === :rhs && _respeciate_rhs!(p, u) && return true

    # φ(bₑ), Leal et al. (2017) Eq. 54: the element amounts carried by the ODE
    # state ARE the constraint of the minimization, and they are handed to the
    # solver as `b`.
    #
    # That is the whole reason the state integrates `bₑ` rather than `nₑ`. Along
    # the way an individual species may want to go negative — the generated
    # dissolution reactions are written in H⁺, and a cement paste contains no
    # acid — and it is the minimizer, not the caller, that redistributes the
    # elements over a feasible set. An earlier version reconstructed `nₑ` from
    # `bₑ` through `pinv(Aₑ)` and clamped the result at `ϵ`; the clamp destroyed
    # the balance the projection had just established, and the solve went on to
    # return amounts of 1e65.
    #
    # The composition below is a starting guess only, and does not have to carry
    # `bₑ`. It is built from the reaction extents, which come free from the
    # kinetic amounts: each reaction carries ν = −1 on its controlling mineral
    # (`_normalise_to_mineral!`), so ξⱼ = nₖⱼ(0) − nₖⱼ.
    be = collect(@view u[1:(p.n_be)])

    nk = @view u[(p.n_be + 1):(p.n_be + p.n_nk)]
    n_eq = p.n_eq_buf

    # WARM START. `p.n_full` holds the speciation left by the previous accepted
    # step, which is an equilibrium for a nearby `bₑ` — by far the best guess
    # available. The stoichiometric reconstruction below is only the cold start,
    # used on the first call before any speciation exists.
    #
    # This matters far more than it looks. The reconstruction places every
    # dissolved element in solution with ZERO hydrates, i.e. a wildly
    # supersaturated composition, which for an aqueous-only system like the
    # calcite reference is harmless and for a cement is nearly the worst possible
    # starting point: more than half the solves failed, no hydrate ever
    # precipitated, and the pore solution came out at pH 6 instead of 12.6.
    if p.eq_warm[]
        for (j, idx) in enumerate(p.idx_equilibrium)
            n_eq[j] = max(_plain(p.n_full[idx]), _EQ_GUESS_FLOOR)
        end
        # The previous speciation was an equilibrium for the PREVIOUS `bₑ`. When
        # an element has since been spent — the sulfate of an OPC once the gypsum
        # is gone — that guess demands more of it than now exists, and the solve
        # starts outside the feasible set. Clipping costs nothing when the guess
        # is already feasible, which is the ordinary case.
        _budget_clip!(n_eq, p.Ae, _plain.(be))
        _restore_feasibility!(n_eq, p.Ae, _plain.(be); maxit = RESTORE_MAXIT[])
        _respeciate_solve!(p, n_eq, be) && return true
        # Infeasible or stalled: fall through to the cold reconstruction rather
        # than carry the bad point into the next step.
    end

    for j in eachindex(n_eq)
        # Start from the composition the specimen was cast with — for a paste,
        # the mixing water and nothing precipitated — and let
        # `_restore_feasibility!` below carry it onto the current `bₑ`.
        #
        # The stoichiometric reconstruction `νₑᵀξ` that used to be added here
        # placed every dissolved element in solution with ZERO hydrates. For an
        # aqueous-only system that is harmless, but for a cement it is close to
        # the worst possible start: it is supersaturated in every phase at once,
        # and its H⁺ entry is strongly negative (−6 per mole of alite) so it is
        # clamped to the floor, losing the acidity that the hydroxides have to
        # balance. Started there, the back-end stops next to its own guess and
        # returned an assemblage demanding 174 % of the sulfate present.
        #
        # Floor strictly inside the box, not at `p.ϵ`: `EquilibriumProblem`
        # raises anything below 1e-16 to exactly its lower bound, and an
        # interior-point method started on its own bound stalls — the calcite
        # reference reported 6 non-converged solves out of 8 steps for that
        # reason alone.
        n_eq[j] = max(p.n_eq_init[j], _EQ_GUESS_FLOOR)
    end
    _budget_clip!(n_eq, p.Ae, _plain.(be))
    _restore_feasibility!(n_eq, p.Ae, _plain.(be); maxit = RESTORE_MAXIT[])

    return _respeciate_solve!(p, n_eq, be; is_reconstruction = true)
end

"""
    _reconstruction_guess!(buf, p, be) -> buf

The independent guess: the composition the specimen was cast with, carried onto
the current element budget. It holds no active set at all, which is exactly what
recommends it where the assemblage is switching.
"""
function _reconstruction_guess!(buf, p, be)
    @inbounds for j in eachindex(buf)
        buf[j] = max(p.n_eq_init[j], _EQ_GUESS_FLOOR)
    end
    _budget_clip!(buf, p.Ae, _plain.(be))
    _restore_feasibility!(buf, p.Ae, _plain.(be); maxit = RESTORE_MAXIT[])
    return buf
end

"""
    _respeciate_solve!(p, n_eq, be) -> Bool

Solve `φ(bₑ)` from the guess `n_eq`, write the result into `p.n_full`, and record
the element-balance residual. Returns `false` if the solve threw.

Two guesses are available and the order between them is chosen, not fixed. The
**warm** one is the previous speciation; the **reconstruction** is the cast
composition carried onto the current budget. Whichever runs first, the other is
tried when the first leaves too much matter unaccounted for, and the better of
the two is kept — so the answer does not depend on the order, only the cost
does.

WHY THE ORDER IS WORTH CHOOSING. The warm start carries the previous ACTIVE SET,
and where the assemblage switches it is the wrong one: an interior-point method
started inside a set of phases that no longer exists does not cross over, it
exhausts its iterations. Measured on a six-hour hydration, 143 respeciations:

| | first solve | then |
|:--|--:|:--|
| the 98 calls the warm guess handled | 1247 ms | — |
| the 45 calls where the assemblage switched | **8087 ms**, rejected | 1467 ms from the reconstruction |

Eight seconds spent learning that a guess is wrong, then one and a half to get
the answer without it. That was 66 % of the whole integration.

An assemblage switch is not a single-step event — a phase takes several steps to
exhaust — so the previous call's outcome says which guess to try first, and
`eq_switching` carries it. This is information the problem already has, not a
tuning parameter: no threshold is introduced, and the tolerance that decides
"too much matter unaccounted for" is the same `_RETRY_ABS_TOL` as before.
"""
function _respeciate_solve!(p, n_eq, be; is_reconstruction::Bool = false)
    # `is_reconstruction` says that `n_eq` IS the reconstruction, because the
    # caller has already fallen back to it. There is then no second guess to
    # try: both would be the same vector, and solving it twice to compare it
    # with itself is the kind of waste that hides in a symmetric-looking branch.
    switching = p.eq_switching[] && !is_reconstruction
    first_guess = switching ? _reconstruction_guess!(p.n_eq_buf2, p, be) : n_eq

    ok, n_e, abs_res = _one_speciation(p, first_guess, be)
    ok || return false

    # The other guess, when the first leaves too much matter unaccounted for.
    # Decided on the ABSOLUTE balance, in moles, because that is the quantity
    # with a meaning: it is the matter the composition fails to account for. The
    # relative measure is for reporting.
    if abs_res > _RETRY_ABS_TOL && is_reconstruction
        # Already on the reconstruction and still short: nothing left to try, but
        # the regime is plainly switching, so the next call should start there.
        p.eq_switching[] = true
    elseif abs_res > _RETRY_ABS_TOL
        other = switching ? n_eq : _reconstruction_guess!(p.n_eq_buf2, p, be)
        ok2, n2, abs2 = _one_speciation(p, other, be)
        if ok2 && abs2 < abs_res
            n_e, abs_res = n2, abs2
        end
        # The regime is switching: prefer the reconstruction next time.
        p.eq_switching[] = true
    else
        # The guess that ran first was enough. Whichever it was, the regime is
        # settled: go back to the warm start, which is the cheap one when it
        # works and is right far more often than not (98 of 143 above).
        p.eq_switching[] = false
    end

    for (j, idx) in enumerate(p.idx_equilibrium)
        p.n_full[idx] = n_e[j]
    end
    p.heat_eq && _heat_sensitivity!(p, n_e, be)

    res = _row_residual(p.Ae, n_e, be)
    abs_res > p.eq_worst_abs[] && (p.eq_worst_abs[] = _plain(abs_res))

    # Separate the trajectory from the other solves, which never enter the
    # result: the first speciation, and, when the right-hand side solves the
    # partition itself (`:rhs`), its stages and rejected steps. Reporting the
    # worst over all of them alarmed about something that does not affect the
    # answer: measured when every evaluation solved the partition, on a full
    # OPC that figure was 1.13 mol while the worst over the 201 accepted steps
    # was 4.3e-4, with a median of 4.8e-9.
    p.on_accepted[] && abs_res > p.eq_worst_abs_acc[] && (p.eq_worst_abs_acc[] = _plain(abs_res))

    # Warm-starting from an INFEASIBLE point locks the error in: the next step
    # starts where this one ended, so a single bad solve poisons every solve
    # after it. Only hand over a speciation that actually satisfies the element
    # balance; otherwise leave `eq_warm` as it was and let the caller retry.
    res <= EQ_RESIDUAL_TOL && (p.eq_warm[] = true)
    res > p.eq_worst_residual[] && (p.eq_worst_residual[] = _plain(res))
    return res <= EQ_RESIDUAL_TOL
end

"""
    _one_speciation(p, guess, be) -> (ok, n_e, abs_res)

One equilibrium solve of `φ(bₑ)` from `guess`. Returns `ok = false` only if the
solve threw; a solve that returns a poor composition still returns `true`, with
its element-balance violation in moles, so the caller can compare attempts.
"""
function _one_speciation(p, guess, be)
    # A run differentiated with respect to its parameters carries dual numbers
    # in its budget, and possibly in its temperature. The partition is then the
    # certified answer, lifted to those duals by the implicit-function theorem
    # (`solve_certified`): the interior point returns no exact zero for an
    # absent phase, so no active set to differentiate at.
    #
    # The values are found first, by the route of a plain run: the interior
    # point, escalated to the certified solve where its balance is poor. The
    # certified solve on the duals then starts from that answer. Started from
    # the guess itself, it is the cold path, which a plain run never takes, and
    # it did not always converge where the plain run did: an accepted step of a
    # differentiated hydration left 3.4e9 mol unaccounted for.
    if eltype(be) <: ForwardDiff.Dual || p.T_q[] isa DynamicQuantities.AbstractQuantity{<:ForwardDiff.Dual}
        p.eq_dual === nothing && throw(
            ArgumentError(
                "differentiating a run under partial equilibrium needs the certified " *
                    "solver of the partition: OptimaSolver loaded, an aqueous phase " *
                    "and `H2O@` in it.",
            ),
        )
        T_v = _plain(ustrip(us"K", p.T_q[])) * u"K"
        P_v = _plain(ustrip(us"Pa", p.P_q[])) * u"Pa"
        ok, n_v, abs_v = _value_speciation(p, guess, _plain.(be), T_v, P_v)
        st0 = ChemicalState(p.eq_system, (ok ? n_v : guess) .* u"mol"; T = p.T_q[], P = p.P_q[])
        eq_c, cert = solve_certified(p.eq_dual, (st0,); b = be, ϵ = p.ϵ)
        n_c = [ustrip(us"mol", x) for x in eq_c.n]
        abs_c = _abs_residual(p.Ae, n_c, be)
        # The rule of a plain run: an uncertified answer that balances worse than
        # the values found is not taken. Those values are then lifted where they
        # are, as `speciated_states` lifts an instant it could not certify, so
        # that the run on dual numbers follows the plain one.
        ok && !cert.optimal && !(abs_c < abs_v) && return true, _lifted_partition(p, n_v, T_v, P_v, be), abs_v
        return true, n_c, abs_c
    end
    return _value_speciation(p, guess, be, p.T_q[], p.P_q[])
end

# The partition `n_v`, found on the values of the budget `be` at `T_v` and `P_v`,
# lifted to the duals of `be` and of the temperature `T`, by default that of the
# run.
function _lifted_partition(p, n_v, T_v, P_v, be; T = p.T_q[])
    at = ChemicalState(p.eq_system, n_v .* u"mol"; T = T, P = p.P_q[])
    eq_d, _ = _lift_equilibrium(
        p.eq_dual, at, ChemicalState(p.eq_system, n_v .* u"mol"; T = T_v, P = P_v), be; ϵ = p.ϵ,
    )
    return [ustrip(us"mol", x) for x in eq_d.n]
end

# `_one_speciation` on plain numbers, at the temperature `T` and pressure `P`.
function _value_speciation(p, guess, be, T, P)
    state_eq = ChemicalState(p.eq_system, guess .* u"mol"; T = T, P = P)

    # `p.eq_solver` is a prebuilt `EquilibriumSolver` over the partition — a
    # solver *object*, not a SciML algorithm — so `solve`, not `equilibrate`.
    local eq_result
    try
        # Polished, if at all, by the certified escalation below: the run keeps
        # its own rule, and its trajectories, for the interior-point partition.
        eq_result = _unpolished(() -> SciMLBase.solve(p.eq_solver, state_eq; ϵ = p.ϵ, b = be))
    catch err
        p.eq_failures[] += 1
        if p.eq_failures[] == 1
            @warn """re-speciation failed; the composition is left frozen for \
            this step. Later failures are counted, not reported.""" exception = err
        end
        return false, Float64[], Inf
    end

    n_e = [ustrip(us"mol", x) for x in eq_result.n]
    abs_res = _abs_residual(p.Ae, n_e, be)

    # ESCALATE, do not certify every time. The interior point is right on most
    # partitions and cheap; where it is not, it is wrong by percent, not by
    # rounding — measured, 3 % on the charge balance of a calcite solution. So
    # the certified route is spent only on the solves that need it.
    #
    # ONE START, NOT TWO, and which one is not a detail. `solve_certified` tries
    # its starts in order and stops at the first that certifies, so a second
    # start is only ever paid for when the first fails — but when it is paid
    # for, it is the whole cost of the step. `state_eq` is the raw guess, with no
    # active set and far from the answer, and a dual solve from there is the
    # COLD path: 45 s against 18 ms from a good start.
    #
    # Measured on a six-hour hydration before this changed: the escalation fired
    # on 37 % of speciations and cost 6.9 s each, which was 99.1 % of the whole
    # integration. The interior point it escalates from cost 19.7 ms.
    #
    # And the second start was redundant besides. `_respeciate_solve!` already
    # retries the entire speciation from an independent guess when the balance
    # is poor — the reconstruction, which carries no active set either and is
    # built to be feasible on the current budget, so it converges in 1.5 s where
    # `state_eq` takes seconds to tens of seconds. Two layers were paying for the
    # same idea and the inner one was the expensive way to have it.
    if abs_res > _RETRY_ABS_TOL && hasproperty(p, :eq_dual) && p.eq_dual !== nothing
        try
            eq_c, cert = solve_certified(
                p.eq_dual, (eq_result,); b = be, ϵ = p.ϵ,
            )
            n_c = [ustrip(us"mol", x) for x in eq_c.n]
            abs_c = _abs_residual(p.Ae, n_c, be)
            if cert.optimal || abs_c < abs_res
                return true, n_c, abs_c
            end
        catch
            # keep the plain answer; the caller judges it on the balance
        end
    end

    return true, n_e, abs_res
end


"""
    _row_residual(Ae, n_e, be) -> Float64

Largest element-balance violation `|Aₑnₑ − bₑ|` relative to that element's own
budget. Unlike a single global scale it cannot hide a small element behind a
large one, which is what let a 0.465 mol sulfur violation report as 1.4e-2.
"""
function _row_residual(Ae, n_e, be)
    # A measure, so on the values: it decides whether a partition is kept, and
    # carries no derivative.
    n_e, be = _plain.(n_e), _plain.(be)
    r = Ae * n_e .- be
    # An element whose total is a millionth of the largest is not tracked
    # meaningfully, and judging it against its own vanishing budget turns a
    # rounding error into an alarm: a full OPC balanced to 1.4e-10 mol at 28 days
    # was reported at 3.2e-2, and the near-empty rows of the first steps saturated
    # the measure at 1.0. Floor every row at a fraction of the system scale.
    scale = 1.0e-6 * maximum(abs, be; init = 0.0)
    worst = 0.0
    for i in eachindex(r)
        flux = zero(eltype(r))
        for j in eachindex(n_e)
            flux += abs(Ae[i, j] * n_e[j])
        end
        worst = max(worst, abs(r[i]) / max(scale, abs(be[i]), flux, 1.0e-30))
    end
    return worst
end

"""
    _abs_residual(Ae, n_e, be) -> Float64

Largest element-balance violation in moles. Reported alongside the relative
measure because it is the one a chemist can judge: 1e-10 mol is machine
precision whatever the system, and 7e-2 mol is not.
"""
_abs_residual(Ae, n_e, be) = maximum(abs, Ae * _plain.(n_e) .- _plain.(be); init = 0.0)

"""
    _budget_clip!(n_eq, Ae, be)

Clip a starting guess to the element budget: no species may exceed what the
totals `bₑ` can supply, `nⱼ ≤ minᵢ bᵢ/Aᵢⱼ` over the rows it consumes.

This does not change the feasible set — it only moves the guess into it. It is
what unblocks the OPC case: once the sulfate is spent the warm start still
carried ettringite at the aluminum budget, three times the sulfur available,
and the interior-point solve could not walk back from there.
"""
function _budget_clip!(n_eq, Ae, be)
    for j in eachindex(n_eq)
        cap = Inf
        for i in axes(Ae, 1)
            a = Ae[i, j]
            # Only rows with a POSITIVE total bound a species from above; a
            # negative total (H⁺ in a hydrating cement) is met by the hydroxides
            # and bounds nothing.
            (a > 0 && be[i] >= 0) && (cap = min(cap, be[i] / a))
        end
        isfinite(cap) && (n_eq[j] = min(n_eq[j], max(cap, _EQ_GUESS_FLOOR)))
    end
    return n_eq
end

"""
    _restore_feasibility!(n_eq, Ae, be; maxit, tol) -> n_eq

Move a starting guess into `{Aₑn = bₑ, n ≥ 0}` by alternating projection: clamp
to the box, then project onto the affine set through the 7×7 system `AₑAₑᵀ`.

The Gibbs minimization is posed with `A n = b` as a hard equality, so a guess
that violates it starts the interior-point method outside its own feasible set.
On a full OPC this was not a detail: the solve stopped on a point demanding
0.732 mol of sulfate against the 0.267 mol available — 174 % over — and no
tolerance, iteration budget, barrier setting or bound changed it, because none
of them addresses an infeasible start. A projected-gradient phase-1 proved the
set is non-empty (residual 3·10⁻¹²); this is the cheap way to land in it.

Ending on the box clamp rather than the affine step matters: the affine
projection alone leaves small negative amounts, which the barrier cannot accept.

`maxit` deserves attention. Alternating projection converges linearly, at a rate
set by the angle between the box and the affine set, and on a cement that angle
can be small: at the six-hour instant of an ordinary Portland cement — where the
iron row carries 0.013 mol across thirteen species — 200 sweeps left a residual
of 8.4e-1, 2000 left 6.7e-2, and 20 000 were needed to reach 6.7e-9.

The default is small on purpose, and measured: it runs at every re-speciation of
a run, and on a full OPC the worst in-run balance is
1.1 mol at 200 sweeps against 8.5 at 2000 and 41 at 100 000. A better guess
producing a worse answer is the back-end's own unpredictability; until that is
understood the in-run budget stays where it measures best. Note that the ranking
depends on the back-end and should be re-measured if it changes.

A replay ([`speciated_states`](@ref)) runs a handful of times and buys accuracy
instead, asking for a much larger budget.
"""
function _restore_feasibility!(n_eq, Ae, be; maxit::Int = 200, tol::Float64 = 1.0e-12)
    G = cholesky(Symmetric(Ae * Ae' + 1.0e-14I))
    sc = max(1.0, maximum(abs, be; init = 1.0))
    for _ in 1:maxit
        @. n_eq = max(n_eq, _EQ_GUESS_FLOOR)
        r = Ae * n_eq .- be
        maximum(abs, r; init = 0.0) <= tol * sc && break
        n_eq .-= Ae' * (G \ r)
    end
    @. n_eq = max(n_eq, _EQ_GUESS_FLOOR)
    return n_eq
end

# ── a right-hand side that reads the partition ──────────────────────────────

"""
    _ReadRecorder(data)

A vector that records which of its entries are read, for
[`_rates_read_speciation`](@ref).
"""
struct _ReadRecorder{T} <: AbstractVector{T}
    data::Vector{T}
    read::BitVector
end
_ReadRecorder(data::AbstractVector) = _ReadRecorder(collect(data), falses(length(data)))
Base.size(r::_ReadRecorder) = size(r.data)
Base.IndexStyle(::Type{<:_ReadRecorder}) = IndexLinear()
Base.getindex(r::_ReadRecorder, i::Int) = (r.read[i] = true; r.data[i])

"""
    _rates_read_speciation(p) -> Bool

Whether a rate law of the run reads the equilibrium partition: an amount or a log
activity of an equilibrium species, or anything the activity model derives from
them. Two probes, on the composition the run starts from: the entries of `n` and
`ln a` each law reads, recorded, and, for a run on plain numbers, the derivative
of each law along the amounts of the partition, through the activity model.

A law that reads the partition makes the right-hand side a function of the
speciation, which `build_kinetics_ode` then solves at every evaluation; one that
does not, a Parrott–Killoh or Avrami law on its own degree of reaction, leaves
it to the step's re-speciation, which is then exact.
"""
function _rates_read_speciation(p)
    p.n_be > 0 || return false
    eq = p.idx_equilibrium
    n = p.n_full
    lna = p.lna_fn(n, p)
    rn, rl = _ReadRecorder(n), _ReadRecorder(lna)
    n0 = StateView(p.n_initial_full, p.species_index)
    t0 = zero(_plain(p.T))
    for kr in p.kin_rxns
        kr.rate_fn(p.T, p.P, t0, StateView(rn, p.species_index), StateView(rl, p.species_index), n0)
    end
    any(i -> rn.read[i] || rl.read[i], eq) && return true
    # Through the activity model: a law reading the activity of a kinetic
    # aqueous species reads the partition through the ionic strength.
    eltype(n) === Float64 && p.T isa Float64 || return false
    D = ForwardDiff.Dual{typeof(ForwardDiff.Tag(_rates_read_speciation, Float64)), Float64, 1}
    eqset = Set(eq)
    nd = [D(n[i], ForwardDiff.Partials((i in eqset ? 1.0 : 0.0,))) for i in eachindex(n)]
    ld = p.lna_fn(nd, p)
    for kr in p.kin_rxns
        r = kr.rate_fn(p.T, p.P, t0, StateView(nd, p.species_index), StateView(ld, p.species_index), n0)
        r isa ForwardDiff.Dual && !iszero(ForwardDiff.partials(r)[1]) && return true
    end
    return false
end

"""
    _rhs_values(p, bv, Tv) -> Union{Nothing, Vector{Float64}}

The partition at the element amounts `bv` and the temperature `Tv`, on plain
numbers, by the certified solve warm-started from the last accepted partition,
then from the cast composition carried onto `bv`; `nothing` when neither
certifies and the better one leaves more than `_RETRY_ABS_TOL` of matter
unaccounted for. The last answer is cached, so that the evaluations of one
point by the integrator and by the step's re-speciation solve it once.
"""
function _rhs_values(p, bv::Vector{Float64}, Tv::Float64)
    c = p.rhs_cache[]
    c !== nothing && c.T == Tv && c.b == bv && return c.n
    P = _plain(ustrip(us"Pa", p.P_q[])) * u"Pa"
    state(n) = ChemicalState(p.eq_system, n .* u"mol"; T = Tv * u"K", P = P)
    warm = Float64[max(_plain(p.n_full[i]), _EQ_GUESS_FLOOR) for i in p.idx_equilibrium]
    eq, cert = _exploring_starts(() -> solve_certified(p.eq_dual, (state(warm),); b = bv, ϵ = p.ϵ))
    if eq === nothing || !cert.optimal
        guess = _reconstruction_guess!(similar(warm), p, bv)
        eq2, cert2 = _exploring_starts(() -> solve_certified(p.eq_dual, (state(guess),); b = bv, ϵ = p.ϵ))
        if eq === nothing
            eq, cert = eq2, cert2
        elseif eq2 !== nothing
            eq, cert = _keep_better(eq, cert, eq2, cert2)
        end
    end
    eq === nothing && return nothing
    n = Float64[ustrip(us"mol", x) for x in eq.n]
    abs_res = _abs_residual(p.Ae, n, bv)
    abs_res > p.eq_worst_abs[] && (p.eq_worst_abs[] = abs_res)
    cert.optimal || abs_res <= _RETRY_ABS_TOL || return nothing
    p.rhs_cache[] = (b = copy(bv), T = Tv, n = n)
    return n
end

"""
    _rhs_partition(p, be, T) -> Union{Nothing, Vector}

The partition at the element amounts `be` and the temperature `T` of a state the
right-hand side is evaluated at, carrying their dual numbers when they carry
some: the Jacobian of a stiff method, a derivative with respect to the
parameters of the run. Solved on the values ([`_rhs_values`](@ref)) and lifted by
the implicit-function theorem at the answer, so that the Jacobian of the
right-hand side holds `∂r/∂nₑ ⋅ ∂nₑ/∂bₑ` exactly.
"""
function _rhs_partition(p, be, T)
    bv = Float64[_plain(x) for x in be]
    Tv = Float64(_plain(T))
    n = _rhs_values(p, bv, Tv)
    n === nothing && return nothing
    (eltype(be) <: ForwardDiff.Dual || T isa ForwardDiff.Dual) || return n
    P = _plain(ustrip(us"Pa", p.P_q[])) * u"Pa"
    return _lifted_partition(p, n, Tv * u"K", P, collect(be); T = T * u"K")
end

# The step's re-speciation when the right-hand side solves the partition: the
# same solve at the accepted state, written where the next one starts from.
function _respeciate_rhs!(p, u)
    be = @view u[1:(p.n_be)]
    Tv = Float64(_plain(ustrip(us"K", p.T_q[])))
    n = _rhs_values(p, Float64[_plain(x) for x in be], Tv)
    n === nothing && return false
    for (j, idx) in enumerate(p.idx_equilibrium)
        p.n_full[idx] = n[j]
    end
    p.eq_warm[] = true
    abs_res = _abs_residual(p.Ae, n, be)
    p.on_accepted[] && abs_res > p.eq_worst_abs_acc[] && (p.eq_worst_abs_acc[] = _plain(abs_res))
    p.heat_eq && _heat_sensitivity!(p, n, collect(be))
    return true
end

# Above this fraction of what the system holds of an element, an amount is
# outside what the chemistry can produce.
const _FEASIBILITY_RTOL = 1.0e-8

"""
    _kinetic_state_infeasible(p, u) -> Bool

Whether the kinetic amounts of the state `u` are outside what the system can
hold: one of them negative beyond rounding, or the kinetic species together
holding more of an element than the system was given. Each amount is judged
against the most of it the element totals allow, so that a trace mineral is held
to its own scale.
"""
function _kinetic_state_infeasible(p, u)
    nk = @view u[(p.n_be + 1):(p.n_be + p.n_nk)]
    E = p.E_kin
    for j in eachindex(nk)
        v = _plain(nk[j])
        isfinite(v) || return true
        cap = minimum((p.B_el[e] / E[e, j] for e in axes(E, 1) if E[e, j] > 0); init = Inf)
        v < -_FEASIBILITY_RTOL * (isfinite(cap) ? cap : 1.0) - 1.0e-14 && return true
    end
    for e in axes(E, 1)
        tot = sum((E[e, j] * _plain(nk[j]) for j in eachindex(nk)); init = 0.0)
        tot > p.B_el[e] * (1 + _FEASIBILITY_RTOL) + 1.0e-14 && return true
    end
    return false
end

# ── build_kinetics_ode ───────────────────────────────────────────────────────

"""
    build_kinetics_ode(kp::KineticsProblem) -> Function

Build the ODE right-hand-side `f!(du, u, p, t)` implementing Leal et al. (2017).

State layout:
  - `u[1:n_be]`                    = bₑ (element amounts in equilibrium partition)
  - `u[n_be+1 : n_be+n_nk]`        = nₖ (moles of kinetic species)
  - `u[n_be+n_nk+1 : n_be+n_nk+M]` = ξ  (extents of the M kinetic reactions)
  - `u[end]`                       = T with a `SemiAdiabaticCalorimeter`,
                                     Q with an `IsothermalCalorimeter`,
                                     absent when there is no calorimeter

ODE equations (Leal 2017, Eq. 66):
  - `dnₖ/dt = νₖᵀ r`
  - `dbₑ/dt = Aₑ νₑᵀ r`
  - `dξ/dt  = r`
  - `dT/dt  = (q̇ − φ(ΔT)) / Cp_total`  (semi-adiabatic)
  - `dQ/dt  = q̇`                       (isothermal)

where `nₑ = φ(bₑ)` is the equilibrium re-speciation constraint.
"""
function build_kinetics_ode(kp::KineticsProblem)
    function f!(du, u, p, t)
        # Promote with `typeof(t)`, not `eltype(u)` alone. Rosenbrock methods
        # (`Rodas5P`, `Rodas4`, …) need a TIME gradient, which they obtain by
        # calling `f!` with a dual `t` and a plain `u`. A rate law that actually
        # depends on `t` — `waller`, and any user law with an explicit time
        # dependence — then returns a `Dual` that cannot be stored in a
        # `Vector{Float64}`, and the solve fails with "First call to automatic
        # differentiation for time gradient failed". Rate laws that ignore `t`,
        # like `parrott_killoh`, never exposed this.
        T_elt = promote_type(eltype(u), typeof(t), eltype(p.n_full))

        # ── 1. Extract state components ──────────────────────────────────
        nk = @view u[(p.n_be + 1):(p.n_be + p.n_nk)]
        T_curr = p.has_T ? u[end] : p.T

        # ── 2. Reconstruct full mole vector ──────────────────────────────
        #
        # When the right-hand side solves the partition (`:rhs`), it works on a
        # buffer of its own: `p.n_full` holds the last accepted partition, the
        # warm start of every solve, and a stage or a rejected step must not
        # move it.
        rhs = p.rhs_mode[] === :rhs && p.n_be > 0
        n_full = T_elt === eltype(p.n_full) ?
            (rhs ? copyto!(p.n_rhs, p.n_full) : p.n_full) : T_elt.(p.n_full)

        # 2a. Kinetic species from nₖ
        for (j, idx) in enumerate(p.idx_kinetic)
            n_full[idx] = max(nk[j], p.ϵ)
        end

        # 2a'. NON-kinetic species from the extents: n = n(0) + νᵀ ξ.
        #
        # Without this the non-kinetic amounts stay pinned to their initial
        # values for the whole run, because the buffer holding them is refreshed
        # only by `respeciate!`. A rate law gating on such a species — sulfate
        # available for ettringite, portlandite available for a pozzolanic
        # reaction — would then never see it move, and would run past depletion
        # while the kinetic mass balance stayed exact and the solver reported
        # success. Skipped when an equilibrium solver is present: there the
        # equilibrium partition is owned by `respeciate!`, which redistributes it
        # in a way the stoichiometry alone cannot reproduce.
        if p.n_be == 0
            ξ = @view u[(p.n_be + p.n_nk + 1):(p.n_be + p.n_nk + p.n_rxn_state)]
            # `νe` holds exactly the equilibrium-partition columns of ν, in the
            # order of `idx_equilibrium`.
            for (k, idx) in enumerate(p.idx_equilibrium)
                acc = zero(T_elt)
                for j in eachindex(ξ)
                    acc += p.νe[j, k] * ξ[j]
                end
                n_full[idx] = max(p.n_initial_full[idx] + acc, p.ϵ)
            end
        end

        # 2b. Equilibrium species: the partition.
        #
        # When no rate law reads it (`:frozen`), the partition is the one the
        # last accepted step left in `p.n_full`, re-speciated once per step by
        # `respeciate!`: the rates, hence the trajectory, do not depend on it,
        # and splitting is exact.
        #
        # When one does (`:rhs`), the partition is solved here, at the `bₑ` and
        # the temperature of the state evaluated, and its derivative with
        # respect to them is lifted into the dual numbers of a Jacobian. Frozen
        # under such a law, the rate is constant within a step, so the method
        # integrates the extent explicitly however implicit it is: on calcite
        # under `r = k(1 − Ω)`, a step longer than the second or so over which
        # `Ω` relaxes overshoots the equilibrium, `Ω` then exceeds one by orders
        # of magnitude, the rate reverses, and the run reached a reaction extent
        # of −457 mol with `Rodas5P` reporting success. A partition that cannot
        # be solved at a trial state makes the right-hand side `NaN`, which
        # rejects the step.
        if rhs
            n_e = _rhs_partition(p, (@view u[1:(p.n_be)]), T_curr)
            if n_e === nothing
                fill!(du, T_elt(NaN))
                return nothing
            end
            for (j, idx) in enumerate(p.idx_equilibrium)
                n_full[idx] = n_e[j]
            end
        end

        # ── 3. Compute log-activities ────────────────────────────────────
        lna = p.lna_fn(n_full, p)

        # ── 4. Build StateViews (O(1) named access) ─────────────────────
        n_sv = StateView(n_full, p.species_index)
        lna_sv = StateView(lna, p.species_index)
        n0_sv = StateView(p.n_initial_full, p.species_index)

        # ── 5. Evaluate kinetic rates r(T, P, t, n, lna, n₀) ────────────
        n_rxn = length(p.kin_rxns)
        rates = Vector{T_elt}(undef, n_rxn)
        for (i, kr) in enumerate(p.kin_rxns)
            rates[i] = kr.rate_fn(T_curr, p.P, t, n_sv, lna_sv, n0_sv)
            if T_elt === eltype(p.rates_buf)
                p.rates_buf[i] = rates[i]
            end
        end

        # ── 6. ODE: dnₖ/dt = νₖᵀ r (Leal Eq. 56) ───────────────────────
        fill!(du, zero(T_elt))
        du_nk = p.νk' * rates
        for j in 1:(p.n_nk)
            du[p.n_be + j] = du_nk[j]
        end

        # ── 6b. ODE: dξ/dt = r, the extents of reaction ────────────────
        for j in 1:(p.n_rxn_state)
            du[p.n_be + p.n_nk + j] = rates[j]
        end

        # ── 7. ODE: dbₑ/dt = Aₑ νₑᵀ r (Leal Eq. 65) ────────────────────
        if p.n_be > 0
            du_be = p.Ae * (p.νe' * rates)
            for j in 1:(p.n_be)
                du[j] = du_be[j]
            end
        end

        # ── 8a. ODE: dQ/dt = q̇ (isothermal) ────────────────────────────
        #
        # The heat of the kinetic reactions when they produce the hydrates. Under
        # partial equilibrium they only dissolve the anhydrous phases into ions,
        # and the heat is then `−dH/dt` over the whole composition, the partition
        # followed through its sensitivity to `bₑ` (`_equilibrium_heat_rate`).
        if p.has_Q
            du[end] = p.heat_eq && p.heat_ready[] ? _equilibrium_heat_rate(p, du, T_curr) :
                heat_rate(p.kin_rxns, rates, T_curr)
        end

        # ── 8b. ODE: dT/dt = (q̇ − φ(ΔT)) / Cp_total (semi-adiabatic) ───
        if p.has_T
            # Heat generation [W]: `−dH/dt` at fixed T under partial
            # equilibrium, `Σᵢ rᵢ (−ΔᵣH⁰ᵢ)` over the kinetic reactions otherwise.
            qdot = p.heat_eq && p.heat_ready[] ? _equilibrium_heat_rate(p, du, T_curr) :
                heat_rate(p.kin_rxns, rates, T_curr)

            # Total heat capacity: Cp_calo + Σᵢ nᵢ Cp°ᵢ(T), and the heat the
            # equilibrium partition takes up as it shifts with T.
            Cp_total = p.Cp_calo
            for (i, cp_fn) in enumerate(p.cp_fns)
                isnothing(cp_fn) && continue
                cp_i = cp_fn(; T = T_curr, unit = false)
                Cp_total = Cp_total + n_full[i] * cp_i
            end
            p.heat_eq && p.heat_ready[] &&
                (Cp_total = Cp_total + _equilibrium_shift_capacity(p, T_curr))

            ΔT = T_curr - p.T_env
            du[end] = (qdot - p.heat_loss_fn(ΔT)) / Cp_total
        end

        return nothing
    end

    return f!
end
