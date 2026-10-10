# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using DynamicQuantities
using LinearAlgebra
using SciMLBase

# ── KineticsProblem ───────────────────────────────────────────────────────────

"""
    struct KineticsProblem{CS, CAL, ES, AM}

Encapsulates a kinetics simulation following [Leal2017](@citet).

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

# Threads

A problem is not meant to be integrated on two threads at once: the runs would
share its equilibrium solver, whose optimizer may keep its last answer as the
next start. Independent trajectories run in parallel with one problem each.

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
    idx_kin = _kinetic_partition(system, kin_rxns)
    idx_eq = setdiff(1:n_sp, idx_kin)
    n_rxn = length(kin_rxns)

    # A species without a rate law is at equilibrium, and the minimization needs
    # its standard Gibbs energy. A glass given a formula for its rate law has
    # none: left without a law, by a mix that does not hold it, it made every
    # re-speciation of the run fail, each step keeping a frozen composition.
    if equilibrium_solver !== nothing
        unpriced = [symbol(system.species[i]) for i in idx_eq if !haskey(system.species[i], :ΔₐG⁰)]
        isempty(unpriced) || throw(
            ArgumentError(
                "$(join(unpriced, ", ")) would be at equilibrium, having no rate law, and " *
                    "$(length(unpriced) == 1 ? "has" : "have") no standard Gibbs energy for the " *
                    "minimization: give $(length(unpriced) == 1 ? "it" : "them") a rate law, or " *
                    "build the system without $(length(unpriced) == 1 ? "it" : "them")."
            )
        )
    end

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
    Ac = _constraint_matrix(_equilibrium_subsystem(system, idx_eq))
    # On plain numbers: the run restores the feasibility of its partition on
    # this matrix, and a site capacity being differentiated (which makes a
    # coupled family's entries dual) is not carried through a kinetic run.
    eltype(Ac) <: ForwardDiff.Dual && throw(
        ArgumentError(
            "KineticsProblem: the conservation matrix of the partition carries dual numbers " *
                "(a site capacity being differentiated); a kinetic run takes it on plain numbers.",
        ),
    )
    Ae = Float64.(Ac)

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

"""
    _kinetic_partition(system, kin_rxns) -> Vector{Int}

The kinetic species of a run: the species each reaction controls and, when one of
them occupies a surface site, every member of its site family.

A family's members share one site budget, so they cannot be split between the
two partitions. A family stays on the equilibrium side when no reaction controls
one of its members, which is the fast surface the partial-equilibrium
formulation assumes. When a reaction controls a member, the adsorption itself is
the slow step: the free site goes with the complexes, and the site row is then
conserved by the stoichiometry of the reactions instead of by the minimization.

Two declarations are refused, each by name:

  - a family whose budget follows its host (`SITES_FOLLOW_HOST`): its sites would
    be created or destroyed by the host's equilibrium while their occupancy is
    integrated, which nothing here accounts for;
  - an electrostatic family sharing its support with a family left at
    equilibrium: the two would see one potential, computed from a charge the
    equilibrium solve does not hold.
"""
function _kinetic_partition(system::ChemicalSystem, kin_rxns)
    idx = unique!(Int[kr.idx_mineral for kr in kin_rxns])
    families = system.site_families
    families === nothing && return idx
    promoted = falses(length(families))
    for (k, grp) in enumerate(system.site_groups)
        any(in(idx), grp) || continue
        f = families[k]
        surface_support(f).coupling === SITES_FOLLOW_HOST && throw(
            ArgumentError(
                "SiteFamily \"$(name(f))\" has a kinetic member, and its site budget " *
                    "follows host \"$(surface_support(f).host)\". A slow surface needs a " *
                    "fixed budget: declare its support with `SITES_FIXED`."
            )
        )
        promoted[k] = true
        for i in grp
            # On the kinetic side an amount moves only through the reactions, so
            # a state no reaction reaches would keep its initial amount forever.
            any(kr -> !iszero(kr.stoich[i]), kin_rxns) || throw(
                ArgumentError(
                    "SiteFamily \"$(name(f))\" is kinetic, but no reaction of the " *
                        "problem changes \"$(symbol(system.species[i]))\": its amount " *
                        "would stay at its initial value. Give each state of a slow " *
                        "family the reaction that forms it."
                )
            )
            i in idx || push!(idx, i)
        end
    end
    groups = support_group(system)
    for (k, f) in enumerate(families)
        promoted[k] || continue
        for g in groups[k]
            promoted[g] && continue
            (is_electrostatic(f.model) || is_electrostatic(families[g].model)) && throw(
                ArgumentError(
                    "SiteFamily \"$(name(f))\" is kinetic and shares its support " *
                        "\"$(surface_support(f).name)\" with \"$(name(families[g]))\", left " *
                        "at equilibrium, under an electrostatic model: both would feel one " *
                        "potential, from a charge the equilibrium solve does not hold. Put " *
                        "the two families on separate supports, or make both kinetic."
                )
            )
        end
    end
    return idx
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

# ── the number type of a run ─────────────────────────────────────────────────
#
# A run is differentiated with respect to whatever carries dual numbers: the
# amounts, temperature or pressure of the initial state, the constants of a
# calorimeter, and the parameters a rate law captures (a rate constant handed in
# as a dual by the function being differentiated). The time span is plain. The state
# of the integrator and every buffer the run writes into its result are of the
# number type that covers them all. Anything narrower either raises or, worse,
# drops the derivative: the ODE interface promotes the state only when it finds
# the duals in the parameter object, and those of a rate law live in closures it
# does not look into.

"""
    _kinetics_number_type(kp) -> Type

The number type of a run of `kp`: `Float64`, or the dual type covering the
initial state, the calorimeter and every rate law.
"""
function _kinetics_number_type(kp::KineticsProblem)
    R = promote_type(Float64, _amount_number_type(kp.initial_state))
    for kr in kp.kinetic_reactions
        R = promote_type(R, _captured_number_type(kr.rate_fn), _captured_number_type(kr.heat_per_mol))
    end
    return promote_type(R, _captured_number_type(kp.calorimeter))
end

"""
    build_u0(kp::KineticsProblem; R = _kinetics_number_type(kp)) -> Vector{R}

Build the initial ODE state vector, in the number type `R` of the run.

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

    # One trailing slot per calorimeter, and only one, `u[end]`. Under partial
    # equilibrium it is the change of the enthalpy of the cell, zero at the
    # start (`_cell_temperature`); in the stoichiometric formulation, the
    # temperature of a semi-adiabatic cell and the accumulated heat of an
    # isothermal one. `p.heat_eq`, `p.has_T` and `p.has_Q` say which.
    if kp.calorimeter isa SemiAdiabaticCalorimeter && isnothing(kp.equilibrium_solver)
        push!(u0, R(safe_ustrip(us"K", kp.calorimeter.T0)))
    elseif kp.calorimeter isa Union{SemiAdiabaticCalorimeter, IsothermalCalorimeter}
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
    # The standard Gibbs energies, when the activities of a phase of the system
    # depend on them (`CompoundEnergyModel`), and their values over RT at the
    # temperature of the run; `nothing` otherwise. See `_lna_params`. A species
    # with none, the glass of a slag, is not a member of such a phase.
    g_fns = _reads_standard_g(kp.system) ?
        _Heterogeneous([haskey(sp, :ΔₐG⁰) ? sp[:ΔₐG⁰] : nothing for sp in kp.system.species]) : nothing
    g_RT = isnothing(g_fns) ? nothing : _standard_g_over_RT(g_fns, T_K, P_Pa)

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

    # A calorimeter under partial equilibrium balances the enthalpy of the whole
    # composition, the partition included (see `_cell_temperature`). It needs
    # the enthalpy of every species, since a species without one would drop out
    # of the sum and take its heat with it.
    heat_eq = (has_T || has_Q) && n_be > 0
    heat_eq && _refuse_missing_enthalpy(kp.system, h_fns)

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
    # The temperature the cell starts at, which a semi-adiabatic run carries
    # from `T0` and an isothermal one holds.
    T0_cell = cal isa SemiAdiabaticCalorimeter ? R(safe_ustrip(us"K", cal.T0)) : T_K
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
        g_fns = g_fns,
        g_RT = g_RT,
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
        # Under partial equilibrium the last state is the change of the
        # enthalpy of the cell, and the temperature and the heat are computed
        # from it (`_cell_temperature`): `H0` is the enthalpy of the paste at the
        # first equilibrium, the reference of that balance, which
        # `_initialize_cell!` sets once the partition at `T0_cell` is known.
        heat_eq = heat_eq,
        H0 = Ref(zero(R)),
        cell_ready = Ref(false),
        T0_cell = T0_cell,
        # The last root of the cell's balance: the plain state it was solved
        # for, the temperature and the heat capacity. The stages, the Jacobian
        # and the re-speciation of one step evaluate a few nearby states, and
        # the passes of a Jacobian the same one.
        cell_root = Ref((Float64[], NaN, NaN)),
        # Equilibrium — Leal et al. (2017) §5. The re-speciation φ(bₑ) is a
        # minimization over the EQUILIBRIUM PARTITION ONLY, at frozen kinetic
        # amounts. Running it over the whole system would let the kinetic
        # minerals equilibrate instantaneously, which is exactly what a kinetic
        # description exists to prevent.
        #
        # In a reference: a `ChemicalSystem` is a vector of species of an
        # abstract type, and an ODE problem whose parameters hold one warns,
        # once per session, that they "can hurt performance", in every page
        # that integrates. The system is only read to build states.
        eq_system = Ref(eq_sys),
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
        T_q = Ref(T0_cell * u"K"),
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
        rhs_cache = Ref{Union{Nothing, _RhsCache}}(nothing),
        # Set while an accessor walks a finished run (`_with_saved_warm_start`):
        # `_rhs_values` then offers the interior point as a last start, which
        # the integration does not need and should not pay for.
        walking = Ref(false),
        n_rhs = similar(n_full),
        # Feasibility of a kinetic state: the element content of the kinetic
        # species and the element totals of the system.
        E_kin = E_all[:, kp.idx_kinetic],
        B_el = B_el,
        # The species whose amount is an amount of substance rather than the
        # total of a component: all but the solutes. Each is judged against the
        # most of it the element totals allow.
        is_solute = [i in kp.system.idx_solutes for i in eachindex(kp.system.species)],
        amount_cap = [
            minimum((B_el[e] / E_all[e, i] for e in axes(E_all, 1) if E_all[e, i] > 0); init = Inf)
                for i in eachindex(kp.system.species)
        ],
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
    # them. A family is fast or slow as a whole (`_kinetic_partition` moves every
    # member of a family whose member a reaction controls), so a split is a
    # declaration error, and it is worth saying which family and why.
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

# The floor of the warm start of a certified solve in the right-hand side: below
# what any amount of a partition means, so that the traces keep their potentials.
const _RHS_GUESS_FLOOR = 1.0e-16

# The last partition `_rhs_values` solved, keyed by the plain element amounts and
# temperature it was solved at.
const _RhsCache = NamedTuple{(:b, :T, :n), Tuple{Vector{Float64}, Float64, Vector{Float64}}}

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
stoichiometry to write down. This is [Lavergne2018; Eqs. 17–21](@citet).
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
    # The slot of each kinetic species in `u`, zero for the others.
    slot = zeros(Int, length(p.h_fns))
    for (j, i) in enumerate(p.idx_kinetic)
        slot[i] = j
    end
    H = zero(promote_type(eltype(u), typeof(T), eltype(p.n_full)))
    @inbounds for (i, h_fn) in enumerate(p.h_fns)
        isnothing(h_fn) && continue
        nᵢ = iszero(slot[i]) ? p.n_full[i] : max(u[p.n_be + slot[i]], p.ϵ)
        H += nᵢ * h_fn(; T = T, unit = false)
    end
    return H
end

# ── the energy balance of a calorimeter under partial equilibrium ───────────
#
# Under partial equilibrium the kinetic reactions only dissolve the anhydrous
# phases into ions, and the hydrates are precipitated by the minimization: the
# heat of the kinetic reactions leaves the precipitation out (on an ordinary
# Portland cement, a semi-adiabatic rise of 207 K). The heat is the fall of the
# enthalpy of the WHOLE composition, `H = Σᵢ nᵢ ΔₐH⁰ᵢ(T)`, whose partition is the
# one the minimization gives at the element amounts and the temperature of the
# state, `nₑ = φ(bₑ, T)`.
#
# The right-hand side is a function of the state alone, as Leal et al. (2015)
# write it for the amounts, and the energy balance is written the same way. The
# state carries the change `ΔH` of the enthalpy of the cell, which only the
# losses through its walls move, and the temperature is the root of
#
#     R(T) = H(φ(bₑ, T), nₖ, T) − H₀ + C_v (T − T₀) − ΔH = 0 ,
#
# solved with the partition at every evaluation. The temperature is not
# integrated: under partial equilibrium the partition moves with it, so its rate
# would hold the derivative of the partition, and a stiff method's Jacobian the
# derivative of that. Written on the enthalpy, every derivative the Jacobian
# holds is a first derivative of `φ`, which the implicit-function theorem gives
# exactly. An isothermal cell is the case where nothing leaves: `ΔH` stays zero,
# and the heat delivered to the bath is `Q = H₀ − H(φ(bₑ, T), nₖ, T)`. The
# equations are those of the theory page *Kinetics under partial equilibrium*.

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
    _paste_enthalpy(p, n_e, nk, T) -> Real

`H = Σᵢ nᵢ ΔₐH⁰ᵢ(T)` [J] of the composition whose equilibrium partition is
`n_e`, in the order of `p.idx_equilibrium`, and whose kinetic amounts are `nk`:
what [`enthalpy`](@ref) gives a state of that composition. Generic in the number
types of all three.
"""
function _paste_enthalpy(p, n_e, nk, T)
    H = zero(promote_type(eltype(n_e), eltype(nk), typeof(T)))
    for (j, i) in enumerate(p.idx_equilibrium)
        H += n_e[j] * p.h_fns[i](; T = T, unit = false)
    end
    for (j, i) in enumerate(p.idx_kinetic)
        H += max(nk[j], p.ϵ) * p.h_fns[i](; T = T, unit = false)
    end
    return H
end

# The temperature of a cell is solved to this step, in kelvin, and in at most
# this many Newton iterations.
const _CELL_T_TOL = 1.0e-9
const _CELL_T_MAXIT = 30

"""
    _cell_residual(p, bv, nkv, n_e, Tv, ΔHv) -> (R, C)

The residual `R(T) = H(φ(bₑ, T), nₖ, T) − H₀ + C_v (T − T₀) − ΔH` of the energy
balance of the cell at the temperature `Tv`, and its derivative along the
equilibrium,

    C = R'(T) = C_v + Σᵢ nᵢ c°ₚ,ᵢ(T) + Σₑ ΔₐH⁰ₑ(T) ∂φₑ/∂T ,

the heat capacity of the cell at equilibrium: the partition `n_e`, solved at
`Tv`, lifted in temperature by the implicit-function theorem at the answer. On
plain numbers.
"""
function _cell_residual(p, bv, nkv, n_e, Tv, ΔHv)
    V = typeof(Tv)
    D = ForwardDiff.Dual{typeof(ForwardDiff.Tag(_cell_residual, V)), V, 1}
    Td = D(Tv, ForwardDiff.Partials((1.0,)))
    P = _plain(ustrip(us"Pa", p.P_q[])) * u"Pa"
    n_d = _lifted_partition(p, n_e, Tv * u"K", P, bv; T = Td * u"K")
    R = _paste_enthalpy(p, n_d, nkv, Td) - _plain(p.H0[]) + _plain(p.Cp_calo) * (Td - _plain(p.T0_cell)) - ΔHv
    return _plain(ForwardDiff.value(R)), _plain(ForwardDiff.partials(R)[1])
end

"""
    _cell_temperature_value(p, bv, nkv, ΔHv) -> Union{Nothing, Tuple{Float64, Float64}}

The temperature of the cell at the plain state `(bv, nkv, ΔHv)` and its heat
capacity at equilibrium there ([`_cell_residual`](@ref)): Newton on `R`, from the
last root found (`p.cell_root`), the nearest state evaluated, or else from the
last accepted temperature, until its step falls below `_CELL_T_TOL`. The
temperature returned is the last one the partition was solved at, which the
evaluation that follows finds in the cache of `_rhs_values`. `nothing` when the
partition cannot be solved, when the capacity is not positive (the balance would
then have no unique root), or after `_CELL_T_MAXIT` iterations.
"""
function _cell_temperature_value(p, bv, nkv, ΔHv)
    last = p.cell_root[][2]
    Tv = isfinite(last) ? last : Float64(_plain(ustrip(us"K", p.T_q[])))
    for _ in 1:_CELL_T_MAXIT
        n_e = _rhs_values(p, bv, Tv)
        n_e === nothing && break
        R, C = _cell_residual(p, bv, nkv, n_e, Tv, ΔHv)
        (isfinite(R) && isfinite(C) && C > 0) || break
        δ = R / C
        abs(δ) <= _CELL_T_TOL && return Tv, C
        Tv -= δ
    end
    # The partition could not be solved, the capacity was not positive, or the
    # iterations ran out.
    return nothing
end

"""
    _cell_temperature(p, u) -> Union{Nothing, Real}

The temperature of a semi-adiabatic cell under partial equilibrium at the state
`u`, the root of its energy balance ([`_cell_temperature_value`](@ref)). On dual
numbers the root is lifted by the implicit-function theorem: the residual is
evaluated on the state's own numbers at the root `T*`, its partition lifted in
`bₑ`, and

    T = T* − (R(u, T*) − R*) / C ,

which carries `∂T/∂u = −(∂R/∂u)/C` exactly and keeps the value `T*` the
partition was solved at. `nothing` when the root cannot be found.
"""
function _cell_temperature(p, u)
    nb, nn = p.n_be, p.n_nk
    uv = Float64[_plain(x) for x in u]
    # The passes of a Jacobian evaluate the same state as the value pass: its
    # root is kept rather than solved again.
    cached = p.cell_root[]
    if cached[1] == uv
        Tv, C = cached[2], cached[3]
    else
        root = _cell_temperature_value(p, uv[1:nb], uv[(nb + 1):(nb + nn)], uv[end])
        root === nothing && return nothing
        Tv, C = root
        p.cell_root[] = (uv, Tv, C)
    end
    eltype(u) === Float64 && return Tv
    n_e = _rhs_partition(p, (@view u[1:nb]), Tv)
    n_e === nothing && return nothing
    R = _paste_enthalpy(p, n_e, (@view u[(nb + 1):(nb + nn)]), Tv) - p.H0[] +
        p.Cp_calo * (Tv - p.T0_cell) - u[end]
    return Tv - (R - _plain(R)) / C
end

"""
    _initialize_cell!(p, u0)

Set the reference of the energy balance of a calorimeter under partial
equilibrium: `H₀`, the enthalpy of the paste at the partition of the initial
element amounts at the initial temperature, from the solve the right-hand side
uses, so that the balance holds at the start with `ΔH = 0` and the first
equilibrium releases nothing. That partition is the warm start of the run.
"""
function _initialize_cell!(p, u0)
    p.heat_eq || return nothing
    p.eq_dual === nothing && throw(
        ArgumentError(
            "a calorimeter under partial equilibrium solves the partition with the certified " *
                "solver at every evaluation: it needs OptimaSolver loaded and an aqueous phase " *
                "with `H2O@` in the partition.",
        ),
    )
    nb, nn = p.n_be, p.n_nk
    T0 = p.T0_cell
    p.T_q[] = _plain(T0) * u"K"
    n_e = _rhs_partition(p, (@view u0[1:nb]), T0)
    n_e === nothing && throw(
        ErrorException(
            "the partition of the initial state cannot be solved at $(_plain(T0)) K: " *
                "the calorimeter has no reference enthalpy.",
        ),
    )
    for (j, i) in enumerate(p.idx_equilibrium)
        p.n_full[i] = n_e[j]
    end
    p.eq_warm[] = true
    p.H0[] = _paste_enthalpy(p, n_e, (@view u0[(nb + 1):(nb + nn)]), T0)
    p.cell_ready[] = true
    return nothing
end

# Run `f`, then restore the warm start of the partition (`p.n_full`, `p.T_q`, the
# cache of `_rhs_values`): an accessor walks a finished run from its start, each
# instant warm-started from the one before, and leaves the run's state as it
# found it.
function _with_saved_warm_start(f, p)
    n, T, c, r, wk = copy(p.n_full), p.T_q[], p.rhs_cache[], p.cell_root[], p.walking[]
    p.walking[] = true
    try
        return f()
    finally
        copyto!(p.n_full, n)
        p.T_q[] = T
        p.rhs_cache[] = c
        p.cell_root[] = r
        p.walking[] = wk
    end
end

"""
    _cell_point(p, u) -> (T, H)

The temperature of the cell and the enthalpy of the paste at the state `u` of a
run under partial equilibrium, on plain numbers, the partition left in `p` as the
warm start of the next instant.
"""
function _cell_point(p, u)
    nb, nn = p.n_be, p.n_nk
    uv = Float64[_plain(x) for x in u]
    T = p.has_T ? _cell_temperature(p, uv) : Float64(_plain(p.T))
    T === nothing && throw(ErrorException("the temperature of the cell cannot be solved at this state."))
    n_e = _rhs_values(p, uv[1:nb], T)
    n_e === nothing && throw(ErrorException("the partition cannot be solved at this state."))
    for (j, i) in enumerate(p.idx_equilibrium)
        p.n_full[i] = n_e[j]
    end
    p.T_q[] = T * u"K"
    return T, Float64(_plain(_paste_enthalpy(p, n_e, (@view uv[(nb + 1):(nb + nn)]), T)))
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
sub-problem [Leal2017; Eq. 54](@cite). `solve` conserves `A·n`, so what has to
be handed to it is a composition whose element totals are exactly `bₑ` — here
the previous speciation, projected onto `bₑ` through the pseudo-inverse of
`Aₑ`. Handing over `p.n_full` unchanged, as an earlier version did, discards
`bₑ` entirely and leaves the element balance to drift.
"""
function respeciate!(p, u)
    p.n_be > 0 || return false

    # The partition of a semi-adiabatic cell is an equilibrium at the cell's
    # temperature, not at the initial one: carried by the state when the
    # formulation is stoichiometric, solved from the cell's enthalpy under
    # partial equilibrium (`_cell_temperature`), once its reference is set.
    if p.has_T
        T = p.heat_eq ? (p.cell_ready[] ? _cell_temperature(p, u) : nothing) : u[end]
        T === nothing || (p.T_q[] = _plain(T) * u"K")
    end

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
        st0 = ChemicalState(p.eq_system[], (ok ? n_v : guess) .* u"mol"; T = p.T_q[], P = p.P_q[])
        eq_c, cert = solve_certified(p.eq_dual, (st0,); b = be, ϵ = p.ϵ, report = false)
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
    # Nested dual numbers, a parameter's derivative under the Jacobian of a stiff
    # method, are lifted one level at a time, as `solve` lifts them: the
    # partition at the inner duals first, then the outer level from it. Handed
    # to the tangent at once, the inner duals met a conversion to `Float64`, so
    # no rate law reading the partition could be differentiated with respect to
    # its parameters under `Rodas5P`.
    Tn = ustrip(us"K", T)
    D = promote_type(typeof(Tn), eltype(be))
    if D <: ForwardDiff.Dual && ForwardDiff.valtype(D) <: ForwardDiff.Dual
        Tg = ForwardDiff.tagtype(D)
        T_in = _strip_tag(Tn, Tg)
        inner = with(_STRIP_TAGS => (_STRIP_TAGS[]..., Tg)) do
            _lifted_partition(p, n_v, T_v, P_v, _strip_tag(collect(be), Tg); T = T_in * u"K")
        end
        at = ChemicalState(p.eq_system[], n_v .* u"mol"; T = T, P = p.P_q[])
        eq_in = ChemicalState(p.eq_system[], inner .* u"mol"; T = T_in * u"K", P = P_v)
        eq_d, _ = _lift_equilibrium(p.eq_dual, at, eq_in, be; ϵ = p.ϵ, strip_tag = Tg)
        return [ustrip(us"mol", x) for x in eq_d.n]
    end
    at = ChemicalState(p.eq_system[], n_v .* u"mol"; T = T, P = p.P_q[])
    eq_d, _ = _lift_equilibrium(
        p.eq_dual, at, ChemicalState(p.eq_system[], n_v .* u"mol"; T = T_v, P = P_v), be; ϵ = p.ϵ,
    )
    return [ustrip(us"mol", x) for x in eq_d.n]
end

# `_one_speciation` on plain numbers, at the temperature `T` and pressure `P`.
function _value_speciation(p, guess, be, T, P)
    state_eq = ChemicalState(p.eq_system[], guess .* u"mol"; T = T, P = P)

    # `p.eq_solver` is a prebuilt `EquilibriumSolver` over the partition — a
    # solver *object*, not a SciML algorithm — so `solve`, not `equilibrate`.
    local eq_result
    try
        # Polished, if at all, by the certified escalation below: the run keeps
        # its own rule, and its trajectories, for the interior-point partition.
        eq_result = SciMLBase.solve(p.eq_solver, state_eq; ϵ = p.ϵ, b = be, polish = false)
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
                p.eq_dual, (eq_result,); b = be, ϵ = p.ϵ, report = false,
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
producing a worse answer is a behavior of the back end not yet understood; until
it is understood the in-run budget stays where it measures best. Note that the
ranking depends on the back-end and should be re-measured if it changes.

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

# ── the parameters of the activity model ────────────────────────────────────

"""
    _lna_params(p, T) -> NamedTuple

The parameters the activity model of a run is called with at the temperature
`T` of the cell: `p`, at that temperature when the cell has one of its own, and
with the standard Gibbs energies over RT at that temperature when a phase of the
system reads them, as the activities of a phase under the compound energy
formalism do ([`CompoundEnergyModel`](@ref)). Without them that model refuses to
compute, since no default would be right.
"""
function _lna_params(p, T)
    q = p.has_T ? merge(p, (T = T,)) : p
    isnothing(p.g_fns) && return q
    g = p.has_T ? _standard_g_over_RT(p.g_fns, T, p.P) : p.g_RT
    return merge(q, (ΔₐG⁰overRT = g,))
end

# `NaN` for a species without a standard Gibbs energy, so that nothing reads it
# unnoticed.
_standard_g_over_RT(g_fns, T, P) =
    _promoted([isnothing(g) ? NaN : g(; T = T, P = P, unit = false) / (R_GAS * T) for g in g_fns])

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
    lna = p.lna_fn(n, _lna_params(p, p.T))
    rn, rl = _ReadRecorder(n), _ReadRecorder(lna)
    n0 = StateView(p.n_initial_full, p.species_index)
    t0 = zero(_plain(p.T))
    for kr in p.kin_rxns
        kr.rate_fn(p.T, p.P, t0, StateView(rn, p.species_index), StateView(rl, p.species_index), n0)
    end
    any(i -> rn.read[i] || rl.read[i], eq) && return true
    # Through the activity model: a law reading the activity of a kinetic
    # aqueous species reads the partition through the ionic strength. Seeded on
    # the values of the amounts, so that a run on dual numbers (a rate constant
    # being differentiated) decides as the run on its values does: skipped
    # there, the probe chose the frozen route where the plain run solved the
    # partition in the right-hand side, and the derivative was that of another
    # trajectory. The seed's dual may end up outside or inside the run's own;
    # `_reads_seed` looks for it in both.
    D = ForwardDiff.Dual{typeof(ForwardDiff.Tag(_rates_read_speciation, Float64)), Float64, 1}
    eqset = Set(eq)
    nd = [D(_plain(n[i]), ForwardDiff.Partials((i in eqset ? 1.0 : 0.0,))) for i in eachindex(n)]
    ld = p.lna_fn(nd, _lna_params(p, p.T))
    for kr in p.kin_rxns
        r = kr.rate_fn(p.T, p.P, t0, StateView(nd, p.species_index), StateView(ld, p.species_index), n0)
        _reads_seed(D, r) && return true
    end
    return false
end

# Whether `r` carries a nonzero derivative along the seed of the dual type `D`,
# however it is nested among the duals of other tags.
_reads_seed(::Type, r) = false
function _reads_seed(::Type{D}, r::ForwardDiff.Dual{S}) where {D <: ForwardDiff.Dual, S}
    S === ForwardDiff.tagtype(D) && return any(c -> !iszero(_plain(c)), ForwardDiff.partials(r))
    return _reads_seed(D, ForwardDiff.value(r)) || any(c -> _reads_seed(D, c), ForwardDiff.partials(r))
end

"""
    _rhs_values(p, bv, Tv) -> Union{Nothing, Vector{Float64}}

The partition at the element amounts `bv` and the temperature `Tv`, on plain
numbers, by the certified solve warm-started from the last accepted partition,
then from the cast composition carried onto `bv`, and, for an accessor walking a
finished run, from the answer of the interior point; `nothing` when none
certifies and the best leaves more than
`_RETRY_ABS_TOL` of matter unaccounted for. The last answer is cached, so that
the evaluations of one point by the integrator and by the step's re-speciation
solve it once.
"""
function _rhs_values(p, bv::Vector{Float64}, Tv::Float64)
    c = p.rhs_cache[]
    c !== nothing && c.T == Tv && c.b == bv && return c.n
    P = _plain(ustrip(us"Pa", p.P_q[])) * u"Pa"
    state(n) = ChemicalState(p.eq_system[], n .* u"mol"; T = Tv * u"K", P = P)
    # The last accepted partition as it is, floored only where it is zero. The
    # floor of the interior point, `_EQ_GUESS_FLOOR`, lifts the traces of a
    # certified partition to 1e-10 mol and their potentials with them; measured
    # on the C100 mortar of Lavergne et al. (2018) at 1.6 h, the dual Newton
    # started there failed to certify where the partition itself, floored at
    # 1e-16, certified at once. The floored start is the next one tried.
    warm = Float64[max(_plain(p.n_full[i]), _RHS_GUESS_FLOOR) for i in p.idx_equilibrium]
    solve_from(n) = _exploring_starts(() -> solve_certified(p.eq_dual, (state(n),); b = bv, ϵ = p.ϵ, report = false))
    eq, cert = solve_from(warm)
    eq, cert = _or_next(eq, cert, () -> solve_from(max.(warm, _EQ_GUESS_FLOOR)))
    # The reconstruction: feasible on the budget, with no active set.
    eq, cert = _or_next(eq, cert, () -> solve_from(_reconstruction_guess!(similar(warm), p, bv)))
    # Where the assemblage switches, the starts above hold the assemblage of the
    # last accepted step, or none, and the dual Newton stalls short of the
    # certificate from either. Measured on a paste of cement c13 of Lavergne et
    # al. (2018) at 23.9 °C and 4.2 h, where hydrogarnet gives way to
    # monosulfate: a KKT error of 1.5e-4 from the warm start, 1.3e-4 from the
    # reconstruction. The interior point crosses to the new assemblage, and the
    # certified solve from its answer proves it at once (7e-15).
    #
    # Offered to an accessor walking a finished run only (`p.walking`). The run
    # passes these states by its step re-speciation, which falls back on the
    # interior point, and a failure inside it only rejects a step; tried at every
    # failure of the integration, this start doubled the cost of a semi-adiabatic
    # run, whose cell temperature is a root over partitions solved at trial
    # temperatures (316 s against 160 s on the C100 mortar of Lavergne et al.
    # 2018). The polish takes the interior point's answer as a start only: its
    # return code is not counted among the solves used anyway, and its answer is
    # kept only when the polish proves it. Kept uncertified, it was accepted on
    # its element balance, which the interior point always meets, and the
    # temperature of that cell then failed where the run had passed.
    if p.walking[] && (eq === nothing || !cert.optimal) && p.eq_solver !== nothing
        cref = Ref{Any}(nothing)
        eq2 = try
            _exploring_starts() do
                SciMLBase.solve(p.eq_solver, state(warm); b = bv, ϵ = p.ϵ, polish = true, certificate = cref)
            end
        catch
            nothing
        end
        cert2 = cref[]
        if eq2 !== nothing && cert2 !== nothing && cert2.optimal
            eq, cert = eq2, cert2
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

# The answer in hand when it is certified; otherwise the better of it and the
# answer `attempt()` gives, which is solved only then: the next start of a
# cascade (`_keep_better`).
function _or_next(eq, cert, attempt)
    (eq !== nothing && cert.optimal) && return eq, cert
    eq2, cert2 = attempt()
    eq === nothing && return eq2, cert2
    eq2 === nothing && return eq, cert
    return _keep_better(eq, cert, eq2, cert2)
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
to its own scale. Without an equilibrium partition, also whether the extents take
a species they consume below zero (`_exhausted_coreactant`, to `rtol`).
"""
_kinetic_state_infeasible(p, u; rtol = _FEASIBILITY_RTOL) =
    _kinetic_amounts_infeasible(p, u) || _exhausted_coreactant(p, u; rtol) !== nothing

# The kinetic amounts alone: negative, or holding more of an element than the
# system was given.
function _kinetic_amounts_infeasible(p, u)
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

"""
    _exhausted_coreactant(p, u; rtol = _FEASIBILITY_RTOL) -> Union{Nothing, Tuple{Int, Float64}}

The first species the extents of the state `u` take below zero, as its index in
the system and its amount, when the species other than the kinetic ones follow
the stoichiometry, `n = n(0) + νᵀξ`, that is without an equilibrium partition;
`nothing` otherwise. Each amount is judged against the most of it the element
totals allow, times `rtol`: the check of a trajectory passes the relative
tolerance of its integration, below which the extents that make the amount are
not known. A phase consumed in the first seconds slightly ahead of the reaction
that forms it, by a few parts in 10⁸ of its calcium, is then not reported.

Only the species whose amount is an amount of substance are judged: a phase, a
gas, the solvent. Without a partition a solute stands for the total of its
component, whose sign carries the acidity: a clinker phase dissolving into the
primaries consumes `H+` below zero, which is hydroxide, not a shortfall.

A rate law that reads only its own phase does not see a co-reactant run out:
under a Parrott–Killoh rate, `C3A + 3 Gp + 26 H2O → ettringite` goes on after
the gypsum is gone. The right-hand side floors the amount it reads at zero, so
the run went on, reported a success, and created the sulfate the extent demanded
beyond what the system held.
"""
function _exhausted_coreactant(p, u; rtol = _FEASIBILITY_RTOL)
    p.n_be == 0 || return nothing
    ξ = @view u[(p.n_nk + 1):(p.n_nk + p.n_rxn_state)]
    for (k, idx) in enumerate(p.idx_equilibrium)
        p.is_solute[idx] && continue
        v = Float64(_plain(p.n_initial_full[idx]))
        for j in eachindex(ξ)
            v += p.νe[j, k] * _plain(ξ[j])
        end
        cap = p.amount_cap[idx]
        v < -rtol * (isfinite(cap) ? cap : 1.0) - 1.0e-14 && return (idx, v)
    end
    return nothing
end

# ── build_kinetics_ode ───────────────────────────────────────────────────────

"""
    build_kinetics_ode(kp::KineticsProblem) -> Function

Build the ODE right-hand-side `f!(du, u, p, t)` implementing [Leal2017](@citet).

State layout:
  - `u[1:n_be]`                    = bₑ (element amounts in equilibrium partition)
  - `u[n_be+1 : n_be+n_nk]`        = nₖ (moles of kinetic species)
  - `u[n_be+n_nk+1 : n_be+n_nk+M]` = ξ  (extents of the M kinetic reactions)
  - `u[end]`                       = with a calorimeter: under partial
                                     equilibrium, the change `ΔH` of the
                                     enthalpy of the cell; otherwise T with a
                                     `SemiAdiabaticCalorimeter`, Q with an
                                     `IsothermalCalorimeter`; absent without one

ODE equations [Leal2017; Eq. 66](@cite):
  - `dnₖ/dt = νₖᵀ r`
  - `dbₑ/dt = Aₑ νₑᵀ r`
  - `dξ/dt  = r`
  - under partial equilibrium, `dΔH/dt = −φ(T − T_env)` (semi-adiabatic) or `0`
    (isothermal), the temperature the root of the cell's energy balance
    (`_cell_temperature`);
  - otherwise `dT/dt = (q̇ − φ(ΔT)) / Cp_total` (semi-adiabatic) or `dQ/dt = q̇`
    (isothermal), `q̇` the heat of the kinetic reactions,

where `nₑ = φ(bₑ, T)` is the equilibrium partition, solved at the element
amounts and the temperature of the state.
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
        # The temperature: the bath's, the state's own in a stoichiometric
        # semi-adiabatic cell, and under partial equilibrium the root of the
        # cell's energy balance, solved with the partition and lifted into the
        # dual numbers of a Jacobian.
        T_curr = if !p.has_T
            p.T
        elseif p.heat_eq
            _cell_temperature(p, u)
        else
            u[end]
        end
        if T_curr === nothing
            fill!(du, T_elt(NaN))
            return nothing
        end

        # ── 2. Reconstruct full mole vector ──────────────────────────────
        #
        # When the right-hand side solves the partition (`:rhs`), it works on a
        # buffer of its own: `p.n_full` holds the last accepted partition, the
        # warm start of every solve, and a stage or a rejected step must not
        # move it.
        rhs = p.n_be > 0 && (p.rhs_mode[] === :rhs || (p.has_T && p.heat_eq))
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
        #
        # At the temperature of the cell: an activity model reads `p.T`, which
        # is the initial temperature.
        lna = p.lna_fn(n_full, _lna_params(p, T_curr))

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

        # ── 8. The calorimeter ──────────────────────────────────────────
        #
        # Under partial equilibrium the state carries the change of the
        # enthalpy of the cell, which only the losses move: zero in an
        # isothermal cell, whose heat is the enthalpy the paste loses
        # (`cumulative_heat`), `−φ(T − T_env)` in a semi-adiabatic one, whose
        # temperature is the root of the balance solved above.
        if p.heat_eq
            p.has_T && (du[end] = -p.heat_loss_fn(T_curr - p.T_env))
        elseif p.has_Q
            # The heat of the kinetic reactions, which produce the hydrates.
            du[end] = heat_rate(p.kin_rxns, rates, T_curr)
        elseif p.has_T
            # dT/dt = (q̇ − φ(ΔT)) / (C_v + Σᵢ nᵢ c°ₚ,ᵢ(T)), the heat of the
            # kinetic reactions over the heat capacity of the vessel and of the
            # composition the extents give.
            qdot = heat_rate(p.kin_rxns, rates, T_curr)
            Cp_total = p.Cp_calo
            for (i, cp_fn) in enumerate(p.cp_fns)
                isnothing(cp_fn) && continue
                Cp_total = Cp_total + n_full[i] * cp_fn(; T = T_curr, unit = false)
            end
            du[end] = (qdot - p.heat_loss_fn(T_curr - p.T_env)) / Cp_total
        end

        return nothing
    end

    return f!
end
