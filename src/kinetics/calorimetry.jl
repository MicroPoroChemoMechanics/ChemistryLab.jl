# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using DynamicQuantities
using LinearAlgebra

# ── AbstractCalorimeter ───────────────────────────────────────────────────────

"""
    abstract type AbstractCalorimeter end

Base type for calorimeter models that can be coupled to a kinetics simulation.

Concrete subtypes:
- [`IsothermalCalorimeter`](@ref): T = constant, the heat Q(t) released to the bath.
- [`SemiAdiabaticCalorimeter`](@ref): variable-T cell [Lavergne2018](@cite).

Under partial equilibrium both balance the enthalpy of the whole composition,
the state carrying the change of the enthalpy of the cell; see the theory page
*Kinetics under partial equilibrium*.
"""
abstract type AbstractCalorimeter end

# _ensure_unit is defined in utils/misc.jl (imported at module level)

# ── Heat-rate from kinetic reactions ─────────────────────────────────────────

"""
    heat_rate(kinetic_reactions, rates, T_K) -> Real

Compute the instantaneous heat generation rate [W = J/s]:

```
q̇ = Σᵢ rᵢ(t) × ΔHᵣ,ᵢ(T)
```

where `rᵢ` [mol/s] is the net rate of the i-th kinetic reaction (positive =
dissolution/forward) and `ΔHᵣ,ᵢ(T)` [J/mol] is the enthalpy of reaction.

AD-compatible: the ΔₐH⁰ callables accept `ForwardDiff.Dual` T inputs.
"""
function heat_rate(
        kinetic_reactions::AbstractVector,
        rates::AbstractVector,
        T_K;
        kwargs...,
    )
    T_val = ustrip(T_K)
    q = zero(promote_type(eltype(rates), typeof(T_val)))
    for (kr, r) in zip(kinetic_reactions, rates)
        ΔHr = _reaction_enthalpy(kr, T_val)
        q = q + r * ΔHr
    end
    return q
end

# The kinetics problem holds its reactions behind `_Heterogeneous`, which is not
# an `AbstractVector` on purpose (see its definition); unwrap it here.
heat_rate(kinetic_reactions::_Heterogeneous, rates::AbstractVector, T_K; kwargs...) =
    heat_rate(kinetic_reactions.fns, rates, T_K; kwargs...)

# ── _reaction_enthalpy dispatch hierarchy ────────────────────────────────────
#
# Priority:
#   1. KineticReaction{R, F, <:Real}   — explicit heat_per_mol, in its number
#                                        type: a heat of reaction being fitted
#                                        to a calorimetric curve carries its
#                                        derivative (`Float64` alone raised a
#                                        `MethodError` on a dual one)
#   2. KineticReaction{R, F, Nothing}  — delegate to reaction stoichiometry
#   3. AbstractReaction                — stoichiometric sum of ΔₐH⁰

function _reaction_enthalpy(kr::KineticReaction{<:Any, <:Any, <:Real}, ::Real)
    return kr.heat_per_mol
end

function _reaction_enthalpy(kr::KineticReaction{<:Any, <:Any, Nothing}, T_K::Real)
    return _reaction_enthalpy(kr.reaction, T_K)
end

function _reaction_enthalpy(reaction::AbstractReaction, T_K::Real)
    # reaction[:ΔᵣH⁰] is a SymbolicFunc or NumericFunc (T-dependent) built lazily
    # by complete_thermo_functions! from species :ΔₐH⁰ properties.
    # Thermodynamic convention: ΔᵣH⁰ < 0 for exothermic.
    # heat_rate requires "heat generated" (positive = exothermic), so we negate.
    if haskey(reaction, :ΔᵣH⁰)
        return -ustrip(reaction[:ΔᵣH⁰](; T = T_K * u"K", unit = true))
    end
    return zero(T_K)
end

# ── IsothermalCalorimeter ─────────────────────────────────────────────────────

"""
    struct IsothermalCalorimeter{T} <: AbstractCalorimeter

Isothermal calorimeter: temperature held constant at `T` [K]; the heat `Q(t)` [J]
the paste releases to the bath. When the kinetic reactions produce the
hydrates, `Q = ∫₀ᵗ q̇ dτ` with `q̇` their heat, integrated as the trailing ODE
state. Under partial equilibrium the trailing state is the change of the
enthalpy of the cell, zero, and `Q = H₀ − H` is the enthalpy the whole
composition has lost — see [`cumulative_heat`](@ref).

# Examples

```julia
cal = IsothermalCalorimeter(298.15u"K")
kp = KineticsProblem(cs, reactions, state0, tspan; calorimeter = cal)
sol = integrate(kp, ks)
t, Q    = cumulative_heat(sol, cal)
t, qdot = heat_flow(sol, cal)
```
"""
struct IsothermalCalorimeter{T} <: AbstractCalorimeter
    T::T
    IsothermalCalorimeter{T}(T_K::T) where {T} = new{T}(T_K)
end

"""
    IsothermalCalorimeter(T) -> IsothermalCalorimeter

Plain `Real` → assumed SI [K]; `Quantity` → converted to K.
"""
function IsothermalCalorimeter(T_K)
    q = _ensure_unit(us"K", T_K)
    return IsothermalCalorimeter{typeof(q)}(q)
end

n_extra_states(::IsothermalCalorimeter) = 1

function extend_u0(u0::AbstractVector, ::IsothermalCalorimeter)
    return vcat(u0, zero(eltype(u0)))
end

function extend_ode!(du, ::Any, p, n_kin::Int, cal::IsothermalCalorimeter)
    T_val = Float64(safe_ustrip(us"K", cal.T))
    qdot = heat_rate(p.kin_rxns, p.rates_buf, T_val)
    du[n_kin + 1] = qdot
    return nothing
end

# ── SemiAdiabaticCalorimeter ──────────────────────────────────────────────────

"""
    struct SemiAdiabaticCalorimeter{C, T, F} <: AbstractCalorimeter

Semi-adiabatic calorimeter following the energy balance of [Lavergne2018](@citet).
When the kinetic reactions produce the hydrates, the temperature is integrated,

```math
\\frac{dT}{dt} = \\frac{\\dot{q}(t) - \\varphi(T - T_{\\rm env})}{C_p + \\sum_i n_i C^\\circ_{p,i}(T)}
```

where:
- `q̇(t)` [W] is the instantaneous heat-generation rate,
- `φ(ΔT)` [W] is the heat-loss function (e.g. linear `L·ΔT` or quadratic `a·ΔT + b·ΔT²`),
- `Cp` [J/K] is the fixed calorimeter heat capacity,
- `Σᵢ nᵢ Cp°ᵢ(T)` is the temperature- and mole-dependent sample heat capacity
  (computed from `p.cp_fns` at every ODE step when available).

Under partial equilibrium the same balance is written on the enthalpy of the
cell: the state carries its change `ΔH`, `dΔH/dt = −φ(T − T_env)`, and the
temperature is the root of `H(φ(bₑ, T), nₖ, T) − H₀ + Cp (T − T₀) = ΔH`, the
enthalpy `H` of the paste taken over its whole composition at the partition
solved at that temperature (see the theory page *Kinetics under partial
equilibrium*). [`temperature_profile`](@ref) returns it.

# Fields

  - `Cp`: heat capacity of everything the ODE does **not** compute for itself —
    the vessel, the flask, an inert filler — in J/K, stored as a `Quantity`.
  - `heat_loss`: callable `φ(ΔT::Real) -> Real [W]`.
  - `T_env`: ambient temperature [K] (stored as `Quantity`).
  - `T0`: initial temperature [K] (stored as `Quantity`).

!!! warning "`Cp` must NOT include the sample"
    The denominator of the balance above is `Cp + Σᵢ nᵢ Cp°ᵢ(T)`, and the second
    term is recomputed from the database at every step. Adding the sample's own
    heat capacity to `Cp` therefore counts it **twice** and understates the
    temperature rise — by a factor of about 1.75 on a cement paste at w/c = 0.4,
    where the paste contributes some 2.5 kJ/K against 0.9 kJ/K for the flask.
    Three scripts shipped with this package made exactly that mistake, and the
    field description above used to invite it by saying "calorimeter + sample".
    Pass the vessel, and let the solver add the paste.

# Examples

```julia
# Linear heat loss — Newton cooling
cal = SemiAdiabaticCalorimeter(; Cp=4000.0u"J/K", T_env=293.15u"K", L=0.5u"W/K", T0=293.15u"K")

# Quadratic heat loss, the form of Lavergne et al. (2018), Eq. (23), with
# illustrative coefficients. `Cp` is the vessel alone; the paste's own heat
# capacity is added by the solver from the database.
cal = SemiAdiabaticCalorimeter(;
    Cp        = 900.0u"J/K",
    T_env     = 293.15u"K",
    heat_loss = ΔT -> 0.3*ΔT + 0.003*ΔT^2,
    T0        = 293.15u"K",
)

kp = KineticsProblem(cs, reactions, state0, tspan; calorimeter = cal)
sol = integrate(kp, ks)
t, T_vec = temperature_profile(sol, cal)
t, qdot  = heat_flow(sol, cal)
```

# References

  - [Lavergne2018](@citet).
"""
struct SemiAdiabaticCalorimeter{C, T, F} <: AbstractCalorimeter
    Cp::C       # heat capacity [J/K]
    heat_loss::F
    T_env::T    # ambient temperature [K]
    T0::T       # initial temperature [K]
end

"""
    SemiAdiabaticCalorimeter(; Cp, T_env, T0, heat_loss=nothing, L=nothing)

Keyword constructor for [`SemiAdiabaticCalorimeter`](@ref).

Exactly one of `heat_loss` or `L` must be provided:
  - `heat_loss`: callable `ΔT -> [W]` (e.g. quadratic `ΔT -> a*ΔT + b*ΔT^2`).
  - `L`: linear Newton cooling coefficient [W/K]. Sets `heat_loss = ΔT -> L * ΔT`.

All scalar fields accept plain `Real` (assumed SI) or `Quantity`:
  - `Cp` → J/K; `T_env` → K; `T0` → K; `L` → W/K.
"""
function SemiAdiabaticCalorimeter(; Cp, T_env, T0, heat_loss = nothing, L = nothing)
    Cp_q = _ensure_unit(us"J/K", Cp)
    # The two temperatures share one number type: an initial temperature being
    # fitted beside a plain ambient one makes both dual.
    Tk = (safe_ustrip(us"K", T_env), safe_ustrip(us"K", T0))
    R = promote_type(map(typeof ∘ float, Tk)...)
    T_env_q = _ensure_unit(us"K", R(Tk[1]))
    T0_q = _ensure_unit(us"K", R(Tk[2]))
    hl = if !isnothing(heat_loss)
        heat_loss
    elseif !isnothing(L)
        L_f = float(safe_ustrip(us"W/K", L))
        ΔT -> L_f * ΔT
    else
        throw(
            ArgumentError(
                "SemiAdiabaticCalorimeter requires either `heat_loss` or `L`",
            ),
        )
    end
    return SemiAdiabaticCalorimeter(Cp_q, hl, T_env_q, T0_q)
end

n_extra_states(::SemiAdiabaticCalorimeter) = 1

function extend_u0(u0::AbstractVector, cal::SemiAdiabaticCalorimeter)
    T0_f = float(safe_ustrip(us"K", cal.T0))
    return vcat(u0, T0_f)
end

"""
    extend_ode!(du, u, p, n_kin, cal::SemiAdiabaticCalorimeter)

Append `dT/dt = (q̇ − φ(ΔT)) / Cp_total(T, n)` to the ODE right-hand side.

`Cp_total = Cp + Σᵢ nᵢ Cp°ᵢ(T)` is recomputed at every ODE step
from `p.cp_fns` and `p.n_full` [Lavergne2018](@cite).
"""
function extend_ode!(du, u, p, n_kin::Int, cal::SemiAdiabaticCalorimeter)
    T_curr = u[n_kin + 1]
    Cp_f = float(safe_ustrip(us"J/K", cal.Cp))
    T_env_f = float(safe_ustrip(us"K", cal.T_env))
    # Variable total heat capacity: Cp_calorimeter + Σᵢ nᵢ Cp°ᵢ(T)
    Cp_total = Cp_f
    for (i, cp_fn) in enumerate(p.cp_fns)
        isnothing(cp_fn) && continue
        cp_i = cp_fn(; T = T_curr, unit = false)   # J/(mol·K)
        Cp_total = Cp_total + p.n_full[i] * cp_i
    end
    qdot = heat_rate(p.kin_rxns, p.rates_buf, T_curr)
    ΔT = T_curr - T_env_f
    du[n_kin + 1] = (qdot - cal.heat_loss(ΔT)) / Cp_total
    return nothing
end

# ── Result extraction ─────────────────────────────────────────────────────────
#
# Under partial equilibrium the state of a calorimeter's run carries the change
# of the enthalpy of the cell, and the temperature and the heat are functions of
# the state (`_cell_temperature`): they are computed here at the instants the
# solution saved, each warm-started from the one before. In the stoichiometric
# formulation they are read from the state, as integrated.

# Whether the run of `sol` carries a calorimeter under partial equilibrium.
_cell_run(sol) = sol.prob.p.heat_eq

# What a read-back computes on the values of a run on dual numbers carries no
# derivative: said, once per read-back, rather than handed back as if it did.
function _warn_values_only(sol, what::AbstractString)
    eltype(first(sol.u)) <: ForwardDiff.Dual || return nothing
    @warn "$what of a run on dual numbers is computed on the values of the run and carries no derivative." maxlog = 1 _id = Symbol(what)
    return nothing
end

# The temperature of the cell and the enthalpy of the paste at each instant the
# solution saved, for a run under partial equilibrium.
function _cell_points(sol)
    p = sol.prob.p
    return _with_saved_warm_start(p) do
        pts = [_cell_point(p, u) for u in sol.u]
        (first.(pts), last.(pts))
    end
end

"""
    heat_flow(sol, cal::IsothermalCalorimeter) -> (t, qdot)

Instantaneous heat-generation rate `q̇(t)` [W] at the instants the solution saved.
In the stoichiometric formulation it is the time derivative of the accumulated
heat the solution carries, read from the solver's own interpolant
(`sol(t, Val{1})`). Under partial equilibrium it is `−dH/dt` of the paste at the
temperature of the bath, the partition lifted by the implicit-function theorem
in the direction the run moves its element amounts.
"""
function heat_flow(sol, cal::IsothermalCalorimeter)
    t = sol.t
    _cell_run(sol) || return t, [sol(ti, Val{1})[end] for ti in t]
    _warn_values_only(sol, "heat_flow")
    p = sol.prob.p
    qdot = _with_saved_warm_start(p) do
        map(eachindex(t)) do i
            u = Float64[_plain(x) for x in sol.u[i]]
            T, _ = _cell_point(p, u)
            n_e = Float64[_plain(p.n_full[k]) for k in p.idx_equilibrium]
            _paste_heat_rate(p, n_e, u, Float64[_plain(x) for x in sol(t[i], Val{1})], T)
        end
    end
    return t, qdot
end

"""
    heat_flow(sol, cal::SemiAdiabaticCalorimeter) -> (t, qdot)

The heat the paste releases, `q̇ = −dH/dt` [W], at the instants the solution
saved: what warms the vessel and what leaves through its walls,
`q̇ = C_v dT/dt + φ(T − T_env)`, `C_v` the heat capacity of the vessel (`cal.Cp`).
`dT/dt` is the time derivative of the solver's interpolant in the stoichiometric
formulation, and under partial equilibrium the derivative of the root of the
cell's energy balance along the rate of the state.
"""
function heat_flow(sol, cal::SemiAdiabaticCalorimeter)
    t = sol.t
    Cp_f = safe_ustrip(us"J/K", cal.Cp)
    T_env_f = safe_ustrip(us"K", cal.T_env)
    if _cell_run(sol)
        _warn_values_only(sol, "heat_flow")
        p = sol.prob.p
        T, _ = _cell_points(sol)
        dTdt = _with_saved_warm_start(p) do
            map(eachindex(t)) do i
                u = Float64[_plain(x) for x in sol.u[i]]
                p.T_q[] = T[i] * u"K"
                _cell_temperature_rate(p, u, Float64[_plain(x) for x in sol(t[i], Val{1})])
            end
        end
        return t, Cp_f .* dTdt .+ cal.heat_loss.(T .- T_env_f)
    end
    n_kin = length(sol.u[1]) - n_extra_states(cal)
    qdot = map(eachindex(t)) do i
        T_i = sol.u[i][n_kin + 1]
        Cp_f * sol(t[i], Val{1})[n_kin + 1] + cal.heat_loss(T_i - T_env_f)
    end
    return t, qdot
end

"""
    cumulative_heat(sol, cal::IsothermalCalorimeter) -> (t, Q)

Cumulative heat `Q(t)` [J] the paste has released to the bath, at the instants
the solution saved.

In the stoichiometric formulation `Q = ∫₀ᵗ q̇ dτ` is the trailing ODE state, `q̇`
being [`heat_rate`](@ref), `Σᵢ rᵢ(−ΔᵣH⁰ᵢ)` over the kinetic reactions, which
produce the hydrates. Under **partial equilibrium** those reactions only
dissolve the anhydrous phases into ions and the hydrates are precipitated by the
Gibbs minimization, so the heat is the enthalpy the paste has lost,
`Q = H₀ − H`, `H = Σᵢ nᵢ ΔₐH⁰ᵢ(T)` over the whole composition, its partition
solved at the state of each instant and `H₀` that of the first equilibrium.
Every species then needs an enthalpy of formation, and a system where one lacks
it is refused. [`heat_release`](@ref) computes the same difference from
certified speciations.
"""
function cumulative_heat(sol, cal::IsothermalCalorimeter)
    if _cell_run(sol)
        _warn_values_only(sol, "cumulative_heat")
        _, H = _cell_points(sol)
        H0 = _plain(sol.prob.p.H0[])
        return sol.t, [H0 - H[i] + _plain(sol.u[i][end]) for i in eachindex(H)]
    end
    n_kin = length(sol.u[1]) - n_extra_states(cal)
    Q = [u[n_kin + 1] for u in sol.u]
    return sol.t, Q
end

"""
    cumulative_heat(sol, cal::SemiAdiabaticCalorimeter) -> (t, Q)

The heat the paste has released [J] at the instants the solution saved. Under
partial equilibrium it is the enthalpy it has lost, which went to warm the
vessel or out through its walls, `Q = C_v (T − T₀) − ΔH`, exactly, `ΔH` the
change of the enthalpy of the cell the state carries. In the stoichiometric
formulation, the integral of [`heat_flow`](@ref) by the rectangle rule.

Under partial equilibrium, as [`heat_flow`](@ref) and
[`temperature_profile`](@ref), it is computed on the values of the run: on a run
carrying dual numbers it carries no derivative, and a warning says so once.
"""
function cumulative_heat(sol, cal::SemiAdiabaticCalorimeter)
    if _cell_run(sol)
        _warn_values_only(sol, "cumulative_heat")
        p = sol.prob.p
        T, _ = _cell_points(sol)
        Cv, T0 = _plain(p.Cp_calo), _plain(p.T0_cell)
        return sol.t, [Cv * (T[i] - T0) - _plain(sol.u[i][end]) for i in eachindex(T)]
    end
    t, qdot = heat_flow(sol, cal)
    Q = similar(qdot)
    Q[1] = zero(eltype(qdot))
    for i in 2:lastindex(t)
        dt = t[i] - t[i - 1]
        Q[i] = Q[i - 1] + qdot[i] * dt
    end
    return t, Q
end

"""
    temperature_profile(sol, cal::SemiAdiabaticCalorimeter; times = sol.t) -> (t, T)

The temperature of the cell [K] at `times`, by default the instants the solution
saved, the state at another instant read from the solver's interpolant: the
state's last entry in the stoichiometric formulation, and under partial
equilibrium the root of the cell's energy balance at that state, each instant
warm-started from the one before.
"""
function temperature_profile(sol, ::SemiAdiabaticCalorimeter; times = sol.t)
    us = times === sol.t ? sol.u : [sol(t) for t in times]
    if _cell_run(sol)
        _warn_values_only(sol, "temperature_profile")
        p = sol.prob.p
        T = _with_saved_warm_start(p) do
            [first(_cell_point(p, u)) for u in us]
        end
        return collect(times), T
    end
    return collect(times), [u[end] for u in us]
end
