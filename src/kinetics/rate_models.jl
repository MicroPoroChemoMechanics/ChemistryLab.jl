# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)
# Portions of this file (the Arrhenius rate constant and the saturation-ratio
# formulation used by the transition-state theory rate model) are Julia ports
# adapted from the Reaktoro C++ library (https://github.com/reaktoro/reaktoro),
# Copyright © 2014-2024 Allan Leal, distributed under the LGPL-2.1-or-later.

using DynamicQuantities
using OrderedCollections

# ── StateView ─────────────────────────────────────────────────────────────────

"""
    StateView{T, I <: AbstractDict}

Thin wrapper giving O(1) named access to a species data vector.

```julia
sv["C3S"] === sv.data[sv.index["C3S"]]
```

The `index` dict is built once at [`KineticsProblem`](@ref) construction;
`data` is a plain vector (mutated in-place or re-wrapped each ODE step) — no
dict allocation in the hot path.

# Examples

```jldoctest
julia> idx = Dict("Ca++" => 1, "C3S" => 2);

julia> sv = StateView([0.5, 1.0], idx);

julia> sv["C3S"]
1.0

julia> haskey(sv, "Ca++")
true
```
"""
struct StateView{T, I <: AbstractDict}
    data::AbstractVector{T}
    index::I
end

Base.getindex(sv::StateView, name::AbstractString) = sv.data[sv.index[name]]
Base.haskey(sv::StateView, name::AbstractString) = haskey(sv.index, name)

# ── KineticFunc ───────────────────────────────────────────────────────────────

"""
    KineticFunc{F, R <: NamedTuple, Q}

Compiled kinetic rate function, analogous to [`NumericFunc`](@ref) for
thermodynamics.

Calling convention (positional, **not** keyword):

```julia
kf(T, P, t, n, lna, n_initial) -> Real   # [mol/s]
```

where:
  - `T` [K], `P` [Pa]: temperature and pressure (plain `Real` or `ForwardDiff.Dual`).
  - `t` [s]: current time.
  - `n::StateView`: moles of all species (named access: `n["C3S"]`).
  - `lna::StateView`: log-activities of all species.
  - `n_initial::StateView`: initial moles (always `Float64`).
  - return: net dissolution rate [mol/s], positive = dissolution.

AD-compatible when the compiled closure is AD-compatible.

# Examples

```jldoctest
julia> idx = Dict("C3S" => 1);

julia> n_sv  = StateView([1.0], idx);

julia> lna_sv = StateView([0.0], idx);

julia> pk = parrott_killoh_avrami(PK84_PARAMS_C3S, "C3S");

julia> pk(293.15, 1e5, 0.0, n_sv, lna_sv, n_sv) > 0
true
```
"""
struct KineticFunc{F, R <: NamedTuple, Q} <: Function
    compiled::F
    vars::NTuple{6, Symbol}   # (T, P, t, n, lna, n_initial) — positional arg names
    refs::R                   # default variable values as Quantity (for documentation/display)
    unit::Q                   # output unit, always u"mol/s"
end

const _KF_VARS = (:T, :P, :t, :n, :lna, :n_initial)
const _KF_DEFAULT_REFS = (T = 298.15u"K", P = 1.0e5u"Pa")

# Convenience constructor — vars defaults to the standard 6-argument names.
KineticFunc(compiled, refs::NamedTuple, unit) = KineticFunc(compiled, _KF_VARS, refs, unit)

# Positional call — hot path for ODE integration (no allocation, no unit handling).
(kf::KineticFunc)(T, P, t, n, lna, n_initial) = kf.compiled(T, P, t, n, lna, n_initial)

# Keyword call — user convenience (REPL, scripts), mirroring NumericFunc/SymbolicFunc.
# T, P, t are ustripped (accept Quantity or plain Real); n, lna, n_initial pass through.
@inline function (kf::KineticFunc)(; kwargs...)
    T_raw = haskey(kwargs, :T) ? kwargs[:T] :
        get(kf.refs, :T, get(_KF_DEFAULT_REFS, :T, nothing))
    P_raw = haskey(kwargs, :P) ? kwargs[:P] :
        get(kf.refs, :P, get(_KF_DEFAULT_REFS, :P, nothing))
    t_raw = haskey(kwargs, :t) ? kwargs[:t] : get(kf.refs, :t, 0.0)
    n_val = get(kwargs, :n, nothing)
    lna_val = get(kwargs, :lna, nothing)
    n0_val = get(kwargs, :n_initial, nothing)
    val = kf.compiled(ustrip(T_raw), ustrip(P_raw), ustrip(t_raw), n_val, lna_val, n0_val)
    return get(kwargs, :unit, false) ? val * kf.unit : val
end

function Base.show(io::IO, kf::KineticFunc)
    print(io, "KineticFunc [", dimension(kf.unit), "]")
    if !isempty(kf.vars)
        print(io, " ◆ vars=(", join(kf.vars, ", "), ")")
    end
    return if !isempty(kf.refs)
        print(io, " ◆ ", join(["$k=$v" for (k, v) in pairs(kf.refs)], ", "))
    end
end

function Base.show(io::IO, ::MIME"text/plain", kf::KineticFunc)
    println(io, "KineticFunc:")
    print(io, "  Unit: [", dimension(kf.unit), "]")
    print(io, "\n  Variables: ", join(kf.vars, ", "))
    return if !isempty(kf.refs)
        print(io, "\n  References: ", join(["$k=$v" for (k, v) in pairs(kf.refs)], ", "))
    end
end

# ── KINETICS_RATE_MODELS / KINETICS_RATE_FACTORIES ────────────────────────────

"""
    KINETICS_RATE_MODELS

Dictionary of raw kinetic rate-constant model expressions, analogous to
[`THERMO_MODELS`](@ref).

Each entry maps a model name (`:arrhenius`, …) to a `Dict` containing:
  - `:k` — symbolic `Expr` for the rate constant as a function of variables.
  - `:vars` — list of variable symbols (e.g. `[:T]`).
  - `:units` — list of `Symbol => Quantity` pairs for parameters and variables.
  - `:output_unit` — `Quantity` representing the output unit.

At package initialization, every entry is compiled into a `ThermoFactory`
stored in [`KINETICS_RATE_FACTORIES`](@ref).

# Example

```julia
# The acid mechanism of calcite dissolution, Palandri & Kharaka (2004), Table 33
calcite = literature_row("PalandriKharaka2004", "carbonate_rates", "calcite")
k_acid = KINETICS_RATE_FACTORIES[:arrhenius](;
    k₀    = 10.0^calcite.acid_log_k,           # mol/(m² s) at T_ref
    Ea    = ustrip(us"J/mol", calcite.acid_E),  # J/mol
    T_ref = 298.15,                            # K
)
k_acid(; T = 310.0)   # → Float64 rate constant
```
"""
const KINETICS_RATE_MODELS = Dict{Symbol, Dict}(
    :arrhenius => Dict(
        :k => :(k₀ * exp(-Ea / R_gas * (1 / T - 1 / T_ref))),
        :vars => [:T],
        :units => [
            :T => u"K",
            :T_ref => u"K",
            :k₀ => u"mol/(m^2*s)",
            :Ea => u"J/mol",
            :R_gas => u"J/(mol*K)",
        ],
        :output_unit => u"mol/(m^2*s)",
    ),
)

"""
    KINETICS_RATE_FACTORIES

Compiled `ThermoFactory` objects for each kinetic rate model.
Populated by `__init__()` from [`KINETICS_RATE_MODELS`](@ref).

Keys are model name symbols (e.g. `:arrhenius`).
Values are `ThermoFactory` callables that return `SymbolicFunc{1}` instances.

# Usage

```julia
factory = KINETICS_RATE_FACTORIES[:arrhenius]
k = factory(; k₀=1e-5, Ea=50000.0, T_ref=298.15, R_gas=R_GAS)
k(; T = 298.15)   # → 1e-5  (rate constant at reference temperature)
```
"""
const KINETICS_RATE_FACTORIES = Dict{Symbol, ThermoFactory}()

"""
    add_kinetics_rate_model(name::Symbol, dict_model::Dict)

Register a new kinetic rate-constant model in [`KINETICS_RATE_MODELS`](@ref)
and compile it into [`KINETICS_RATE_FACTORIES`](@ref).

`dict_model` must contain at minimum `:k` (expression), `:vars` (variable list),
`:units` (parameter units), and `:output_unit`.

# Example

```julia
add_kinetics_rate_model(:power_law, Dict(
    :k    => :(k₀ * (T / T_ref)^n),
    :vars => [:T],
    :units => [:T => u"K", :T_ref => u"K", :k₀ => u"mol/(m^2*s)", :n => u"1"],
    :output_unit => u"mol/(m^2*s)",
))
```
"""
function add_kinetics_rate_model(name::Symbol, dict_model::Dict)
    KINETICS_RATE_MODELS[name] = dict_model
    KINETICS_RATE_FACTORIES[name] = _build_kinetics_rate_factory(dict_model)
    return nothing
end

# Internal: compile one KINETICS_RATE_MODELS entry → ThermoFactory
function _build_kinetics_rate_factory(d::Dict)
    return ThermoFactory(
        d[:k],
        get(d, :vars, [:T]);
        units = get(d, :units, nothing),
        output_unit = get(d, :output_unit, u"1"),
    )
end

# ── arrhenius_rate_constant ────────────────────────────────────────────────────

"""
    arrhenius_rate_constant(k₀, Ea; T_ref=298.15, R_gas=R_GAS) -> NumericFunc

Build a temperature-dependent Arrhenius rate constant as a [`NumericFunc`](@ref):

```
k(T) = k₀ × exp(-Eₐ / R × (1/T - 1/T_ref))
```

The returned object is callable as `k(; T=...)` and fully AD-compatible
(ForwardDiff-safe: the closure captures `k₀`, `Ea`, `T_ref`, `R_gas` directly,
so dual numbers propagate correctly through all parameters).

Arithmetic between `SymbolicFunc`/`NumericFunc` objects is supported, so rate
constants can be composed with activity or surface-area functions.

# Arguments

  - `k₀`: pre-exponential factor at `T_ref`. Plain `Real` → SI [mol/(m² s)];
    `Quantity` → automatically converted (e.g. `5e-4u"mol/(m^2*s)"`).
  - `Ea`: activation energy. Plain `Real` → SI [J/mol]; `Quantity` → converted
    (e.g. `62.0u"kJ/mol"`).
  - `T_ref`: reference temperature. Plain `Real` → SI [K]; `Quantity` → converted
    (e.g. `298.15u"K"`). Default `298.15`.
  - `R_gas`: gas constant [J/(mol K)] (plain `Real` only; default [`R_GAS`](@ref),
    the CODATA value taken from `DynamicQuantities.Constants`).

# Returns

A `NumericFunc` with variable `T` (in K) and `refs = (T = T_ref * u"K",)`.

# Examples

```jldoctest
julia> k = arrhenius_rate_constant(5.0e-4, 62000.0);

julia> isapprox(k(; T = 298.15), 5.0e-4; rtol = 1e-10)
true

julia> k(; T = 350.0) > k(; T = 298.15)   # higher T → higher k
true
```

Unit-aware: `k₀` in mmol/(m²·s), `Ea` in kJ/mol, `T_ref` in K — all converted to SI:
```julia
k = arrhenius_rate_constant(0.5u"mmol/(m^2*s)", 62.0u"kJ/mol"; T_ref = 298.15u"K")
```

AD-compatible through all parameters:
```julia
ForwardDiff.derivative(T  -> arrhenius_rate_constant(5e-4, 62000.0)(; T = T),  298.15)
ForwardDiff.derivative(Ea -> arrhenius_rate_constant(5e-4, Ea)(; T = 350.0),   62000.0)
ForwardDiff.derivative(k₀ -> arrhenius_rate_constant(k₀,   62000.0)(; T = 298.15), 5e-4)
```
"""
function arrhenius_rate_constant(
        k₀,
        Ea;
        T_ref = 298.15,
        R_gas::Real = R_GAS,
    )
    k₀_si = safe_ustrip(us"mol/(m^2*s)", k₀)
    Ea_si = safe_ustrip(us"J/mol", Ea)
    T_ref_si = safe_ustrip(us"K", T_ref)
    # Closure captures SI values; no Float64 cast → ForwardDiff.Dual propagates correctly
    # through k₀, Ea, or T_ref when differentiating through construction.
    f = (T) -> k₀_si * exp(-Ea_si / R_gas * (1 / T - 1 / T_ref_si))
    # refs is metadata for default call values — always stored as plain Float64
    refs = (T = Float64(_primal(T_ref_si)) * u"K",)
    return NumericFunc(f, (:T,), refs, u"mol/(m^2*s)")
end

# ── Saturation ratio ───────────────────────────────────────────────────────────

"""
    saturation_ratio(stoich::AbstractVector, lna::AbstractVector,
                     ΔₐG⁰overRT::AbstractVector; ϵ=1e-16) -> Real

Compute the saturation ratio Ω = IAP / K for a kinetic reaction.

```
ln Ω = Σᵢ νᵢ ln aᵢ − ln K
     = Σᵢ νᵢ ln aᵢ + ΔᵣG⁰/(RT)   (note: ln K = −ΔᵣG⁰/RT = −Σᵢ νᵢ ΔₐG⁰ᵢ/RT)
```

where `stoich[i]` is the stoichiometric coefficient (positive for products,
negative for reactants), `lna[i]` is the log-activity of species `i`,
and `ΔₐG⁰overRT[i]` is the dimensionless standard Gibbs energy of formation
`ΔₐG⁰ᵢ / RT` for species `i`.

# Arguments

  - `stoich`: stoichiometric coefficient vector for this reaction (length = number of species).
  - `lna`: log-activity vector (same indexing as species in system).
  - `ΔₐG⁰overRT`: dimensionless standard Gibbs energies `ΔₐG⁰ᵢ/RT`.
  - `ϵ`: floor to avoid `exp` overflow when Ω → ∞.

# Returns

`Ω = exp(ln_IAP - ln_K)` where `ln_K = -ΔᵣG⁰/RT`.

AD-compatible (ForwardDiff-safe).
"""
function saturation_ratio(
        stoich::AbstractVector,
        lna::AbstractVector,
        ΔₐG⁰overRT::AbstractVector;
        ϵ::Real = _AMOUNT_FLOOR,
    )
    # ln IAP = Σᵢ νᵢ ln aᵢ
    ln_iap = sum(stoich[i] * lna[i] for i in eachindex(stoich))
    # ln K = -ΔᵣG⁰/RT = -Σᵢ νᵢ ΔₐG⁰ᵢ/RT
    ln_K = -sum(stoich[i] * ΔₐG⁰overRT[i] for i in eachindex(stoich))
    return exp(ln_iap - ln_K)
end

# ── RateModelCatalyst ─────────────────────────────────────────────────────────

"""
    struct RateModelCatalyst{T<:Real}

Describes the contribution of a catalyst species to a reaction mechanism rate.

The catalyst multiplies the base rate by `exp(n * ln aᵢ) = aᵢ^n`, where
`aᵢ` is the activity of the catalyst species.

# Fields

  - `species`: PHREEQC-format [ParkhurstAppelo2013](@cite) formula string of the
    catalyst species (e.g. `"H+"`, `"OH-"`).
  - `n`: power exponent (dimensionless).

# Examples

```julia
acid_catalyst   = RateModelCatalyst("H+",  0.5)    # ∝ a(H+)^0.5
base_catalyst   = RateModelCatalyst("OH-", 0.5)    # ∝ a(OH-)^0.5
co2_catalyst    = RateModelCatalyst("CO2", 1.0)    # ∝ a(CO2)
```
"""
struct RateModelCatalyst{T <: Real}
    species::String
    n::T
end

# ── RateModelInhibitor ────────────────────────────────────────────────────────

"""
    RateModelInhibitor(species, K; m = 1)

A species that slows a mechanism down: the mechanism's rate is multiplied by

```math
\\left(1 + K\\,a_i\\right)^{-m},
```

``a_i`` the activity of `species`. The factor is one where the species is
absent, and falls as ``a_i^{-m}`` where ``K a_i \\gg 1``: unlike a catalyst
with a negative order, ``a_i^{n}`` with ``n < 0``, which grows without bound as
the species vanishes, it leaves the rate measured without the species unchanged.

# Fields

  - `species`: the species, by its symbol in the system.
  - `K`: the inverse of the activity at which the factor is ``2^{-m}``.
  - `m`: the order where ``K a_i \\gg 1``.

# Examples

```julia
inh = RateModelInhibitor("Ca+2", 2000.0)       # halves the rate at a(Ca²⁺) = 5e-4
```
"""
struct RateModelInhibitor{T <: Real}
    species::String
    K::T
    m::T
end
function RateModelInhibitor(species::AbstractString, K::Real; m::Real = 1)
    K >= 0 || throw(ArgumentError("RateModelInhibitor($species): K must not be negative; got $K."))
    m >= 0 || throw(ArgumentError("RateModelInhibitor($species): m must not be negative; got $m."))
    Kp, mp = promote(float(K), float(m))
    return RateModelInhibitor{typeof(Kp)}(String(species), Kp, mp)
end

# ── RateMechanism ─────────────────────────────────────────────────────────────

"""
    struct RateMechanism{F<:AbstractFunc, T<:Real}

A single kinetic mechanism (acid/neutral/base/…) contributing to the overall
mineral dissolution or precipitation rate.

The mechanism rate is:
```
r_mech = k(T) × [Π_catalysts aᵢ^nᵢ] × [Π_inhibitors (1 + Kⱼ aⱼ)^(-mⱼ)] × sign(1 - Ω) × |1 - Ω^p|^q
```

# Fields

  - `k`: rate constant as `AbstractFunc` (typically `SymbolicFunc{1}` from
    [`arrhenius_rate_constant`](@ref)). Called as `k(; T=...)`.
  - `p`: saturation exponent `p` in `(1 - Ω^p)^q`. Default 1.0.
  - `q`: outer exponent `q`. Default 1.0.
  - `catalysts`: vector of [`RateModelCatalyst`](@ref) (may be empty).
  - `inhibitors`: vector of [`RateModelInhibitor`](@ref) (may be empty, which
    it is unless given).

# Examples

```julia
calcite = literature_row("PalandriKharaka2004", "carbonate_rates", "calcite")
k_acid = arrhenius_rate_constant(10.0^calcite.acid_log_k, calcite.acid_E)
mech   = RateMechanism(k_acid, 1.0, 1.0, [RateModelCatalyst("H+", calcite.acid_n_H)])
```
"""
struct RateMechanism{F <: AbstractFunc, T <: Real}
    k::F
    p::T
    q::T
    catalysts::Vector{RateModelCatalyst{T}}
    inhibitors::Vector{RateModelInhibitor{T}}
end
RateMechanism{F, T}(k::F, p, q, catalysts) where {F <: AbstractFunc, T <: Real} =
    RateMechanism{F, T}(k, T(p), T(q), catalysts, RateModelInhibitor{T}[])

"""
    RateMechanism(k::AbstractFunc, p::Real, q::Real) -> RateMechanism

Construct a [`RateMechanism`](@ref) with no catalyst contributions.
"""
function RateMechanism(k::AbstractFunc, p::Real, q::Real)
    T = typeof(promote(p, q)[1])
    return RateMechanism{typeof(k), T}(k, T(p), T(q), RateModelCatalyst{T}[])
end

"""
    RateMechanism(k::AbstractFunc, p::Real, q::Real, catalysts,
                  inhibitors = RateModelInhibitor[]) -> RateMechanism

Construct a [`RateMechanism`](@ref) with catalysts and, optionally, inhibitors,
their numbers promoted to one type.
"""
function RateMechanism(
        k::AbstractFunc, p::Real, q::Real, catalysts::AbstractVector{<:RateModelCatalyst},
        inhibitors::AbstractVector{<:RateModelInhibitor} = RateModelInhibitor{Float64}[],
    )
    T = promote_type(
        typeof(float(p)), typeof(float(q)), (typeof(c.n) for c in catalysts)...,
        (typeof(i.K) for i in inhibitors)...,
    )
    return RateMechanism{typeof(k), T}(
        k, T(p), T(q), [RateModelCatalyst{T}(c.species, T(c.n)) for c in catalysts],
        [RateModelInhibitor{T}(i.species, T(i.K), T(i.m)) for i in inhibitors],
    )
end

# ── Palandri & Kharaka (2004): the mechanisms of a mineral ───────────────────

# The tables of the report as transcribed, by mineral: those whose name ends in
# `_rates` hold one row per mineral and one group of columns per mechanism.
function _pk04_row(mineral::AbstractString)
    key = "PalandriKharaka2004"
    for t in sort!(collect(keys(literature(key).tables)))
        endswith(t, "_rates") || continue
        tb = literature_table(key, t)
        k = findfirst(==(mineral), tb.mineral)
        k === nothing || return literature_row(key, t, mineral)
    end
    throw(
        ArgumentError(
            "palandri_kharaka: no mineral \"$mineral\" in the transcribed tables; " *
                "`palandri_kharaka_minerals()` lists them."
        )
    )
end

"""
    palandri_kharaka_minerals() -> Vector{String}

The minerals whose dissolution mechanisms [`palandri_kharaka`](@ref) builds: the
rows of the tables of [PalandriKharaka2004](@citet) transcribed in
`data/literature/PalandriKharaka2004.json`.
"""
function palandri_kharaka_minerals()
    key = "PalandriKharaka2004"
    tables = sort!([t for t in keys(literature(key).tables) if endswith(t, "_rates")])
    return [String(m) for t in tables for m in literature_table(key, t).mineral]
end

const _PK04_MECHANISMS = (:acid, :neutral, :base, :carbonate)

"""
    palandri_kharaka(mineral; mechanisms = (:acid, :neutral, :base), pco2 = nothing,
                     assume_Ea = nothing) -> Vector{RateMechanism}

The dissolution mechanisms of `mineral` as [PalandriKharaka2004](@citet) tabulate
them, ready for [`transition_state`](@ref): of the mechanisms named in
`mechanisms`, each one the report gives a rate constant for, with its Arrhenius
rate constant (`log k` at 25 °C and the activation energy, through
[`arrhenius_rate_constant`](@ref)) and its catalyst.

| mechanism | rate per unit area, far from equilibrium | catalyst |
|:--|:--|:--|
| `:acid` | ``k_\\text{acid}(T)\\,a_{\\text{H}^+}^{\\,n}`` | `"H+"`, `n > 0` |
| `:neutral` | ``k_\\text{neutral}(T)`` | none |
| `:base` | ``k_\\text{base}(T)\\,a_{\\text{H}^+}^{\\,n}`` | `"H+"`, `n < 0` |
| `:carbonate` | ``k_\\text{carbonate}(T)\\,P_{\\text{CO}_2}^{\\,n}`` | the species `pco2` |

The report writes the base mechanism as an order in `H⁺` that is negative, which
is how it is applied here rather than converted to an order in `OH⁻`. Each
mechanism takes the saturation term `1 − Ω` (`p = q = 1`): the report gives
other exponents for a few minerals in its text, not in the tables.

`mineral` is a row of the transcribed tables (`"calcite"`, `"quartz, BET surface
area"`, `"gibbsite"`, `"anhydrite"`, …), listed by
[`palandri_kharaka_minerals`](@ref). Where the report gives several sets for one
mineral (quartz normalized by its geometric or by its BET area, amorphous
silica from two data sets), the row says which: the rate per unit area has to be
used with the area it was normalized by.

The carbonate mechanism reads the partial pressure of CO₂, which a rate law sees
as the activity of a gas species, `P/P°` ([`P_STANDARD`](@ref)): asking for
`:carbonate` requires `pco2`, the symbol of that species in the system.

An activation energy the report does not give (that of the neutral mechanism of
gypsum, which its text says the data could not determine) is refused, unless
`assume_Ea` supplies one; the rate constant is then that value's, an assumption
of the caller, and exact at 25 °C whatever it is.

# Examples

```julia
mech = palandri_kharaka("calcite")             # acid and neutral
tst = transition_state(mech, cs, rxn, BETSurfaceArea(90.0))
```

See also: [`palandri_kharaka_minerals`](@ref), [`RateMechanism`](@ref).
"""
function palandri_kharaka(
        mineral::AbstractString; mechanisms = (:acid, :neutral, :base),
        pco2::Union{Nothing, AbstractString} = nothing, assume_Ea = nothing,
    )
    for m in mechanisms
        m in _PK04_MECHANISMS || throw(ArgumentError("palandri_kharaka: no mechanism :$m; the report's are $(_PK04_MECHANISMS)."))
    end
    (:carbonate in mechanisms && pco2 === nothing) && throw(
        ArgumentError("palandri_kharaka: the carbonate mechanism reads P(CO₂); name the gas species with `pco2`.")
    )
    row = _pk04_row(mineral)
    out = RateMechanism[]
    for m in mechanisms
        # A mechanism the report gives no rate constant for, a dash in its
        # table, is not one of this mineral's.
        logk = get(row, Symbol(m, "_log_k"), missing)
        ismissing(logk) && continue
        E = get(row, Symbol(m, "_E"), missing)
        if ismissing(E)
            assume_Ea === nothing && throw(
                ArgumentError(
                    "palandri_kharaka: the report gives no activation energy for the $m mechanism of " *
                        "$mineral; pass `assume_Ea` to use one of your own."
                )
            )
            E = assume_Ea
        end
        k = arrhenius_rate_constant(10.0^ustrip(logk), E)
        cat = if m === :acid || m === :base
            [RateModelCatalyst("H+", float(ustrip(row[Symbol(m, "_n_H")])))]
        elseif m === :carbonate
            [RateModelCatalyst(String(pco2), float(ustrip(row.carbonate_n_PCO2)))]
        else
            RateModelCatalyst{Float64}[]
        end
        push!(out, RateMechanism(k, 1.0, 1.0, cat))
    end
    isempty(out) && throw(ArgumentError("palandri_kharaka: the report gives $mineral none of the mechanisms $(mechanisms)."))
    return out
end

# ── Snellings (2013): the dissolution of a glass at pH 13 ────────────────────

# The abscissa of Fig. 8 of Snellings (2013), with which its regression holds:
# the moles of the modifiers' cations over twice the moles of aluminum plus the
# moles of silicon (the figure counts aluminum twice; see
# data/literature/Snellings2013.json, whose notes give the check).
function _snellings2013_abscissa(oxides, modifiers)
    T = promote_type(Float64, (typeof(float(v)) for v in values(oxides))...)
    moles = Dict{Symbol, T}()
    for (ox, f) in oxides
        el, nu = _oxide_cation(ox)
        moles[el] = get(moles, el, zero(T)) + float(f) * nu / _oxide_molar_mass(ox)
    end
    num = zero(T)
    for m in modifiers
        num += get(moles, first(_oxide_cation(m)), zero(T))
    end
    den = 2 * get(moles, :Al, zero(T)) + get(moles, :Si, zero(T))
    den > 0 || throw(ArgumentError("snellings2013_glass: the analysis holds neither Al2O3 nor SiO2."))
    return num / den
end

"""
    snellings2013_glass(oxides; Ea, T_ref = 293.15, modifiers = ("CaO",),
                        inhibitors = RateModelInhibitor{Float64}[],
                        extrapolate = false) -> RateMechanism

The dissolution of a calcium aluminosilicate glass far from equilibrium at
pH 13, as [Snellings2013](@citet) measured it on six synthetic glasses from the
composition of a slag to silica: per unit BET area and per mole of the glass's
cations, at 20 °C,

```math
\\log_{10} k_{20} = a\\,x + b ,
\\qquad
x = \\frac{\\sum_{\\text{modifiers}} n_M}{2\\,n_{\\text{Al}} + n_{\\text{Si}}} ,
```

``n`` the moles of each cation in the analysis `oxides` (mass fractions over
the molar masses of the library), and ``a = 2.74``, ``b = -8.47`` the regression
printed on [Snellings2013; Fig. 8](@cite), read from
`data/literature/Snellings2013.json`. The abscissa is that of the figure, the
ratio above with aluminum counted twice, as the compositions of the six glasses
of [Snellings2013; Table I](@cite) give it; the regression holds with that
ratio.

What [Snellings2013](@citet) covers, and what it does not:

  - **Composition.** The glasses hold CaO, Al₂O₃ and SiO₂ only, so only the
    oxides named in `modifiers` enter the numerator: CaO by default; the MgO,
    Na₂O or K₂O of a real slag are left out unless named, as an assumption of
    the caller. ``x`` runs from 0 (silica) to 0.556 (the most calcic glass);
    outside, the constant is refused unless `extrapolate = true`.
  - **Temperature.** The rates were measured at 20 °C only: the activation
    energy `Ea` (J/mol, or a quantity) has no default, and the constant is
    [`arrhenius_rate_constant`](@ref) referred to `T_ref` (K), 20 °C.
  - **pH.** At 13 only; [Snellings2013](@citet) gives no dependence on pH.
  - **Solution.** In NaOH alone. Dissolved calcium slows every glass and
    aluminum the glasses whose calcium only balances their aluminum
    ([Snellings2013; Table II](@cite)); the paper gives the measurements, not a
    law, so `inhibitors` ([`RateModelInhibitor`](@ref)) are the caller's.
    Dissolved silicon up to 11 mM, and the affinity above ``15\\,RT``, change
    nothing measurable.

The rate is per mole of cations: the glass species is
[`glass_species`](@ref) with `M = `[`cation_molar_mass`](@ref)`(oxides)`, and its
area a [`BETSurfaceArea`](@ref) with the material's BET surface.

# Examples

```julia
slag = Dict("CaO" => 0.43, "Al2O3" => 0.19, "SiO2" => 0.38)    # glass G1
mech = snellings2013_glass(slag; Ea = 60.0e3)                 # Ea: the caller's
log10(mech.k(; T = 293.15))                                   # 2.74 x - 8.47
```

See also: [`transition_state`](@ref), [`glass_species`](@ref).
"""
function snellings2013_glass(
        oxides::AbstractDict{<:AbstractString, <:Real}; Ea,
        T_ref = 273.15 + ustrip(literature_value("Snellings2013", "temperature_C")),
        modifiers = ("CaO",), inhibitors::AbstractVector{<:RateModelInhibitor} = RateModelInhibitor{Float64}[],
        extrapolate::Bool = false,
    )
    key = "Snellings2013"
    a = ustrip(literature_value(key, "fig8_slope"))
    b = ustrip(literature_value(key, "fig8_intercept"))
    x = _snellings2013_abscissa(oxides, modifiers)
    # The domain is that of the glasses of Table I, their abscissas computed
    # the same way.
    g = literature_table(key, "glasses")
    x_max = maximum(
        _snellings2013_abscissa(
            Dict("CaO" => ustrip(g.CaO_percent[i]) / 100, "Al2O3" => ustrip(g.Al2O3_percent[i]) / 100, "SiO2" => ustrip(g.SiO2_percent[i]) / 100),
            ("CaO",),
        ) for i in eachindex(g.glass)
    )
    # To rounding: the most calcic glass itself sits on the bound.
    (extrapolate || 0 <= x <= x_max + 1.0e-12) || throw(
        DomainError(
            x,
            "snellings2013_glass: the composition gives x = $(round(_primal(x); digits = 3)), outside the " *
                "glasses measured (0 to $(round(x_max; digits = 3))); pass `extrapolate = true` to use the regression there."
        )
    )
    k = arrhenius_rate_constant(10.0^(a * x + b), Ea; T_ref)
    return RateMechanism(k, 1.0, 1.0, RateModelCatalyst{Float64}[], inhibitors)
end

# ── parrott_killoh factory ──────────────────────────────────────────────────────

"""
    parrott_killoh(params::NamedTuple, mineral_name::AbstractString; α_max=1.0) -> KineticFunc

Build a smoothed three-mechanism clinker hydration rate as a
[`KineticFunc`](@ref).

!!! danger "Deprecated, and no longer attributed to Parrott & Killoh"
    Use [`parrott_killoh_avrami`](@ref) with [`PK84_PARAMS_C3S`](@ref) and
    siblings instead. This function is kept so that existing scripts keep
    running, and it warns once per session.

    **Why the attribution is withdrawn.** The formulas below are not those of
    [ParrottKilloh1984](@citet): the nucleation–growth term carries no Avrami logarithm,
    and `K₃` — a shell-formation coefficient — sits in the *diffusion*
    expression where the canonical formulation uses `K₂`. Nor do the shipped
    parameters match any published set: `N₁ = 3.3` is the canonical `n₃`, and
    the canonical `k₃ = 1.1` has no counterpart at all. The primary source,
    [ParrottKilloh1984](@citet), is a conference proceedings without a DOI that
    could not be consulted, so the attribution is retracted rather than repaired
    by an invented calibration.

    **Why it matters.** With `PK_PARAMS_*` the diffusion branch takes over very
    early — measured at `α/α_max` = 0.003 for C₂S, 0.013 for C₃S, 0.057 for C₃A,
    and over the whole range for C₄AF — and governs throughout the interval a
    seven-day run traverses, because its prefactor `3K₃/N₃ = 0.0018 d⁻¹` is 28
    times smaller than the canonical `k₂ = 0.05 d⁻¹`. The rate then integrates
    in closed form,

        α(t) = α_max · [1 − (1 − √(2·K₃·t / N₃))³]

    which is **independent of the phase**, because `K₃ = 0.0024 d⁻¹` and
    `N₃ = 4` are identical in all four parameter sets. Measured on a CEM I at
    w/c = 0.40 over seven days, C₃S, C₂S and C₃A all land on `α ≈ 0.239`, and
    C₄AF lower still at 0.193 — there its own nucleation-growth branch is
    slower than diffusion and limits instead. The weighted mean comes to 0.234. The
    signature is
    unmistakable: `K₁` spans a factor of 18 across the four phases and changes
    almost nothing.

`params` must be a `NamedTuple` with keys `K₁`, `N₁`, `K₂`, `N₂`, `K₃`, `N₃`,
`B`, `Ea`, `T_ref`. All dimensional values accept plain `Real` (SI) or
`DynamicQuantities.Quantity`.

`mineral_name` is the PHREEQC formula string (e.g. `"C3S"`) used to look up the
mineral moles in the `n` and `n_initial` [`StateView`](@ref)s.

Three competing mechanisms determine the rate [ParrottKilloh1984](@cite):

| Mechanism | Formula |
|-----------|---------|
| Nucleation–growth | `r_NG = (K₁/N₁)(1-ξ)^N₁ / (1 + B·ξ^N₃)` |
| Interaction | `r_I = K₂(1-ξ)^N₂` |
| Diffusion | `r_D = 3K₃(1-ξ)^(2/3) / (N₃·(1-(1-ξ)^(1/3)))` |

The rate [mol/s] is `n_initial × Aₜ × min(max(r_NG, r_I), r_D)` where
`ξ = α / α_max` is the normalized degree of hydration and
`Aₜ = exp(-Ea/R × (1/T - 1/T_ref))` is the Arrhenius factor.

`α_max` can be set to apply the water/cement ratio limit of [Powers1948](@citet):
`α_max = powers_alpha_max(w_c)`.

# Returns

A [`KineticFunc`](@ref) — callable as
`pk(T, P, t, n::StateView, lna::StateView, n_initial::StateView) -> Real [mol/s]`.
AD-compatible (ForwardDiff-safe): no `Float64` casts in the evaluation path.

# Examples

Shown rather than run: the constructor warns by design, and a deprecation
warning on stderr is not something a doctest should have to reproduce verbatim.
The numerical behavior of this variant is pinned by the test suite instead — see
the branch oracle in `test/kinetics/test_rate_models.jl`, which asserts that
`PK_PARAMS_*` puts C3S, C2S and C3A on the diffusion branch within a few percent
of hydration and that C4AF is limited by its own nucleation branch.

```julia
pk = parrott_killoh(PK_PARAMS_C3S, "C3S")
idx = Dict("C3S" => 1)
n0 = StateView([1.0], idx)
lna = StateView([0.0], idx)
pk(293.15, 1e5, 0.0, n0, lna, n0) > 0      # true
```

See also: [`PK_PARAMS_C3S`](@ref), [`PK_PARAMS_C2S`](@ref),
[`PK_PARAMS_C3A`](@ref), [`PK_PARAMS_C4AF`](@ref).
"""
function parrott_killoh(params::NamedTuple, mineral_name::AbstractString; α_max::Real = 1.0)
    K₁ = safe_ustrip(us"1/s", params.K₁)
    N₁ = float(params.N₁)
    K₂ = safe_ustrip(us"1/s", params.K₂)
    N₂ = float(params.N₂)
    K₃ = safe_ustrip(us"1/s", params.K₃)
    N₃ = float(params.N₃)
    B = float(params.B)
    Ea = safe_ustrip(us"J/mol", params.Ea)
    T_ref = safe_ustrip(us"K", params.T_ref)
    @warn """`parrott_killoh` is deprecated: its formulas and parameters are not those of \
    Parrott & Killoh. With `PK_PARAMS_*` the diffusion branch takes over within the first \
    few percent of hydration (α ≈ 0.003 for C2S, 0.013 for C3S, 0.057 for C3A), so those \
    three phases reach α(7 d) = 0.2386 whatever their K₁; C4AF is limited by its own \
    nucleation branch instead and reaches only 0.193. Use `parrott_killoh_avrami` with \
    `PK84_PARAMS_*`.""" maxlog = 1

    α_max_f = float(α_max)

    f = (T, _P, _t, n, _lna, n_initial) -> begin
        n_m = n[mineral_name]
        n_init = max(n_initial[mineral_name], oneunit(n_m) * 1.0e-30)
        # degree of hydration α ∈ [0, α_max)
        α = min(max(one(T) - n_m / n_init, zero(T)), α_max_f - oftype(T, 1.0e-10))
        ξ = α / α_max_f
        # Arrhenius temperature correction
        Aₜ = exp(-Ea / R_GAS * (one(T) / T - one(T) / T_ref))
        one_m_ξ = one(ξ) - ξ
        # r_NG: nucleation–growth [s⁻¹]
        r_NG = (K₁ / N₁) * one_m_ξ^N₁ / (one(ξ) + B * ξ^N₃)
        # r_I: interaction [s⁻¹]
        r_I = K₂ * one_m_ξ^N₂
        # r_D: diffusion [s⁻¹] (denominator clamped to avoid 0/0 at α=0)
        denom_D = max(one(ξ) - one_m_ξ^(one(ξ) / 3), oftype(ξ, 1.0e-10))
        r_D = 3 * K₃ * one_m_ξ^(2 * one(ξ) / 3) / (N₃ * denom_D)
        return n_init * Aₜ * min(max(r_NG, r_I), r_D)
    end

    refs = (T = Float64(_primal(T_ref)) * u"K", P = 1.0e5u"Pa")
    return KineticFunc(f, refs, u"mol/s")
end

# ── Parameters of the deprecated smoothed variant ────────────────────────────
#
# Read from `data/literature/ParrottKilloh1984.json`, table `smoothed_variant`,
# where their provenance is recorded as unstated: neither the rate constants nor
# the activation energies (attributed to Schindler & Folliard, 2005) could be
# traced to a source.

function _pk_smoothed_params(phase::AbstractString)
    t = literature_table("ParrottKilloh1984", "smoothed_variant")
    i = findfirst(==(phase), t.phase)
    return (
        K₁ = t.K1[i], N₁ = t.N1[i], K₂ = t.K2[i], N₂ = t.N2[i],
        K₃ = t.K3[i], N₃ = t.N3[i], B = t.B[i], Ea = t.Ea[i],
        T_ref = literature_value("ParrottKilloh1984", "T_ref"),
    )
end

"""
    PK_PARAMS_C3S :: NamedTuple

Parameters of the deprecated smoothed variant [`parrott_killoh`](@ref) for alite
(C₃S = Ca₃SiO₅), with keys `K₁`, `N₁`, `K₂`, `N₂`, `K₃`, `N₃`, `B`, `Ea`,
`T_ref`.

Their provenance is unestablished. They are read from the table
`smoothed_variant` of `data/literature/ParrottKilloh1984.json`, which records
them as such.

```julia
pk = parrott_killoh(PK_PARAMS_C3S, "C3S")
# or with the water availability cap of Powers (1948):
pk = parrott_killoh(PK_PARAMS_C3S, "C3S"; α_max = powers_alpha_max(w_c))
```
"""
const PK_PARAMS_C3S = _pk_smoothed_params("C3S")

"""
    PK_PARAMS_C2S :: NamedTuple

Parameters of the deprecated smoothed variant for belite (C₂S = Ca₂SiO₄). See
[`PK_PARAMS_C3S`](@ref).
"""
const PK_PARAMS_C2S = _pk_smoothed_params("C2S")

"""
    PK_PARAMS_C3A :: NamedTuple

Parameters of the deprecated smoothed variant for tricalcium aluminate
(C₃A = Ca₃Al₂O₆) in the presence of sulfate (gypsum), corresponding to ettringite
formation. See [`PK_PARAMS_C3S`](@ref).
"""
const PK_PARAMS_C3A = _pk_smoothed_params("C3A")

"""
    PK_PARAMS_C4AF :: NamedTuple

Parameters of the deprecated smoothed variant for tetracalcium aluminoferrite
(C₄AF = Ca₄Al₂Fe₂O₁₀). See [`PK_PARAMS_C3S`](@ref).
"""
const PK_PARAMS_C4AF = _pk_smoothed_params("C4AF")

# ── parrott_killoh_avrami — the canonical 1984 formulation ────────────────────

"""
    PK_AVRAMI_SEED

Lower bound imposed on the normalized degree of hydration inside the Avrami
branch of [`parrott_killoh_avrami`](@ref), so that the rate is strictly positive
at `α = 0` and the ODE leaves its degenerate initial point. See the discussion
in [`parrott_killoh_avrami`](@ref).

One visible consequence: a phase whose rate is governed by the Avrami branch
near `α = 0` (C₃S, C₃A, C₄AF) starts more slowly than one governed by the power
law (C₂S), so alite only overtakes belite after a few minutes. The crossover
falls inside the induction period, which this model does not describe anyway.
"""
const PK_AVRAMI_SEED = 1.0e-6

"""
    parrott_killoh_avrami(params::NamedTuple, mineral_name::AbstractString;
                         α_max = 1.0, blaine = nothing, humidity = nothing,
                         w_c = nothing, H = nothing) -> KineticFunc

Build the clinker hydration rate of [ParrottKilloh1984](@citet) in its **canonical
formulation**, as reported by [Lothenbach2008](@citet) and used by
[Lavergne2018](@citet).

`params` must be a `NamedTuple` with keys `k₁`, `n₁`, `k₂`, `k₃`, `n₃`, `Ea`,
`T_ref` — see [`PK84_PARAMS_C3S`](@ref) and siblings. Dimensional values accept
plain `Real` (SI) or `DynamicQuantities.Quantity`.

Three competing mechanisms limit the rate, and the **slowest one wins**:

| Mechanism | Formula |
|-----------|---------|
| Nucleation–growth (Avrami) | `α̇₁ = (k₁/n₁)(1-ξ)(-ln(1-ξ))^(1-n₁)` |
| Diffusion (Jander) | `α̇₂ = k₂(1-ξ)^(2/3) / (1-(1-ξ)^(1/3))` |
| Shell formation (power law) | `α̇₃ = k₃(1-ξ)^n₃` |

so that `α̇ = min(α̇₁, α̇₂, α̇₃)`, with `ξ = α/α_max` the normalized degree of
hydration. The returned rate [mol/s] is `n_initial × Aₜ × β_B × β_h × f × α̇`,
where `Aₜ = exp(-Ea/R × (1/T - 1/T_ref))` is the Arrhenius factor, `β_B` the
Blaine fineness factor ([`blaine_factor`](@ref)), `β_h` the relative-humidity
reduction ([`humidity_factor`](@ref)) and `f` the water/cement factor
([`pk_wc_factor`](@ref)). Each is 1 when its keyword is `nothing`.

!!! note "Two Parrott–Killoh variants ship with ChemistryLab"
    [`parrott_killoh`](@ref) implements a *different*, smoothed variant
    (`min(max(r_NG, r_I), r_D)` with a `B`-damped nucleation term) together with
    the parameter set of [`PK_PARAMS_C3S`](@ref) and siblings. The two are not
    interchangeable: their parameters are **not** transferable, and only
    `parrott_killoh_avrami` with [`PK84_PARAMS_C3S`](@ref) reproduces the α(t)
    curves published in the cement literature cited above.

With the canonical parameters, C₂S has no nucleation–growth stage and C₃S has no
diffusion-controlled stage — an artifact of the 1984 fit that
[ParrottKilloh1984](@citet) acknowledged, and a convenient signature to check an
implementation against.

# Keyword arguments

  - `α_max`: the water availability cap of [Powers1948](@citet) — see [`powers_alpha_max`](@ref).
  - `blaine`: Blaine fineness of the binder, as a `Quantity`, a plain `Real` in
    m²/kg, a [`BlaineSurfaceArea`](@ref) — all three frozen for the whole
    integration — or a [`ShrinkingCoreArea`](@ref), which makes the factor
    follow the grains as they are consumed. `nothing` (default) means no
    correction. See the warning below before using the last one.
  - `humidity`: internal relative humidity, either a constant in `[0, 1]` or a
    callable `t -> h(t)`. `nothing` (default) means no correction.
  - `w_c`: the water/cement ratio, which slows the hydration once the degree
    exceeds a critical fraction of it ([`pk_wc_factor`](@ref)); `H` the critical
    degree of the phase in place of that of [ParrottKilloh1984](@citet), as
    [Lothenbach2008](@citet) fit one per clinker phase. `nothing` (default) means no correction.

!!! warning "An evolving fineness is not a free improvement"
    The constants of [ParrottKilloh1984](@citet) were fitted with `β_B` **constant**, so
    passing a [`ShrinkingCoreArea`](@ref) leaves this law outside the
    calibration it came with and it has to be recalibrated — the machinery is
    `scripts/hydration_calibration.jl`.

    Worse, the extra freedom largely **already exists** in the law. With
    `α_max = 1` the remaining fraction is `n/n₀ = 1 - ξ`, so multiplying by
    `(1-ξ)^p` turns the shell-formation branch `k₃(1-ξ)^{n₃}` into
    `k₃(1-ξ)^{n₃+p}`: wherever that branch is the active one, `p` and `n₃` are
    the same parameter written twice, and fitting both is fitting a sum. The
    Jander branch has no such exponent, so `p` is a genuine degree of freedom
    only while diffusion controls. Which of the two holds over a given dataset
    is a measurement, and [`identifiability`](@ref) is what makes it — do not
    report `p` and `n₃` from one curve without it.

# Returns

A [`KineticFunc`](@ref) — callable as
`pk(T, P, t, n::StateView, lna::StateView, n_initial::StateView) -> Real [mol/s]`.
AD-compatible (ForwardDiff-safe): no `Float64` casts in the evaluation path.

# Examples

```jldoctest
julia> pk = parrott_killoh_avrami(PK84_PARAMS_C3S, "C3S"; blaine = 380u"m^2/kg");

julia> idx = Dict("C3S" => 1);

julia> n0 = StateView([1.0], idx);

julia> lna = StateView([0.0], idx);

julia> pk(293.15, 1e5, 3600.0, StateView([0.9], idx), lna, n0) > 0
true
```

See also: [`PK84_PARAMS_C3S`](@ref), [`waller`](@ref), [`blaine_factor`](@ref),
[`humidity_factor`](@ref), [`powers_alpha_max`](@ref).
"""
function parrott_killoh_avrami(
        params::NamedTuple, mineral_name::AbstractString;
        α_max::Real = 1.0, blaine = nothing, humidity = nothing,
        w_c = nothing, H = nothing,
    )
    k₁ = safe_ustrip(us"1/s", params.k₁)
    n₁ = float(params.n₁)
    k₂ = safe_ustrip(us"1/s", params.k₂)
    k₃ = safe_ustrip(us"1/s", params.k₃)
    n₃ = float(params.n₃)
    Ea = safe_ustrip(us"J/mol", params.Ea)
    T_ref = safe_ustrip(us"K", params.T_ref)
    α_max_f = float(α_max)
    β_B0, β_Bp = _fineness_parts(blaine, PK_BLAINE_REF)

    f = (T, _P, t, n, _lna, n_initial) -> begin
        n_m = n[mineral_name]
        n_init = max(n_initial[mineral_name], oneunit(n_m) * 1.0e-30)
        α = min(max(one(T) - n_m / n_init, zero(T)), α_max_f - oftype(T, 1.0e-10))
        ξ = α / α_max_f
        Aₜ = exp(-Ea / R_GAS * (one(T) / T - one(T) / T_ref))
        β_B = _fineness_at(β_B0, β_Bp, max(n_m / n_init, zero(ξ)))
        β_h = humidity === nothing ? one(ξ) : humidity_factor(_humidity_at(humidity, t, n))
        one_m_ξ = max(one(ξ) - ξ, oftype(ξ, 1.0e-12))
        # α̇₁ — Avrami nucleation and growth. For n₁ < 1 the (-ln(1-ξ))^(1-n₁)
        # factor VANISHES at ξ = 0, so α̇ = 0 and α ≡ 0 solves the ODE: hydration
        # would never start. Parrott & Killoh's own discrete scheme escapes this by
        # evaluating the integrated Avrami law over the first finite time step;
        # a continuous solver cannot, so the argument is floored at ξ_seed. The
        # seed only sets how fast the solution leaves the degenerate point — by
        # ξ ≈ 1e-3 the Avrami branch is already three orders of magnitude above it.
        ξ_avrami = max(ξ, oftype(ξ, PK_AVRAMI_SEED))
        r₁ = (k₁ / n₁) * one_m_ξ * (-log(one(ξ) - ξ_avrami))^(one(ξ) - n₁)
        # α̇₂ — Jander diffusion through the hydrate layer. Denominator → 0 as ξ → 0.
        denom = max(one(ξ) - one_m_ξ^(one(ξ) / 3), oftype(ξ, 1.0e-12))
        r₂ = k₂ * one_m_ξ^(2 * one(ξ) / 3) / denom
        # α̇₃ — power law, thick shell around the unreacted grain.
        r₃ = k₃ * one_m_ξ^n₃
        f_wc = w_c === nothing ? one(ξ) : pk_wc_factor(α, w_c; H)
        return n_init * Aₜ * β_B * β_h * f_wc * min(r₁, r₂, r₃)
    end

    refs = (T = Float64(_primal(T_ref)) * u"K", P = 1.0e5u"Pa")
    return KineticFunc(f, refs, u"mol/s")
end

# The constants of the water/cement factor, read from the papers that print
# them: Parrott and Killoh's form as Lothenbach and Winnefeld (2006, Section 4.1)
# give it, and the slope of the form with a critical degree per phase as
# Lothenbach et al. (2008, Section 3.2) write it.
const _PK_WC = (
    intercept = ustrip(literature_value("LothenbachWinnefeld2006", "pk_wc_intercept")),
    slope = ustrip(literature_value("LothenbachWinnefeld2006", "pk_wc_slope")),
    critical = ustrip(literature_value("LothenbachWinnefeld2006", "pk_wc_critical_degree")),
    slope_H = ustrip(literature_value("LothenbachLeSaout2008", "pk_wc_slope")),
)

"""
    pk_wc_factor(α, w_c; H = nothing) -> Real

The water/cement factor of the law of [ParrottKilloh1984](@citet): 1 while the degree
of hydration `α` of the phase is below a critical fraction of `w_c`, and beyond it

```math
f = (1 + 4.444\\,w/c - 3.333\\,\\alpha)^4 \\quad (\\alpha > 1.333\\,w/c),
```

as [LothenbachWinnefeld2006; Section 4.1](@citet) state it, which stops the phase at
``\\alpha = (1 + 4.444\\,w/c)/3.333``; the printed constants being rounded
(``3.333 \\times 1.333 = 4.443``), it steps up by ``(1 + 0.0011\\,w/c)^4`` at the
threshold. With `H`, the critical degree of the phase
that [Lothenbach2008; Section 3.2](@citet) fit for each clinker phase,

```math
f = (1 + 3.333\\,(H\\,w/c - \\alpha))^4 \\quad (\\alpha > H\\,w/c).
```

The constants are read from `data/literature/`. The water/cement ratio is a
proxy for the room left to grow in, so the factor is the long-term ceiling of the
law; [`powers_alpha_max`](@ref) is another statement of the same limit.
"""
function pk_wc_factor(α, w_c; H = nothing)
    if H === nothing
        α <= _PK_WC.critical * w_c && return one(α)
        return max(1 + _PK_WC.intercept * w_c - _PK_WC.slope * α, zero(α))^4
    end
    α <= H * w_c && return one(α)
    return max(1 + _PK_WC.slope_H * (H * w_c - α), zero(α))^4
end

# Internal: a humidity keyword is a constant, a function of TIME, or a
# `PoreHumidity`, which reads the current COMPOSITION. The composition is passed
# to all three so the last one is reachable; the first two ignore it, so nothing
# a caller wrote before changes.
@inline _humidity_at(h::Real, _t, _n) = h
@inline _humidity_at(h, t, _n) = h(t)
# The `PoreHumidity` method is further down, where that type is defined.

# ── Canonical Parrott & Killoh (1984) parameters ──────────────────────────────
#
# Table 3 of Lavergne et al. (2018), themselves quoting Parrott & Killoh (1984)
# as reported by Lothenbach et al. (2008); activation energies from Table 4
# (Maekawa et al.), which differ markedly from the uniform values often assumed:
# the minerals that hydrate later have the *lower* apparent activation energy,
# E_C3A > E_C3S > E_C4AF > E_C2S. Both tables are read from
# `data/literature/Lavergne2018.json`.

function _pk84_params(phase::AbstractString)
    t = literature_table("Lavergne2018", "parrott_killoh_1984")
    e = literature_table("Lavergne2018", "activation_energies")
    i = findfirst(==(phase), t.phase)
    j = findfirst(==(phase), e.phase)
    return (
        k₁ = t.k1[i], n₁ = t.n1[i], k₂ = t.k2[i], k₃ = t.k3[i], n₃ = t.n3[i],
        Ea = e.Ea[j], T_ref = literature_value("Lavergne2018", "T_ref"),
    )
end

"""
    PK84_PARAMS_C3S :: NamedTuple

Canonical parameters of [ParrottKilloh1984](@citet) for alite (C₃S = Ca₃SiO₅), with keys
`k₁`, `n₁`, `k₂`, `k₃`, `n₃`, `Ea`, `T_ref`. They are valid for the Blaine
fineness [`PK_BLAINE_REF`](@ref) and a reference temperature of 20 °C.

The rate constants and exponents are those of Table 3 of [Lavergne2018](@citet),
the activation energy that of their Table 4, both read from
`data/literature/Lavergne2018.json`.

Pass to [`parrott_killoh_avrami`](@ref), **not** to [`parrott_killoh`](@ref) —
the two use different functional forms and their parameters are not transferable.
"""
const PK84_PARAMS_C3S = _pk84_params("C3S")

"""
    PK84_PARAMS_C2S :: NamedTuple

Canonical parameters of [ParrottKilloh1984](@citet) for belite (C₂S = Ca₂SiO₄).

With `n₁ = 1` the Avrami branch reduces to `k₁(1-ξ)`, which never limits the
rate: belite hydration is governed by the power law throughout.

See [`PK84_PARAMS_C3S`](@ref).
"""
const PK84_PARAMS_C2S = _pk84_params("C2S")

"""
    PK84_PARAMS_C3A :: NamedTuple

Canonical parameters of [ParrottKilloh1984](@citet) for tricalcium aluminate
(C₃A = Ca₃Al₂O₆). See [`PK84_PARAMS_C3S`](@ref).
"""
const PK84_PARAMS_C3A = _pk84_params("C3A")

"""
    PK84_PARAMS_C4AF :: NamedTuple

Canonical parameters of [ParrottKilloh1984](@citet) for tetracalcium aluminoferrite
(C₄AF = Ca₄Al₂Fe₂O₁₀). See [`PK84_PARAMS_C3S`](@ref).
"""
const PK84_PARAMS_C4AF = _pk84_params("C4AF")

# ── waller — pozzolanic and latent-hydraulic additions ───────────────────────

"""
    waller(params::NamedTuple, mineral_name::AbstractString;
           α_max = 1.0, blaine = nothing, humidity = nothing) -> KineticFunc

Build the reaction rate of a pozzolanic or latent-hydraulic addition (fly ash,
silica fume, ground granulated slag) after [Waller1999](@citet), as a
[`KineticFunc`](@ref).

The degree of reaction follows a sigmoid in log-time,

```math
α(t) = \\frac{1}{1 + (τ/t)^n}
```

whose rate, written as a function of the current degree so that it composes with
temperature, fineness and humidity corrections, is

```math
α̇ = \\frac{n}{τ} (1 - α)^{1 + 1/n} α^{1 - 1/n}.
```

`params` must be a `NamedTuple` with keys `τ`, `n`, `Ea`, `T_ref` and, optionally,
`blaine_ref` — see [`WALLER_PARAMS_FLY_ASH`](@ref).

Pozzolanic reactions are markedly more temperature-sensitive than the hydraulic
reactions of clinker: the shipped activation energy is higher than that of any
clinker phase in [`PK84_PARAMS_C3S`](@ref) and its siblings.

Silica fume is far finer than the cement, by an order of magnitude when its
surface is measured by BET. Its reactivity is nonetheless represented here
through an effective *Blaine* fineness ([`WALLER_PARAMS_SILICA_FUME`](@ref)) —
the two measurements probe different physical phenomena and are not
interchangeable.

# Keyword arguments

Identical to [`parrott_killoh_avrami`](@ref). The Blaine correction is taken
relative to `params.blaine_ref`, the fineness of the fly ash the kinetics were
adjusted to, not to the clinker reference of [`PK_BLAINE_REF`](@ref).

!!! note "Here an evolving area is a real degree of freedom"
    The degeneracy warned about on [`parrott_killoh_avrami`](@ref) does **not**
    carry over in the same form. The sigmoid rate is
    `(n/τ)(1-ξ)^{1+1/n} ξ^{1-1/n}`, and `n` sets both exponents at once, in
    opposite directions; a [`ShrinkingCoreArea`](@ref) exponent `p` shifts only
    the `(1-ξ)` one. So `p` is not a rewriting of `n` — it is a separate shape,
    correlated with `n` but distinguishable. The calibration caveat still
    applies in full: these constants were fitted with `β_B` frozen.

# Examples

```jldoctest
julia> fa = waller(WALLER_PARAMS_FLY_ASH, "FlyAsh");

julia> idx = Dict("FlyAsh" => 1);

julia> n0 = StateView([1.0], idx);

julia> lna = StateView([0.0], idx);

julia> fa(293.15, 1e5, 86400.0, StateView([0.9], idx), lna, n0) > 0
true
```

No parameter set is shipped for a slag: give its characteristic time from a
source of your own, `waller(merge(WALLER_PARAMS_FLY_ASH, (τ = τ_slag,)), "GGBS")`.

See also: [`WALLER_PARAMS_FLY_ASH`](@ref), [`WALLER_PARAMS_SILICA_FUME`](@ref),
[`parrott_killoh_avrami`](@ref).
"""
function waller(
        params::NamedTuple, mineral_name::AbstractString;
        α_max::Real = 1.0, blaine = nothing, humidity = nothing
    )
    τ = safe_ustrip(us"s", params.τ)
    n_w = float(params.n)
    Ea = safe_ustrip(us"J/mol", params.Ea)
    T_ref = safe_ustrip(us"K", params.T_ref)
    α_max_f = float(α_max)
    blaine_ref = hasproperty(params, :blaine_ref) ? params.blaine_ref :
        WALLER_PARAMS_FLY_ASH.blaine_ref
    β_B0, β_Bp = _fineness_parts(blaine, blaine_ref)

    f = (T, _P, t, n, _lna, n_initial) -> begin
        n_m = n[mineral_name]
        n_init = max(n_initial[mineral_name], oneunit(n_m) * 1.0e-30)
        α = min(max(one(T) - n_m / n_init, zero(T)), α_max_f - oftype(T, 1.0e-10))
        ξ = α / α_max_f
        Aₜ = exp(-Ea / R_GAS * (one(T) / T - one(T) / T_ref))
        β_B = _fineness_at(β_B0, β_Bp, max(n_m / n_init, zero(ξ)))
        β_h = humidity === nothing ? one(ξ) : humidity_factor(_humidity_at(humidity, t, n))
        # At ξ = 0 the closed form α̇(α) is singular (α^(1-1/n) → ∞ for n < 1).
        # Fall back to the explicit α̇(t) of the sigmoid, which is finite for t > 0
        # and vanishes as t → 0 — the induction period the sigmoid encodes.
        if _primal(ξ) < 1.0e-12
            t_pos = max(t, oftype(ξ, 1.0e-6))
            x = (τ / t_pos)^n_w
            r = (n_w / t_pos) * x / (one(ξ) + x)^2
        else
            one_m_ξ = max(one(ξ) - ξ, oftype(ξ, 1.0e-12))
            r = (n_w / τ) * one_m_ξ^(one(ξ) + one(ξ) / n_w) * ξ^(one(ξ) - one(ξ) / n_w)
        end
        return n_init * Aₜ * β_B * β_h * r
    end

    refs = (T = Float64(_primal(T_ref)) * u"K", P = 1.0e5u"Pa")
    return KineticFunc(f, refs, u"mol/s")
end

# The fly-ash parameters are those of p. 42 of Lavergne et al. (2018), adjusted
# to the results of Waller (1999), read from `data/literature/`. No slag set:
# the 100-day time this package attributed to the thesis is not in it (see
# `data/literature/Waller1999.json`).

_waller_params(τ) = (
    τ = τ, n = literature_value("Lavergne2018", "waller_n"),
    blaine_ref = literature_value("Lavergne2018", "waller_blaine_ref"),
    Ea = literature_value("Lavergne2018", "waller_Ea"),
    T_ref = literature_value("Lavergne2018", "T_ref"),
)

"""
    WALLER_PARAMS_FLY_ASH :: NamedTuple

Parameters of [Waller1999](@citet) for class-F fly ash, with keys `τ`, `n`, `blaine_ref`,
`Ea`, `T_ref`, as adjusted by [Lavergne2018; p. 42](@citet) to SEM image
analyses of a fly ash of assumed pozzolanic activity. The values are read from
`data/literature/Lavergne2018.json`.

Pass to [`waller`](@ref).
"""
const WALLER_PARAMS_FLY_ASH = _waller_params(literature_value("Lavergne2018", "waller_tau_fly_ash"))

"""
    WALLER_PARAMS_SILICA_FUME :: NamedTuple

Parameters of [Waller1999](@citet) applied to silica fume — identical kinetics to
[`WALLER_PARAMS_FLY_ASH`](@ref), the higher reactivity being carried by the
fineness. Pass to [`waller`](@ref) the effective Blaine fineness recommended by
[Lavergne2018](@citet),
`blaine = literature_value("Lavergne2018", "blaine_silica_fume")`; the BET
surface of silica fume is **not** a Blaine fineness and must not be used here.
"""
const WALLER_PARAMS_SILICA_FUME = WALLER_PARAMS_FLY_ASH


# ── Correction factors ───────────────────────────────────────────────────────

"""
    PK_BLAINE_REF :: Quantity

The Blaine fineness the constants of [ParrottKilloh1984](@citet) were calibrated
at, B₀ of [Lavergne2018; p. 39](@citet), read from `data/literature/Lavergne2018.json`.

It is the default reference of [`blaine_factor`](@ref) and the one
[`parrott_killoh_avrami`](@ref) applies, written once so the two cannot drift
apart — a rate constant and the fineness it was measured against are one datum
in two places.
"""
const PK_BLAINE_REF = literature_value("Lavergne2018", "blaine_ref_clinker")

# A bare number or quantity is a Blaine fineness, which is what every caller of
# `blaine_factor` has always meant. A typed area passes through, and that is the
# whole guard: a `BETSurfaceArea` reaches `area_ratio` as a BET area and is
# refused there rather than divided.
#
# Defined here rather than just above the function: a comment between a
# docstring's closing quotes and the definition detaches the docstring, silently.
_as_blaine(x::AbstractSpecificArea) = x
_as_blaine(x) = BlaineSurfaceArea(x)

"""
    _fineness_parts(blaine, blaine_ref) -> (β₀, p)

The fineness factor split into the part that can be computed once and the part
that has to be read from the state at every evaluation.

`p === nothing` means the factor is **frozen**, which is what a bare fineness has
always meant and what every rate law shipped in this package was calibrated
with. A [`ShrinkingCoreArea`](@ref) instead returns its exponent, and the caller
multiplies `β₀` by `g(n/n₀)` at each evaluation.

The split exists so that the frozen path keeps costing nothing: a constant
fineness is still one `area_ratio` for the whole integration, evaluated here.
"""
_fineness_parts(::Nothing, _ref) = (1.0, nothing)
_fineness_parts(m::ShrinkingCoreArea, ref) =
    (area_ratio(m.initial, _as_blaine(ref)), m.exponent)
_fineness_parts(m, ref) = (blaine_factor(m; blaine_ref = ref), nothing)

"""
    _fineness_at(β₀, p, fraction) -> Real

The fineness factor at a remaining fraction `n/n₀`. `p === nothing` returns `β₀`
**identically** — same object, no arithmetic — so a law with a frozen factor
computes exactly the number it computed before this existed.
"""
@inline _fineness_at(β₀, ::Nothing, _fraction) = β₀
@inline _fineness_at(β₀, p, fraction) = β₀ * _shrink_fraction(fraction, p)

"""
    blaine_factor(blaine; blaine_ref = PK_BLAINE_REF) -> Real

Fineness correction of the hydration rate: the parameters of
[ParrottKilloh1984](@citet) were
adjusted for a cement of Blaine fineness `blaine_ref`, and the rate scales as
`blaine / blaine_ref`.

Both arguments accept a `DynamicQuantities.Quantity`, a plain `Real` in m²/kg,
or a [`BlaineSurfaceArea`](@ref). The default reference is [`PK_BLAINE_REF`](@ref)
for clinker phases; pass `blaine_ref = WALLER_PARAMS_FLY_ASH.blaine_ref` for the
Waller kinetics of additions.

!!! warning "A BET area is not a Blaine fineness"
    Passing a [`BETSurfaceArea`](@ref) raises instead of returning a number.
    Silica fume is about 20 000 m²/kg by BET and about 2 000 m²/kg by effective
    Blaine, so the substitution would multiply its hydration rate by ten. A bare
    number is still read as a Blaine fineness, which is the historical contract;
    type the argument to get the check.

# Examples

```jldoctest
julia> blaine_factor(PK_BLAINE_REF)
1.0

julia> round(blaine_factor(1.2 * PK_BLAINE_REF); digits = 4)
1.2
```
"""
function blaine_factor(blaine; blaine_ref = PK_BLAINE_REF)
    return area_ratio(_as_blaine(blaine), _as_blaine(blaine_ref))
end

"""
    humidity_factor(h) -> Real

Reduction coefficient `β_h` applied to the hydration rate at internal relative
humidity `h ∈ [0, 1]`, as [Lavergne2018; Eq. 10](@citet) write it:

```math
β_h = \\left(\\frac{h - 0.55}{0.45}\\right)^4 \\ \\text{if } h > 0.80,
\\qquad β_h = 0 \\ \\text{otherwise.}
```

Hydration is taken to stop below 80 % relative humidity, on thermodynamic
grounds. The cut is a genuine discontinuity: the one-sided limit from above is
`β_h(0.80⁺) ≈ 0.0953` while the value at and below 0.80 is exactly 0. The jump
is mild in practice, the rate having already fallen by an order of magnitude
from `β_h(0.99) ≈ 0.914`.

# Examples

```jldoctest
julia> humidity_factor(0.75)
0.0

julia> humidity_factor(0.80)
0.0

julia> round(humidity_factor(0.801); digits = 4)
0.0968
```
"""
function humidity_factor(h::Real)
    return _primal(h) > 0.8 ? ((h - oftype(h, 0.55)) / oftype(h, 0.45))^4 : zero(h)
end

"""
    PoreHumidity(retention, system; reference, T = temperature(reference))

The internal relative humidity of a sealed paste, computed from the composition
it currently has.

Pass it as the `humidity` keyword of [`parrott_killoh_avrami`](@ref) or
[`waller`](@ref) and the rate law stops reading a humidity imposed from outside
and starts reading the one the material makes for itself. That is what closes the
loop: hydration consumes water, the pore space empties, the humidity falls, and
[`humidity_factor`](@ref) throttles the reaction — self-desiccation, which is
what the `α_max = w/c / 0.42` of [Powers1948](@citet) describes empirically and
what [`powers_alpha_max`](@ref) otherwise supplies as an input.

The humidity is the water activity the `retention` law returns at the current
degree of saturation of the pore space,
`S = V_liquid / (V_ref − V_solid)`, both volumes recomputed from the composition
the rate law is handed. `V_ref` is the fresh paste's total volume — the same
reference the two-argument [`porosity`](@ref) uses, and the same sealed-curing
convention: the volume the reactions destroy stays inside as empty porosity.

!!! warning "This is where the arrest comes from, and it is kinetic"
    The Kelvin term is far too small to arrest hydration thermodynamically.
    Measured on a CEM I paste, imposing a water activity anywhere from 0.95 down
    to 0.05 leaves the equilibrium assemblage unchanged: the shift is
    `RT ln a_w = −553 J/mol` of water at `a_w = 0.80`, worth about 1.8 kJ per mole
    of alite against a hydration Gibbs energy of order −100 kJ/mol. Nulling that
    would need `a_w ≈ 5e-6`, a Kelvin radius smaller than a water molecule.

    A real paste stops at 75–80 % RH because transport and nucleation stop, not
    because the reaction has become unfavorable. So the humidity belongs in the
    **rate law**, through `humidity_factor`, and [`CapillaryWater`](@ref) is what
    makes the water activity of the equilibrium state mean the same thing.

# Examples

```julia
h = PoreHumidity(retention, cs; reference = fresh)
rxn[:rate] = parrott_killoh_avrami(PK84_PARAMS_C3S, "C3S"; humidity = h)
```

See also: [`humidity_factor`](@ref), [`WaterRetention`](@ref),
[`CapillaryWater`](@ref), [`powers_alpha_max`](@ref).
"""
struct PoreHumidity{R <: WaterRetention, T <: Real}
    retention::R
    V̄::Vector{T}            # standard molar volume per species, m³/mol
    idx_liquid::Vector{Int}
    idx_solid::Vector{Int}
    V_ref::T                # m³
    V_m_w::T                # the solvent's molar volume, for the Kelvin conversion
    T_K::T
end

function PoreHumidity(
        retention::WaterRetention, system::ChemicalSystem;
        reference::ChemicalState, T = temperature(reference),
    )
    missing_V = missing_molar_volumes(reference)
    isempty(missing_V) || throw(
        ArgumentError(
            "PoreHumidity needs a volume balance it can trust, and these species " *
                "are present in the reference with no standard molar volume: " *
                join(missing_V, ", ") * ". They contribute zero to the pore volume " *
                "in silence, so the saturation — and with it the humidity — would " *
                "be wrong."
        )
    )
    isempty(system.idx_solvent) && throw(
        ArgumentError("PoreHumidity needs an aqueous solvent to compute a humidity for.")
    )
    P = pressure(reference)
    V̄ = Float64[
        _has_molar_volume(sp) ? ustrip(us"m^3/mol", _molar_volume(sp)(T = T, P = P; unit = true)) : 0.0
            for sp in system.species
    ]
    V_ref = ustrip(us"m^3", volume(reference).total)
    V_ref > 0 || throw(ArgumentError("PoreHumidity: the reference volume is zero."))
    return PoreHumidity(
        retention, V̄, collect(system.idx_aqueous), collect(system.idx_crystal),
        V_ref, V̄[only(system.idx_solvent)], ustrip(us"K", T),
    )
end

"""
    pore_saturation(h::PoreHumidity, n) -> Real

Degree of saturation of the pore space at composition `n`, `V_liquid / V_pore`
with `V_pore = V_ref − V_solid`. Clamped to `[0, 1]`.
"""
function pore_saturation(h::PoreHumidity, n::AbstractVector)
    V_liq = zero(eltype(n))
    for i in h.idx_liquid
        V_liq += h.V̄[i] * n[i]
    end
    V_sol = zero(eltype(n))
    for i in h.idx_solid
        V_sol += h.V̄[i] * n[i]
    end
    V_pore = h.V_ref - V_sol
    V_pore > 0 || return one(eltype(n))
    return clamp(V_liq / V_pore, zero(eltype(n)), one(eltype(n)))
end

function (h::PoreHumidity)(n::AbstractVector)
    S = pore_saturation(h, n)
    # A saturated pore is the singular point of every retention law of the van
    # Genuchten family: `dp_c/dS` is **unbounded** as `S → 1`, so a rate law
    # differentiated there — which is what an implicit ODE solver does at the
    # very first step, a fresh paste being saturated — inherits an infinite
    # Jacobian entry and the integration cannot start at all. Measured: the
    # gradient of this function comes back with an infinite norm at `S = 1`.
    #
    # The value at saturation is not in doubt: water held at zero suction has
    # `a_w = 1`. So it is returned as a constant, branching on the primal the
    # way `_hkf_sigma` and `_pitzer_g` do, which sets the derivative to zero
    # there instead of to infinity.
    #
    # That is a **regularization and not an identity** — the true derivative is
    # unbounded, not zero. It applies only within `1e-10` of full saturation,
    # which in practice is the initial condition alone: the rate at `S = 1` is
    # unhindered, so water is consumed and the solver leaves the singular point
    # on its first successful step, after which the real derivative is used. It
    # is large just below saturation (of order `(1-S)^(-m)`) and that is a
    # stiffness, which is what the implicit solver is for.
    _primal(S) >= 1 - 1.0e-10 && return one(eltype(n))
    return _retention_activity(h.retention, S, h.V_m_w, h.T_K)
end

# The third `_humidity_at` method lives here rather than beside the other two:
# Julia needs `PoreHumidity` to exist when the method is defined, and the other
# two are declared long before this type is.
@inline _humidity_at(h::PoreHumidity, _t, n) = h(n.data)

# The two water/cement ratios at which Powers (1948) has a paste hydrate
# completely, sealed and with curing water supplied from outside. Their
# difference is the chemical shrinkage -- the volume the reaction loses because
# the hydrates are denser than the reagents. Sealed, that volume empties into the pore space and the paste
# desiccates itself; immersed, it is refilled from the bath, so the same paste
# reaches full hydration from a lower mixing water content.
#
# The two ratios are read from `data/literature/Powers1948.json`, where they are
# kept with their source, rather than typed here.
const POWERS_W_SEALED = literature_value("Powers1948", "w_c_sealed")
const POWERS_W_SATURATED = literature_value("Powers1948", "w_c_saturated")

"""
    powers_alpha_max(w_c; curing = :sealed) -> Real

The upper bound of [Powers1948](@citet) on the degree of hydration set by the availability of
water, `α_max = min(1, w/c / k)`: a paste below `w/c = k` cannot hydrate
completely, `k` being $(POWERS_W_SEALED) sealed or $(POWERS_W_SATURATED) water-cured, according to `curing`.

The 0.42 is **not** a stoichiometric demand, and reading it as one leads to the
wrong conclusion about what a Gibbs minimization should return. It is about
0.23 g of *non-evaporable* water per gram of cement — the water written into the
hydrate formulae, which is a mass balance — plus about 0.19 g of **gel water**
held in the C-S-H gel pores, which is physically present and chemically
unavailable. In a sealed paste hydration stops by self-desiccation with water
still in the specimen, so this bound is a statement about transport and access,
not about thermodynamics: an equilibrium calculation on the same mix consumes all
the clinker well below 0.42, and only runs out of water near the stoichiometric
demand.

# The two curing conventions

`curing = :sealed` is a specimen that exchanges nothing with its surroundings —
the convention of [`porosity`](@ref) and of the [w/c example](@ref sec-wc-ratio).
`curing = :saturated` is a specimen kept under water after setting, free to draw
in what the chemical shrinkage empties; the bound is then 0.36 and a mix that
would arrest sealed can go on reacting. **Neither is a property of the cement**:
they are two boundary conditions on the same paste, and which one applies is the
caller's to state.

Where water is abundant — `w/c` above the coefficient, or a cure that keeps
supplying it — the bound is 1 and this function stops doing anything, which is
the correct answer rather than a degenerate case: nothing about water is then
limiting the reaction.

Pass the result as the `α_max` keyword of [`parrott_killoh`](@ref),
[`parrott_killoh_avrami`](@ref) or [`waller`](@ref).

!!! note "It is a ceiling, not a schedule"
    `α_max` says how far the reaction can go, never how far it has got. At an
    early age the degree of reaction is set by the kinetics and is far below this
    bound; the bound binds only at long times, and only for the constituents
    whose kinetics would otherwise have taken them past it. For a constituent
    that reacts slowly — a slag, and a fly ash still more — the binding limit at
    28 days is its own dissolution rate, not the water. Take the **smaller** of
    the two.

# Examples

```jldoctest
julia> powers_alpha_max(0.5)
1.0

julia> round(powers_alpha_max(0.32); digits = 4)
0.7619

julia> round(powers_alpha_max(0.32; curing = :saturated); digits = 4)
0.8889
```
"""
function powers_alpha_max(w_c::Real; curing::Symbol = :sealed)
    k = if curing === :sealed
        POWERS_W_SEALED
    elseif curing === :saturated
        POWERS_W_SATURATED
    else
        throw(ArgumentError("curing must be :sealed or :saturated, got :$curing"))
    end
    return min(one(w_c), w_c / oftype(w_c, k))
end

# The author of the 1984 model is L. J. Parrott; the functions were first
# exported under a misspelling of the name, which keeps working.
@deprecate parrot_killoh(args...; kwargs...) parrott_killoh(args...; kwargs...)
@deprecate parrot_killoh_avrami(args...; kwargs...) parrott_killoh_avrami(args...; kwargs...)
