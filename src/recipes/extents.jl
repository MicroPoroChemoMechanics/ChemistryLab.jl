# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# ── How much of a constituent has reacted ─────────────────────────────────────

"""
    abstract type AbstractExtent end

The degree of reaction of a material or of one of its constituents, as a
function of time: `extent(e, t)` is a number between 0 and 1, `t` a time
(a `Quantity`, or a number of days). A [`Recipe`](@ref) multiplies the extent of
a material by the extent of each of its constituents; what has not reacted is
kept aside as the residue, with its mass and its volume.

Concrete extents: [`ConstantExtent`](@ref), [`TabulatedExtent`](@ref),
[`LogisticExtent`](@ref), [`ParrottKillohExtent`](@ref) and
[`CappedExtent`](@ref).
"""
abstract type AbstractExtent end

# Times are read in days; a bare number is a number of days.
_days(t::Real) = float(t)
_days(t::DynamicQuantities.AbstractQuantity) = ustrip(us"d", t)

"""
    ConstantExtent(value)

The same degree of reaction at every time; `ConstantExtent(1)` is complete
reaction, the extent of water and of salts.
"""
struct ConstantExtent{T <: Real} <: AbstractExtent
    value::T
    function ConstantExtent(value::Real)
        0 <= value <= 1 || throw(ArgumentError("ConstantExtent: a degree of reaction is in [0, 1]; got $value."))
        return new{typeof(value)}(value)
    end
end
extent(e::ConstantExtent, _ = nothing) = e.value

"""
    TabulatedExtent(times, values)

A measured degree of reaction, interpolated linearly in time between the
tabulated points and held at the first and last values outside them. `times`
are in days (or `Quantity`s), `values` between 0 and 1. The usual source is a
table of `data/literature`, for instance
`literature_table("Durdzinski2017", "degree_of_reaction"; material = "S1", technique = "SEM-IA")`.
"""
struct TabulatedExtent{T <: Real} <: AbstractExtent
    days::Vector{T}
    values::Vector{T}
    function TabulatedExtent(times::AbstractVector, values::AbstractVector)
        length(times) == length(values) >= 1 || throw(
            ArgumentError("TabulatedExtent: $(length(times)) times for $(length(values)) values.")
        )
        d = [_days(t) for t in times]
        v = float.(values)
        p = sortperm(d)
        d, v = d[p], v[p]
        all(x -> 0 <= x <= 1, v) || throw(ArgumentError("TabulatedExtent: degrees of reaction are in [0, 1]."))
        allunique(d) || throw(ArgumentError("TabulatedExtent: two values at one time."))
        T = promote_type(eltype(d), eltype(v))
        return new{T}(Vector{T}(d), Vector{T}(v))
    end
end
function extent(e::TabulatedExtent, t)
    d = _days(t)
    d <= first(e.days) && return first(e.values)
    d >= last(e.days) && return last(e.values)
    i = searchsortedlast(e.days, d)
    w = (d - e.days[i]) / (e.days[i + 1] - e.days[i])
    return (1 - w) * e.values[i] + w * e.values[i + 1]
end

"""
    LogisticExtent(; final, half_time, slope, initial = 0)

A logistic curve in the logarithm of time,

```math
\\alpha(t) = \\alpha_0 + \\frac{\\alpha_\\infty - \\alpha_0}{1 + (t_{1/2}/t)^{s}} ,
```

a generic shape to fit to measured degrees of reaction when a table is too
sparse to interpolate. `half_time` is the time at which half of the final gain
is reached (days, or a `Quantity`); nothing here is fitted by the package.
"""
struct LogisticExtent{T <: Real} <: AbstractExtent
    initial::T
    final::T
    half_time::T
    slope::T
end
function LogisticExtent(; final::Real, half_time, slope::Real, initial::Real = 0.0)
    (0 <= initial <= final <= 1) || throw(ArgumentError("LogisticExtent: need 0 ≤ initial ≤ final ≤ 1."))
    th = _days(half_time)
    (th > 0 && slope > 0) || throw(ArgumentError("LogisticExtent: half_time and slope must be positive."))
    return LogisticExtent(promote(float(initial), float(final), th, float(slope))...)
end
function extent(e::LogisticExtent, t)
    d = _days(t)
    d <= 0 && return e.initial
    return e.initial + (e.final - e.initial) / (1 + (e.half_time / d)^e.slope)
end

"""
    CappedExtent(inner, cap)

`inner`, never above `cap`. The cap of a clinker is usually Powers' water
limit, `CappedExtent(inner, powers_alpha_max(w_c))`.
"""
struct CappedExtent{E <: AbstractExtent, T <: Real} <: AbstractExtent
    inner::E
    cap::T
    function CappedExtent(inner::AbstractExtent, cap::Real)
        0 <= cap <= 1 || throw(ArgumentError("CappedExtent: the cap is a degree of reaction, in [0, 1]."))
        return new{typeof(inner), typeof(cap)}(inner, cap)
    end
end
extent(e::CappedExtent, t) = min(extent(e.inner, t), e.cap)

"""
    ParrottKillohExtent(phase; T = 293.15u"K", α_max = 1.0, blaine = nothing, w_c = nothing)

The degree of hydration of the clinker phase `phase` ("C3S", "C2S", "C3A" or
"C4AF") under the rate law of Parrott and Killoh (1984) in the form and with the
parameters [`parrott_killoh_avrami`](@ref) uses, at the constant temperature
`T`: the ordinary differential equation of that law integrated from zero, on a
logarithmic grid of time. `α_max` is the ceiling of the law (Powers' water limit,
for instance), `blaine` the fineness correction.

`w_c` applies instead the water/cement factor of Parrott and Killoh, as
Lothenbach and Winnefeld (2006, Section 4.1) state it: the rate is multiplied by

```math
f = \\begin{cases} 1 & \\alpha \\le 1.333\\, w/c \\\\
(1 + 4.444\\, w/c - 3.333\\, \\alpha)^4 & \\alpha > 1.333\\, w/c \\end{cases}
```

which is continuous at ``\\alpha = 1.333\\, w/c`` and stops the hydration at
``\\alpha = (1 + 4.444\\, w/c)/3.333``.
"""
struct ParrottKillohExtent <: AbstractExtent
    phase::String
    days::Vector{Float64}
    values::Vector{Float64}
end
function ParrottKillohExtent(
        phase::AbstractString; T = 293.15u"K", α_max::Real = 1.0, blaine = nothing,
        w_c = nothing, horizon_days::Real = 3650.0,
    )
    rate = parrott_killoh_avrami(_pk84_params(String(phase)), String(phase); α_max = α_max, blaine = blaine)
    TK = _days_free_temperature(T)
    # dα/dt, in 1/s, from the rate on one mole of the phase (n = 1 − α); the
    # positional call of a `KineticFunc` takes and returns bare SI numbers.
    ph = String(phase)
    fwc(α) = (w_c === nothing || α <= 1.333 * w_c) ? 1.0 : max(1 + 4.444 * w_c - 3.333 * α, 0.0)^4
    dα(α) = rate(TK, 1.0e5, 0.0, Dict(ph => 1 - α), nothing, Dict(ph => 1.0)) * fwc(α)
    grid = exp.(range(log(1.0e-4), log(horizon_days); length = 4001))   # days
    α = 0.0
    vals = zeros(length(grid))
    tprev = 0.0
    for (k, d) in enumerate(grid)
        h = (d - tprev) * 86400.0
        # Four stages of Runge-Kutta, the rate being smooth once the Avrami seed
        # has been left; α is held below the ceiling.
        k1 = dα(α)
        k2 = dα(min(α + h * k1 / 2, α_max - 1.0e-12))
        k3 = dα(min(α + h * k2 / 2, α_max - 1.0e-12))
        k4 = dα(min(α + h * k3, α_max - 1.0e-12))
        α = clamp(α + h * (k1 + 2k2 + 2k3 + k4) / 6, 0.0, α_max)
        vals[k] = α
        tprev = d
    end
    return ParrottKillohExtent(String(phase), grid, vals)
end
_days_free_temperature(T::Real) = float(T)
_days_free_temperature(T::DynamicQuantities.AbstractQuantity) = ustrip(us"K", T)
function extent(e::ParrottKillohExtent, t)
    d = _days(t)
    d <= first(e.days) && return first(e.values) * d / first(e.days)
    d >= last(e.days) && return last(e.values)
    i = searchsortedlast(e.days, d)
    w = (d - e.days[i]) / (e.days[i + 1] - e.days[i])
    return (1 - w) * e.values[i] + w * e.values[i + 1]
end

# A number or a function of time given where an extent is expected.
_as_extent(e::AbstractExtent) = e
_as_extent(x::Real) = ConstantExtent(x)
_as_extent(f::Function) = _FunctionExtent(f)
struct _FunctionExtent{F} <: AbstractExtent
    f::F
end
extent(e::_FunctionExtent, t) = e.f(t)
