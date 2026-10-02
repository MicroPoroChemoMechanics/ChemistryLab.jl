# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using DynamicQuantities

# ── From a total mass loss to a thermogram ───────────────────────────────────
#
# `ignition_loss` gives what an assemblage loses. A thermogravimetric curve says
# WHEN it loses it, and that is the part that identifies phases: C-S-H, AFt and
# AFm all release below 200 °C and are told apart by the shape of the release,
# not by its total.
#
# The windows are not consequences of the formulas. They are either taken from a
# publication or **identified from a measured thermogram**, and this file exists
# so that both are possible and neither is silent about which it is: every
# parameter is a `Traced`, so a window standing in for one nobody measured
# prints as a placeholder for as long as it is one.

"""
    struct DecompositionWindow{T<:Real}

When a phase releases what it releases, as a step in temperature, in one of two
forms.

**A logistic step**, the form a peak fitted to a thermogram takes:

```math
f(T) = \\frac{1}{1 + \\exp\\!\\left(-\\dfrac{T - T_{1/2}}{w}\\right)}
```

A logistic rather than a step because a decomposition is not instantaneous, and
rather than something with more shape parameters because two — a midpoint and a
width — are already at the edge of what one peak in a thermogram determines.

**A temperature interval**, the form a thermogravimetric reading takes when it
attributes the mass lost between two temperatures to one phase ("the weight
loss between 350 and 500 °C" for portlandite): the whole release happens between
`T₁` and `T₂`, and none outside,

```math
f(T) = 3t^2 - 2t^3, \\qquad t = \\operatorname{clamp}\\!\\left(\\frac{T - T_1}{T_2 - T_1},\\ 0,\\ 1\\right),
```

a smooth step whose rate vanishes at both ends, so that the curve and its
derivative are continuous and the loss between `T₁` and `T₂` is exactly the
phase's content. The logistic has tails: a phase centered in an interval loses
part of its release outside it.

In either form `f` is the **fraction already released** at temperature `T`, so
the phase's contribution to a thermogram is `m_i f(T)` and to its derivative
`m_i f'(T)`.

# Fields

  - `phase`: the species symbol the window belongs to, as the system spells it.
  - `midpoint`: `T₁/₂` in **kelvin**, where half of the release has happened
    (for an interval, its center).
  - `width`: `w` in kelvin. The logistic release runs from roughly `T₁/₂ − 3w`
    to `T₁/₂ + 3w`, so a peak 100 K wide has `w ≈ 17 K`; for an interval, its
    half-span, `(T₂ − T₁)/2`.
  - `releases`: `:water` or `:carbon_dioxide`.
  - `fraction`: how much of that phase's release this window accounts for, `1`
    by default. **A phase can go in stages** — gypsum loses its two waters in
    two steps, `CaSO₄·2H₂O → CaSO₄·½H₂O → CaSO₄` — and one window per stage with
    fractions summing to one is how that is written. The sum is checked, in
    [`thermogram`](@ref), against each phase and product it is given for.
  - `shape`: `:logistic` or `:interval`.

Both temperature parameters are [`Traced`](@ref), which is the point rather than a
decoration: a window read from a paper and a window guessed to get a picture on
the screen are the same two numbers and are not the same claim.

See also: [`thermogram`](@ref), [`window_parameters`](@ref), [`window_interval`](@ref).
"""
struct DecompositionWindow{T <: Real}
    phase::String
    midpoint::Traced{T}
    width::Traced{T}
    releases::Symbol
    fraction::T
    shape::Symbol
end

"""
    DecompositionWindow(phase, midpoint, width; releases = :water, fraction = 1,
                        kind = PROV_UNSTATED, source = "") -> DecompositionWindow
    DecompositionWindow(phase; between = (T₁, T₂), releases = :water, fraction = 1,
                        kind = PROV_UNSTATED, source = "") -> DecompositionWindow

Build a window: a logistic step from its `midpoint` and `width`, or a
temperature interval from its two ends `between`. Temperatures are in kelvin, as
plain numbers, as quantities, or as [`Traced`](@ref) values that keep their own
provenance (an interval read with [`literature_table`](@ref) carries its
source).

Passing bare numbers with no `kind` leaves them `PROV_UNSTATED`, which is the
weakest claim there is — deliberately, so a window nobody sourced never
strengthens a result.
"""
function DecompositionWindow(
        phase::AbstractString; between,
        releases::Symbol = :water, fraction::Real = 1,
        kind::ProvenanceKind = PROV_UNSTATED, source::AbstractString = "",
    )
    length(between) == 2 || throw(
        ArgumentError("`between` takes the two ends of the interval, (T₁, T₂); got $(length(between)) values."),
    )
    lo, hi = (_traced_kelvin(t, kind, source) for t in between)
    value(lo) < value(hi) || throw(
        ArgumentError("an interval needs T₁ < T₂; got ($(value(lo)) K, $(value(hi)) K)."),
    )
    # The center and the half-span, with the weaker of the two ends' standing.
    prov = weakest(lo, hi)
    src = lo.source == hi.source ? lo.source : join(filter(!isempty, [lo.source, hi.source]), "; ")
    mid = Traced((value(lo) + value(hi)) / 2, prov, src)
    half = Traced((value(hi) - value(lo)) / 2, prov, src)
    return _decomposition_window(phase, mid, half, releases, fraction, :interval)
end

# A temperature as a `Traced` value in kelvin: a number is kelvin, a quantity is
# converted, a `Traced` keeps its standing.
_traced_kelvin(t::Traced, kind, source) = Traced(float(_kelvin(value(t))), provenance(t), t.source)
_traced_kelvin(t, kind, source) = Traced(float(_kelvin(t)), kind, source)

function DecompositionWindow(
        phase::AbstractString, midpoint, width;
        releases::Symbol = :water, fraction::Real = 1,
        kind::ProvenanceKind = PROV_UNSTATED, source::AbstractString = "",
    )
    m = midpoint isa Traced ? midpoint : Traced(float(midpoint), kind, source)
    w = width isa Traced ? width : Traced(float(width), kind, source)
    return _decomposition_window(phase, m, w, releases, fraction, :logistic)
end

function _decomposition_window(phase, m, w, releases, fraction, shape)
    releases in (:water, :carbon_dioxide) || throw(
        ArgumentError(
            "a window releases :water or :carbon_dioxide; got :$releases. " *
                "Those are the two `ignition_loss` accounts for."
        ),
    )
    value(w) > 0 || throw(
        ArgumentError("a window's width must be positive; got $(value(w)) K."),
    )
    0 < fraction <= 1 || throw(
        ArgumentError(
            "a window accounts for a fraction in (0, 1] of its phase's release; " *
                "got $fraction. Use several windows to split a release in stages."
        ),
    )
    v = promote(value(m), value(w), float(fraction))
    # `m.source` and not `source(m)`: the keyword argument of this constructor is
    # itself named `source`, so inside here that name is a `String` and calling
    # it is the error it sounds like. The field is unambiguous.
    return DecompositionWindow{eltype(v)}(
        String(phase),
        Traced(convert(eltype(v), value(m)), provenance(m), m.source),
        Traced(convert(eltype(v), value(w)), provenance(w), w.source),
        releases,
        convert(eltype(v), fraction),
        shape,
    )
end

"""
    window_interval(w::DecompositionWindow) -> (T₁, T₂)

The two ends, in kelvin, of an interval window; for a logistic one, the
temperatures at which 1 % and 99 % of its release have happened.
"""
function window_interval(w::DecompositionWindow)
    m, h = value(w.midpoint), value(w.width)
    w.shape === :interval && return (m - h, m + h)
    return (m - h * log(99), m + h * log(99))
end

function Base.show(io::IO, w::DecompositionWindow)
    if w.shape === :interval
        T1, T2 = window_interval(w)
        print(io, "DecompositionWindow(", w.phase, ", between ", T1, " K and ", T2, " K, ", w.releases)
    else
        print(
            io, "DecompositionWindow(", w.phase, ", ", value(w.midpoint), " K ± ",
            value(w.width), " K, ", w.releases,
        )
    end
    isone(w.fraction) || print(io, ", ", round(100 * w.fraction; digits = 1), " %")
    return print(io, ", ", _prov_label(weakest(w.midpoint, w.width)), ")")
end

"""
    released_fraction(w::DecompositionWindow, T) -> Real

The fraction of `w`'s phase already released at temperature `T` in kelvin.
"""
function released_fraction(w::DecompositionWindow, T)
    w.shape === :interval || return 1 / (1 + exp(-(T - value(w.midpoint)) / value(w.width)))
    t = _interval_position(w, T)
    return t * t * (3 - 2t)
end

# Where `T` lies in an interval window, from 0 at `T₁` to 1 at `T₂`.
_interval_position(w, T) = clamp((T - value(w.midpoint)) / (2 * value(w.width)) + 1 // 2, 0, 1)

"""
    released_rate(w::DecompositionWindow, T) -> Real

`df/dT` — the shape one peak of a DTG curve has, per kelvin.
"""
function released_rate(w::DecompositionWindow, T)
    if w.shape === :interval
        t = _interval_position(w, T)
        return 6 * t * (1 - t) / (2 * value(w.width))
    end
    f = released_fraction(w, T)
    return f * (1 - f) / value(w.width)
end

"""
    thermogram(state, windows; temperatures, relative_to = :initial) -> NamedTuple

A thermogravimetric curve for `state` under `windows`, as
`(; temperature, mass, loss, dtg, by_phase, reference_mass, mass_percent, loss_percent)`:

  - `temperature`: the grid, in kelvin, as given.
  - `mass`: the sample mass remaining, in kilograms — the **solid** mass of
    `state`, residue included, because that is what sits on the pan.
  - `loss`: what has been released, in kilograms, counted from nothing rather
    than from the first point of the grid — so `mass = mass₀ - loss` with
    `mass₀` the solid mass, and `loss[1]` is the tail already gone before the
    grid begins rather than zero by definition.
  - `dtg`: `-dm/dT` in kilograms per kelvin, which is the curve a
    thermogravimetric analysis actually resolves peaks in.
  - `by_phase`: the same loss, per phase, so a peak can be attributed.
  - `reference_mass`, `mass_percent`, `loss_percent`: the curve as published
    thermograms give it, the mass and the mass lost since the first temperature
    of the grid, in percent of a reference mass, which `relative_to` names:
      - `:initial`, the sample at the first temperature of the grid;
      - a temperature (kelvin, or a quantity), the sample at that temperature: a
        dry mass ([Scholer2015](@cite) take the weight at 500 °C,
        [Shi2016](@cite) at 800 °C) or an ignited one ([ShiLothenbach2020](@cite)
        give bound water in percent of the sample ignited at 980 °C);
      - `:ignited`, the sample once every window has released.
    See [the thermogravimetry page](@ref sec-example-tga) for the conventions.

A window whose phase `state` has nothing to release from contributes nothing and
raises nothing — [`windows_without_phases`](@ref) is how a typo in a phase name
is caught rather than read off a missing peak.

Every phase of `state` that carries hydrogen or carbon and has **no window**
contributes its mass to the starting point and never leaves — which is the
honest behavior, because a phase nobody said when to decompose has not been
modeled. [`phases_without_windows`](@ref) is how to find out that this happened
rather than to discover it from a curve that integrates to the wrong total.

!!! warning "A logistic has infinite tails"
    The curve does not start at exactly zero: a window centered at 400 K with a
    width of 15 K is already 1.3 ‰ through at 300 K, and that mass is counted as
    released before the grid begins. It is the model rather than an error — no
    renormalization happens here, because silently rescaling a curve to make it
    start at zero would put the discrepancy somewhere a reader cannot see.

    Start the grid at least `6w` below the lowest midpoint and the leak is under
    0.3 %; `10w` puts it under 5 × 10⁻⁵. `thermogram(...).loss[1]` is what it
    actually is.

# What this is

A **deconvolution model**, not a kinetic one: it says what fraction of a phase
has gone at a temperature, not how fast it goes at a heating rate. That is the
form published TGA interpretations take, and it is what makes the windows
identifiable from one curve — a rate model would need several heating rates
before its parameters meant anything.
"""
function thermogram(
        state::ChemicalState, windows::AbstractVector{<:DecompositionWindow};
        temperatures, relative_to = :initial,
    )
    _check_fractions(windows)          # fail before doing any of the work

    per_phase = bound_water_per_phase(state)
    water = Dict(p.first => p.second for p in per_phase)
    co2 = _co2_per_phase(state)

    T = collect(float.(temperatures))
    # The number type of everything the curves are computed from: the
    # temperatures, the state (a composition being differentiated), and the
    # midpoints, widths and fractions of the windows (fitted to a curve).
    ET = promote_type(
        eltype(T), Float64, _realtype(eltype(state.n)),
        (typeof(value(w.midpoint)) for w in windows)...,
        (typeof(value(w.width)) for w in windows)...,
        (typeof(float(w.fraction)) for w in windows)...,
    )
    # The sample mass is the SOLID mass, residue included — a thermogram plots
    # what is on the pan, not only the part that will leave it.
    m0 = ustrip(us"kg", mass(state).solid)

    by_phase = Dict{String, Vector{ET}}()
    loss = zeros(ET, length(T))
    dtg = zeros(ET, length(T))
    released = Tuple{DecompositionWindow, ET}[]        # each window and its mass
    for w in windows
        pool = w.releases === :water ? water : co2
        haskey(pool, w.phase) || continue
        m = ustrip(us"kg", pool[w.phase]) * w.fraction
        push!(released, (w, m))
        curve = [m * released_fraction(w, t) for t in T]
        by_phase[w.phase] = get(by_phase, w.phase, zeros(ET, length(T))) .+ curve
        loss .+= curve
        dtg .+= [m * released_rate(w, t) for t in T]
    end
    mass_curve = m0 .- loss
    # The reference mass, exactly at its temperature rather than read off the grid.
    ref = if relative_to === :initial
        isempty(T) ? m0 : mass_curve[1]
    elseif relative_to === :ignited
        m0 - sum((m for (_, m) in released); init = zero(ET))
    elseif relative_to isa Union{Real, DynamicQuantities.AbstractQuantity}
        Tr = _kelvin(relative_to)
        m0 - sum((m * released_fraction(w, Tr) for (w, m) in released); init = zero(ET))
    else
        throw(
            ArgumentError(
                "`relative_to` is :initial, :ignited or a temperature; got $(repr(relative_to))."
            ),
        )
    end
    start = isempty(T) ? m0 : mass_curve[1]
    return (;
        temperature = T,
        mass = mass_curve,
        loss = loss,
        dtg = dtg,
        by_phase = by_phase,
        reference_mass = ref,
        mass_percent = 100 .* mass_curve ./ ref,
        loss_percent = 100 .* (start .- mass_curve) ./ ref,
    )
end

"""
    _check_fractions(windows)

Refuse a set of windows whose fractions do not sum to one for some phase and
product.

Two windows on one phase with the default fraction would release its mass
**twice**, and the only symptom is a curve that integrates to more than
[`ignition_loss`](@ref) — a silent doubling rather than an error. The check is
here rather than at construction because a window does not know what it will be
used with.
"""
function _check_fractions(windows::AbstractVector{<:DecompositionWindow})
    # A check, so on the values: a fraction being differentiated is still a
    # fraction, and the derivatives of a sum that must be one say nothing here.
    sums = Dict{Tuple{String, Symbol}, Float64}()
    for w in windows
        k = (w.phase, w.releases)
        sums[k] = get(sums, k, 0.0) + _plain(float(w.fraction))
    end
    for (k, total) in sums
        isapprox(total, 1.0; atol = 1.0e-8) && continue
        throw(
            ArgumentError(
                "the windows for $(k[1]) releasing $(k[2]) have fractions summing " *
                    "to $total, not 1. A phase that goes in stages needs one " *
                    "window per stage with the fractions splitting its release; " *
                    "two windows both accounting for all of it would release its " *
                    "mass twice."
            ),
        )
    end
    return nothing
end

"""
    _co2_per_phase(state) -> Dict{String, Quantity}

The carbon dioxide each solid species would release, by symbol.
"""
function _co2_per_phase(state::ChemicalState)
    system = state.system
    mc = _ignition_molar_mass(system, "CO2@", "CO2", _CO2_ATOMS)
    out = Dict{String, typeof(uconvert(us"kg", 1.0u"mol" * mc))}()
    for i in _solid_indices(system)
        sp = system.species[i]
        c = get(atoms(sp), :C, 0)
        iszero(c) && continue
        m = uconvert(us"kg", c * state.n[i] * mc)
        ustrip(us"kg", m) > 0 || continue
        out[symbol(sp)] = m
    end
    return out
end

"""
    phases_without_windows(state, windows) -> Vector{Pair{String,Symbol}}

What `state` would release on heating and has no window saying when, as
`phase => product` pairs.

A thermogram computed with one of these present integrates to less than
[`ignition_loss`](@ref) says, and the difference is silent. This is how to see
it before reading a curve.

# Why the product is in the answer

**A phase can need two windows.** A carbonated hydrate carries hydrogen and
carbon, releases water and carbon dioxide, and does so at different
temperatures — a hemicarboaluminate is not an exotic case in a cement. Coverage
counted per phase rather than per `phase => product` reports such a phase as
covered when only half of it is, and the missing half is precisely the silent
shortfall this exists to catch.
"""
function phases_without_windows(
        state::ChemicalState, windows::AbstractVector{<:DecompositionWindow},
    )
    covered = Set((w.phase, w.releases) for w in windows)
    out = Pair{String, Symbol}[]
    for p in bound_water_per_phase(state)
        (p.first, :water) in covered || push!(out, p.first => :water)
    end
    for k in sort!(collect(keys(_co2_per_phase(state))))
        (k, :carbon_dioxide) in covered || push!(out, k => :carbon_dioxide)
    end
    return out
end

"""
    windows_without_phases(state, windows) -> Vector{String}

Which windows name a phase `state` has nothing to release from.

The mirror of [`phases_without_windows`](@ref), and the one that catches a
**typo**: a window on `"Portlandit"` contributes nothing, raises nothing, and
leaves a curve that is simply missing a peak. Both directions are silent by
themselves, which is why both are reported rather than one.

A window on a phase that is genuinely absent — a paste with no calcite — is not
an error and shows up here too; the two are indistinguishable from inside, and
naming them is all this can honestly do.
"""
function windows_without_phases(
        state::ChemicalState, windows::AbstractVector{<:DecompositionWindow},
    )
    water = Set(p.first for p in bound_water_per_phase(state))
    co2 = Set(keys(_co2_per_phase(state)))
    return [
        w.phase for w in windows
            if !(w.phase in (w.releases === :water ? water : co2))
    ]
end

"""
    window_parameters(windows) -> (θ, names)

The windows' parameters as one vector, `[T₁/₂, w]` per logistic window and
`[T₁, T₂]` per interval, with names — what an optimizer and
[`identifiability`](@ref) take.

The two are returned together because a parameter vector whose entries are not
named is a parameter vector nobody can report.
"""
function window_parameters(windows::AbstractVector{<:DecompositionWindow})
    θ = Any[]
    names = String[]
    for w in windows
        if w.shape === :interval
            T1, T2 = window_interval(w)
            push!(θ, T1); push!(names, "T₁($(w.phase))")
            push!(θ, T2); push!(names, "T₂($(w.phase))")
        else
            push!(θ, value(w.midpoint)); push!(names, "T½($(w.phase))")
            push!(θ, value(w.width)); push!(names, "w($(w.phase))")
        end
    end
    # In the number type of the windows, which a nested fit makes dual.
    return _promoted(θ), names
end

"""
    with_window_parameters(windows, θ; kind = PROV_FITTED, source = "") -> Vector

`windows` with their parameters replaced by `θ`, in the order
[`window_parameters`](@ref) produced.

The replacements are marked `PROV_FITTED` by default, because that is what a
number coming out of an optimizer is. Pass `kind` explicitly when it is not —
stepping a parameter to differentiate against it, for instance, does not make
the result a fit.
"""
function with_window_parameters(
        windows::AbstractVector{<:DecompositionWindow}, θ;
        kind::ProvenanceKind = PROV_FITTED, source::AbstractString = "",
    )
    length(θ) == 2 * length(windows) || throw(
        ArgumentError(
            "expected $(2 * length(windows)) parameters for $(length(windows)) " *
                "windows (two each); got $(length(θ))."
        ),
    )
    return [
        w.shape === :interval ?
            DecompositionWindow(
                w.phase; between = (θ[2i - 1], θ[2i]),
                releases = w.releases, fraction = w.fraction, kind = kind, source = source,
            ) :
            DecompositionWindow(
                w.phase, θ[2i - 1], θ[2i];
                releases = w.releases, fraction = w.fraction,
                kind = kind, source = source,
            ) for (i, w) in enumerate(windows)
    ]
end
