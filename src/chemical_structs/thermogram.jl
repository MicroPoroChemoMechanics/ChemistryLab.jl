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

When a phase releases what it releases, as a logistic step in temperature.

```math
f(T) = \\frac{1}{1 + \\exp\\!\\left(-\\dfrac{T - T_{1/2}}{w}\\right)}
```

`f` is the **fraction already released** at temperature `T`, so the phase's
contribution to a thermogram is `m_i f(T)` and to its derivative `m_i f'(T)`.
A logistic rather than a step because a decomposition is not instantaneous, and
rather than something with more shape parameters because two — a midpoint and a
width — are already at the edge of what one peak in a thermogram determines.

# Fields

  - `phase`: the species symbol the window belongs to, as the system spells it.
  - `midpoint`: `T₁/₂` in **kelvin**, where half of the release has happened.
  - `width`: `w` in kelvin. The release runs from roughly `T₁/₂ − 3w` to
    `T₁/₂ + 3w`, so a peak 100 K wide has `w ≈ 17 K`.
  - `releases`: `:water` or `:carbon_dioxide`.

Both parameters are [`Traced`](@ref), which is the point rather than a
decoration: a window read from a paper and a window guessed to get a picture on
the screen are the same two numbers and are not the same claim.

See also: [`thermogram`](@ref), [`window_parameters`](@ref).
"""
struct DecompositionWindow{T <: Real}
    phase::String
    midpoint::Traced{T}
    width::Traced{T}
    releases::Symbol
end

"""
    DecompositionWindow(phase, midpoint, width; releases = :water,
                        kind = PROV_UNSTATED, source = "") -> DecompositionWindow

Build a window. `midpoint` and `width` are in kelvin, as plain numbers or as
[`Traced`](@ref) values that keep their own provenance.

Passing bare numbers with no `kind` leaves them `PROV_UNSTATED`, which is the
weakest claim there is — deliberately, so a window nobody sourced never
strengthens a result.
"""
function DecompositionWindow(
        phase::AbstractString, midpoint, width;
        releases::Symbol = :water, kind::ProvenanceKind = PROV_UNSTATED,
        source::AbstractString = "",
    )
    releases in (:water, :carbon_dioxide) || throw(
        ArgumentError(
            "a window releases :water or :carbon_dioxide; got :$releases. " *
                "Those are the two `ignition_loss` accounts for."
        ),
    )
    m = midpoint isa Traced ? midpoint : Traced(float(midpoint), kind, source)
    w = width isa Traced ? width : Traced(float(width), kind, source)
    value(w) > 0 || throw(
        ArgumentError("a window's width must be positive; got $(value(w)) K."),
    )
    v = promote(value(m), value(w))
    # `m.source` and not `source(m)`: the keyword argument of this constructor is
    # itself named `source`, so inside here that name is a `String` and calling
    # it is the error it sounds like. The field is unambiguous.
    return DecompositionWindow{eltype(v)}(
        String(phase),
        Traced(convert(eltype(v), value(m)), provenance(m), m.source),
        Traced(convert(eltype(v), value(w)), provenance(w), w.source),
        releases,
    )
end

function Base.show(io::IO, w::DecompositionWindow)
    return print(
        io, "DecompositionWindow(", w.phase, ", ", value(w.midpoint), " K ± ",
        value(w.width), " K, ", w.releases, ", ",
        _prov_label(weakest(w.midpoint, w.width)), ")",
    )
end

"""
    released_fraction(w::DecompositionWindow, T) -> Real

The fraction of `w`'s phase already released at temperature `T` in kelvin.
"""
released_fraction(w::DecompositionWindow, T) =
    1 / (1 + exp(-(T - value(w.midpoint)) / value(w.width)))

"""
    released_rate(w::DecompositionWindow, T) -> Real

`df/dT` — the shape one peak of a DTG curve has, per kelvin.
"""
function released_rate(w::DecompositionWindow, T)
    f = released_fraction(w, T)
    return f * (1 - f) / value(w.width)
end

"""
    thermogram(state, windows; temperatures) -> NamedTuple

A thermogravimetric curve for `state` under `windows`, as
`(; temperature, mass, loss, dtg, by_phase)`:

  - `temperature`: the grid, in kelvin, as given.
  - `mass`: the sample mass remaining, in kilograms — the **solid** mass of
    `state`, residue included, because that is what sits on the pan.
  - `loss`: what has been released, in kilograms — `mass[1] - mass`.
  - `dtg`: `-dm/dT` in kilograms per kelvin, which is the curve a
    thermogravimetric analysis actually resolves peaks in.
  - `by_phase`: the same loss, per phase, so a peak can be attributed.

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
        temperatures,
    )
    system = state.system
    per_phase = bound_water_per_phase(state)
    water = Dict(p.first => p.second for p in per_phase)
    co2 = _co2_per_phase(state)

    T = collect(float.(temperatures))
    ET = promote_type(
        eltype(T), Float64,
        (typeof(value(w.midpoint)) for w in windows)...,
    )
    # The sample mass is the SOLID mass, residue included — a thermogram plots
    # what is on the pan, not only the part that will leave it.
    m0 = ustrip(us"kg", mass(state).solid)

    by_phase = Dict{String, Vector{ET}}()
    loss = zeros(ET, length(T))
    dtg = zeros(ET, length(T))
    for w in windows
        pool = w.releases === :water ? water : co2
        haskey(pool, w.phase) || continue
        m = ustrip(us"kg", pool[w.phase])
        curve = [m * released_fraction(w, t) for t in T]
        by_phase[w.phase] = get(by_phase, w.phase, zeros(ET, length(T))) .+ curve
        loss .+= curve
        dtg .+= [m * released_rate(w, t) for t in T]
    end
    return (;
        temperature = T,
        mass = m0 .- loss,
        loss = loss,
        dtg = dtg,
        by_phase = by_phase,
    )
end

"""
    _co2_per_phase(state) -> Dict{String, Quantity}

The carbon dioxide each solid species would release, by symbol.
"""
function _co2_per_phase(state::ChemicalState)
    system = state.system
    mc = _ignition_molar_mass(system, "CO2@", "CO2", 0.0440095u"kg/mol")
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
    phases_without_windows(state, windows) -> Vector{String}

Which phases of `state` would release something on heating and have no window
saying when.

A thermogram computed with one of these present integrates to less than
[`ignition_loss`](@ref) says, and the difference is silent. This is how to see
it before reading a curve.
"""
function phases_without_windows(
        state::ChemicalState, windows::AbstractVector{<:DecompositionWindow},
    )
    covered = Set(w.phase for w in windows)
    out = String[]
    for p in bound_water_per_phase(state)
        p.first in covered || push!(out, p.first)
    end
    for k in keys(_co2_per_phase(state))
        (k in covered || k in out) || push!(out, k)
    end
    return out
end

"""
    window_parameters(windows) -> (θ, names)

The windows' parameters as one vector, `[T₁/₂, w]` per window, with names —
what an optimizer and [`identifiability`](@ref) take.

The two are returned together because a parameter vector whose entries are not
named is a parameter vector nobody can report.
"""
function window_parameters(windows::AbstractVector{<:DecompositionWindow})
    θ = Float64[]
    names = String[]
    for w in windows
        push!(θ, value(w.midpoint)); push!(names, "T½($(w.phase))")
        push!(θ, value(w.width)); push!(names, "w($(w.phase))")
    end
    return θ, names
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
                "windows (a midpoint and a width each); got $(length(θ))."
        ),
    )
    return [
        DecompositionWindow(
            w.phase, θ[2i - 1], θ[2i];
            releases = w.releases, kind = kind, source = source,
        ) for (i, w) in enumerate(windows)
    ]
end
