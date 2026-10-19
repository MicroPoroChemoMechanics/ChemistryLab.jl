# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2026 Jean-François Barthélémy (Cerema, UMR MCD)

# The windows of the thermogram: the interval over which each release of the
# solids (the water of a phase, or its carbon dioxide) is spread. Each release
# takes the interval in which the cement literature reads it (ChemistryLab's
# `literature_windows`), unless the user sets it an interval or a logistic peak
# of their own. A release is keyed by its phase and its product
# (shared/releases.jl).

"""
    WindowChoice

The window the user sets for one release of the solids.
"""
abstract type WindowChoice <: Choice end

"The whole release between `from` and `to` °C, along a smooth step."
struct IntervalWindowChoice <: WindowChoice
    from::Float64
    to::Float64
end
kind(::Type{IntervalWindowChoice}) = "window_interval"
fields(::Type{IntervalWindowChoice}) = [
    Field("from", :number; default = 30.0, unit = "°C", min = 20, max = 1100),
    Field("to", :number; default = 550.0, unit = "°C", min = 20, max = 1100),
]

"A logistic peak: half of the release at `midpoint` °C, 98 % of it within 4.6 `width` on either side."
struct LogisticWindowChoice <: WindowChoice
    midpoint::Float64
    width::Float64
end
kind(::Type{LogisticWindowChoice}) = "window_logistic"
fields(::Type{LogisticWindowChoice}) = [
    Field("midpoint", :number; default = 450.0, unit = "°C", min = 20, max = 1100),
    Field("width", :number; default = 15.0, unit = "K", min = 0.1, max = 300),
]

members(::Type{WindowChoice}) = (IntervalWindowChoice, LogisticWindowChoice)

"""
    TgaChoice

The windows of the thermogram of a calculation.
"""
abstract type TgaChoice <: Choice end

"The windows the user sets, one per release; every other release keeps the literature's."
struct TgaWindowsChoice <: TgaChoice
    windows::Vector{Pair{String, WindowChoice}}
end
kind(::Type{TgaWindowsChoice}) = "tga_windows"
fields(::Type{TgaWindowsChoice}) = [Field("windows", :keyed; default = Any[], source = "window")]

members(::Type{TgaChoice}) = (TgaWindowsChoice,)

"What keeps the windows `t` from being drawn."
function tga_problems(t::TgaWindowsChoice)
    out = Problem[]
    for (key, w) in t.windows
        if release_of(key) === nothing
            push!(
                out, Problem(
                    "window_key", "\"$key\" names no release of the solids: a window is set for \"<phase> (water)\" or \"<phase> (CO2)\".",
                    Dict{String, Any}("key" => key),
                )
            )
        end
        append!(out, window_problems(key, w))
    end
    return out
end

window_problems(_, ::LogisticWindowChoice) = Problem[]
function window_problems(key, w::IntervalWindowChoice)
    w.from < w.to && return Problem[]
    return [
        Problem(
            "window_order", "The window of $key ends at $(w.to) °C, before it starts at $(w.from) °C.",
            Dict{String, Any}("key" => key, "from" => w.from, "to" => w.to),
        ),
    ]
end

# A temperature in °C as the script writes it, in kelvin.
kelvin_code(c::Real) = "T_ZERO_CELSIUS + $(julia_number(c))"

releases_keyword(r::Symbol) = r === :water ? "" : "; releases = :carbon_dioxide"

"The `DecompositionWindow` of the release `key`, as Julia source."
function window_code(key::AbstractString, w::IntervalWindowChoice)
    phase, r = release_of(key)
    return "DecompositionWindow($(repr(phase)); between = ($(kelvin_code(w.from)), $(kelvin_code(w.to)))$(releases_keyword(r)))"
end
function window_code(key::AbstractString, w::LogisticWindowChoice)
    phase, r = release_of(key)
    return "DecompositionWindow($(repr(phase)), $(kelvin_code(w.midpoint)), $(julia_number(w.width))$(releases_keyword(r)))"
end

"""
    tga_lines(t) -> Vector{String}

The windows of the thermogram, after the calculation: those the user set, which
the thermogram completes with the literature's (`literature_windows`).
"""
function tga_lines(t::TgaWindowsChoice)
    head = [
        "# The thermogram: each release of the solids (the water of a phase, its carbon",
        "# dioxide) is spread over the interval in which the cement literature reads it",
        isempty(t.windows) ? "# (`literature_windows`)." :
            "# (`literature_windows`), except those set here (temperatures in kelvin).",
    ]
    body = isempty(t.windows) ? ["windows = DecompositionWindow[]"] :
        vcat(["windows = ["], ["    " * window_code(k, w) * "," for (k, w) in t.windows], ["]"])
    return vcat(head, body, ["study = merge(study, (; windows))"])
end
