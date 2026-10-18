# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2026 Jean-François Barthélémy (Cerema, UMR MCD)

# How much of a constituent has reacted: the degree of reaction, between 0 and
# 1, constant or a law of the age. Each law is one of ChemistryLab's extents
# (src/recipes/extents.jl), written in the script with the parameters the page
# gives. A constituent with no law reacts completely, as in the templates.

"""
    ExtentChoice

A law of the degree of reaction of a constituent (`ConstantExtent`,
`LogisticExtent`, `TabulatedExtent`, `ParrottKillohExtent`, `CappedExtent`).
"""
abstract type ExtentChoice <: Choice end

"A degree of reaction that does not change with the age."
struct ConstantExtentChoice <: ExtentChoice
    value::Float64
end
kind(::Type{ConstantExtentChoice}) = "constant"
fields(::Type{ConstantExtentChoice}) = [Field("value", :number; default = 1.0, min = 0, max = 1)]
extent_code(x::ConstantExtentChoice, _) = julia_number(x.value)

"A logistic law in the logarithm of the age, fitted to measured degrees of reaction."
struct LogisticExtentChoice <: ExtentChoice
    final::Float64
    half_time::Float64
    slope::Float64
    initial::Float64
end
kind(::Type{LogisticExtentChoice}) = "logistic"
fields(::Type{LogisticExtentChoice}) = [
    Field("final", :number; default = 0.8, min = 0, max = 1),
    Field("half_time", :number; default = 3.0, unit = "d", min = 1.0e-6),
    Field("slope", :number; default = 1.5, min = 1.0e-6),
    Field("initial", :number; default = 0.0, min = 0, max = 1, level = :expert),
]
extent_code(x::LogisticExtentChoice, _) =
    "LogisticExtent(; final = $(julia_number(x.final)), half_time = $(julia_number(x.half_time)), " *
    "slope = $(julia_number(x.slope)), initial = $(julia_number(x.initial)))"

"Measured degrees of reaction at given ages, interpolated linearly between them."
struct TabulatedExtentChoice <: ExtentChoice
    days::Vector{Float64}
    values::Vector{Float64}
end
kind(::Type{TabulatedExtentChoice}) = "tabulated"
fields(::Type{TabulatedExtentChoice}) = [
    Field("days", :numbers; default = [1.0, 7.0, 28.0], unit = "d", min = 0),
    Field("values", :numbers; default = [0.3, 0.6, 0.8], min = 0, max = 1),
]
extent_code(x::TabulatedExtentChoice, _) =
    "TabulatedExtent([$(join(julia_number.(x.days), ", "))], [$(join(julia_number.(x.values), ", "))])"

"""
The law of Parrott and Killoh for a clinker phase, at the temperature and
water-to-cement ratio of the mix, and for a Blaine fineness; capped, when asked,
at the degree the water allows (Powers).
"""
struct ParrottKillohExtentChoice <: ExtentChoice
    blaine::Float64
    powers_cap::Bool
end
kind(::Type{ParrottKillohExtentChoice}) = "parrott_killoh"
fields(::Type{ParrottKillohExtentChoice}) = [
    Field("blaine", :number; default = 385.0, unit = "m²/kg", min = 1),
    Field("powers_cap", :bool; default = true),
]
const CLINKER_PHASES = ("C3S", "C2S", "C3A", "C4AF")
function extent_code(x::ParrottKillohExtentChoice, phase::AbstractString)
    cap = x.powers_cap ? ", α_max = powers_alpha_max(w_b)" : ""
    return "ParrottKillohExtent($(repr(phase)); T, w_c = w_b, blaine = $(julia_number(x.blaine))$cap)"
end

"Another law, never above a cap."
struct CappedExtentChoice <: ExtentChoice
    inner::ExtentChoice
    cap::Float64
end
kind(::Type{CappedExtentChoice}) = "capped"
fields(::Type{CappedExtentChoice}) = [
    Field("inner", :nested; default = Dict{String, Any}("kind" => "logistic"), source = "extent"),
    Field("cap", :number; default = 0.9, min = 0, max = 1),
]
extent_code(x::CappedExtentChoice, phase) = "CappedExtent($(extent_code(x.inner, phase)), $(julia_number(x.cap)))"

members(::Type{ExtentChoice}) =
    (ConstantExtentChoice, LogisticExtentChoice, TabulatedExtentChoice, ParrottKillohExtentChoice, CappedExtentChoice)

"Whether a law depends on the age, so that a calculation needs one."
depends_on_age(::ConstantExtentChoice) = false
depends_on_age(::ExtentChoice) = true
depends_on_age(x::CappedExtentChoice) = depends_on_age(x.inner)

"The problems of a law that its bounds do not say, for the constituent `key` of `where_`."
extent_problems(::ExtentChoice, key, where_) = Problem[]
function extent_problems(x::TabulatedExtentChoice, key, where_)
    length(x.days) == length(x.values) && length(x.days) >= 2 && issorted(x.days; lt = <=) && return Problem[]
    return [
        Problem(
            "tabulated_extent", "The table of $key in $where_ needs as many degrees as ages, at least two, the ages increasing.",
            Dict{String, Any}("constituent" => key, "material" => where_),
        ),
    ]
end
function extent_problems(::ParrottKillohExtentChoice, key, where_)
    key in CLINKER_PHASES && return Problem[]
    return [
        Problem(
            "parrott_killoh_phase", "The law of Parrott and Killoh applies to C3S, C2S, C3A and C4AF, not to $key ($where_).",
            Dict{String, Any}("constituent" => key, "material" => where_),
        ),
    ]
end
extent_problems(x::CappedExtentChoice, key, where_) = extent_problems(x.inner, key, where_)
