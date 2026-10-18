# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2026 Jean-François Barthélémy (Cerema, UMR MCD)

# What the page edits: the choices of a calculation, kept in the server so that
# restarting the compute process loses nothing, and embedded in every script
# the studio writes so that the script reopens exactly. The model holds choices
# and the user's own numbers only; published data stay in ChemistryLab's files
# and are read by name.

"""
    MaterialChoice(template, fraction, extents, material_extent)

A material of the mix: a ChemistryLab template (`material_templates()`), its
mass fraction of the binder, the degree of reaction of each constituent named
(`"C3S" => 0.8`), and the degree of the material as a whole, `nothing` to leave
the template's.
"""
mutable struct MaterialChoice
    template::String
    fraction::Float64
    extents::Vector{Pair{String, Float64}}
    material_extent::Union{Nothing, Float64}
end

"""
    StudioModel

The choices of a calculation: its title, the database file, the published phase
list the system is built from, the activity model, the materials of the mix,
the water-to-binder ratio, the temperature (°C) and the calculation mode.
"""
mutable struct StudioModel
    title::String
    database::String
    phase_list::String
    activity_model::String
    materials::Vector{MaterialChoice}
    w_b::Float64
    temperature_c::Float64
    mode::String
end

"The activity models the page offers, by key."
const ACTIVITY_MODELS = ("cemdata18_KOH", "cemdata18_NaOH")

"The calculation modes implemented so far."
const MODES = ("equilibrium",)

const MODEL_VERSION = 1

"""
    default_model() -> StudioModel

The calculation the studio opens with, and compiles at start: a paste of the
CEM I 52.5 N of Gruyaert et al. (2010) at their water-to-binder ratio, 0.5, in
the phase list of Lothenbach and Winnefeld (2006). Every part of this cement
that is left unreacted has a sourced density, so that its volumes and porosity
are known. The degrees of reaction are illustrative.
"""
function default_model()
    cement = MaterialChoice(
        "CEM I 52.5 N (Gruyaert 2010), Bogue", 1.0,
        ["C3S" => 0.8, "C2S" => 0.45, "C3A" => 0.9, "C4AF" => 0.6], nothing,
    )
    return StudioModel(
        "CEM I 52.5 N paste", "cemdata18-thermofun.json",
        "Portland paste (Lothenbach and Winnefeld 2006)", "cemdata18_KOH",
        [cement], 0.5, 20.0, "equilibrium",
    )
end

"The model as the JSON object the page and the embedded block carry."
function to_dict(m::StudioModel)
    return Dict{String, Any}(
        "version" => MODEL_VERSION,
        "title" => m.title, "database" => m.database, "phase_list" => m.phase_list,
        "activity_model" => m.activity_model,
        "materials" => [
            Dict{String, Any}(
                "template" => c.template, "fraction" => c.fraction,
                "extents" => [Dict{String, Any}("constituent" => k, "extent" => v) for (k, v) in c.extents],
                "material_extent" => c.material_extent,
            ) for c in m.materials
        ],
        "w_b" => m.w_b, "temperature_c" => m.temperature_c, "mode" => m.mode,
    )
end

# Field readers that name the field when the page sends a wrong type: the
# page is ours, so a mismatch is a bug to report, not a value to guess.
function _field(d::AbstractDict, key::AbstractString, ::Type{String}, default)
    v = get(d, key, default)
    v isa AbstractString || throw(ArgumentError("model field \"$key\" must be text, got $(repr(v))"))
    return String(v)
end
function _field(d::AbstractDict, key::AbstractString, ::Type{Float64}, default)
    v = get(d, key, default)
    (v isa Real && !(v isa Bool)) || throw(ArgumentError("model field \"$key\" must be a number, got $(repr(v))"))
    return Float64(v)
end
function _optional_number(d::AbstractDict, key::AbstractString)
    v = get(d, key, nothing)
    v === nothing && return nothing
    return _field(d, key, Float64, nothing)
end

function _material(d)
    d isa AbstractDict || throw(ArgumentError("each material must be an object, got $(repr(d))"))
    raw = get(d, "extents", Any[])
    raw isa AbstractVector || throw(ArgumentError("model field \"extents\" must be a list"))
    extents = Pair{String, Float64}[
        _field(e, "constituent", String, nothing) => _field(e, "extent", Float64, nothing) for e in raw
    ]
    return MaterialChoice(
        _field(d, "template", String, nothing), _field(d, "fraction", Float64, nothing),
        extents, _optional_number(d, "material_extent"),
    )
end

"""
    model_from_dict(d) -> StudioModel

The model a JSON object describes, the defaults filling what it omits. A field
of the wrong type is refused by name.
"""
function model_from_dict(d::AbstractDict)
    version = get(d, "version", MODEL_VERSION)
    version == MODEL_VERSION || throw(ArgumentError("model version $version is not one this studio reads ($MODEL_VERSION)"))
    base = default_model()
    raw = get(d, "materials", nothing)
    materials = raw === nothing ? base.materials : raw isa AbstractVector ? MaterialChoice[_material(x) for x in raw] :
        throw(ArgumentError("model field \"materials\" must be a list"))
    return StudioModel(
        _field(d, "title", String, base.title), _field(d, "database", String, base.database),
        _field(d, "phase_list", String, base.phase_list), _field(d, "activity_model", String, base.activity_model),
        materials, _field(d, "w_b", Float64, base.w_b), _field(d, "temperature_c", Float64, base.temperature_c),
        _field(d, "mode", String, base.mode),
    )
end

"""
    Problem(code, text, values)

Something that keeps a model from running: a code the page translates, the
English text, and the values the text quotes.
"""
struct Problem
    code::String
    text::String
    values::Dict{String, Any}
end

Problem(code, text) = Problem(code, text, Dict{String, Any}())

to_dict(p::Problem) = Dict{String, Any}("code" => p.code, "text" => p.text, "values" => p.values)

"""
    problems(m) -> Vector{Problem}

What keeps the model from running, in the user's terms. Empty when it can run.
"""
function problems(m::StudioModel)
    out = Problem[]
    isempty(m.materials) && push!(out, Problem("no_material", "The mix has no material."))
    for c in m.materials
        c.fraction > 0 || push!(
            out, Problem(
                "fraction_not_positive", "The fraction of $(c.template) must be positive.",
                Dict{String, Any}("material" => c.template),
            )
        )
        for (k, v) in c.extents
            0 <= v <= 1 || push!(
                out, Problem(
                    "extent_out_of_range",
                    "The degree of reaction of $k in $(c.template) must lie between 0 and 1.",
                    Dict{String, Any}("constituent" => k, "material" => c.template),
                )
            )
        end
        e = c.material_extent
        (e === nothing || 0 <= e <= 1) || push!(
            out, Problem(
                "extent_out_of_range",
                "The degree of reaction of $(c.template) must lie between 0 and 1.",
                Dict{String, Any}("constituent" => c.template, "material" => c.template),
            )
        )
    end
    total = sum(c.fraction for c in m.materials; init = 0.0)
    (isempty(m.materials) || isapprox(total, 1; atol = 1.0e-9)) || push!(
        out, Problem(
            "fractions_sum", "The binder fractions add up to $(round(total; digits = 4)), not 1.",
            Dict{String, Any}("total" => round(total; digits = 4)),
        )
    )
    m.w_b > 0 || push!(out, Problem("w_b_not_positive", "The water-to-binder ratio must be positive."))
    0 <= m.temperature_c <= 99 || push!(
        out, Problem("temperature_out_of_range", "The temperature must lie between 0 and 99 °C.")
    )
    m.activity_model in ACTIVITY_MODELS || push!(
        out, Problem(
            "unknown_activity_model", "Unknown activity model $(m.activity_model).",
            Dict{String, Any}("model" => m.activity_model),
        )
    )
    m.mode in MODES || push!(
        out, Problem("unknown_mode", "Unknown calculation mode $(m.mode).", Dict{String, Any}("mode" => m.mode))
    )
    return out
end
