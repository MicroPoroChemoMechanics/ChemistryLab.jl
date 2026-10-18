# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2026 Jean-François Barthélémy (Cerema, UMR MCD)

# What the page edits: the choices of a calculation, kept in the server so that
# restarting the compute process loses nothing, and embedded in every script
# the studio writes so that the script reopens exactly. The model holds choices
# and the user's own numbers only; published data stay in ChemistryLab's files
# and are read by name. Each part is a family of choices (src/schema.jl,
# src/features/), read, written and checked by the same declarations.

"""
    StudioModel

The choices of a calculation: its title, the database, the chemical system, the
activity model, the materials of the binder, the mix, the calculation and the
constraint of the equilibrium.
"""
mutable struct StudioModel
    title::String
    database::DatabaseChoice
    system::SystemChoice
    activity::ActivityChoice
    materials::Vector{MaterialChoice}
    mix::MixChoice
    calculation::CalculationChoice
    constraint::ConstraintChoice
end

const MODEL_VERSION = 2

"""
    default_model() -> StudioModel

The calculation the studio opens with, and compiles at start: a paste of the
CEM I 52.5 N of Gruyaert et al. (2010) at their water-to-binder ratio, 0.5, in
the phase list of Lothenbach and Winnefeld (2006). Every part of this cement
that is left unreacted has a sourced density, so that its volumes and porosity
are known. The degrees of reaction are illustrative.
"""
function default_model()
    extents = Pair{String, ExtentChoice}[
        "C3S" => ConstantExtentChoice(0.8), "C2S" => ConstantExtentChoice(0.45),
        "C3A" => ConstantExtentChoice(0.9), "C4AF" => ConstantExtentChoice(0.6),
    ]
    return StudioModel(
        "CEM I 52.5 N paste",
        BuiltinDatabase("cemdata18-thermofun.json"),
        PhaseListSystem("Portland paste (Lothenbach and Winnefeld 2006)", String[], String[], "CSHQ", String[]),
        Cemdata18Activity("KOH", true),
        MaterialChoice[TemplateMaterial("CEM I 52.5 N (Gruyaert 2010), Bogue", 1.0, extents, nothing)],
        Mix(0.5, 20.0, 100.0, 1.0),
        EquilibriumCalc(28.0),
        FixedTPChoice(),
    )
end

register_family!("database", DatabaseChoice)
register_family!("system", SystemChoice)
register_family!("activity", ActivityChoice)
register_family!("material", MaterialChoice)
register_family!("extent", ExtentChoice)
register_family!("mix", MixChoice)
register_family!("calculation", CalculationChoice)
register_family!("constraint", ConstraintChoice)
register_family!("rate", RateChoice)
register_family!("calorimeter", CalorimeterChoice)

"The model as the JSON object the page and the embedded block carry."
function to_dict(m::StudioModel)
    return Dict{String, Any}(
        "version" => MODEL_VERSION, "title" => m.title,
        "database" => to_dict(m.database), "system" => to_dict(m.system), "activity" => to_dict(m.activity),
        "materials" => to_dict.(m.materials), "mix" => to_dict(m.mix),
        "calculation" => to_dict(m.calculation), "constraint" => to_dict(m.constraint),
    )
end

function _part(d::AbstractDict, key::AbstractString, ::Type{F}, default::F) where {F <: Choice}
    haskey(d, key) || return default
    return choice_from_dict(F, d[key])
end

"""
    model_from_dict(d) -> StudioModel

The model a JSON object describes, the defaults filling the parts it omits. A
field of the wrong type is refused by name. A model of version 1 (the first
studio) is read as the same calculation.
"""
function model_from_dict(d::AbstractDict)
    version = get(d, "version", MODEL_VERSION)
    version == 1 && return model_from_dict(migrate_v1(d))
    version == MODEL_VERSION || throw(ArgumentError("model version $version is not one this studio reads ($MODEL_VERSION)"))
    base = default_model()
    title = get(d, "title", base.title)
    title isa AbstractString || throw(ArgumentError("model field \"title\" must be text, got $(repr(title))"))
    raw = get(d, "materials", nothing)
    materials = raw === nothing ? base.materials : raw isa AbstractVector ?
        MaterialChoice[choice_from_dict(MaterialChoice, x) for x in raw] :
        throw(ArgumentError("model field \"materials\" must be a list"))
    return StudioModel(
        String(title),
        _part(d, "database", DatabaseChoice, base.database), _part(d, "system", SystemChoice, base.system),
        _part(d, "activity", ActivityChoice, base.activity), materials, _part(d, "mix", MixChoice, base.mix),
        _part(d, "calculation", CalculationChoice, base.calculation), _part(d, "constraint", ConstraintChoice, base.constraint),
    )
end

"""
    migrate_v1(d) -> Dict

A model of the first studio (one equilibrium of template materials, constant
degrees of reaction) in the form of version 2.
"""
function migrate_v1(d::AbstractDict)
    materials = Any[
        Dict{String, Any}(
            "kind" => "template", "template" => get(x, "template", ""), "fraction" => get(x, "fraction", 1.0),
            "extents" => Any[
                Dict{String, Any}("key" => e["constituent"], "value" => Dict{String, Any}("kind" => "constant", "value" => e["extent"]))
                    for e in get(x, "extents", Any[])
            ],
            "material_extent" => get(x, "material_extent", nothing) === nothing ? nothing :
                Dict{String, Any}("kind" => "constant", "value" => x["material_extent"]),
        ) for x in get(d, "materials", Any[])
    ]
    am = get(d, "activity_model", "cemdata18_KOH")
    return Dict{String, Any}(
        "version" => MODEL_VERSION, "title" => get(d, "title", "study"),
        "database" => Dict{String, Any}("kind" => "builtin", "name" => get(d, "database", "cemdata18-thermofun.json")),
        "system" => Dict{String, Any}("kind" => "phase_list", "phase_list" => get(d, "phase_list", "Portland paste (Lothenbach and Winnefeld 2006)")),
        "activity" => Dict{String, Any}("kind" => "cemdata18", "electrolyte" => am == "cemdata18_NaOH" ? "NaOH" : "KOH"),
        "materials" => materials,
        "mix" => Dict{String, Any}("kind" => "mix", "w_b" => get(d, "w_b", 0.5), "temperature" => get(d, "temperature_c", 20.0)),
        "calculation" => Dict{String, Any}("kind" => "equilibrium"),
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

"Whether some law of the degrees of reaction depends on the age."
depends_on_age(m::StudioModel) = any(depends_on_age, m.materials)

"""
    problems(m) -> Vector{Problem}

What keeps the model from running, in the user's terms. Empty when it can run.
"""
function problems(m::StudioModel)
    out = Problem[]
    append!(out, database_problems(m.database))
    append!(out, bound_problems(m.system, "the system"))
    append!(out, system_problems(m.system))
    append!(out, bound_problems(m.activity, "the activity model"))
    isempty(m.materials) && push!(out, Problem("no_material", "The mix has no material."))
    for x in m.materials
        x.fraction > 0 || push!(
            out, Problem(
                "fraction_not_positive", "The fraction of $(label(x)) must be positive.",
                Dict{String, Any}("material" => label(x)),
            )
        )
        append!(out, material_problems(x))
    end
    total = sum(x.fraction for x in m.materials; init = 0.0)
    (isempty(m.materials) || isapprox(total, 1; atol = 1.0e-9)) || push!(
        out, Problem(
            "fractions_sum", "The binder fractions add up to $(round(total; digits = 4)), not 1.",
            Dict{String, Any}("total" => round(total; digits = 4)),
        )
    )
    append!(out, bound_problems(m.mix, "the mix"))
    append!(out, calculation_problems(m.calculation))
    if m.calculation isa HydrationPath && !depends_on_age(m)
        push!(out, Problem("no_law_of_age", "A path in time needs a degree of reaction that depends on the age (a law, not a constant)."))
    end
    append!(out, constraint_problems(m.constraint))
    if m.calculation isa KineticsCalc && !(m.constraint isa FixedTPChoice)
        push!(out, Problem("kinetics_constraint", "A kinetic run holds the temperature and the pressure only: choose no constraint."))
    end
    return out
end
