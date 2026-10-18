# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2026 Jean-François Barthélémy (Cerema, UMR MCD)

# The materials of a mix: a published one by its template, or the user's own
# from an oxide analysis or a mineral composition entered in a table. Each is a
# ChemistryLab `Material` in the script, with a law of the degree of reaction
# for the constituents that need one.

"""
    MaterialChoice

A material of the binder, with its mass fraction in the binder.
"""
abstract type MaterialChoice <: Choice end

const MATERIAL_KINDS = ("cement", "scm", "other")

"A published material, by the name of its ChemistryLab template."
struct TemplateMaterial <: MaterialChoice
    template::String
    fraction::Float64
    extents::Vector{Pair{String, ExtentChoice}}
    material_extent::Union{Nothing, ExtentChoice}
end
kind(::Type{TemplateMaterial}) = "template"
fields(::Type{TemplateMaterial}) = [
    Field("template", :choice; default = "CEM I 52.5 N (Gruyaert 2010), Bogue", source = "templates"),
    Field("fraction", :number; default = 1.0, min = 0, max = 1),
    Field("extents", :keyed; default = Any[], source = "extent"),
    Field("material_extent", :nested; default = nothing, source = "extent", level = :expert),
]

"A material of the user's, from its oxide analysis: by Bogue's calculation, or as its oxides."
struct OxideMaterial <: MaterialChoice
    name::String
    material_kind::String
    fraction::Float64
    oxides::Vector{AmountRow}
    route::String
    sulfate::String
    density::Float64
    extents::Vector{Pair{String, ExtentChoice}}
    material_extent::Union{Nothing, ExtentChoice}
end
kind(::Type{OxideMaterial}) = "oxides"
fields(::Type{OxideMaterial}) = [
    Field("name", :text; default = "my cement"),
    Field("material_kind", :choice; default = "cement", options = collect(MATERIAL_KINDS)),
    Field("fraction", :number; default = 1.0, min = 0, max = 1),
    Field("oxides", :rows; default = Any[]),
    Field("route", :choice; default = "bogue", options = ["bogue", "oxides"]),
    Field("sulfate", :choice; default = "Gp", options = ["Gp", "Anh", "hemihydrate"]),
    Field("density", :number; default = 0.0, unit = "g/cm³", min = 0, level = :expert),
    Field("extents", :keyed; default = Any[], source = "extent"),
    Field("material_extent", :nested; default = nothing, source = "extent", level = :expert),
]

"A material of the user's, from its mineral composition: species of the database by mass."
struct MineralMaterial <: MaterialChoice
    name::String
    material_kind::String
    fraction::Float64
    minerals::Vector{AmountRow}
    extents::Vector{Pair{String, ExtentChoice}}
    material_extent::Union{Nothing, ExtentChoice}
end
kind(::Type{MineralMaterial}) = "minerals"
fields(::Type{MineralMaterial}) = [
    Field("name", :text; default = "my clinker"),
    Field("material_kind", :choice; default = "cement", options = collect(MATERIAL_KINDS)),
    Field("fraction", :number; default = 1.0, min = 0, max = 1),
    Field("minerals", :rows; default = Any[]),
    Field("extents", :keyed; default = Any[], source = "extent"),
    Field("material_extent", :nested; default = nothing, source = "extent", level = :expert),
]

members(::Type{MaterialChoice}) = (TemplateMaterial, OxideMaterial, MineralMaterial)

"The name the page shows for a material."
label(m::TemplateMaterial) = m.template
label(m::Union{OxideMaterial, MineralMaterial}) = m.name

"Whether a material's laws depend on the age."
depends_on_age(m::MaterialChoice) =
    any(depends_on_age(last(p)) for p in m.extents) || (m.material_extent !== nothing && depends_on_age(m.material_extent))

"""
    material_problems(m) -> Vector{Problem}

The problems of a material its fields' bounds do not say: the laws, and the
table of an analysis.
"""
function material_problems(m::MaterialChoice)
    out = bound_problems(m, label(m))
    for (k, e) in m.extents
        append!(out, extent_problems(e, k, label(m)))
    end
    m.material_extent === nothing || append!(out, extent_problems(m.material_extent, label(m), label(m)))
    append!(out, table_problems(m))
    return out
end

table_problems(::TemplateMaterial) = Problem[]
function table_problems(m::OxideMaterial)
    isempty(m.oxides) && return [Problem("empty_table", "The oxide analysis of $(m.name) is empty.", Dict{String, Any}("material" => m.name))]
    any(r -> r.unit != "%", m.oxides) &&
        return [Problem("oxide_unit", "An oxide analysis is entered in mass %; $(m.name) has other units.", Dict{String, Any}("material" => m.name))]
    total = sum(r.amount for r in m.oxides)
    total <= 105 && return Problem[]
    return [
        Problem(
            "oxide_total", "The oxides of $(m.name) add up to $(round(total; digits = 2)) %, more than 105 %.",
            Dict{String, Any}("material" => m.name, "total" => round(total; digits = 2)),
        ),
    ]
end
function table_problems(m::MineralMaterial)
    isempty(m.minerals) && return [Problem("empty_table", "The composition of $(m.name) is empty.", Dict{String, Any}("material" => m.name))]
    any(r -> r.unit != "%", m.minerals) &&
        return [Problem("mineral_unit", "A mineral composition is entered in mass %; $(m.name) has other units.", Dict{String, Any}("material" => m.name))]
    total = sum(r.amount for r in m.minerals)
    total <= 100 + 1.0e-9 && return Problem[]
    return [
        Problem(
            "mineral_total", "The minerals of $(m.name) add up to $(round(total; digits = 2)) %, more than 100 %.",
            Dict{String, Any}("material" => m.name, "total" => round(total; digits = 2)),
        ),
    ]
end

# ── script ────────────────────────────────────────────────────────────────────

"The `with_extents` line of a material, if it has laws."
function extents_lines(var::AbstractString, m::MaterialChoice)
    isempty(m.extents) && m.material_extent === nothing && return String[]
    pairs = join(("$(repr(k)) => $(extent_code(e, k))" for (k, e) in m.extents), ", ")
    call = "$var = with_extents($var, Dict{String, Any}($pairs)"
    m.material_extent === nothing || (call *= "; material_extent = $(extent_code(m.material_extent, label(m)))")
    return [call * ")"]
end

"Mass fractions as Julia source, from rows in percent."
fractions_code(rows) = "Dict(" * join(("$(repr(r.formula)) => $(julia_number(r.amount / 100))" for r in rows), ", ") * ")"

"""
    material_lines(var, m) -> Vector{String}

The script lines that build material `m` under the name `var`.
"""
function material_lines(var::AbstractString, m::TemplateMaterial)
    return vcat(["$var = material_template($(repr(m.template)), cs)"], extents_lines(var, m))
end

function material_lines(var::AbstractString, m::OxideMaterial)
    lines = [
        "# The oxide analysis of $(m.name), in mass fractions.",
        "$(var)_oxides = $(fractions_code(m.oxides))",
    ]
    if m.route == "oxides"
        density = m.density > 0 ? "; density = $(julia_number(m.density))" : ""
        sep = isempty(density) ? "; " : ", "
        push!(lines, "$var = oxide_material($(repr(m.name)), $(var)_oxides$density$(sep)kind = :$(m.material_kind))")
    else
        append!(
            lines, [
                "# Bogue's calculation: the clinker phases the oxides make, the sulfate as $(m.sulfate);",
                "# what no phase takes stays as one oxide constituent.",
                "$(var)_bogue = bogue($(var)_oxides, cs; sulfate = $(repr(m.sulfate)))",
                "$(var)_rest = sum(values($(var)_bogue.remainder); init = 0.0)",
                "$var = Material(",
                "    $(repr(m.name)), :$(m.material_kind); constituents = vcat(",
                "        AbstractConstituent[MineralConstituent(cs[p]; mass_fraction = f) for (p, f) in $(var)_bogue.phases if f > 0],",
                "        $(var)_rest > 0 ? [OxideConstituent(\"minor oxides\", Dict(k => v / $(var)_rest for (k, v) in $(var)_bogue.remainder); mass_fraction = $(var)_rest)] : AbstractConstituent[],",
                "    ),",
                ")",
            ]
        )
    end
    return vcat(lines, extents_lines(var, m))
end

function material_lines(var::AbstractString, m::MineralMaterial)
    constituents = ["        MineralConstituent(cs[$(repr(r.formula))]; mass_fraction = $(julia_number(r.amount / 100)))," for r in m.minerals]
    lines = vcat(
        [
            "# $(m.name), by its minerals (mass fractions); they must be species of the system.",
            "$var = Material(",
            "    $(repr(m.name)), :$(m.material_kind); constituents = AbstractConstituent[",
        ],
        constituents, ["    ],", ")"],
    )
    return vcat(lines, extents_lines(var, m))
end
