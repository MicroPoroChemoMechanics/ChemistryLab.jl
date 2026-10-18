# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2026 Jean-François Barthélémy (Cerema, UMR MCD)

# The script a model stands for. It is written for a reader who knows little
# Julia and less thermochemistry: every section says what it does and why, and
# the script runs on its own in the studio's environment. The studio runs the
# text of sections 1 to 5 itself, so what the page shows is what the script
# computes; sections 6 and 7 (tables and figures) serve the script alone.

"""
    GeneratedScript(source, compute)

A generated script: `source`, the whole file, model block included, and
`compute`, sections 1 to 5, which the studio runs.
"""
struct GeneratedScript
    source::String
    compute::String
end

const RULE_WIDTH = 78

"A section heading as a comment line: `## 1. Packages ----…`."
function heading(n::Integer, title::AbstractString)
    head = "## $n. $title "
    return head * repeat("-", max(RULE_WIDTH - textwidth(head), 3))
end

"A number as Julia source, keeping it a float: `0.45`, `1.0`, `20.0`."
julia_number(x::Real) = repr(Float64(x))

"""
    material_name(template, taken) -> String

A Julia variable name for a material, from its template's name: the words
before the parenthesis, lowercased (`"S1 slag (Durdzinski 2017)"` gives
`s1_slag`), unique among `taken` and clear of the script's other names.
"""
function material_name(template::AbstractString, taken)
    words = lowercase(strip(first(split(template, '('))))
    base = strip(replace(words, r"[^a-z0-9]+" => "_"), '_')
    (isempty(base) || !isletter(first(base))) && (base = "material_" * base)
    base in RESERVED_NAMES && (base *= "_material")
    name, k = base, 2
    while name in taken
        name = "$(base)_$k"
        k += 1
    end
    return String(name)
end

const RESERVED_NAMES = Set(["substances", "cs", "recipe", "model", "state", "certificate", "study", "masses", "width"])

function activity_model_code(key::AbstractString)
    key == "cemdata18_KOH" && return "cemdata18_activity_model(:KOH)"
    key == "cemdata18_NaOH" && return "cemdata18_activity_model(:NaOH)"
    throw(ArgumentError("unknown activity model $key"))
end

function activity_model_comment(key::AbstractString)
    salt = key == "cemdata18_NaOH" ? "NaOH" : "KOH"
    return "# How the dissolved species interact: extended Debye–Hückel, with the\n" *
        "# parameters Cemdata18 recommends for a pore solution dominated by $salt."
end

function extents_code(name::AbstractString, c::MaterialChoice)
    isempty(c.extents) && c.material_extent === nothing && return String[]
    pairs = join(("$(repr(k)) => $(julia_number(v))" for (k, v) in c.extents), ", ")
    call = "$name = with_extents($name, Dict{String, Any}($pairs)"
    c.material_extent === nothing || (call *= "; material_extent = $(julia_number(c.material_extent))")
    return [call * ")"]
end

function materials_section(m::StudioModel)
    lines = String[
        heading(4, "Materials and mix"),
        "# Each material comes from a published composition shipped with ChemistryLab",
        "# (`material_templates()` lists them). A degree of reaction says which",
        "# fraction of a constituent has reacted, between 0 and 1.",
    ]
    names = String[]
    for c in m.materials
        name = material_name(c.template, names)
        push!(names, name)
        push!(lines, "$name = material_template($(repr(c.template)), cs)")
        append!(lines, extents_code(name, c))
    end
    fractions = join(("$n => $(julia_number(c.fraction))" for (n, c) in zip(names, m.materials)), ", ")
    t = julia_number(m.temperature_c)
    append!(
        lines, [
            "# The mix: the mass fraction of each material in the binder, the",
            "# water-to-binder ratio w/b, and the temperature ($(t) °C, written in kelvins).",
            "recipe = Recipe($fractions; w_b = $(julia_number(m.w_b)), T = T_ZERO_CELSIUS_Q + $(t)u\"K\")",
        ]
    )
    return lines
end

function calculation_section(m::StudioModel)
    return String[
        heading(5, "Calculation"),
        activity_model_comment(m.activity_model),
        "model = $(activity_model_code(m.activity_model))",
        "# The equilibrium minimizes the Gibbs energy of the paste. It comes with a",
        "# certificate saying what was proven (`certificate.optimal`, `certificate.scope`).",
        "state, certificate = equilibrate_certified(recipe, cs; model)",
        "# Everything computed above, gathered for the tables and figures below.",
        "study = (; recipe, cs, model, output = state, certificate)",
    ]
end

const TABLES_SECTION = String[
    heading(6, "Tables"),
    "println(\"Certificate: optimal = \", certificate.optimal, \", scope = \", certificate.scope)",
    "println(\"pH of the pore solution: \", round(pore_solution(state).pH; digits = 2))",
    "println(\"Bound water: \", round(bound_water(state); digits = 3), \" g per g of binder\")",
    "println(\"Phase masses, in g for the binder mass of the recipe (100 g by default):\")",
    "masses = phase_masses(state)",
    "width = maximum(length, keys(masses)) + 2",
    "for (phase, mass) in masses",
    "    println(\"  \", rpad(phase, width), round(mass; digits = 2))",
    "end",
]

function figures_section(stem::AbstractString)
    file = repr(stem * "_phases.png")
    return String[
        heading(7, "Figures"),
        "# Plots draws the figure, saved next to this script.",
        "using Plots",
        "bar(collect(keys(masses)), collect(values(masses)); legend = false, xrotation = 45,",
        "    ylabel = \"mass (g)\", title = \"Phase assemblage\")",
        "savefig(joinpath(@__DIR__, $file))",
    ]
end

"""
    script_stem(title) -> String

A file name stem from a title: lowercase words joined by underscores.
"""
function script_stem(title::AbstractString)
    s = strip(replace(lowercase(title), r"[^a-z0-9]+" => "_"), '_')
    return isempty(s) ? "chemistrylab_study" : String(s)
end

"""
    generate(m; version, date) -> GeneratedScript

The script of model `m`, written by the studio of version `version` on `date`.
"""
function generate(m::StudioModel; version::VersionNumber = STUDIO_VERSION, date::Dates.Date = Dates.today())
    stem = script_stem(m.title)
    header = String[
        "# " * repeat("=", RULE_WIDTH - 2),
        "# $(m.title)",
        "# Written by ChemistryLab Studio $version on $date.",
        "#",
        "# This script runs on its own, outside the studio:",
        "#     julia --project=@chemistrylab-studio $stem.jl",
        "# Sections 1 to 5 compute; section 6 prints tables; section 7 draws figures.",
        "# " * repeat("=", RULE_WIDTH - 2),
    ]
    compute = vcat(
        header, [""],
        [
            heading(1, "Packages"),
            "# ChemistryLab holds the chemistry, OptimaSolver solves the equilibrium,",
            "# and DynamicQuantities gives the units (u\"g\", u\"K\").",
            "using ChemistryLab, OptimaSolver, DynamicQuantities",
            "",
            heading(2, "Thermodynamic database"),
            "# The species and their thermodynamic data, downloaded from the publisher of",
            "# the database at the first use and kept afterwards.",
            "substances = build_species(datapath($(repr(m.database))); verbose = false)",
            "",
            heading(3, "Chemical system"),
            "# The phases the paste may form, from a published phase list",
            "# (`phase_lists()` lists them).",
            "cs = phase_list_system($(repr(m.phase_list)), substances)",
            "",
        ],
        materials_section(m), [""],
        calculation_section(m),
    )
    rest = vcat([""], TABLES_SECTION, [""], figures_section(stem))
    compute_text = join(compute, '\n') * '\n'
    body = compute_text * join(rest, '\n') * '\n'
    return GeneratedScript(embed_model(body, m), compute_text)
end
