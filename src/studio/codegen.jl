# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2026 Jean-François Barthélémy (Cerema, UMR MCD)

# The script a model stands for. It is written for a reader who knows little
# Julia and less thermochemistry: every section says what it does and why, and
# the script runs on its own in the studio's environment. The studio runs the
# text of sections 1 to 5 itself, so what the page shows is what the script
# computes; sections 6 and 7 (tables and figures) serve the script alone. Each
# part of the model writes its own lines (src/features/).

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

"The names the script gives itself, which a material may not take."
const RESERVED_NAMES = Set(
    [
        "substances", "cs", "recipe", "model", "state", "certificate", "study", "masses", "width", "w_b", "T",
        "path", "paste", "states", "temperatures", "ratios", "fractions", "scm", "ages", "column", "solution",
        "V", "M_water", "fresh", "db", "species", "members", "solid_solutions", "table", "phases", "x",
        "rates", "horizon", "problem", "solution", "w_c", "calorimeter", "degrees", "curves", "times", "key",
        "glasses", "surface",
    ]
)

"""
    material_name(label, taken) -> String

A Julia variable name for a material, from its name: the words before the
parenthesis, lowercased (`"S1 slag (Durdzinski 2017)"` gives `s1_slag`), unique
among `taken` and clear of the script's other names.
"""
function material_name(label::AbstractString, taken)
    words = lowercase(strip(first(split(label, '('))))
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

"""
    script_stem(title) -> String

A file name stem from a title: lowercase words joined by underscores.
"""
function script_stem(title::AbstractString)
    s = strip(replace(lowercase(title), r"[^a-z0-9]+" => "_"), '_')
    return isempty(s) ? "chemistrylab_study" : String(s)
end

function materials_section(m::StudioModel)
    lines = vcat([heading(4, "Materials and mix")], mix_lines(m.mix))
    push!(
        lines,
        "# Each material comes from a published composition shipped with ChemistryLab",
        "# (`material_templates()` lists them) or from the user's own analysis. A degree",
        "# of reaction says which fraction of a constituent has reacted, between 0 and 1.",
    )
    names = String[]
    for x in m.materials
        name = material_name(label(x), names)
        push!(names, name)
        append!(lines, material_lines(name, x))
    end
    fractions = join(("$n => $(julia_number(x.fraction))" for (n, x) in zip(names, m.materials)), ", ")
    push!(
        lines,
        "# The mix: the mass fraction of each material in the binder, and the conditions.",
        "recipe = Recipe($fractions; $(recipe_keywords(m.mix)))",
    )
    return lines
end

"""
    diagnose_call(m) -> String

The call to `diagnose_system` for model `m`: the recipe, the model, the age the
laws are read at, and the reactants of the phase list.
"""
function diagnose_call(m::StudioModel)
    kws = ["recipe", "model"]
    a = diagnosis_age(m.calculation)
    depends_on_age(m) && a !== nothing && push!(kws, "t = $(julia_number(a))")
    m.system isa PhaseListSystem && push!(kws, "reactants = phase_list($(repr(m.system.phase_list))).reactants")
    return "diagnose_system(cs; $(join(kws, ", ")))"
end

function calculation_section(m::StudioModel)
    clines, ckw = constraint_code(m.constraint)
    return vcat(
        [heading(5, "Calculation")], activity_lines(m.activity),
        [
            "# What deserves a look before computing: a constituent the system cannot",
            "# hold, data read beyond their range, a phase counted twice…",
            "foreach(println, $(diagnose_call(m)))",
        ],
        calculation_lines(m.calculation, depends_on_age(m), clines, ckw),
        surface_lines(m.surface, m.calculation),
    )
end

"""
    diagnostic_source(m) -> String

What the studio runs to diagnose a model before it is computed: the database,
the system, the materials and the activity model, then `diagnostics`.
"""
function diagnostic_source(m::StudioModel)
    lines = vcat(
        ["using ChemistryLab, OptimaSolver, DynamicQuantities"], database_lines(m.database), system_lines(m.system),
        materials_section(m), activity_lines(m.activity), ["diagnostics = $(diagnose_call(m))"],
    )
    return join(lines, '\n') * '\n'
end

"The tables of section 6, for a single equilibrium, a path or a kinetic run."
tables_section(m::StudioModel) = vcat([heading(6, "Tables")], table_lines(m.calculation))

function table_lines(::EquilibriumCalc)
    return [

        "println(\"Certificate: optimal = \", certificate.optimal, \", scope = \", certificate.scope)",
        "println(\"pH of the pore solution: \", round(pore_solution(state).pH; digits = 2))",
        "println(\"Bound water: \", round(bound_water(state); digits = 3), \" g per g of binder\")",
        "println(\"Phase masses, in g for the binder mass of the recipe:\")",
        "masses = phase_masses(state)",
        "width = maximum(length, keys(masses)) + 2",
        "for (phase, mass) in masses",
        "    println(\"  \", rpad(phase, width), round(mass; digits = 2))",
        "end",
    ]
end

function table_lines(::CalculationChoice)
    return [
        "# One row per point of the path: the parameter, whether the equilibrium is",
        "# certified, the pH, the porosity, the bound water and the mass of each phase.",
        "table = process_table(path)",
        "show(stdout, MIME\"text/plain\"(), table; allrows = true, allcols = true)",
        "println()",
    ]
end

function table_lines(::KineticsCalc)
    return [
        "# The degree of reaction of each constituent with a rate law, at the ages.",
        "degrees = degrees_of_hydration(solution, problem; times = ages .* 86400)",
        "println(\"Degree of reaction at \", join(ages, \", \"), \" days:\")",
        "for (name, α) in degrees",
        "    println(\"  \", rpad(name, 10), join((rpad(round(a; digits = 3), 8) for a in α)))",
        "end",
        "println(\"pH of the pore solution: \", [round(pore_solution(s, model).pH; digits = 2) for s in states])",
    ]
end

"The figure of section 7, for a single equilibrium, a path or a kinetic run."
function figures_section(m::StudioModel, stem::AbstractString)
    head = [heading(7, "Figures"), "# Plots draws the figure, saved next to this script.", "using Plots"]
    return vcat(head, figure_lines(m.calculation, stem))
end

function figure_lines(::EquilibriumCalc, stem)
    return [
        "bar(collect(keys(masses)), collect(values(masses)); legend = false, xrotation = 45,",
        "    ylabel = \"mass (g)\", title = \"Phase assemblage\")",
        "savefig(joinpath(@__DIR__, $(repr(stem * "_phases.png"))))",
    ]
end

function figure_lines(::KineticsCalc, stem)
    return [
        "# The degree of reaction of each constituent in time, from the integration's own steps.",
        "times = filter(>(0), solution.t)",
        "curves = degrees_of_hydration(solution, problem; times)",
        "plot(times ./ 86400, reduce(hcat, collect(values(curves))); labels = permutedims(collect(keys(curves))),",
        "    xscale = :log10, xlabel = \"age (d)\", ylabel = \"degree of reaction\", legend = :topleft)",
        "savefig(joinpath(@__DIR__, $(repr(stem * "_degrees.png"))))",
    ]
end

function figure_lines(::CalculationChoice, stem)
    return [
        "# The mass of each phase along the path, stacked.",
        "masses = [phase_masses(rs) for rs in path]",
        "phases = unique(k for d in masses for k in keys(d) if !startswith(k, \"unreacted\"))",
        "x = Float64.(path.parameter)",
        "areaplot(x, reduce(hcat, [[get(d, p, 0.0) for d in masses] for p in phases]);",
        "    labels = permutedims(phases), xlabel = string(path.label), ylabel = \"mass (g)\", legend = :outerright)",
        "savefig(joinpath(@__DIR__, $(repr(stem * "_phases.png"))))",
    ]
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
        "# Written by ChemistryLab Studio (ChemistryLab $version) on $date.",
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
        ],
        database_lines(m.database), [""],
        [heading(3, "Chemical system")], system_lines(m.system), [""],
        materials_section(m), [""],
        calculation_section(m),
    )
    rest = vcat([""], tables_section(m), [""], figures_section(m, stem))
    compute_text = join(compute, '\n') * '\n'
    body = compute_text * join(rest, '\n') * '\n'
    return GeneratedScript(embed_model(body, m), compute_text)
end
