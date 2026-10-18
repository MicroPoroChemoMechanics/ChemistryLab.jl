# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2026 Jean-François Barthélémy (Cerema, UMR MCD)

# What to compute: one equilibrium of the paste, or a path of equilibria (in
# time, along a replacement, a carbonation, a salt or another titrant, renewals
# of a leaching or attacking water, a range of temperatures or of water). The
# paths are ChemistryLab's processes (src/recipes/processes.jl); the sweeps are
# written the way its documentation writes them, and tabulated alike.

"""
    CalculationChoice

The calculation of a study; the `mode` of its kind is `"equilibrium"`,
`"extents"` (a path of equilibria) or `"kinetics"`, the button of the page that
offers it.
"""
abstract type CalculationChoice <: Choice end

mode(::Type{<:CalculationChoice}) = "extents"

"One equilibrium of the paste, at the age that the laws of the degrees of reaction read."
struct EquilibriumCalc <: CalculationChoice
    age::Float64
end
kind(::Type{EquilibriumCalc}) = "equilibrium"
fields(::Type{EquilibriumCalc}) = [Field("age", :number; default = 28.0, unit = "d", min = 0)]
mode(::Type{EquilibriumCalc}) = "equilibrium"

"The paste at each of a list of ages, the degrees of reaction read at each."
struct HydrationPath <: CalculationChoice
    ages::Vector{Float64}
end
kind(::Type{HydrationPath}) = "hydration"
fields(::Type{HydrationPath}) = [Field("ages", :numbers; default = [0.5, 1.0, 2.0, 7.0, 28.0, 90.0], unit = "d", min = 0)]

"Part of the binder replaced by a supplementary material, at increasing fractions."
struct BlendPath <: CalculationChoice
    scm::String
    scm_extent::Float64
    fractions::Vector{Float64}
    age::Float64
end
kind(::Type{BlendPath}) = "blend"
fields(::Type{BlendPath}) = [
    Field("scm", :choice; default = "blast-furnace slag (Gruyaert 2010)", source = "templates"),
    Field("scm_extent", :number; default = 0.5, min = 0, max = 1),
    Field("fractions", :numbers; default = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5], min = 0, max = 1),
    Field("age", :number; default = 28.0, unit = "d", min = 0),
]

"The hydrated paste carbonated by increasing amounts of CO₂."
struct CarbonationPath <: CalculationChoice
    amounts::Vector{Float64}
    age::Float64
end
kind(::Type{CarbonationPath}) = "carbonation"
fields(::Type{CarbonationPath}) = [
    Field("amounts", :numbers; default = [0.0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0], unit = "mol", min = 0),
    Field("age", :number; default = 28.0, unit = "d", min = 0),
]

"""
A salt added to the hydrated paste, in increasing amounts. The system must hold
its elements: for a chloride, a phase that binds chlorine (Friedel's salt,
`C4AClH10`), added in the System tab.
"""
struct SaltPath <: CalculationChoice
    salt::String
    amounts::Vector{Float64}
    age::Float64
end
kind(::Type{SaltPath}) = "salt"
fields(::Type{SaltPath}) = [
    Field("salt", :text; default = "NaCl"),
    Field("amounts", :numbers; default = [0.0, 0.01, 0.02, 0.04, 0.06], unit = "mol", min = 0),
    Field("age", :number; default = 28.0, unit = "d", min = 0),
]

"Any species or formula the system can express, added to the hydrated paste."
struct TitrationPath <: CalculationChoice
    species::String
    amounts::Vector{Float64}
    age::Float64
end
kind(::Type{TitrationPath}) = "titration"
fields(::Type{TitrationPath}) = [
    Field("species", :text; default = "HCl"),
    Field("amounts", :numbers; default = [0.0, 0.1, 0.2, 0.4], unit = "mol", min = 0),
    Field("age", :number; default = 28.0, unit = "d", min = 0, level = :expert),
]

"The pore solution removed and renewed with pure water, step after step."
struct LeachingPath <: CalculationChoice
    steps::Int
    renewal::Float64
    age::Float64
end
kind(::Type{LeachingPath}) = "leaching"
fields(::Type{LeachingPath}) = [
    Field("steps", :integer; default = 12, min = 1, max = 1000),
    Field("renewal", :number; default = 5000.0, unit = "g", min = 1.0e-9),
    Field("age", :number; default = 28.0, unit = "d", min = 0),
]

"""
The pore solution removed and renewed with a solution of the user's (seawater, a
sulfate solution), given by its solutes in a table, step after step.
"""
struct AttackPath <: CalculationChoice
    solution::Vector{AmountRow}
    volume::Float64
    steps::Int
    age::Float64
end
kind(::Type{AttackPath}) = "attack"
fields(::Type{AttackPath}) = [
    Field("solution", :rows; default = Any[Dict("formula" => "Na2SO4", "amount" => 50.0, "unit" => "mmol/L", "notation" => "mendeleev")]),
    Field("volume", :number; default = 500.0, unit = "mL", min = 1.0e-9),
    Field("steps", :integer; default = 10, min = 1, max = 1000),
    Field("age", :number; default = 28.0, unit = "d", min = 0),
]

"The same recipe at a range of temperatures."
struct TemperatureSweep <: CalculationChoice
    temperatures::Vector{Float64}
    age::Float64
end
kind(::Type{TemperatureSweep}) = "temperature_sweep"
fields(::Type{TemperatureSweep}) = [
    Field("temperatures", :numbers; default = [5.0, 20.0, 40.0, 60.0], unit = "°C", min = 0, max = 99),
    Field("age", :number; default = 28.0, unit = "d", min = 0),
]

"The same recipe at a range of water-to-binder ratios."
struct WaterSweep <: CalculationChoice
    ratios::Vector{Float64}
    age::Float64
end
kind(::Type{WaterSweep}) = "water_sweep"
fields(::Type{WaterSweep}) = [
    Field("ratios", :numbers; default = [0.3, 0.4, 0.5, 0.6], min = 0.01, max = 100),
    Field("age", :number; default = 28.0, unit = "d", min = 0),
]

"The age a calculation reads the laws at; a hydration reads them at each of its ages."
age(c::CalculationChoice) = c.age

"The age the diagnosis reads the laws at: the calculation's, the first of a hydration."
diagnosis_age(c::CalculationChoice) = age(c)
diagnosis_age(c::HydrationPath) = isempty(c.ages) ? nothing : first(c.ages)

calculation_problems(c::CalculationChoice) = bound_problems(c, "the calculation")
function calculation_problems(c::Union{HydrationPath, BlendPath, CarbonationPath, SaltPath, TitrationPath, TemperatureSweep, WaterSweep})
    out = bound_problems(c, "the calculation")
    values = parameter_values(c)
    isempty(values) && push!(out, Problem("empty_path", "The list of values of the calculation is empty."))
    return out
end
function calculation_problems(c::AttackPath)
    out = bound_problems(c, "the calculation")
    isempty(c.solution) && push!(out, Problem("empty_solution", "The attacking solution has no solute."))
    any(r -> r.unit ∉ ("mmol/L", "mol/L"), c.solution) &&
        push!(out, Problem("solution_unit", "The solutes of the attacking solution are given in mol/L or mmol/L."))
    return out
end

parameter_values(c::HydrationPath) = c.ages
parameter_values(c::BlendPath) = c.fractions
parameter_values(c::Union{CarbonationPath, SaltPath, TitrationPath}) = c.amounts
parameter_values(c::TemperatureSweep) = c.temperatures
parameter_values(c::WaterSweep) = c.ratios

vector_code(v) = "[" * join(julia_number.(v), ", ") * "]"

"`t = 28.0, ` when the laws read an age, nothing otherwise."
age_keyword(c::CalculationChoice, timed::Bool) = timed ? "t = $(julia_number(age(c))), " : ""

"""
    calculation_lines(c, timed, constraint_lines, constraint_kw) -> Vector{String}

The lines that compute, ending with `study`, whose `output` is a `RecipeState`
or a `ProcessResult`: `timed` when the laws read an age, after the lines and
with the keywords of the constraint.
"""
function calculation_lines(c::EquilibriumCalc, timed::Bool, constraint_lines, constraint_kw)
    kw = isempty(constraint_kw) ? "" : ", " * constraint_kw
    return vcat(
        constraint_lines, [
            "# The equilibrium minimizes the Gibbs energy of the paste. It comes with a",
            "# certificate saying what was proven (`certificate.optimal`, `certificate.scope`).",
            "state, certificate = equilibrate_certified(recipe, cs; $(age_keyword(c, timed))model$kw)",
            "# Everything computed above, gathered for the tables and figures below.",
            "study = (; recipe, cs, model, output = state, certificate)",
        ]
    )
end

# The hydrated paste the processes start from.
function base_lines(c::CalculationChoice, timed::Bool, constraint_lines, constraint_kw)
    kw = isempty(constraint_kw) ? "" : ", " * constraint_kw
    return vcat(
        constraint_lines, [
            "# The hydrated paste the process starts from.",
            "paste, certificate = equilibrate_certified(recipe, cs; $(age_keyword(c, timed))model$kw)",
        ]
    )
end

const STUDY_LINE = "study = (; recipe, cs, model, output = path, certificate = nothing)"

function calculation_lines(c::HydrationPath, timed::Bool, constraint_lines, constraint_kw)
    return vcat(
        constraint_lines, [
            "# The paste at each age, each equilibrium started from the previous one.",
            "ages = $(vector_code(c.ages))   # days",
            "path = hydrate(recipe, cs, ages; model$(isempty(constraint_kw) ? "" : ", " * constraint_kw))",
            STUDY_LINE,
        ]
    )
end

function calculation_lines(c::BlendPath, timed::Bool, constraint_lines, constraint_kw)
    return [
        "# The supplementary material, and the part of it that has reacted.",
        "scm = with_extents(material_template($(repr(c.scm)), cs), Dict{String, Any}(); material_extent = $(julia_number(c.scm_extent)))",
        "# The binder with a growing fraction of it, the other materials scaled down.",
        "fractions = $(vector_code(c.fractions))",
        "path = blend(recipe, scm, fractions, cs; $(age_keyword(c, timed))model)",
        STUDY_LINE,
    ]
end

function calculation_lines(c::CarbonationPath, timed::Bool, constraint_lines, constraint_kw)
    return vcat(
        base_lines(c, timed, constraint_lines, constraint_kw), [
            "# CO₂ added in increasing amounts (mol for the binder mass of the recipe).",
            "path = carbonate(paste, $(vector_code(c.amounts)))",
            STUDY_LINE,
        ]
    )
end

function calculation_lines(c::SaltPath, timed::Bool, constraint_lines, constraint_kw)
    return vcat(
        base_lines(c, timed, constraint_lines, constraint_kw), [
            "# $(c.salt) added in increasing amounts (mol for the binder mass of the recipe).",
            "path = add_salt(paste, $(repr(c.salt)), $(vector_code(c.amounts)))",
            STUDY_LINE,
        ]
    )
end

function calculation_lines(c::TitrationPath, timed::Bool, constraint_lines, constraint_kw)
    return vcat(
        base_lines(c, timed, constraint_lines, constraint_kw), [
            "# $(c.species) added in increasing amounts (mol for the binder mass of the recipe).",
            "path = titrate(paste, $(repr(c.species)), $(vector_code(c.amounts)))",
            STUDY_LINE,
        ]
    )
end

function calculation_lines(c::LeachingPath, timed::Bool, constraint_lines, constraint_kw)
    return vcat(
        base_lines(c, timed, constraint_lines, constraint_kw), [
            "# At each step the pore solution is removed and replaced by pure water.",
            "path = leach(paste, $(c.steps); renewal = $(julia_number(c.renewal))u\"g\")",
            STUDY_LINE,
        ]
    )
end

function calculation_lines(c::AttackPath, timed::Bool, constraint_lines, constraint_kw)
    terms = join(
        ("$(julia_number(r.unit == "mmol/L" ? r.amount / 1000 : r.amount)) .* column($(repr(r.formula)))" for r in c.solution),
        " .+\n    ",
    )
    return vcat(
        base_lines(c, timed, constraint_lines, constraint_kw), [
            "# The attacking solution as a budget over the primaries of the system: each",
            "# solute (mol/L) and the water it is dissolved in, its density taken as one.",
            "column(f) = (j = findfirst(==(f), symbol.(cs.species));",
            "    j === nothing ? primary_decomposition(Species(f), cs.SM.primaries) : Float64.(cs.SM.A[:, j]))",
            "V = $(julia_number(c.volume))   # mL of solution per renewal",
            "M_water = ustrip(us\"g/mol\", cs.dict_species[\"H2O@\"][:M])",
            "solution = V / M_water .* column(\"H2O@\") .+ V / 1000 .* (",
            "    $terms",
            ")",
            "# At each step the pore solution is removed and replaced by the solution.",
            "path = leach(paste, $(c.steps); solution)",
            STUDY_LINE,
        ]
    )
end

function calculation_lines(c::TemperatureSweep, timed::Bool, constraint_lines, constraint_kw)
    kw = isempty(constraint_kw) ? "" : ", " * constraint_kw
    return [
        "# The same recipe at each temperature (°C).",
        "temperatures = $(vector_code(c.temperatures))",
        "states = [",
        "    first(equilibrate_certified(Recipe(recipe.binder...; w_b, T = T_ZERO_CELSIUS_Q + θ * u\"K\"), cs; $(age_keyword(c, timed))model$kw))",
        "    for θ in temperatures",
        "]",
        "path = ProcessResult(:temperature, temperatures, states)",
        STUDY_LINE,
    ]
end

function calculation_lines(c::WaterSweep, timed::Bool, constraint_lines, constraint_kw)
    kw = isempty(constraint_kw) ? "" : ", " * constraint_kw
    return [
        "# The same binder at each water-to-binder ratio.",
        "ratios = $(vector_code(c.ratios))",
        "states = [",
        "    first(equilibrate_certified(Recipe(recipe.binder...; w_b = r, T), cs; $(age_keyword(c, timed))model$kw))",
        "    for r in ratios",
        "]",
        "path = ProcessResult(:w_b, ratios, states)",
        STUDY_LINE,
    ]
end

members(::Type{CalculationChoice}) = (
    EquilibriumCalc, HydrationPath, BlendPath, CarbonationPath, SaltPath, TitrationPath,
    LeachingPath, AttackPath, TemperatureSweep, WaterSweep,
)
member_traits(::Type{T}) where {T <: CalculationChoice} = Dict{String, Any}("mode" => mode(T))
