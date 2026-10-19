# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2026 Jean-François Barthélémy (Cerema, UMR MCD)

# Rate equations integrated in time: a constituent given a rate law reacts as
# the law says, the others are at equilibrium from the start (to their degree of
# reaction), and the equilibrium of the products is solved at each step. This
# is ChemistryLab's `KineticsProblem` built from the recipe
# (src/recipes/recipe_kinetics.jl), integrated by OrdinaryDiffEq, then read at
# the ages to report as certified equilibria (`speciated_states`). The laws are
# those of ChemistryLab's library, or one the user writes.

"""
    RateChoice

The rate law of a constituent that reacts in time (`parrott_killoh_avrami`,
`waller`, a first-order law, or a law of the user's).
"""
abstract type RateChoice <: Choice end

"The law of Parrott and Killoh for a clinker phase, with its published parameters."
struct PKAvramiRate <: RateChoice
    blaine::Float64
    water_factor::Bool
    powers_cap::Bool
end
kind(::Type{PKAvramiRate}) = "parrott_killoh_avrami"
fields(::Type{PKAvramiRate}) = [
    Field("blaine", :number; default = 385.0, unit = "m²/kg", min = 0),
    Field("water_factor", :bool; default = true),
    Field("powers_cap", :bool; default = true, level = :expert),
]

"The law of Waller for the reaction of a fly ash or a silica fume."
struct WallerRate <: RateChoice
    preset::String
    tau::Float64
    blaine::Float64
end
kind(::Type{WallerRate}) = "waller"
fields(::Type{WallerRate}) = [
    Field("preset", :choice; default = "fly_ash", options = ["fly_ash", "silica_fume"]),
    Field("tau", :number; default = 0.0, unit = "d", min = 0),
    Field("blaine", :number; default = 0.0, unit = "m²/kg", min = 0, level = :expert),
]

"A rate proportional to what is left of the constituent."
struct FirstOrderRate <: RateChoice
    k::Float64
end
kind(::Type{FirstOrderRate}) = "first_order"
fields(::Type{FirstOrderRate}) = [Field("k", :number; default = 0.1, unit = "1/d", min = 0)]

"A rate law of the user's, as a Julia expression of `T`, `P`, `t`, `n`, `lna`, `n0` and `key`."
struct CustomRate <: RateChoice
    expression::String
end
kind(::Type{CustomRate}) = "custom_rate"
fields(::Type{CustomRate}) = [Field("expression", :text; default = "1.0e-7 * n[key]")]

members(::Type{RateChoice}) = (PKAvramiRate, WallerRate, FirstOrderRate, CustomRate)

"The clinker phases with a published set of the law of Parrott and Killoh."
const PK_PHASES = ("C3S", "C2S", "C3A", "C4AF")

"""
    rate_problems(r, key) -> Vector{Problem}

What keeps the law `r` of constituent `key` from running besides its bounds
(which the calculation's own check reads): a phase the law of Parrott and
Killoh has no parameters for, an expression that does not read as Julia.
"""
rate_problems(::RateChoice, key) = Problem[]
function rate_problems(r::PKAvramiRate, key)
    out = Problem[]
    key in PK_PHASES || push!(
        out, Problem(
            "pk_rate_phase", "The law of Parrott and Killoh has parameters for C3S, C2S, C3A and C4AF, not for $key.",
            Dict{String, Any}("constituent" => key),
        )
    )
    return out
end
function rate_problems(r::CustomRate, key)
    out = Problem[]
    ex = Meta.parse(r.expression; raise = false)
    if isempty(strip(r.expression)) || (ex isa Expr && ex.head in (:error, :incomplete))
        push!(
            out, Problem(
                "rate_syntax", "The rate law of $key does not read as a Julia expression.",
                Dict{String, Any}("constituent" => key),
            )
        )
    end
    return out
end

"""
    rate_code(r, key) -> String

The rate law of constituent `key` as Julia source; `w_c`, the water over the
cement, is defined before when a law reads it.
"""
function rate_code(r::PKAvramiRate, key)
    kws = String[]
    r.blaine > 0 && push!(kws, "blaine = $(julia_number(r.blaine))")
    r.water_factor && push!(kws, "w_c")
    r.powers_cap && push!(kws, "α_max = powers_alpha_max(w_c)")
    return "parrott_killoh_avrami(PK84_PARAMS_$key, $(repr(key))$(isempty(kws) ? "" : "; " * join(kws, ", ")))"
end
function rate_code(r::WallerRate, key)
    params = r.preset == "silica_fume" ? "WALLER_PARAMS_SILICA_FUME" : "WALLER_PARAMS_FLY_ASH"
    r.tau > 0 && (params = "merge($params, (τ = $(julia_number(r.tau))u\"d\",))")
    return "waller($params, $(repr(key))$(r.blaine > 0 ? "; blaine = $(julia_number(r.blaine))" : ""))"
end
rate_code(r::FirstOrderRate, key) =
    "KineticFunc((T, P, t, n, lna, n0) -> $(julia_number(r.k)) / 86400 * n[$(repr(key))], NamedTuple(), u\"mol/s\")"
rate_code(r::CustomRate, key) =
    "let key = $(repr(key))\n        KineticFunc((T, P, t, n, lna, n0) -> $(strip(r.expression)), NamedTuple(), u\"mol/s\")\n    end"

"Whether a law reads the water over the cement."
reads_w_c(r::PKAvramiRate) = r.water_factor || r.powers_cap
reads_w_c(::RateChoice) = false

"""
    CalorimeterChoice

The thermal conditions of a kinetic run: isothermal, or a semi-adiabatic cell
whose temperature follows the heat the paste releases.
"""
abstract type CalorimeterChoice <: Choice end

"The temperature of the mix, held."
struct IsothermalChoice <: CalorimeterChoice end
kind(::Type{IsothermalChoice}) = "isothermal"
fields(::Type{IsothermalChoice}) = Field[]

"A semi-adiabatic cell: its heat capacity, its loss to the room (linear and quadratic in the difference of temperature), the room's temperature."
struct SemiAdiabaticChoice <: CalorimeterChoice
    heat_capacity::Float64
    loss::Float64
    ambient::Float64
    loss_quadratic::Float64
end
kind(::Type{SemiAdiabaticChoice}) = "semi_adiabatic"
fields(::Type{SemiAdiabaticChoice}) = [
    Field("heat_capacity", :number; default = 900.0, unit = "J/K", min = 0),
    Field("loss", :number; default = 0.05, unit = "W/K", min = 0),
    Field("ambient", :number; default = 20.0, unit = "°C", min = -20, max = 99),
    Field("loss_quadratic", :number; default = 0.0, unit = "W/K²", min = 0, level = :expert),
]

members(::Type{CalorimeterChoice}) = (IsothermalChoice, SemiAdiabaticChoice)

calorimeter_lines(::IsothermalChoice) = String[], ""
function calorimeter_lines(c::SemiAdiabaticChoice)
    # the loss in watts, linear in the difference of temperature or with the
    # quadratic term of a calibration (NF EN 196-9)
    loss = c.loss_quadratic > 0 ?
        "heat_loss = ΔT -> $(julia_number(c.loss)) * ΔT + $(julia_number(c.loss_quadratic)) * ΔT^2," :
        "L = $(julia_number(c.loss))u\"W/K\","
    return [
            "# A semi-adiabatic cell: the heat capacity of the vessel and the inert parts (the",
            "# paste's own is added), its loss to the room (W), the room's temperature.",
            "calorimeter = SemiAdiabaticCalorimeter(;",
            "    Cp = $(julia_number(c.heat_capacity))u\"J/K\", $loss",
            "    T_env = T_ZERO_CELSIUS_Q + $(julia_number(c.ambient))u\"K\", T0 = T,",
            ")",
        ], ", calorimeter"
end

"""
    KineticsCalc

A kinetic run of the paste up to `horizon` days, reported at `ages`: the
constituents named in `rates` react by their laws, the others are at
equilibrium to their degree of reaction from the start.
"""
struct KineticsCalc <: CalculationChoice
    horizon::Float64
    ages::Vector{Float64}
    rates::Vector{Pair{String, RateChoice}}
    calorimeter::CalorimeterChoice
    reltol::Float64
end
kind(::Type{KineticsCalc}) = "kinetics"
fields(::Type{KineticsCalc}) = [
    Field("horizon", :number; default = 28.0, unit = "d", min = 1.0e-3, max = 36500),
    Field("ages", :numbers; default = [1.0, 2.0, 7.0, 28.0], unit = "d", min = 0),
    Field("rates", :keyed; default = Any[], source = "rate"),
    Field("calorimeter", :nested; default = Dict{String, Any}("kind" => "isothermal"), source = "calorimeter"),
    Field("reltol", :number; default = 1.0e-6, min = 1.0e-12, max = 1.0e-2, level = :expert),
]
mode(::Type{KineticsCalc}) = "kinetics"
parameter_values(c::KineticsCalc) = c.ages

"The start of the run, where the laws of the degrees of reaction are read."
diagnosis_age(::KineticsCalc) = 0.0

function calculation_problems(c::KineticsCalc)
    out = bound_problems(c, "the calculation")
    isempty(c.ages) && push!(out, Problem("empty_path", "The list of values of the calculation is empty."))
    isempty(c.rates) && push!(out, Problem("no_rate", "A kinetic run needs a rate law for at least one constituent."))
    any(>(c.horizon), c.ages) && push!(
        out, Problem(
            "ages_beyond_horizon", "The ages to report go beyond the $(julia_number(c.horizon)) days of the run.",
            Dict{String, Any}("horizon" => c.horizon),
        )
    )
    for (key, r) in c.rates
        append!(out, rate_problems(r, key))
    end
    return out
end

function calculation_lines(c::KineticsCalc, timed::Bool, constraint_lines, constraint_kw)
    cal, cal_kw = calorimeter_lines(c.calorimeter)
    water = any(reads_w_c ∘ last, c.rates) ? [
            "# The water over the cement alone, which the clinker laws read: the other",
            "# materials of the binder leave the cement more water.",
            "w_c = w_b / sum(f for (m, f) in recipe.binder if m.kind === :cement; init = 0.0)",
        ] : String[]
    rates = ["    $(repr(key)) => $(rate_code(r, key))," for (key, r) in c.rates]
    glass = [
        "# A glass (a constituent known by its oxides) has no formula to dissolve: each",
        "# one given a rate law takes, in its material and in the system, a pseudo-species",
        "# written from its own oxides, with an enthalpy estimated from them. With no",
        "# glass among the laws, the recipe and the system are left as they are.",
        "glasses = with_glass_species(recipe, cs, keys(rates))",
        "recipe, cs = glasses.recipe, glasses.system",
    ]
    return vcat(
        [
            "# Kinetics is integrated in time by OrdinaryDiffEq.",
            "using OrdinaryDiffEq",
        ],
        water,
        [
            "# The rate law of each constituent that reacts in time (mol/s); the others",
            "# are at equilibrium from the start, to their degree of reaction.",
            "rates = Dict{String, Any}(",
        ],
        rates,
        [")"],
        glass,
        cal,
        [
            "horizon = $(julia_number(c.horizon))   # days",
            "problem = KineticsProblem(",
            "    recipe, cs, rates, (0.0, horizon * 86400); activity_model = model,",
            "    equilibrium_solver = EquilibriumSolver(cs, model, OptimaOptimizer())$cal_kw,",
            ")",
            "solution = integrate(problem, KineticsSolver(; ode_solver = Rodas5P(), reltol = $(julia_number(c.reltol)), abstol = 1.0e-10))",
            "# The paste at the ages to report that the run reached (an integration that",
            "# fails stops short of its horizon), each a certified equilibrium of what has",
            "# reacted by then.",
            "ages = filter(a -> a * 86400 <= last(solution.t), $(vector_code(c.ages)))   # days",
            "states = speciated_states(solution, problem; times = ages .* 86400)",
            "study = (; recipe, cs, model, output = (; problem, solution, ages, states), certificate = nothing)",
        ],
    )
end
