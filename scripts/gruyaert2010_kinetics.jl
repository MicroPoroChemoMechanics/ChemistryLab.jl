# =============================================================================
#  gruyaert2010_kinetics.jl — the pastes of Gruyaert et al. (2010), in time
#
#  The CEM I 52.5 N and slag pastes of scripts/gruyaert2010.jl, no longer at the
#  degrees of hydration the image analysis measured, but integrated from the
#  mixing: the four clinker phases dissolve under the Parrott–Killoh law, the
#  slag glass under the Waller law, and the equilibrium is solved along the way.
#  The page docs/src/examples/blended_slag_kinetics.md compares the run with the
#  degrees of hydration (Table 5), the heat (Tables 2 and 6) and the bound water
#  (Fig. 7) of the article, which the run was not fitted to, except where said.
#
#  The materials, the system and the activity model are those of
#  scripts/gruyaert2010.jl, so that the two routes cannot drift apart.
# =============================================================================

# Once per module: the test suite includes the other file as well.
isdefined(@__MODULE__, :gruyaert_system) || include(joinpath(@__DIR__, "gruyaert2010.jl"))

using OrdinaryDiffEq

const G10K_CLINKER_PARAMS = Dict(
    "C3S" => PK84_PARAMS_C3S, "C2S" => PK84_PARAMS_C2S,
    "C3A" => PK84_PARAMS_C3A, "C4AF" => PK84_PARAMS_C4AF,
)

"""
    gruyaert_slag_species(batch = "CAL") -> Species

The glass of the slag of Table 1 as a pseudo-species: one formula unit carries
the reported oxides of 100 g of slag ([`glass_species`](@ref)). The sulfur is
taken as sulfate, as in `gruyaert_budget`.
"""
function gruyaert_slag_species(batch = "CAL")
    r = gruyaert_oxides("BFS-" * batch)
    ox = Dict(
        o => ustrip(getproperty(r, Symbol(o))) / 100
            for o in ("CaO", "SiO2", "Al2O3", "Fe2O3", "MgO", "SO3")
    )
    return glass_species(ox; symbol = "BFS", M = 100.0u"g/mol")
end

"""
    gruyaert_kinetic_system(slag_species) -> ChemicalSystem

The system of `gruyaert_system`, with the slag glass as one more species.
"""
function gruyaert_kinetic_system(slag_species)
    cs = gruyaert_system()
    return ChemicalSystem(
        vcat(collect(cs.species), [slag_species]), CEMDATA_PRIMARIES;
        solid_solutions = collect(cs.solid_solutions),
    )
end

"""
    gruyaert_recipe(; slag, slag_species = nothing, batch = "CAL") -> Recipe

100 g of binder at w/b = 0.5 with the fraction `slag` of slag: the cement as its
Bogue phases with the gypsum and the calcite of its analysis, the slag as the
glass `slag_species`. Nothing has reacted: the rates decide.
"""
function gruyaert_recipe(; slag, slag_species = nothing, batch = "CAL")
    c = gruyaert_bogue("OPC-" * batch)
    cement = Material(
        "CEM I 52.5 N", :cement; constituents = vcat(
            [MineralConstituent(G10_DB[ph]; mass_fraction = f) for (ph, f) in c.clinker],
            [
                MineralConstituent(G10_DB["Gp"]; mass_fraction = c.gypsum),
                MineralConstituent(G10_DB["Cal"]; mass_fraction = c.calcite),
            ],
        ),
    )
    wb = ustrip(literature_value("Gruyaert2010", "water_binder_ratio"))
    slag == 0 && return Recipe(cement => 1.0; w_b = wb)
    glass = Material("slag", :scm; constituents = [MineralConstituent(slag_species; mass_fraction = 1.0)])
    return Recipe(cement => 1.0 - slag, glass => slag; w_b = wb)
end

"""
    gruyaert_slag_time(; n = WALLER_PARAMS_FLY_ASH.n) -> Quantity

The characteristic time of the Waller law for this slag, **calibrated** on one
measurement: the 72 % of slag the image analysis found at 28 months in the paste
with 50 % slag (Table 5). At 293.15 K, the reference temperature of the law, and
at the slag's Blaine fineness, equal to the reference one, the law integrates to
`α = 1/(1 + (τ/t)ⁿ)`, which is solved for `τ`. The exponent `n` is the one
Lavergne et al. (2018) fitted for a fly ash: ASSUMED for a slag.
"""
function gruyaert_slag_time(; n = WALLER_PARAMS_FLY_ASH.n)
    row = literature_table("Gruyaert2010", "hydration_degree_slag"; age_days = 852u"d", slag_to_binder = 0.5)
    α = only(row.alpha_slag) / 100
    t = 852.0
    return t * (1 / α - 1)^(1 / n) * u"d"
end

"""
    gruyaert_rates(; slag = 0.0, w_c = nothing, τ_slag = nothing, batch = "CAL") -> Dict

The rate laws: Parrott–Killoh for the four clinker phases at the cement's Blaine
fineness, slowed by the water/cement factor at `w_c` ([`pk_wc_factor`](@ref)),
and, given `τ_slag`, the Waller law for the slag glass with that characteristic
time, the apparent activation energy the article fits for this slag at the
paste's cement-to-binder ratio (Eq. 3), and the exponent of the fly-ash set.
"""
function gruyaert_rates(; slag = 0.0, w_c = nothing, τ_slag = nothing, batch = "CAL")
    blaine(m) = ustrip(gruyaert_oxides(m * "-" * batch).blaine) * u"m^2/kg"
    rates = Dict{String, Any}(
        ph => parrott_killoh_avrami(p, ph; blaine = blaine("OPC"), w_c) for (ph, p) in G10K_CLINKER_PARAMS
    )
    if τ_slag !== nothing
        G10q(k) = literature_value("Gruyaert2010", k)
        E_S = G10q("activation_energy_slag_slope") * (1 - slag) + G10q("activation_energy_slag_intercept")
        rates["BFS"] = waller(merge(WALLER_PARAMS_FLY_ASH, (τ = τ_slag, Ea = E_S)), "BFS"; blaine = blaine("BFS"))
    end
    return rates
end

"""
    gruyaert_run(cs; slag, τ_slag = nothing, slag_species = nothing, days = 1019)
        -> (; kp, sol, recipe)

The paste with the fraction `slag` of slag, integrated over `days` at 293.15 K.
The water/cement ratio of the clinker law is the water over the cement alone,
`w/b / (1 − slag)`: the slag leaves the cement more water to hydrate in. The
heat is read afterwards off the certified replay ([`heat_release`](@ref)), and
only for the paste without slag: a glass has no enthalpy of formation, so a
blend runs in the system of `gruyaert_kinetic_system` and its heat cannot be
computed.
"""
function gruyaert_run(cs; slag, τ_slag = nothing, slag_species = nothing, days = 1019)
    recipe = gruyaert_recipe(; slag, slag_species)
    wb = ustrip(literature_value("Gruyaert2010", "water_binder_ratio"))
    rates = gruyaert_rates(; slag, w_c = wb / (1 - slag), τ_slag = slag == 0 ? nothing : τ_slag)
    model = HKFActivityModel(å = 0.0, Ḃ = G10_BDOT, Kₙ = 0.0)
    kp = KineticsProblem(
        recipe, cs, rates, (0.0, days * 86400.0);
        activity_model = model,
        equilibrium_solver = EquilibriumSolver(cs, model, OptimaOptimizer()),
    )
    sol = integrate(kp, KineticsSolver(; ode_solver = Rodas5P(), reltol = 1.0e-6, abstol = 1.0e-10))
    return (; kp, sol, recipe)
end

"""
    gruyaert_degree(run, names, t) -> Float64

The degree of reaction at the time `t` (s) of the kinetic species `names`
together, by mass: one minus what is left of them over what there was.
"""
function gruyaert_degree(run, names, t)
    kp, sol = run.kp, run.sol
    u = sol(t)
    nbe = size(kp.Ae, 1)                 # the element amounts come first
    left, there = 0.0, 0.0
    for nm in names
        i = findfirst(s -> symbol(s) == nm, kp.system.species)
        k = findfirst(==(i), kp.idx_kinetic)
        M = ustrip(us"g/mol", kp.system.species[i][:M])
        n0 = ustrip(us"mol", kp.initial_state.n[i])
        left += u[nbe + k] * M
        there += n0 * M
    end
    return 1 - left / there
end
