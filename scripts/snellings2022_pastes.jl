# =============================================================================
#  snellings2022_pastes.jl — the slag-limestone cement integrated at 5, 20, 40 °C
#
#  The ternary cement of Snellings et al. (2022), 50 % CEM I 52.5 R, 40 %
#  slag and 10 % limestone, at w/b 0.5, integrated from the mixing at each of
#  the three curing temperatures: the four clinker phases under Parrott–Killoh
#  with the published constants and activation energies, the glass of the slag
#  under the Waller law fitted on its degree of reaction at the three
#  temperatures (scripts/snellings2022_kinetics.jl), everything else at
#  equilibrium, the limestone included. The page
#  docs/src/examples/slag_temperature_pastes.md compares the runs with the
#  bound water and the portlandite of Fig. 8 and the hydrates of Fig. 10, none
#  of which was fitted on.
#
#  The materials are the templates of data/recipe_templates.toml, built from
#  Table 1 in data/literature/Snellings2022.json. Each assumption is stated
#  where it is made.
# =============================================================================

isdefined(@__MODULE__, :sn22_slag_fit) || include(joinpath(@__DIR__, "snellings2022_kinetics.jl"))
using OptimaSolver
using OrdinaryDiffEq

const SN22_DB = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-thermofun.json"); verbose = false))

# The equilibrium phases: those of the quaternary pastes of Schöler et al.
# (2015), a CEM I 52.5 R with slag and limestone too, CSHQ with its alkali end
# members and the siliceous hydrogarnet that Cemdata18 gives the iron.
const SN22_PURE = split(
    "Gp Anh hemihydrate syngenite K2SO4 Cal Portlandite ettringite monosulphate12 " *
        "monocarbonate hemicarbonate C4AH13 C3AH6 C3FH6 straetlingite " *
        "hydrotalcite Mg2AlC0.5OH Brc FeOOHmic AlOHmic Amor-Sl"
)
const SN22_SS = [
    "CSHQ" => ["CSHQ-JenD", "CSHQ-JenH", "CSHQ-TobD", "CSHQ-TobH", "KSiOH", "NaSiOH"],
    "C3(AF)S0.84H" => ["C3AFS0.84H4.32", "C3FS0.84H4.32"],
]

"""
    sn22p_setup() -> (; cs, pc, slag, ls, glass)

The system and the three materials, built together as for the pastes of
Schöler et al. (2015): the glass of the slag becomes the pseudo-species
`"SLAG"`, which a rate law can dissolve, once the primaries of the system are
known, and the system then holds it.

The quartz of the three materials is inert. The calcite of the slag and of the
limestone, the dolomite, the sulfates and the minor oxides of the cement are at
equilibrium from the mixing. ASSUMED: the oxides of the glass and of the minor
oxides the system has no element for (TiO₂, MnO, P₂O₅) stay in their mass and
out of the formula; the slag's sulfur, which its analysis does not report, is
left out.
"""
function sn22p_setup()
    members = reduce(vcat, last.(SN22_SS))
    sp = speciation(
        collect(values(SN22_DB)), vcat(SN22_PURE, collect(keys(SN22_PK)), members);
        aggregate_state = [AS_AQUEOUS], exclude_species = split("H2@ O2@ CH4@"),
    )
    ss = [SolidSolutionPhase(n, [SN22_DB[m] for m in ms]) for (n, ms) in SN22_SS]
    base = ChemicalSystem(sp, CEMDATA_PRIMARIES; solid_solutions = ss)
    inert_quartz(m) = with_extents(m, Dict("Quartz" => 0.0))
    pc = inert_quartz(material_template("PC (Snellings 2022)", SN22_DB))
    slag = inert_quartz(material_template("slag (Snellings 2022)", SN22_DB))
    ls = inert_quartz(material_template("limestone (Snellings 2022)", SN22_DB))
    glass = glass_species(only(c for c in slag.constituents if c.name == "glass"), base; symbol = "SLAG")
    slag = with_species(slag, Dict("glass" => glass))
    cs = ChemicalSystem(vcat(sp, [glass]), CEMDATA_PRIMARIES; solid_solutions = ss)
    return (; cs, pc, slag, ls, glass)
end

"""
    sn22p_recipe(setup, T_C; w_b = 0.5) -> Recipe

100 g of the ternary cement, 50:40:10, at `w_b` and `T_C` °C. Nothing has
reacted yet except what is at equilibrium from the mixing.
"""
sn22p_recipe(setup, T_C; w_b = 0.5) = Recipe(
    setup.pc => sn22_value("pc_percent") / 100, setup.slag => sn22_value("slag_percent") / 100,
    setup.ls => sn22_value("limestone_percent") / 100; w_b, T = _sn22_kelvin(T_C),
)

"""
    sn22p_rates(setup, θ; w_b = 0.5) -> Dict

The rate laws, by constituent: Parrott–Killoh for the four clinker phases of
the cement, at its Blaine fineness and at the water/cement ratio of the paste
(`w_b` over the share of the cement), with the published constants and
activation energies; the Waller law of the slag glass with the constants `θ`,
one entry of the `θ` of [`sn22_slag_fit`](@ref), at the slag's own fineness, on
which it was fitted.
"""
function sn22p_rates(setup, θ; w_b = 0.5)
    blaine = literature_value(SN22, "blaine_pc")
    w_c = w_b / (sn22_value("pc_percent") / 100)
    rates = Dict{String, Any}()
    for c in setup.pc.constituents
        c isa MineralConstituent || continue
        s = symbol(c.species)
        haskey(SN22_PK, s) && (rates[c.name] = parrott_killoh_avrami(SN22_PK[s], s; blaine, w_c))
    end
    rates["SLAG"] = sn22_slag_law(θ; name = "SLAG")
    return rates
end

"""
    sn22p_run(setup, T_C, θ; w_b = 0.5, days = 180) -> (; kp, sol, recipe, T_C)

The paste integrated over `days` at the constant temperature `T_C` °C, in
Cemdata18's extended Debye–Hückel model for a KOH solution.
"""
function sn22p_run(setup, T_C, θ; w_b = 0.5, days = 180)
    recipe = sn22p_recipe(setup, T_C; w_b)
    model = cemdata18_activity_model(:KOH)
    kp = KineticsProblem(
        recipe, setup.cs, sn22p_rates(setup, θ; w_b), (0.0, days * 86400.0);
        activity_model = model,
        equilibrium_solver = EquilibriumSolver(setup.cs, model, OptimaOptimizer()),
    )
    sol = integrate(kp, KineticsSolver(; ode_solver = Rodas5P(), reltol = 1.0e-6, abstol = 1.0e-10))
    return (; kp, sol, recipe, T_C)
end

# What Figs. 8 and 10 measure, by the species of the system that each counts.
const SN22_XRD = (
    portlandite = ("Portlandite",), ettringite = ("ettringite",),
    hydrotalcite = ("hydrotalcite", "Mg2AlC0.5OH"), hemicarbonate = ("hemicarbonate",),
    monocarbonate = ("monocarbonate",),
)

"""
    sn22p_observables(run, days) -> Vector{NamedTuple}

At each of `days`, on the certified replay of `run`: the bound water and the
portlandite as the thermogravimetry of Fig. 8 reports them, over the mass at
550 °C, and the hydrates as the diffraction of Fig. 10 reports them, in g per
100 g of binder.

ASSUMED: every hydrate has lost all its water by 550 °C and nothing else has
been lost, the carbonates keeping their CO₂; the mass at 550 °C is every solid,
the unreacted part of the binder included, less that water.
"""
function sn22p_observables(run, days)
    inert = sum(x.mass for x in budget(run.recipe, run.kp.system).residual; init = 0.0)
    g(q) = ustrip(uconvert(us"g", q))
    states = speciated_states(run.sol, run.kp; times = collect(float.(days)) .* 86400.0)
    sd = run.kp.system.dict_species
    return map(states) do st
        water = g(bound_water(st))
        dry = g(mass(st).solid) + inert - water
        m(names) = sum(g(mass(st, sd[n])) for n in names if haskey(sd, n); init = 0.0)
        xrd = map(m, SN22_XRD)
        (; bound_water = 100 * water / dry, portlandite_tga = 100 * xrd.portlandite / dry, xrd...)
    end
end

"""
    sn22p_measured(phase, technique, T_C) -> Vector

The measured values of Figs. 8 and 10 for the paste at w/b 0.5 cured at `T_C`
°C, at the six ages, in percent.
"""
function sn22p_measured(phase, technique, T_C)
    t = literature_table(SN22, "hydrates_wb05"; phase, technique)
    k = findall(x -> ustrip(x) == T_C, t.temperature_C)
    p = sortperm(ustrip.(us"d", t.age[k]))
    return Float64.(ustrip.(t.percent[k][p]))
end
