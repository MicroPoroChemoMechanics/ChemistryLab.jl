# =============================================================================
#  scholer2015_kinetics.jl — the quaternary cements of Schöler et al. (2015), in time
#
#  Ten pastes of 50 % OPC with blast-furnace slag, siliceous fly ash and
#  limestone (their Table 4), at w/b = 0.45 and 20 °C, integrated from the
#  mixing: the four clinker phases under the Parrott–Killoh law, the glass of
#  the slag and that of the fly ash under the Waller law, everything else at
#  equilibrium, the limestone included. The page
#  docs/src/examples/quaternary_kinetics.md compares the runs with the bound
#  water and the portlandite the authors measured by thermogravimetry (their
#  Table 8), which nothing here was fitted to.
#
#  The materials are the templates of data/recipe_templates.toml, built from
#  the analyses of Tables 1 to 3 in data/literature/Scholer2015.json. Each
#  assumption is stated where it is made.
# =============================================================================

using ChemistryLab
using DynamicQuantities
using OptimaSolver
using OrdinaryDiffEq

const S15 = "Scholer2015"
const S15_DB = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-thermofun.json"); verbose = false))

s15_value(k) = literature_value(S15, k)

# The equilibrium phases: those of the Portland pages, the sulfates and the
# alkali sulfate of the cement, and CSHQ with its alkali end members, since the
# cement and both additions carry potassium and sodium.
const S15_PURE = split(
    "Gp Anh hemihydrate syngenite K2SO4 Cal Portlandite ettringite monosulphate12 " *
        "monocarbonate hemicarbonate C4AH13 C3AH6 C3FH6 straetlingite " *
        "hydrotalcite Mg2AlC0.5OH Brc FeOOHmic AlOHmic Amor-Sl"
)
const S15_SS = [
    "CSHQ" => ["CSHQ-JenD", "CSHQ-JenH", "CSHQ-TobD", "CSHQ-TobH", "KSiOH", "NaSiOH"],
    "C3(AF)S0.84H" => ["C3AFS0.84H4.32", "C3FS0.84H4.32"],
]
const S15_CLINKER = ("C3S", "C2S", "C3A", "C4AF")

# The database of the CASH+ models, read the first time a setup asks for it.
const _S15_CASHPLUS = Ref{Any}(nothing)
function _s15_cashplus()
    _S15_CASHPLUS[] === nothing &&
        (_S15_CASHPLUS[] = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-cashplus.json"); verbose = false)))
    return _S15_CASHPLUS[]
end

# The aqueous ion pairs Miron et al. (2022a) left out when fitting CASH+NK
# (their Sections 3.2 and 7.4), left out with it.
const S15_CASHPLUS_EXCLUDED = ["NaOH@", "KOH@", "NaHSiO3@", "KHSiO3@"]

"""
    s15_setup(; gel = "CSHQ") -> (; cs, mats)

The system and the materials, built together because each needs the other: the
glass of a material becomes a species only once the primaries of the system are
known, and the system then has to hold that species.

`mats` holds the four materials of Table 1 from their templates and the
anhydrite the authors add to bring every mix to 3 % SO₃ (`opc`, `bfs`, `fa`,
`ls`, `anhydrite`), the glass of the slag and that of the fly ash replaced by
their pseudo-species, `"BFS"` and `"FA"`, which a rate law can dissolve
(`glasses` holds the two species). `cs` is the system: the aqueous species of
Cemdata18 for the elements of the materials, the equilibrium phases, the four
clinker phases and the two glasses.

As in the authors' own calculations (Section 2.2), only the glass of the two
additions reacts: their crystals, quartz and mullite among them, are inert, and
so is the quartz of the cement. Everything else in the cement and the limestone
is at equilibrium from the mixing: the sulfates, the alkali sulfates, the
periclase, and the calcite, which the minimization turns into carboaluminates as
far as the phases competing for the aluminum leave it any.

ASSUMED: the sulfur of the slag, sulfide in the glass, is taken as the sulfate
the analysis reports it as; the oxides of the glasses the system has no element
for (TiO₂, MnO, P₂O₅) stay in their mass and out of their formula.

`gel` is the model of the C-S-H: `"CSHQ"`, or `"CNASH_ss"` (Myers et al. 2014)
or `"CASH+NK"` (Miron et al. 2022a, on `cemdata18-cashplus.json`, without the
aqueous ion pairs its authors left out when fitting it), in the form
`data/solid_solutions.toml` ships them.
"""
function s15_setup(; gel = "CSHQ")
    gel in ("CSHQ", "CNASH_ss", "CASH+NK") ||
        throw(ArgumentError("s15_setup: the gel is \"CSHQ\", \"CNASH_ss\" or \"CASH+NK\"; got \"$gel\""))
    db = gel == "CASH+NK" ? _s15_cashplus() : S15_DB
    ss = [SolidSolutionPhase(n, [db[m] for m in ms]) for (n, ms) in S15_SS]
    if gel != "CSHQ"
        shipped = only(p for p in build_solid_solutions(datapath("solid_solutions.toml"), db) if name(p) == gel)
        ss = [shipped; ss[2:end]]
    end
    members = [symbol(m) for p in ss for m in p.end_members]
    excluded = vcat(split("H2@ O2@ CH4@"), gel == "CASH+NK" ? S15_CASHPLUS_EXCLUDED : String[])
    sp = speciation(
        collect(values(db)), vcat(S15_PURE, collect(S15_CLINKER), members);
        aggregate_state = [AS_AQUEOUS], exclude_species = excluded,
    )
    base = ChemicalSystem(sp, CEMDATA_PRIMARIES; solid_solutions = ss)

    inert(m) = with_extents(m, Dict(c.name => 0.0 for c in m.constituents if c.name != "glass"))
    opc = with_extents(material_template("OPC (Schöler 2015)", S15_DB), Dict("Quartz" => 0.0))
    bfs = inert(material_template("blast-furnace slag (Schöler 2015)", S15_DB))
    fa = inert(material_template("siliceous fly ash (Schöler 2015)", S15_DB))
    ls = material_template("limestone (Schöler 2015)", S15_DB)
    glass(m, sym) = glass_species(only(c for c in m.constituents if c.name == "glass"), base; symbol = sym)
    g_bfs, g_fa = glass(bfs, "BFS"), glass(fa, "FA")
    anhydrite = Material("anhydrite", :other; constituents = [MineralConstituent(S15_DB["Anh"]; mass_fraction = 1.0)])
    mats = (;
        opc, ls, anhydrite,
        bfs = with_species(bfs, Dict("glass" => g_bfs)),
        fa = with_species(fa, Dict("glass" => g_fa)),
        glasses = (g_bfs, g_fa),
    )
    cs = ChemicalSystem(vcat(sp, [g_bfs, g_fa]), CEMDATA_PRIMARIES; solid_solutions = ss)
    return (; cs, mats)
end

"""
    s15_mix(name) -> NamedTuple

The row of Table 4: the mass percent of OPC, slag, fly ash and limestone.
"""
function s15_mix(name)
    t = literature_table(S15, "mixes")
    k = findfirst(==(name), t.mix)
    k === nothing && throw(ArgumentError("Table 4 has no mix $name; it has " * join(t.mix, ", ")))
    return (; opc = ustrip(t.OPC[k]) / 100, bfs = ustrip(t.BFS[k]) / 100, fa = ustrip(t.FA[k]) / 100, ls = ustrip(t.LS[k]) / 100)
end

"""
    s15_anhydrite(mix) -> Float64

The anhydrite (g per 100 g of the binder of Table 4) that brings the SO₃ of the
whole to the 3 % the authors set (Section 2.1), from the SO₃ of each analysis
and that of anhydrite, from its formula.

ASSUMED: the anhydrite is added to the binder of Table 4, whose four materials
make up 100 %, rather than counted in it, and the water/binder ratio is that of
those 100 %. The article states the target and not the dose.
"""
function s15_anhydrite(mix)
    so3(m) = get(literature_oxides(S15, "oxides", m), "SO3", 0.0)
    s = 100 * (mix.opc * so3("OPC") + mix.bfs * so3("BFS") + mix.fa * so3("FA") + mix.ls * so3("LS"))
    target = s15_value("so3_of_the_binder") / 100
    f = oxide_content(S15_DB["Anh"], ["SO3"])["SO3"]
    return max(target * 100 - s, 0.0) / (f - target)
end

"""
    s15_recipe(mats, mix) -> Recipe

100 g of the binder of `mix` at the water/binder ratio of the article, with its
anhydrite as an addition. Nothing has reacted yet except what is at
equilibrium from the mixing.
"""
function s15_recipe(mats, mix)
    parts = Pair{Material, Float64}[mats.opc => mix.opc]
    mix.bfs > 0 && push!(parts, mats.bfs => mix.bfs)
    mix.fa > 0 && push!(parts, mats.fa => mix.fa)
    mix.ls > 0 && push!(parts, mats.ls => mix.ls)
    return Recipe(
        parts...; w_b = ustrip(s15_value("water_binder_ratio")),
        additions = [mats.anhydrite => s15_anhydrite(mix) * u"g"],
    )
end

"""
    s15_glass_time(which; n = WALLER_PARAMS_FLY_ASH.n) -> Quantity

The characteristic time of the Waller law for the glass of the slag
(`:bfs`) or of the fly ash (`:fa`), **calibrated** on the long-term degree of
reaction the authors assume for it in their own calculations, 71.1 % and
43.6 % of the glass after one year (Section 2.2): the closed form of the law
at its reference temperature, `α = 1/(1 + (τ/t)ⁿ)`, solved for `τ`. Those
degrees are not measurements of these pastes, and the comparison is made on
other observables, at other ages. The exponent `n` is the one Lavergne et al.
(2018) fitted for a fly ash: ASSUMED for both glasses.
"""
function s15_glass_time(which; n = WALLER_PARAMS_FLY_ASH.n)
    key = which === :bfs ? "assumed_reaction_slag_glass" : which === :fa ? "assumed_reaction_fly_ash_glass" :
        throw(ArgumentError("s15_glass_time: :bfs or :fa, got $which"))
    α = ustrip(s15_value(key)) / 100
    t = ustrip(us"d", s15_value("assumed_reaction_age"))
    return t * (1 / α - 1)^(1 / n) * u"d"
end

"""
    s15_rates(mix) -> Dict

Parrott–Killoh for the four clinker phases at the cement's Blaine fineness and
at the water/cement ratio of the paste, `w/b` over the cement's share; the
Waller law for each glass present, with its calibrated time and no fineness
correction, the calibration being made for that material.
"""
function s15_rates(mix)
    blaine = literature_row(S15, "oxides", "OPC").blaine
    w_c = ustrip(s15_value("water_binder_ratio")) / mix.opc
    pk = Dict("C3S" => PK84_PARAMS_C3S, "C2S" => PK84_PARAMS_C2S, "C3A" => PK84_PARAMS_C3A, "C4AF" => PK84_PARAMS_C4AF)
    rates = Dict{String, Any}(ph => parrott_killoh_avrami(p, ph; blaine, w_c) for (ph, p) in pk)
    mix.bfs > 0 && (rates["BFS"] = waller(merge(WALLER_PARAMS_FLY_ASH, (τ = s15_glass_time(:bfs),)), "BFS"))
    mix.fa > 0 && (rates["FA"] = waller(merge(WALLER_PARAMS_FLY_ASH, (τ = s15_glass_time(:fa),)), "FA"))
    return rates
end

"""
    s15_system(cs, mix) -> ChemicalSystem

The system of `s15_setup` without the glass of an addition the mix does not
hold: that glass has no rate law in the mix, and a species without one is at
equilibrium, which a glass, known by its formula and not by its standard Gibbs
energy, cannot be.
"""
function s15_system(cs, mix)
    absent = [g for (g, x) in (("BFS", mix.bfs), ("FA", mix.fa)) if x == 0]
    isempty(absent) && return cs
    return ChemicalSystem(
        [s for s in cs.species if !(symbol(s) in absent)], CEMDATA_PRIMARIES;
        solid_solutions = cs.solid_solutions,
    )
end

"""
    s15_run(cs, mats, name; days = 182) -> (; kp, sol, recipe, mix)

The paste `name` of Table 4 integrated over `days` at 20 °C, in the activity
model the authors used, Cemdata18's extended Debye–Hückel equation for a KOH
solution ([`cemdata18_activity_model`](@ref)).
"""
function s15_run(cs, mats, name; days = 182)
    mix = s15_mix(name)
    recipe = s15_recipe(mats, mix)
    cs = s15_system(cs, mix)
    model = cemdata18_activity_model(:KOH)
    kp = KineticsProblem(
        recipe, cs, s15_rates(mix), (0.0, days * 86400.0);
        activity_model = model,
        equilibrium_solver = EquilibriumSolver(cs, model, OptimaOptimizer()),
    )
    sol = integrate(kp, KineticsSolver(; ode_solver = Rodas5P(), reltol = 1.0e-6, abstol = 1.0e-10))
    return (; kp, sol, recipe, mix)
end

"""
    s15_tga(run, days) -> Vector{(; bound_water, portlandite)}

What Table 8 reports, computed at each of `days`: the bound water without that
of portlandite and the portlandite, in percent of the sample dried at 500 °C
(their Eqs. 3 and 4). The sample is every solid of the replayed state and the
inert crystals of the recipe; the dry sample is the same less all its bound
water.

ASSUMED: every hydrate has lost all its water at 500 °C and nothing else has
been lost. The thermobalance sees less: the sample was dried at 40 °C after a
solvent exchange, which takes some of the water of the C-S-H and the AFm
before the run starts, and hydrotalcite keeps part of its hydroxyls past
500 °C.
"""
function s15_tga(run, days)
    inert = sum(x.mass for x in budget(run.recipe, run.kp.system).residual; init = 0.0)
    g(q) = ustrip(uconvert(us"g", q))
    states = speciated_states(run.sol, run.kp; times = collect(float.(days)) .* 86400.0)
    return map(states) do st
        water = g(bound_water(st))
        per = Dict(bound_water_per_phase(st))
        ch_water = g(get(per, "Portlandite", 0.0u"g"))
        dry = g(mass(st).solid) + inert - water
        m_ch = g(mass(st, st.system.dict_species["Portlandite"]))
        (; bound_water = 100 * (water - ch_water) / dry, portlandite = 100 * m_ch / dry)
    end
end

"""
    s15_measured(name) -> (; days, bound_water, portlandite)

Table 8 for the mix `name`, in percent of the sample dried at 500 °C.
"""
function s15_measured(name)
    t = literature_table(S15, "tga"; mix = name)
    return (; days = ustrip.(us"d", t.age_days), bound_water = ustrip.(t.bound_water), portlandite = ustrip.(t.portlandite))
end
