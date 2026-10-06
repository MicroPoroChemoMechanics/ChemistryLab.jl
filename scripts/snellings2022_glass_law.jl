# =============================================================================
#  snellings2022_glass_law.jl — the slag of a ternary paste dissolved by the law
#  of its glass measured in dilute solutions
#
#  The pastes of scripts/snellings2022_pastes.jl (50 % CEM I 52.5 R, 40 % slag,
#  10 % limestone, w/b 0.5, at 5, 20 and 40 °C), with the glass of the slag
#  dissolved by the far-from-equilibrium rate of Snellings (2013) instead of
#  the Waller law fitted on these pastes: the rate per unit area from the
#  composition of the glass (`snellings2013_glass`), slowed by the calcium of
#  the pore solution with the factor fitted on the dilute solutions of that
#  paper (scripts/snellings2013_glass.jl), over the BET surface the authors
#  measured on the slag, shrinking as the glass dissolves. Nothing is fitted
#  on the pastes: the page docs/src/examples/glass_dissolution.md compares the
#  degree of reaction with the one Snellings et al. (2022) measured.
#
#  ASSUMED, each where it is made: the BET surface of the slag is that of its
#  glass, and the grains shrink as spheres (exponent 2/3); the activation
#  energy is the apparent one the authors fitted on the degree of reaction of
#  the slag at w/b 0.5 (their Table 2), the dilute experiments being at 20 °C
#  only; the glass is percalcic (CaO/Al2O3 far above one, as G1 and G2), so
#  aluminum does not slow it; its MgO is left out of the composition law,
#  which was measured on glasses without it, and the law is used beyond the
#  most calcic glass measured.
# =============================================================================

isdefined(@__MODULE__, :sn22p_setup) || include(joinpath(@__DIR__, "snellings2022_pastes.jl"))
isdefined(@__MODULE__, :sn13_rows) || include(joinpath(@__DIR__, "snellings2013_glass.jl"))

"""
    sn22g_setup() -> (; cs, pc, slag, ls, glass, oxides)

The system and the materials of [`sn22p_setup`](@ref), with the glass of the
slag counted per mole of its cations ([`cation_molar_mass`](@ref)), the unit of
the rate of Snellings (2013); `oxides` is the part of the analysis of the glass
its formula holds.
"""
function sn22g_setup()
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
    c = only(c for c in slag.constituents if c.name == "glass")
    oxides = Dict(k => v for (k, v) in c.oxides if ChemistryLab._representable(k, base.SM.primaries))
    glass = glass_species(c, base; symbol = "SLAG", M = cation_molar_mass(oxides))
    slag = with_species(slag, Dict("glass" => glass))
    cs = ChemicalSystem(vcat(sp, [glass]), CEMDATA_PRIMARIES; solid_solutions = ss)
    return (; cs, pc, slag, ls, glass, oxides)
end

"""
    sn22g_slag_law(setup, fit; w_b = 0.5, modifiers = ("CaO",), area = :bet) -> KineticFunc

The rate of the glass of the slag: [`snellings2013_glass`](@ref) of its
analysis, with the activation energy of Table 2 of Snellings et al. (2022) at
`w_b` and the calcium inhibitor of `fit` ([`sn13_fit`](@ref)), over the surface
of the slag shrinking as spheres: its BET surface (`area = :bet`, the
measurement the rates of the glasses were normalized by) or its Blaine surface
(`area = :blaine`, the area of a smooth grain, about five times smaller).
"""
function sn22g_slag_law(setup, fit; w_b = 0.5, modifiers = ("CaO",), area = :bet)
    Ea = only(literature_table(SN22, "activation_energies"; constituent = "slag", w_b).Ea)
    mech = snellings2013_glass(setup.oxides; Ea, modifiers, inhibitors = sn13_inhibitors(fit, :percalcic), extrapolate = true)
    cs = setup.cs
    A = Float64.(cs.SM.A)
    j = findfirst(s -> symbol(s) == "SLAG", cs.species)
    products = [cs.dict_species[symbol(p)] for (q, p) in enumerate(cs.SM.primaries) if !iszero(A[q, j])]
    rxn = Reaction([cs.species[j]], products; symbol = "SLAG dissolution")
    s = area === :bet ? BETSurfaceArea(literature_value(SN22, "bet_slag")) :
        area === :blaine ? BlaineSurfaceArea(literature_value(SN22, "blaine_slag")) :
        throw(ArgumentError("sn22g_slag_law: area is :bet or :blaine; got $area."))
    return transition_state([mech], cs, rxn, ShrinkingCoreArea(s; exponent = 2 // 3))
end

"""
    sn22g_rates(setup, fit; w_b = 0.5) -> Dict

The rate laws of [`sn22p_rates`](@ref) for the clinker, and
[`sn22g_slag_law`](@ref) for the glass of the slag.
"""
function sn22g_rates(setup, fit; w_b = 0.5, area = :bet)
    blaine = literature_value(SN22, "blaine_pc")
    w_c = w_b / (sn22_value("pc_percent") / 100)
    rates = Dict{String, Any}()
    for c in setup.pc.constituents
        c isa MineralConstituent || continue
        s = symbol(c.species)
        haskey(SN22_PK, s) && (rates[c.name] = parrott_killoh_avrami(SN22_PK[s], s; blaine, w_c))
    end
    rates["SLAG"] = sn22g_slag_law(setup, fit; w_b, area)
    return rates
end

"""
    sn22g_run(setup, fit, T_C; w_b = 0.5, days = 180) -> (; kp, sol, recipe, T_C)

The paste of [`sn22p_run`](@ref) at `T_C` °C with the slag under
[`sn22g_slag_law`](@ref) and the clinker as there.
"""
function sn22g_run(setup, fit, T_C; w_b = 0.5, days = 180, area = :bet)
    recipe = sn22p_recipe(setup, T_C; w_b)
    model = cemdata18_activity_model(:KOH)
    kp = KineticsProblem(
        recipe, setup.cs, sn22g_rates(setup, fit; w_b, area), (0.0, days * 86400.0);
        activity_model = model,
        equilibrium_solver = EquilibriumSolver(setup.cs, model, OptimaOptimizer()),
    )
    sol = integrate(kp, KineticsSolver(; ode_solver = Rodas5P(), reltol = 1.0e-6, abstol = 1.0e-10))
    return (; kp, sol, recipe, T_C)
end

"""
    sn22g_slag_degree(run, days) -> Vector

The degree of reaction of the glass of the slag (percent) at `days`, from the
amount of the kinetic species in the integrated state.
"""
function sn22g_slag_degree(run, days)
    p = run.sol.prob.p
    i = findfirst(s -> symbol(s) == "SLAG", run.kp.system.species)
    j = findfirst(==(i), p.idx_kinetic)
    n0 = ustrip(us"mol", run.kp.initial_state.n[i])
    return [100 * (1 - run.sol(d * 86400.0)[p.n_be + j] / n0) for d in days]
end
