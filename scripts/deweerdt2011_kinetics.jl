# =============================================================================
#  deweerdt2011_kinetics.jl — the pastes of De Weerdt et al. (2011), in time
#
#  The four pastes of scripts/de_weerdt_2011.jl (a clinker interground with
#  gypsum, blended with limestone powder, siliceous fly ash or both, at w/b 0.5
#  and 20 °C) integrated from the mixing instead of computed at measured
#  extents: the four clinker phases under the law of Parrott and Killoh, with the
#  parameters Lothenbach et al. (2008) report and the fineness of Table 1; the
#  glass of the fly ash at the degree of reaction the authors measured, through
#  the fit printed on their Fig. 7; everything else at equilibrium. The page
#  docs/src/examples/ternary_kinetics.md compares the runs with what Table 7
#  measured over six months: the clinker phases left, the portlandite and the
#  ettringite, none of which was fitted.
# =============================================================================

isdefined(@__MODULE__, :dw11_recipe) || include(joinpath(@__DIR__, "de_weerdt_2011.jl"))
using OrdinaryDiffEq
using ForwardDiff: ForwardDiff
using LinearAlgebra: Diagonal, diag

const DW11_CLINKER = ("C3S", "C2S", "C3A", "C4AF")
const DW11_PK = Dict(
    "C3S" => PK84_PARAMS_C3S, "C2S" => PK84_PARAMS_C2S,
    "C3A" => PK84_PARAMS_C3A, "C4AF" => PK84_PARAMS_C4AF,
)

# The database of the CASH+ models, read the first time a run asks for it.
const _DW11K_CASHPLUS = Ref{Any}(nothing)
function _dw11k_cashplus()
    _DW11K_CASHPLUS[] === nothing &&
        (_DW11K_CASHPLUS[] = build_species(datapath("cemdata18-cashplus.json"); verbose = false))
    return _DW11K_CASHPLUS[]
end

"""
    dw11k_setup(; gel = "CSHQ") -> (; cs, fly_ash, glass)

The system of the phase list of `dw11_system`, which holds the clinker phases,
with the glass of the fly ash added as a species ([`glass_species`](@ref)), and
the fly ash whose glass is that species ([`with_species`](@ref)), its crystals
inert. `gel` is the model of the C-S-H: `"CSHQ"`, or one of the two the
validation page also computes these pastes with, `"CNASH_ss"` (on Cemdata18)
and `"CASH+NK"` (on `cemdata18-cashplus.json`).
"""
function dw11k_setup(; gel = "CSHQ")
    base = if gel == "CSHQ"
        dw11_system()
    elseif gel == "CNASH_ss"
        phase_list_system(DW11_PHASES, DW11_SUBSTANCES; replace = Dict("CSHQ" => "CNASH_ss"))
    elseif gel == "CASH+NK"
        phase_list_system(DW11_PHASES, _dw11k_cashplus(); replace = Dict("CSHQ" => "CASH+NK"))
    else
        throw(ArgumentError("dw11k_setup: the gel is \"CSHQ\", \"CNASH_ss\" or \"CASH+NK\"; got \"$gel\""))
    end
    fa = material_template("siliceous fly ash (De Weerdt 2011)", DW11_DB)
    glass = glass_species(only(c for c in fa.constituents if c.name == "glass"), base; symbol = "FA")
    fa = with_extents(fa, Dict(c.name => 0.0 for c in fa.constituents if c.name != "glass"))
    cs = ChemicalSystem(vcat(base.species, [glass]), CEMDATA_PRIMARIES; solid_solutions = base.solid_solutions)
    return (; cs, fly_ash = with_species(fa, Dict("glass" => glass)), glass)
end

"""
    dw11k_system(setup, mix) -> ChemicalSystem

The system of `mix`: that of `dw11k_setup`, without the glass of the fly ash in
the two pastes that hold none. A glass has no standard Gibbs energy, and left
without a rate law it would be at equilibrium, which `KineticsProblem` refuses.
"""
function dw11k_system(setup, mix)
    has_fa = ustrip(dw11_table("mixes").fly_ash[findfirst(==(mix), dw11_table("mixes").mix)]) > 0
    has_fa && return setup.cs
    return ChemicalSystem(
        [s for s in setup.cs.species if symbol(s) != "FA"], CEMDATA_PRIMARIES;
        solid_solutions = setup.cs.solid_solutions,
    )
end

"""
    dw11k_recipe(setup, mix) -> Recipe

The paste of `mix` as `dw11_recipe` builds it, but with every clinker phase
whole at the mixing, for its law to dissolve, and the fly ash of `setup`.

ASSUMED: the minor oxides of the clinker (free lime, the alkalis, magnesia, its
sulfate) are in the equilibrium from the mixing. The page at measured extents
releases them with the clinker as a whole; here a constituent without a rate law
enters at its extent at the start, so they enter whole, which puts the alkalis
in the pore solution early.
"""
function dw11k_recipe(setup, mix)
    m = dw11_table("mixes")
    i = findfirst(==(mix), m.mix)
    opc, fa, ls = ustrip(m.opc[i]) / 100, ustrip(m.fly_ash[i]) / 100, ustrip(m.limestone[i]) / 100
    g = dw11_value("gypsum_percent_of_opc") / 100
    # The clinker phases enter whole whatever their extent, `KineticsProblem`
    # taking back what `budget` would have reacted; their extent is left as the
    # template has it, so that the residue the recipe sets aside holds none of them.
    clinker = with_extents(material_template("clinker (De Weerdt 2011)", DW11_DB), Dict("minor oxides" => 1.0))
    parts = Pair{Material, Float64}[clinker => opc * (1 - g), dw11_gypsum() => opc * g]
    fa > 0 && push!(parts, setup.fly_ash => fa)
    ls > 0 && push!(parts, dw11_limestone() => ls)
    return Recipe(parts...; w_b = dw11_value("water_binder_ratio"), T = dw11_value("curing_temperature") * u"K")
end

"""
    dw11k_fly_ash_law(share) -> KineticFunc

The glass of the fly ash reacting as the authors measured it: the fit printed on
their Fig. 7, `y = a + b ln(t + c)` percent of the fly ash reacted at `t` days,
written on the state rather than on the time, as the Waller law is. With `α` the
degree of reaction of the glass, `share` the glass's mass fraction of the fly
ash (the authors take the glass alone to react), `y = 100 share α`, and

    dα/dt = b / (100 share) · exp(-(100 share α - a) / b)   per day,

since `t + c = exp((y - a) / b)`. At 20 °C only, the temperature of the fit.
"""
function dw11k_fly_ash_law(share)
    a, b = dw11_value("fly_ash_fit_a"), dw11_value("fly_ash_fit_b")
    day = 86400.0
    f = (_T, _P, _t, n, _lna, n_initial) -> begin
        n_init = max(n_initial["FA"], oneunit(n["FA"]) * 1.0e-30)
        α = max(one(n_init) - n["FA"] / n_init, zero(n_init))
        y = 100 * share * α
        return n["FA"] > 0 ? n_init * b / (100 * share) * exp(-(y - a) / b) / day : zero(y)
    end
    return KineticFunc(f, (T = dw11_value("curing_temperature") * u"K", P = 1.0e5u"Pa"), u"mol/s")
end

# The names Lothenbach et al. (2008) give the clinker phases in their Table 3.
const DW11_PHASE_NAMES_2008 = Dict("C3S" => "alite", "C2S" => "belite", "C3A" => "aluminate", "C4AF" => "ferrite")

"""
    dw11k_clinker_law(phase; blaine, w_c, pk = :lavergne2018) -> KineticFunc

The Parrott–Killoh law of a clinker phase, with one of two published parameter
sets:

  - `:lavergne2018`, the constants Lavergne et al. (2018) tabulate
    (`PK84_PARAMS_*`, those of Lothenbach and Winnefeld 2006), and the
    water/cement factor that slows every phase past `1.333 w/c`;
  - `:lothenbach2008`, the constants of Lothenbach et al. (2008, Table 3), whose
    belite has larger `K₂` and `K₃`, and the critical degree `H` of each phase
    in the same factor (`pk_wc_factor`), which lets the alite past `0.9` at
    w/c = 0.5 where the other set stops it at `0.67`.

The activation energies are those of `PK84_PARAMS_*` in both; at 20 °C, the
temperature of these pastes, they do not act.
"""
function dw11k_clinker_law(phase; blaine, w_c, pk = :lavergne2018, alite = nothing)
    if alite !== nothing && phase == "C3S"
        # The alite of `dw11k_alite_fit`: its constants over the published ones.
        return parrott_killoh_avrami(_dw11k_alite_params(alite), phase; blaine, w_c, H = get(alite, :H, nothing))
    end
    pk === :lavergne2018 && return parrott_killoh_avrami(DW11_PK[phase], phase; blaine, w_c)
    pk === :lothenbach2008 || throw(ArgumentError("pk is :lavergne2018 or :lothenbach2008; got :$pk"))
    r = literature_row("LothenbachLeSaout2008", "parrott_killoh", DW11_PHASE_NAMES_2008[phase])
    per_day(x) = ustrip(x) * u"1/d"
    params = merge(
        DW11_PK[phase],
        (k₁ = per_day(r.K1), n₁ = ustrip(r.N1), k₂ = per_day(r.K2), k₃ = per_day(r.K3), n₃ = ustrip(r.N3)),
    )
    return parrott_killoh_avrami(params, phase; blaine, w_c, H = ustrip(r.H))
end

"""
    dw11k_rates(setup, mix; pk = :lavergne2018) -> Dict

Parrott–Killoh for the four clinker phases at the fineness of the OPC (Table 1)
and the water/clinker ratio of the paste, `w/b` over the OPC's share, with the
parameter set `pk` (`dw11k_clinker_law`); the law of the fly-ash glass in the
pastes that hold it.
"""
function dw11k_rates(setup, mix; pk = :lavergne2018, alite = nothing)
    m = dw11_table("mixes")
    i = findfirst(==(mix), m.mix)
    opc, fa = ustrip(m.opc[i]) / 100, ustrip(m.fly_ash[i]) / 100
    blaine = literature_value(DW11, "blaine_opc")
    w_c = dw11_value("water_binder_ratio") / opc
    rates = Dict{String, Any}(p => dw11k_clinker_law(p; blaine, w_c, pk, alite) for p in DW11_CLINKER)
    if fa > 0
        # `with_species` names the constituent after its species.
        share = only(c.mass_fraction for c in setup.fly_ash.constituents if c.name == "FA")
        rates["FA"] = dw11k_fly_ash_law(share)
    end
    return rates
end

"""
    dw11k_run(setup, mix; days = 180, pk = :lavergne2018, alite = nothing) -> (; kp, sol, recipe)

The paste `mix` integrated over `days` at 20 °C, in the activity model of the
page at measured extents, Cemdata18's for a KOH solution, the clinker under the
parameter set `pk` (`dw11k_clinker_law`), the alite under the constants
`alite` of `dw11k_alite_fit` when given.
"""
function dw11k_run(setup, mix; days = 180, pk = :lavergne2018, alite = nothing)
    cs = dw11k_system(setup, mix)
    recipe = dw11k_recipe(setup, mix)
    model = cemdata18_activity_model(:KOH)
    kp = KineticsProblem(
        recipe, cs, dw11k_rates(setup, mix; pk, alite), (0.0, days * 86400.0);
        activity_model = model,
        equilibrium_solver = EquilibriumSolver(cs, model, OptimaOptimizer()),
    )
    sol = integrate(kp, KineticsSolver(; ode_solver = Rodas5P(), reltol = 1.0e-6, abstol = 1.0e-10))
    return (; kp, sol, recipe)
end

"""
    dw11k_replay(run, days) -> Vector{ChemicalState}

The certified replay of `run` at each of `days` ([`speciated_states`](@ref)).
"""
dw11k_replay(run, days) = speciated_states(run.sol, run.kp; times = collect(float.(days)) .* 86400.0)

"""
    dw11k_contents(run, states) -> Vector{Dict{String, Float64}}

What Table 7 reports, computed on `states`, the replay of `run`: the content of
portlandite, ettringite and the four clinker phases in wt.% of the solids of the
paste (every solid with its water, and the unreacted part of the binder set
aside by the recipe), and `opc_reacted`, the table's own measure: one minus the
sum of the contents of the four clinker phases over that sum at the mixing, when
the solids are the binder. Like the table's, it counts the bound water that
dilutes the clinker as clinker reacted.
"""
function dw11k_contents(run, states)
    inert = sum(x.mass for x in budget(run.recipe, run.kp.system).residual; init = 0.0)
    g(q) = ustrip(uconvert(us"g", q))
    sp = run.kp.system.dict_species
    s0 = run.kp.initial_state
    clinker0 = 100 * sum(g(mass(s0, sp[p])) for p in DW11_CLINKER) / run.recipe.binder_mass   # grams
    labels = vcat(["portlandite" => "Portlandite", "ettringite" => "ettringite"], [p => p for p in DW11_CLINKER])
    return map(states) do st
        solids = g(mass(st).solid) + inert
        c = Dict(q => 100 * g(mass(st, sp[s])) / solids for (q, s) in labels)
        c["opc_reacted"] = 100 * (1 - sum(c[p] for p in DW11_CLINKER) / clinker0)
        c
    end
end

"""
    dw11k_amount(state, name) -> Float64

The moles of the species `name` in `state`, those of its later instances
(`name#2`, …) included.
"""
dw11k_amount(state, name) = sum(
    (
        ustrip(us"mol", state.n[i]) for (i, s) in enumerate(state.system.species)
            if symbol(s) == name || startswith(symbol(s), name * "#")
    );
    init = 0.0,
)

"""
    dw11k_gel(state, gel) -> (; Ca_Si, Al_Si)

The calcium-to-silicon and aluminum-to-silicon ratios of the C-S-H `gel`, the
name of its solid solution, in `state`.
"""
function dw11k_gel(state, gel)
    e = solid_solution_totals(state, gel).elements
    return (; Ca_Si = e[:Ca] / e[:Si], Al_Si = get(e, :Al, 0.0) / e[:Si])
end

# ── The alite calibrated on the plain cement ─────────────────────────────────

# The published constants of the alite, rate constants per day, and the critical
# degree 1.333 of the water/cement factor the published set uses.
_dw11k_alite_published() = (
    k₁ = ustrip(us"1/d", PK84_PARAMS_C3S.k₁), n₁ = PK84_PARAMS_C3S.n₁, k₂ = ustrip(us"1/d", PK84_PARAMS_C3S.k₂),
    k₃ = ustrip(us"1/d", PK84_PARAMS_C3S.k₃), n₃ = PK84_PARAMS_C3S.n₃, H = ustrip(literature_value("LothenbachWinnefeld2006", "pk_wc_critical_degree")),
)

# The constants of the law from a set of them (rate constants per day), H aside.
function _dw11k_alite_params(θ)
    q = merge(_dw11k_alite_published(), θ)
    return merge(
        PK84_PARAMS_C3S,
        (k₁ = q.k₁ * u"1/d", n₁ = q.n₁, k₂ = q.k₂ * u"1/d", k₃ = q.k₃ * u"1/d", n₃ = q.n₃),
    )
end

"""
    dw11k_alite_fit(run; spec = (:k₂, :k₃, :n₃, :H)) -> NamedTuple

The Parrott–Killoh law of the alite fitted on the alite of the plain cement,
Table 7 at its five ages, the constants of `spec` free and the others published:
Levenberg–Marquardt on their logarithms, with the exact Jacobian
(`ForwardDiff`). The law is integrated on its own (`ParrottKillohExtent`), its
degree of reaction turned into a content with the solids of `run`, the plain
cement integrated with the published set, which the alite changes little.

Returns `θ` (the fitted constants, rate constants per day), `published`, the
root-mean-square misfit in wt.% before and after (`rms_published`, `rms`), the
contents at the five ages and the measured ones, and the
[`identifiability`](@ref) of `spec` at the fit, for a measurement good to
±1 wt.%.
"""
function dw11k_alite_fit(run; spec = (:k₂, :k₃, :n₃, :H))
    days = [1, 7, 28, 90, 180]
    states = dw11k_replay(run, days)
    inert = sum(x.mass for x in budget(run.recipe, run.kp.system).residual; init = 0.0)
    g(q) = ustrip(uconvert(us"g", q))
    solids = [g(mass(st).solid) + inert for st in states]
    m0 = g(mass(run.kp.initial_state, run.kp.system.dict_species["C3S"]))
    blaine = literature_value(DW11, "blaine_opc")
    w_c = dw11_value("water_binder_ratio")
    measured = Float64[dw11_phase_content("OPC", d, "C3S") for d in days]
    published = _dw11k_alite_published()
    function contents(θ)
        q = NamedTuple{spec}(Tuple(θ))
        e = ParrottKillohExtent(
            "C3S"; blaine, w_c, H = get(q, :H, published.H),
            parameters = (; (k => getproperty(merge(published, q), k) for k in (:k₁, :n₁, :k₂, :k₃, :n₃))...),
        )
        return [100 * m0 * (1 - extent(e, d)) / solids[k] for (k, d) in enumerate(days)]
    end
    residual(z) = contents(exp.(z)) .- measured
    z = log.([getproperty(published, k) for k in spec])
    r = residual(z)
    f = sum(abs2, r)
    f0 = f
    λ = 1.0e-2
    for _ in 1:80
        J = ForwardDiff.jacobian(residual, z)
        A, gr = J' * J, J' * r
        accepted = false
        for _ in 1:25
            zn = z .- (A + λ * Diagonal(diag(A) .+ 1.0e-12)) \ gr
            rn = residual(zn)
            fn = sum(abs2, rn)
            if fn < f
                z, r, f, accepted = zn, rn, fn, true
                λ = max(λ / 3, 1.0e-9)
                break
            end
            λ *= 4
        end
        accepted || break
    end
    θ = NamedTuple{spec}(Tuple(exp.(z)))
    id = identifiability(contents, collect(values(θ)); observed = measured, noise = 1.0, names = string.(collect(spec)))
    return (;
        θ, published = NamedTuple{spec}(Tuple(getproperty(published, k) for k in spec)),
        rms_published = sqrt(f0 / length(days)), rms = sqrt(f / length(days)),
        contents = contents(collect(values(θ))), measured, days, identifiability = id,
    )
end
