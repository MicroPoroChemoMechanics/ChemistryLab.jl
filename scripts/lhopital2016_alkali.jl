# =============================================================================
#  lhopital2016_alkali.jl — alkali uptake by C-A-S-H, L'Hôpital et al. (2016)
#
#  L'Hôpital, Lothenbach, Scrivener and Kulik (2016) synthesized C-S-H and
#  C-A-S-H (Al/Si 0.05) at a Ca/Si from 0.6 to 1.6, 2 g of lime, silica fume and
#  monocalcium aluminate in 90 mL of water or of a KOH or NaOH solution from
#  0.01 to 0.5 mol/L, equilibrated each at 20 °C for 91, 182 or 364 days, and
#  reported the solutions (Appendix B) and the C-S-H, its alkali taken from the
#  fall of the dissolved alkali (Appendix C). This file computes the same
#  batches at equilibrium with two models of the gel, CSHQ with its alkali end
#  members (Cemdata18) and CASH+NK (Miron et al. 2022a), the aluminum free to
#  form the hydrates of Cemdata18 beside either. The page
#  docs/src/tutorials/validation_alkali_uptake.md compares them.
#
#  Read the circularity first. The alkali end members of CSHQ were fitted on
#  the isotherms of Hong and Glasser (1999), not on these: for CSHQ the
#  comparison is a prediction. CASH+NK was fitted on these very data, among
#  others (Miron et al. 2022a, Table 4): for CASH+NK it checks that the model
#  is applied as its authors fitted it.
# =============================================================================

using ChemistryLab
using DynamicQuantities
using OptimaSolver

const LH16 = "LHopital2016b"
lh16_value(name) = ustrip(literature_value(LH16, name))
lh16_table(name; where...) = literature_table(LH16, name; where...)

const LH16_DB = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-thermofun.json"); verbose = false))
# The database of the CASH+ models, read the first time a gel asks for it.
const _LH16_CASHPLUS = Ref{Any}(nothing)
function _lh16_cashplus()
    _LH16_CASHPLUS[] === nothing &&
        (_LH16_CASHPLUS[] = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-cashplus.json"); verbose = false)))
    return _LH16_CASHPLUS[]
end

# The phases: the three reactants, which dissolve, portlandite and amorphous
# silica, and the calcium aluminate hydrates of Cemdata18, hydrogarnet and its
# siliceous member, strätlingite and microcrystalline gibbsite among them, for
# the aluminum the gel does not take.
const LH16_PURE = split("Lim Amor-Sl CA Portlandite C3AH6 C3AS0.84H4.32 straetlingite C4AH13 CAH10 AlOHmic")
const LH16_CSHQ = ["CSHQ-JenD", "CSHQ-JenH", "CSHQ-TobD", "CSHQ-TobH"]
const LH16_ALKALI = Dict("KOH" => ("K+", "KSiOH"), "NaOH" => ("Na+", "NaSiOH"))
# The aqueous ion pairs Miron et al. (2022a) left out when fitting CASH+NK
# (their Sections 3.2 and 7.4); of them the database holds NaOH@ and KOH@.
const LH16_CASHPLUS_EXCLUDED = ["NaOH@", "KOH@", "NaHSiO3@", "KHSiO3@"]

"""
    lh16_system(gel, alkali) -> ChemicalSystem

The system of a batch in `alkali` (`"KOH"`, `"NaOH"` or `"none"`) with the gel
`gel`: `"CSHQ"`, with the alkali end member of that hydroxide, or `"CASH+NK"`,
on `cemdata18-cashplus.json`, whose twelve end members carry both alkalis.

For `CASH+NK` the aqueous ion pairs of the alkalis are left out
(`LH16_CASHPLUS_EXCLUDED`), and `Ca(OH)2@` is kept, as Miron et al. (2022a)
fitted the model (their Sections 3.2 and 7.4). For `CSHQ` a system holds the alkali of its batch only: a budget without
potassium would leave `KSiOH` and the potassium species at the floor.
"""
function lh16_system(gel, alkali)
    if gel == "CSHQ"
        db = LH16_DB
        members = alkali == "none" ? LH16_CSHQ : vcat(LH16_CSHQ, [last(LH16_ALKALI[alkali])])
        ss = [SolidSolutionPhase("CSHQ", [db[m] for m in members])]
        sp = speciation(
            collect(values(db)), vcat(LH16_PURE, members);
            aggregate_state = [AS_AQUEOUS], exclude_species = split("H2@ O2@ CH4@"),
        )
        return ChemicalSystem(sp, CEMDATA_PRIMARIES; solid_solutions = ss)
    elseif gel == "CASH+NK"
        db = _lh16_cashplus()
        nk = only(p for p in build_solid_solutions(datapath("solid_solutions.toml"), db) if name(p) == "CASH+NK")
        members = [symbol(m) for m in nk.end_members]
        sp = speciation(
            collect(values(db)), vcat(LH16_PURE, members);
            aggregate_state = [AS_AQUEOUS], exclude_species = vcat(split("H2@ O2@ CH4@"), LH16_CASHPLUS_EXCLUDED),
        )
        return ChemicalSystem(sp, CEMDATA_PRIMARIES; solid_solutions = [nk])
    end
    throw(ArgumentError("lh16_system: the gel is \"CSHQ\" or \"CASH+NK\"; got \"$gel\""))
end

"""
    lh16_state(cs, Al_Si, Ca_Si, alkali, c) -> ChemicalState

The batch of Appendix A at the target `Al_Si` and `Ca_Si`, its lime, silica fume
and monocalcium aluminate, in 90 mL of `alkali` at `c` mol/L, at 20 °C.

ASSUMED: the 90 mL of solution are 90 g of water, the hydroxide dissolved in
them; the solids are the oxides they are written as, the silica fume amorphous
silica.
"""
function lh16_state(cs, Al_Si, Ca_Si, alkali, c)
    mix = lh16_table("mixing_proportions"; Al_Si_target = Al_Si, Ca_Si_target = Ca_Si)
    g(q) = ustrip(us"g", only(q))
    st = ChemicalState(cs; T = literature_value(LH16, "temperature"))
    V_mL = ustrip(u"mL", literature_value(LH16, "solution_volume"))
    set_quantity!(st, "H2O@", V_mL * u"g")
    set_quantity!(st, "Lim", g(mix.CaO) * u"g")
    set_quantity!(st, "Amor-Sl", g(mix.SiO2) * u"g")
    g(mix.CaO_Al2O3) > 0 && set_quantity!(st, "CA", g(mix.CaO_Al2O3) * u"g")
    if alkali != "none"
        n = c * V_mL / 1000
        set_quantity!(st, first(LH16_ALKALI[alkali]), n * u"mol")
        set_quantity!(st, "OH-", n * u"mol")
    end
    return st
end

lh16_model(alkali) = cemdata18_activity_model(alkali == "NaOH" ? :NaOH : :KOH)

"""
    lh16_observables(eq, gel, alkali) -> NamedTuple

What Appendices B and C report, on the equilibrium `eq`: the dissolved Si, Ca,
Al and alkali in mmol/L, the pH, and the Ca/Si and alkali/Si of the gel.
"""
function lh16_observables(eq, gel, alkali)
    cs = eq.system
    n = ustrip.(us"mol", eq.n)
    V_L = ustrip(uconvert(us"L", volume(eq).liquid))
    dissolved(e) = 1000 * sum(n[k] * Float64(get(atoms_charge(cs.species[k]), e, 0)) for k in cs.idx_aqueous) / V_L
    el = alkali == "NaOH" ? :Na : :K
    e = solid_solution_totals(eq, gel).elements
    return (;
        Si = dissolved(:Si), Ca = dissolved(:Ca), Al = dissolved(:Al), alkali = alkali == "none" ? 0.0 : dissolved(el),
        pH = pH(eq, lh16_model(alkali)), Ca_Si = e[:Ca] / e[:Si], alkali_Si = get(e, el, 0.0) / e[:Si],
    )
end

"""
    lh16_compute(gel) -> Dict

Every batch of Appendices B and C computed with `gel`, keyed by its target
Al/Si, Ca/Si, hydroxide and concentration (the times are separate batches of
one composition, which an equilibrium does not tell apart). Each series of one
composition is walked down from its most concentrated solution, each solve
starting from the last ([`equilibrate_path`](@ref)).
"""
function lh16_compute(gel)
    rows = lh16_table("pH")
    batches = unique(zip(rows.Al_Si_target, rows.Ca_Si_target, rows.alkali, ustrip.(us"mol/L", rows.alkali_concentration)))
    out = Dict{NTuple{4, Any}, Any}()
    for alkali in ("none", "KOH", "NaOH")
        cs = lh16_system(gel, alkali)
        A = Float64.(cs.SM.A)
        model = lh16_model(alkali)
        for (al, ca) in unique((b[1], b[2]) for b in batches if b[3] == alkali)
            series = sort([b for b in batches if b[1] == al && b[2] == ca && b[3] == alkali]; by = b -> -b[4])
            states = [lh16_state(cs, b[1], b[2], b[3], b[4]) for b in series]
            eqs, certs = equilibrate_path(first(states), [A * ustrip.(us"mol", st.n) for st in states]; model)
            for (b, eq, cert) in zip(series, eqs, certs)
                out[b] = (; lh16_observables(eq, gel, alkali)..., certified = cert.optimal, state = eq)
            end
        end
    end
    return out
end

# The batch of a row of the tables: target Al/Si and Ca/Si, hydroxide, and its
# concentration in mol/L.
_lh16_batch(t, i) = (t.Al_Si_target[i], t.Ca_Si_target[i], t.alkali[i], ustrip(us"mol/L", t.alkali_concentration[i]))

"""
    lh16_measured(batch, column) -> Vector{Float64}

The values of `column` of Appendix C (`:alkali_Si`, `:Ca_Si`, …) measured on the
samples of `batch`, one per time, those not measured left out.
"""
function lh16_measured(batch, column)
    t = lh16_table("solid_composition")
    v = getproperty(t, column)
    return Float64[v[i] for i in eachindex(v) if _lh16_batch(t, i) == batch && !ismissing(v[i])]
end

"""
    lh16_measured_solution(batch, element) -> Vector{Float64}

The concentrations of `element` measured in the solutions of `batch`, one per
time, in mmol/L, those below the detection limit left out; `"pH"` for the pH.
"""
function lh16_measured_solution(batch, element)
    if element == "pH"
        t = lh16_table("pH")
        return Float64[t.pH[i] for i in eachindex(t.pH) if _lh16_batch(t, i) == batch]
    end
    t = lh16_table("pore_solution")
    return Float64[
        ustrip(us"mmol/L", t.concentration[i]) for i in eachindex(t.element)
            if _lh16_batch(t, i) == batch && t.element[i] == element && t.qualifier[i] == "measured"
    ]
end

_lh16_mean(x) = sum(x) / length(x)

"""
    lh16_uptake_ratios(results, Ca_Si) -> Vector{Float64}

The alkali over silicon of the gel computed in `results` (one entry of
[`lh16_compute`](@ref)), over the mean of the measured ones, in the batches at
the target `Ca_Si` with an alkali whose every measurement is positive.
"""
function lh16_uptake_ratios(results, Ca_Si)
    r = Float64[]
    for b in sort(collect(keys(results)); by = b -> (b[3], b[1], b[4]))
        b[2] == Ca_Si && b[3] != "none" || continue
        m = lh16_measured(b, :alkali_Si)
        !isempty(m) && all(>(0), m) && push!(r, results[b].alkali_Si / _lh16_mean(m))
    end
    return r
end

"""
    lh16_solution_ratios(results; low) -> (; alkali, Si, Ca, ΔpH)

The dissolved alkali, silicon and calcium computed in `results`, over the mean
of the measured ones, and the absolute difference of pH, in the batches below a
Ca/Si of 1.1 (`low = true`) or above.
"""
function lh16_solution_ratios(results; low)
    out = (alkali = Float64[], Si = Float64[], Ca = Float64[], ΔpH = Float64[])
    for b in sort(collect(keys(results)); by = b -> (b[3], b[2], b[1], b[4]))
        (b[2] < 1.1) == low || continue
        v = results[b]
        if b[3] != "none"
            m = lh16_measured_solution(b, b[3] == "NaOH" ? "Na" : "K")
            isempty(m) || push!(out.alkali, v.alkali / _lh16_mean(m))
        end
        for (e, q) in (("Si", :Si), ("Ca", :Ca))
            m = lh16_measured_solution(b, e)
            isempty(m) || push!(getproperty(out, q), getproperty(v, q) / _lh16_mean(m))
        end
        m = lh16_measured_solution(b, "pH")
        isempty(m) || push!(out.ΔpH, abs(v.pH - _lh16_mean(m)))
    end
    return out
end
