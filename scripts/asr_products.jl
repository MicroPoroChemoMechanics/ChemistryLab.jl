# =============================================================================
#  asr_products.jl — the syntheses of the products of the alkali-silica
#  reaction at 80 °C
#
#  Shi and Lothenbach (2019) synthesized at 80 °C, from silica fume, lime and
#  KOH or NaOH at an alkali/Si of 0.5 and a Ca/Si from 0 to 0.5, the
#  crystalline products of the alkali-silica reaction: K- and Na-shlykovite and
#  the potassium calcium silicate hydrate they call ASR-P1. They identified the
#  solids by XRD after 90 days, analyzed the solutions (their Tables 1 and 6), and
#  derived from them the solubility products of the three solids at 80 °C (their
#  Table 5). Jin et al. (2023) then gave the two shlykovites standard properties
#  at 25 °C, with a pH treatment of their own, which `cemdata18-asr.json` carries.
#  This file computes the 17 syntheses with Cemdata18, its C-S-H (CSHQ with the
#  alkali end member), amorphous silica and portlandite, and the products in
#  one of two sets:
#
#    :jin  the shlykovites of `cemdata18-asr.json`, at any temperature, and
#          ASR-P1 from Table 5 of Shi and Lothenbach, at 80 °C only;
#    :shi  all three from Table 5 of Shi and Lothenbach, at 80 °C only.
#
#  The page docs/src/examples/asr_products.md compares phases and solutions with
#  the measurements.
#
#  ASSUMED: the nominal alkali/Si of 0.5, the hydroxide pellets taken as pure
#  (the paper gives a lower bound on their purity); a solid of Table 5 held at
#  the Gibbs energy its solubility product gives at 80 °C, valid at that
#  temperature and no other.
# =============================================================================

using ChemistryLab
using DynamicQuantities
using OptimaSolver

const SL19 = "ShiLothenbach2019"
const SL19_T = 353.15
const SL19_DB = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-asr.json"); verbose = false))
const SL19_CSHQ = ["CSHQ-JenD", "CSHQ-JenH", "CSHQ-TobD", "CSHQ-TobH"]
const SL19_MODEL = Dict("K" => cemdata18_activity_model(:KOH), "Na" => cemdata18_activity_model(:NaOH))

"""
    sl19_reaction(phase) -> Vector{Pair{String, Float64}}

The dissolution products of `phase` as Table 5 of Shi and Lothenbach (2019)
writes its solubility product, with their coefficients.
"""
function sl19_reaction(phase)
    t = literature_table(SL19, "asr_product_solubility")
    i = only(findall(==(phase), t.phase))
    return map(split(strip(last(split(t.reaction[i], "="))), " + ")) do term
        parts = split(strip(term))
        length(parts) == 2 ? String(parts[2]) => parse(Float64, parts[1]) : String(parts[1]) => 1.0
    end
end

"""
    sl19_table5(phase) -> Species

A product as Table 5 of Shi and Lothenbach (2019) defines it: the Gibbs energy
its solubility product and the Cemdata18 energies of its dissolution products
give it at 80 °C, held constant, so that the record holds at that temperature
only.
"""
function sl19_table5(phase)
    t = literature_table(SL19, "asr_product_solubility")
    i = only(findall(==(phase), t.phase))
    G = sum(ν * SL19_DB[s][:ΔₐG⁰](T = SL19_T) for (s, ν) in sl19_reaction(phase))
    G += R_GAS * SL19_T * log(10) * Float64(ustrip(t.log_Ks0[i]))
    sp = Species(String(t.formula[i]); name = phase, symbol = phase, aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)
    sp[:ΔₐG⁰] = SymbolicFunc(G * u"J/mol")
    return sp
end

"""
    sl19_products(set) -> Vector{Species}

K-shlykovite, Na-shlykovite and ASR-P1 in the set `:jin` or `:shi` (see the
header of this file).
"""
function sl19_products(set)
    set in (:jin, :shi) || throw(ArgumentError("sl19_products: the set is :jin or :shi; got $set."))
    shly = set === :jin ? [SL19_DB["K-shlykovite"], SL19_DB["Na-shlykovite"]] : sl19_table5.(["K-shlykovite", "Na-shlykovite"])
    return vcat(shly, [sl19_table5("ASR-P1")])
end

"""
    sl19_system(alkali, set) -> ChemicalSystem

The aqueous species of Cemdata18 for Ca, Si and the alkali (`"K"` or `"Na"`),
amorphous silica, lime and portlandite, CSHQ with the end member of that alkali,
and the products of that alkali in `set`.
"""
function sl19_system(alkali, set)
    member = alkali == "K" ? "KSiOH" : "NaSiOH"
    gel = [SL19_DB[m] for m in vcat(SL19_CSHQ, [member])]
    sp = speciation(
        collect(values(SL19_DB)), vcat(["Amor-Sl", "Lim", "Portlandite"], SL19_CSHQ, [member]);
        aggregate_state = [AS_AQUEOUS], exclude_species = split("H2@ O2@ CH4@"),
    )
    other = alkali == "K" ? :Na : :K
    products = [p for p in sl19_products(set) if !haskey(atoms(p), other)]
    return ChemicalSystem(vcat(sp, products), CEMDATA_PRIMARIES; solid_solutions = [SolidSolutionPhase("CSHQ", gel)])
end

"""
    sl19_state(cs, i) -> ChemicalState

The mix of row `i` of Table 1 at 80 °C: its silica fume and lime as weighed,
its water, and the hydroxide at the nominal alkali/Si.
"""
function sl19_state(cs, i)
    mixes = literature_table(SL19, "mixes")
    alkali = startswith(mixes.series[i], "K") ? "K" : "Na"
    n_Si = ustrip(us"g", mixes.SiO2_am[i]) / ustrip(us"g/mol", SL19_DB["Amor-Sl"][:M])
    st = ChemicalState(cs; T = SL19_T * u"K")
    set_quantity!(st, "H2O@", mixes.H2O[i])
    set_quantity!(st, "Amor-Sl", n_Si * u"mol")
    iszero(ustrip(mixes.CaO[i])) || set_quantity!(st, "Lim", mixes.CaO[i])
    n_alkali = Float64(ustrip(mixes.alkali_Si[i])) * n_Si
    set_quantity!(st, alkali * "+", n_alkali * u"mol")
    set_quantity!(st, "OH-", n_alkali * u"mol")
    return st
end

# The solids as the XRD names them: amorphous silica, the products, C-S-H for
# any member of CSHQ, portlandite.
function _sl19_solids(symbols)
    names = String[]
    "Amor-Sl" in symbols && push!(names, "SiO2 (am)")
    for p in ("K-shlykovite", "Na-shlykovite", "ASR-P1")
        p in symbols && push!(names, p)
    end
    any(s -> startswith(s, "CSHQ") || s in ("KSiOH", "NaSiOH"), symbols) && push!(names, "C-S-H")
    "Portlandite" in symbols && push!(names, "portlandite")
    return join(names, ", ")
end

"""
    sl19_syntheses(set) -> Vector{NamedTuple}

Each mix of Table 1 at equilibrium at 80 °C with the products of `set`: its
sample and series, whether it certified, the solids present, the Si, alkali and
Ca of the solution (mmol/L) and the measured ones with the phases of the XRD
(Table 6, `NaN` and an empty string where the paper gives none), the ionic
strength (mol/kg) and the equilibrium.
"""
function sl19_syntheses(set)
    mixes = literature_table(SL19, "mixes")
    sol = literature_table(SL19, "solutions")
    systems = Dict{String, Any}()
    return map(eachindex(mixes.sample)) do i
        alkali = startswith(mixes.series[i], "K") ? "K" : "Na"
        cs = get!(() -> sl19_system(alkali, set), systems, alkali)
        eq, cert = equilibrate_certified(sl19_state(cs, i); model = SL19_MODEL[alkali])
        n = ustrip.(us"mol", eq.n)
        V_L = ustrip(uconvert(us"L", volume(eq).liquid))
        dissolved(e) = 1000 * sum(n[k] * Float64(get(atoms(cs.species[k]), e, 0)) for k in cs.idx_aqueous) / V_L
        solids = [symbol(cs.species[k]) for k in eachindex(n) if !(k in cs.idx_aqueous) && n[k] > 1.0e-7]
        j = findfirst(k -> sol.sample[k] == mixes.sample[i] && sol.series[k] == mixes.series[i], eachindex(sol.sample))
        meas(c) = j === nothing ? NaN : (x = getproperty(sol, c)[j]; x === nothing || ismissing(x) ? NaN : Float64(ustrip(us"mmol/L", x)))
        (;
            sample = mixes.sample[i], series = mixes.series[i], alkali, certified = cert.optimal,
            solids = _sl19_solids(solids), Si = dissolved(:Si), A = dissolved(Symbol(alkali)), Ca = dissolved(:Ca),
            Si_meas = meas(:Si), A_meas = meas(Symbol(alkali)), Ca_meas = meas(:Ca),
            xrd = j === nothing ? "" : String(sol.phases_xrd[j]), I = ionic_strength(eq), state = eq,
        )
    end
end

"""
    sl19_ionic_strength_equation(eq, alkali) -> Function

`F(s) = ln(½ Σ zᵢ² mᵢ(eˢ)) − s` at the potentials of the equilibrium `eq`: each
ion keeps the activity it has there, and its molality at an ionic strength `eˢ`
is that activity over the activity coefficient the model gives it at `eˢ`. Its
zeros are the ionic strengths the solve can choose between; `eq` sits on one.
Uses the activity coefficients the solver itself inverts (`_aqueous_form`, an
internal function of ChemistryLab).
"""
function sl19_ionic_strength_equation(eq, alkali)
    cs, model = eq.system, SL19_MODEL[alkali]
    ions = [i for i in cs.idx_solutes if !iszero(charge(cs.species[i]))]
    form = ChemistryLab._aqueous_form(model, cs, ions)
    A, B = form.AB((; T = SL19_T, P = 1.0e5))
    lna = log_activities(eq, model)
    lnγ(t, I) = log(10) * form.log10γ(t, form.z[t], I, sqrt(I), A, B)
    return s -> log(sum(form.z[t]^2 / 2 * exp(lna[symbol(cs.species[ions[t]])] - lnγ(t, exp(s))) for t in eachindex(ions))) - s
end
