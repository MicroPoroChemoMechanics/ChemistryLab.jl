# =============================================================================
#  snellings2013_glass.jl — the dissolution of the glasses of supplementary
#  cementitious materials at pH 13, and how dissolved calcium and aluminum
#  slow it down
#
#  Snellings (2013) dissolved six synthetic calcium aluminosilicate glasses,
#  from the composition of a slag to silica, in NaOH solutions at pH 13 and
#  20 °C, alone and with Al, Ca or Si added, and gave the initial rates in its
#  Table II. The paper relates the rate to the composition of the glass
#  (`snellings2013_glass`) but gives no law for the effect of the solution. This
#  file computes the activities of Ca²⁺ and of AlO₂⁻ (aluminate) in each of the
#  solutions of Table II and fits to them one inhibitor of each, the form
#  `(1 + K a)⁻¹` of `RateModelInhibitor`. The page
#  docs/src/examples/glass_dissolution.md shows the fit against the 51 rates.
#
#  ASSUMED: the concentration of the NaOH is not given; each solution is
#  computed at the pH of 13 the paper sets, with sodium hydroxide as the reagent
#  that holds it, at 20 °C, a liter of solution taken as a kilogram of water.
#  Al is added as Al(NO3)3, Ca as Ca(OH)2 and Si as SiO2 (Section II(3)); no
#  solid forms, as none did (Section IV(2)).
# =============================================================================

using ChemistryLab
using DynamicQuantities
using OptimaSolver
using ForwardDiff

const SN13 = "Snellings2013"
const SN13_DB = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-thermofun.json"); verbose = false))
const SN13_MODEL = cemdata18_activity_model(:NaOH)

"""
    sn13_system() -> ChemicalSystem

The aqueous species of Cemdata18 holding Ca, Al, Si, Na, N, H and O, nitrate
the only nitrogen species (no redox of the added nitrate), no solid.
"""
function sn13_system()
    sp = speciation(
        collect(values(SN13_DB)), [:Ca, :Al, :Si, :Na, :N, :H, :O, :Zz];
        aggregate_state = [AS_AQUEOUS], exclude_species = split("H2@ O2@ NH4+ NH3@ N2@"),
    )
    return ChemicalSystem(sp, CEMDATA_PRIMARIES)
end

"""
    sn13_solution(cs, al, ca, si; T = 293.15) -> NamedTuple

The solution of NaOH at pH 13 with `al` mol/L of Al(NO3)3, `ca` of Ca(OH)2 and
`si` of SiO2: the activities of Ca²⁺ and AlO₂⁻, the ionic strength, the NaOH it
took, and whether the equilibrium certified.
"""
function sn13_solution(cs, al, ca, si; T = 293.15)
    st = ChemicalState(cs; T = T * u"K")
    set_quantity!(st, "H2O@", 1.0u"kg")
    set_quantity!(st, "NaOH@", 0.1u"mol")       # a start; the pH constraint adjusts it
    al > 0 && (set_quantity!(st, "Al+3", al * u"mol"); set_quantity!(st, "NO3-", 3al * u"mol"))
    ca > 0 && (set_quantity!(st, "Ca+2", ca * u"mol"); set_quantity!(st, "OH-", 2ca * u"mol"))
    si > 0 && set_quantity!(st, "SiO2@", si * u"mol")
    b = budget(st)
    q = Ref{Vector{Float64}}()
    eq, cert = equilibrate_certified(st; model = SN13_MODEL, b, constraint = FixedpH(13.0; titrant = "NaOH@"), parameters = q)
    lna = log_activities(eq, SN13_MODEL)
    return (;
        a_Ca = exp(lna["Ca+2"]), a_Al = exp(lna["AlO2-"]), I = ionic_strength(eq), NaOH = 0.1 + q[][1],
        certified = cert.optimal,
    )
end

"""
    sn13_rows(; cs = sn13_system()) -> Vector{NamedTuple}

The 51 experiments of Table II with the activities of their solutions: glass,
added Al, Ca and Si (mol/L), the measured log10 of the rate (mol/m²/s), and
the log10 of the rate of the same glass in NaOH alone.
"""
function sn13_rows(; cs = sn13_system())
    t = literature_table(SN13, "initial_rates")
    mol(x) = ustrip(us"mol/L", x)
    base = Dict(
        t.glass[k] => ustrip(t.log_rate[k]) for k in eachindex(t.glass)
            if all(iszero ∘ mol, (t.Al_initial[k], t.Ca_initial[k], t.Si_initial[k]))
    )
    out = NamedTuple[]
    for k in eachindex(t.glass)
        al, ca, si = mol(t.Al_initial[k]), mol(t.Ca_initial[k]), mol(t.Si_initial[k])
        s = sn13_solution(cs, al, ca, si)
        push!(out, (; glass = t.glass[k], al, ca, si, log_rate = ustrip(t.log_rate[k]), log_base = base[t.glass[k]], s...))
    end
    return out
end

# The glasses whose calcium only balances their aluminum, along CaO/Al2O3 = 1,
# and silica: the "tectosilicate" glasses of the paper, against the percalcic
# G1 and G2 (Section IV(1)).
const SN13_TECTOSILICATE = ("G3", "G4", "G5", "G6")

"""
    sn13_fit_K(pairs) -> NamedTuple

The `K` of one inhibitor of order one fitted, in log10 of the rate, to `pairs`
of (activity, measured change of log10 of the rate against NaOH alone): Newton
on ln K from the best of a grid, its derivatives by automatic differentiation.
"""
function sn13_fit_K(pairs)
    sse(u) = sum((d + log10(1 + exp(u) * a))^2 for (a, d) in pairs)
    u = argmin(sse, range(log(1.0), log(1.0e7); length = 400))
    for _ in 1:50
        g = ForwardDiff.derivative(sse, u)
        h = ForwardDiff.derivative(v -> ForwardDiff.derivative(sse, v), u)
        h > 0 || break
        step = g / h
        u -= step
        abs(step) < 1.0e-12 && break
    end
    res = [d + log10(1 + exp(u) * a) for (a, d) in pairs]
    return (; K = exp(u), rms = sqrt(sum(abs2, res) / length(res)), worst = maximum(abs, res), n = length(res))
end

"""
    sn13_fit(rows) -> NamedTuple

The inhibitors fitted on Table II: calcium on every glass; aluminum on the
tectosilicate glasses, and, to report what the paper calls within the error,
on the percalcic ones apart.
"""
function sn13_fit(rows)
    d(r) = r.log_rate - r.log_base
    ca = sn13_fit_K([(r.a_Ca, d(r)) for r in rows if r.ca > 0])
    al = sn13_fit_K([(r.a_Al, d(r)) for r in rows if r.al > 0 && r.glass in SN13_TECTOSILICATE])
    al_percalcic = sn13_fit_K([(r.a_Al, d(r)) for r in rows if r.al > 0 && !(r.glass in SN13_TECTOSILICATE)])
    si = [d(r) for r in rows if r.si > 0]
    return (; ca, al, al_percalcic, si_rms = sqrt(sum(abs2, si) / length(si)), si_worst = maximum(abs, si))
end

"""
    sn13_inhibitors(fit, glass) -> Vector{RateModelInhibitor}

The inhibitors of `fit` for a glass: calcium always; aluminum for a
tectosilicate glass only (`glass` names one of the paper's, or is `:tectosilicate`
or `:percalcic`).
"""
function sn13_inhibitors(fit, glass)
    tecto = glass === :tectosilicate || glass in SN13_TECTOSILICATE
    inh = [RateModelInhibitor("Ca+2", fit.ca.K)]
    tecto && push!(inh, RateModelInhibitor("AlO2-", fit.al.K))
    return inh
end
