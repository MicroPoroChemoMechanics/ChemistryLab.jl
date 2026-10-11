# =============================================================================
#  georget2022_hemicarbonate.jl — hemicarbonate and monocarbonate as the
#  carbonate of a paste grows
#
#  Georget et al. (2022) hydrated 10 g of C3A with portlandite and calcite in
#  the proportions of the AFm phases, for 28 days at 20 °C, replacing
#  portlandite by calcite from one sample to the next: the carbonate taken up
#  by a paste at fixed Ca/Al, a carbonation in the AFm. They measured the
#  phases and the solution, and published with their data set (Georget et al.
#  2021) the same series computed with GEMS and Cemdata18 in 101 steps. This
#  file computes those 101 steps with this package; the page
#  docs/src/examples/hemicarbonate.md compares both with the measurement.
#
#  ASSUMED: the water of the GEMS run is not published; each step here has 2.7
#  times the mass of its solids, the water/solid ratio of the experiments. The
#  1e-7 g of calcite the data set gives its last step is taken as none.
# =============================================================================

using ChemistryLab
using DynamicQuantities
using OptimaSolver

const GE22_DB = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-thermofun.json"); verbose = false))
const GE22_SOLIDS = split("C3A Lim Portlandite Cal hemicarbonate monocarbonate C3AH6 C4AH13 C2AH7.5 CAH10 AlOHmic")

"""
    ge22_system() -> ChemicalSystem

C3A, lime, portlandite, calcite, the carbonate and hydroxide AFm phases,
katoite, the other calcium aluminate hydrates and gibbsite, in water.
"""
function ge22_system()
    sp = speciation(
        collect(values(GE22_DB)), GE22_SOLIDS;
        aggregate_state = [AS_AQUEOUS], exclude_species = split("H2@ O2@ CH4@"),
    )
    return ChemicalSystem(sp, CEMDATA_PRIMARIES)
end

const GE22_MODEL = cemdata18_activity_model(:KOH)

"""
    ge22_steps() -> NamedTuple

The 101 steps of the data set: the portlandite and calcite added to 10 g of
C3A (g), the ζ of each, and the phases GEMS computed (g).
"""
ge22_steps() = literature_table("Georget2021data", "carbonate_series_cemdata18")

"""
    ge22_series(; T = 293.15) -> Vector{NamedTuple}

Each step of the series at equilibrium, from the previous one: the C3A of the
samples (97 % of the 10 g, the rest free lime), the portlandite and calcite of
the step, and 2.7 times their mass of water. The masses of katoite,
hemicarbonate, monocarbonate, calcite and portlandite (g), the pH, and whether
each step certified.
"""
function ge22_series(; T = 293.15)
    cs = ge22_system()
    t = ge22_steps()
    M(s) = ustrip(us"g/mol", GE22_DB[s][:M])
    out = NamedTuple[]
    prev = nothing
    for k in eachindex(t.step)
        ch, cc = ustrip(us"g", t.CaOH2_added[k]), ustrip(us"g", t.CaCO3_added[k])
        # The 1e-7 g of calcite of the last step, a placeholder for none, would
        # leave a carbon budget of 1e-9 mol, below what a balance resolves.
        cc < 1.0e-6 && (cc = 0.0)
        st = ChemicalState(cs; T = T * u"K")
        set_quantity!(st, "C3A", 9.7u"g")
        set_quantity!(st, "Lim", 0.3u"g")
        set_quantity!(st, "Portlandite", ch * u"g")
        set_quantity!(st, "Cal", cc * u"g")
        set_quantity!(st, "H2O@", 2.7 * (10 + ch + cc) * u"g")
        b = budget(st)
        eq, cert = equilibrate_certified(prev === nothing ? st : prev; model = GE22_MODEL, b)
        n = ustrip.(us"mol", eq.n)
        mass(s) = n[findfirst(x -> symbol(x) == s, cs.species)] * M(s)
        push!(
            out, (;
                zeta = t.zeta_CO3[k], certified = cert.optimal, katoite = mass("C3AH6"),
                hemicarbonate = mass("hemicarbonate"), monocarbonate = mass("monocarbonate"),
                calcite = mass("Cal"), portlandite = mass("Portlandite"), pH = pH(eq, GE22_MODEL),
            ),
        )
        prev = eq
    end
    return out
end

"""
    ge22_breakpoints(rows) -> NamedTuple

The ζ at which katoite is gone and at which hemicarbonate is gone, along a
series of increasing ζ (the data set runs it decreasing).
"""
function ge22_breakpoints(rows)
    r = sort(rows; by = x -> x.zeta)
    i = findfirst(x -> x.katoite < 1.0e-6, r)
    j = findfirst(x -> x.hemicarbonate < 1.0e-6 && x.monocarbonate > 1.0e-6, r)
    return (; katoite_gone = r[i].zeta, hemicarbonate_gone = r[j].zeta)
end

"""
    ge22_phases(rows, zeta) -> Vector{String}

The phases present at the step of the series nearest `zeta`, in the
abbreviations of Table 7 of Georget et al. (2022): Hc, Mc, katoite, CH, CC.
"""
function ge22_phases(rows, zeta)
    r = rows[argmin([abs(x.zeta - zeta) for x in rows])]
    out = String[]
    for (f, lab) in ((:hemicarbonate, "Hc"), (:monocarbonate, "Mc"), (:katoite, "katoite"), (:portlandite, "CH"), (:calcite, "CC"))
        getproperty(r, f) > 1.0e-6 && push!(out, lab)
    end
    return out
end
