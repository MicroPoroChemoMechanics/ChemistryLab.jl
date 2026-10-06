# =============================================================================
#  berner1992_leaching.jl — the leaching of a cement paste, in zero dimensions
#
#  Water that renews itself around a paste takes its calcium: portlandite
#  dissolves first, then the AFm and the ettringite, and the C-S-H loses
#  calcium incongruently, its Ca/Si falling with the calcium of the solution.
#  Two sources measure it:
#
#    - Berner (1992), Appendix A, Tables 8 to 13: six sets of solubility data
#      of synthetic C-S-H, the Ca/Si of the solid against the calcium and the
#      silicon of the solution in equilibrium with it, from 17 to 30 °C, as he
#      compiled them;
#    - Adenot and Buil (1992): the zones of a paste leached by deionized water,
#      the portlandite dissolved first, then the monosulfate, then the
#      ettringite.
#
#  This file computes the first with the C-S-H of Cemdata18, CSHQ, and the
#  second with a paste leached by successive renewals of its solution. The page
#  docs/src/examples/leaching.md compares them. Nothing is fitted.
# =============================================================================

using ChemistryLab
using DynamicQuantities
using OptimaSolver

# The paste: the CEM I 42.5 N of Lothenbach and Winnefeld (2006), whose
# materials, recipe and phases are those of its own validation page.
isdefined(@__MODULE__, :lw06_recipe) || include(joinpath(@__DIR__, "lothenbach_winnefeld_2006.jl"))

const BE92_DB = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-thermofun.json"); verbose = false))
const BE92_CSHQ = ["CSHQ-JenD", "CSHQ-JenH", "CSHQ-TobD", "CSHQ-TobH"]
const BE92_PURE = split("Lim Amor-Sl Portlandite")

"""
    be92_system() -> ChemicalSystem

Lime, amorphous silica, portlandite and the C-S-H of Cemdata18 (CSHQ, its four
alkali-free end members) in water, with the aqueous species of Ca and Si.
"""
function be92_system()
    sp = speciation(
        collect(values(BE92_DB)), vcat(BE92_PURE, BE92_CSHQ);
        aggregate_state = [AS_AQUEOUS], exclude_species = split("H2@ O2@ CH4@"),
    )
    ss = [SolidSolutionPhase("CSHQ", [BE92_DB[m] for m in BE92_CSHQ])]
    return ChemicalSystem(sp, CEMDATA_PRIMARIES; solid_solutions = ss)
end

"""
    be92_state(cs, Ca_Si; silica = 0.05, water = 1000.0, T = 298.15) -> ChemicalState

`silica` mol of amorphous silica and `Ca_Si` times as much lime in `water` g of
water at `T` (K): a synthesis of C-S-H at a bulk Ca/Si, as the compiled
experiments made it, from lime and silica.
"""
function be92_state(cs, Ca_Si; silica = 0.05, water = 1000.0, T = 298.15)
    st = ChemicalState(cs; T = T * u"K")
    set_quantity!(st, "H2O@", water * u"g")
    set_quantity!(st, "Amor-Sl", silica * u"mol")
    set_quantity!(st, "Lim", Ca_Si * silica * u"mol")
    return st
end

const BE92_MODEL = cemdata18_activity_model(:KOH)

"""
    be92_sweep(ratios; T = 298.15) -> Vector{NamedTuple}

The C-S-H made at each bulk Ca/Si of `ratios`, in increasing order, each
equilibrium started from the previous one: the Ca/Si of the CSHQ formed, the
calcium and silicon of the solution (mmol/L), the pH, whether portlandite or
amorphous silica is present beside the gel, and whether the equilibrium
certified.
"""
function be92_sweep(ratios; T = 298.15)
    cs = be92_system()
    A = Float64.(cs.SM.A)
    states = [be92_state(cs, r; T) for r in ratios]
    eqs, certs = equilibrate_path(first(states), [A * ustrip.(us"mol", st.n) for st in states]; model = BE92_MODEL)
    return map(zip(ratios, eqs, certs)) do (r, eq, cert)
        n = ustrip.(us"mol", eq.n)
        V_L = ustrip(uconvert(us"L", volume(eq).liquid))
        dissolved(e) = 1000 * sum(n[k] * Float64(get(atoms_charge(cs.species[k]), e, 0)) for k in cs.idx_aqueous) / V_L
        gel = solid_solution_totals(eq, "CSHQ").elements
        amount(s) = n[findfirst(x -> symbol(x) == s, cs.species)]
        (;
            bulk = r, Ca_Si = gel[:Ca] / gel[:Si], Ca = dissolved(:Ca), Si = dissolved(:Si),
            pH = pH(eq, BE92_MODEL), portlandite = amount("Portlandite") > 1.0e-9,
            silica = amount("Amor-Sl") > 1.0e-9, certified = cert.optimal,
        )
    end
end

"""
    be92_measured() -> NamedTuple

The rows of Berner's compilation (Appendix A, Tables 8 to 13): their source,
temperature (K), the Ca/Si of the solid, the calcium and silicon of the solution
(mmol/L, `missing` where not reported) and the pH.
"""
function be92_measured()
    t = literature_table("Berner1992", "csh_solubility")
    v(x) = ismissing(x) || x === nothing ? missing : Float64(ustrip(x))
    return (;
        source = t.original_source, T = v.(t.T), Ca_Si = v.(t.CaSi_solid),
        Ca = [ismissing(x) || x === nothing ? missing : ustrip(us"mmol/L", x) for x in t.Ca_total],
        Si = [ismissing(x) || x === nothing ? missing : ustrip(us"mmol/L", x) for x in t.Si_total],
        pH = v.(t.pH),
    )
end

"""
    be92_interpolate(sweep, field, Ca_Si) -> Float64

`field` of the computed C-S-H at the Ca/Si `Ca_Si` of a measured solid, by
linear interpolation in the Ca/Si of the gel, on the logarithm for the
concentrations; `NaN` outside the range the gel spans.
"""
function be92_interpolate(sweep, field, Ca_Si)
    pts = sort([(p.Ca_Si, getproperty(p, field)) for p in sweep if !p.portlandite && !p.silica]; by = first)
    x = first.(pts)
    (Ca_Si < x[1] || Ca_Si > x[end]) && return NaN
    k = clamp(searchsortedlast(x, Ca_Si), 1, length(x) - 1)
    w = (Ca_Si - x[k]) / (x[k + 1] - x[k])
    y0, y1 = pts[k][2], pts[k + 1][2]
    field === :pH && return y0 + w * (y1 - y0)
    return exp(log(y0) + w * (log(y1) - log(y0)))
end

"""
    be92_leaching(; t = 28, stages = ((40, 500.0), (40, 5000.0))) -> Vector{NamedTuple}

The CEM I paste of Lothenbach and Winnefeld (2006) at `t` days, its solution
then removed and replaced by pure water again and again ([`leach`](@ref)): for
each `(steps, renewal)` of `stages`, `steps` renewals of `renewal` g of water
per 100 g of cement. At each step: the water that has passed (g per 100 g of
cement), the calcium, silicon, sulfur and aluminum of the solution (mmol/L),
the pH, the mass of portlandite, of the AFm phases, of ettringite and of the
C-S-H (g per 100 g of cement), the Ca/Si of the C-S-H, and whether the
equilibrium certified.
"""
function be92_leaching(; t = 28, stages = ((40, 500.0), (40, 5000.0)))
    cs = lw06_system()
    model = cemdata18_activity_model(:KOH)
    rs, cert = equilibrate_certified(lw06_recipe(), cs; t, model)
    cert.optimal || error("be92_leaching: the paste at $t days did not certify.")
    rows = NamedTuple[]
    water = 0.0
    push!(rows, _be92_row(rs, water))
    for (steps, w) in stages
        p = leach(rs, steps; renewal = w * rs.recipe.binder_mass / 100)
        for s in p.states
            water += w
            push!(rows, _be92_row(s, water))
        end
        rs = p.states[end]
    end
    return rows
end

function _be92_row(rs, water)
    m = phase_masses(rs)
    c = pore_solution_mmol(rs.state)
    gel = solid_solution_totals(rs.state, "CSHQ")
    afm = sum(v for (k, v) in m if startswith(k, "AFm") || k in ("monocarbonate", "hemicarbonate"); init = 0.0)
    return (;
        water, Ca = get(c, "Ca", 0.0), Si = get(c, "Si", 0.0), S = get(c, "S", 0.0), Al = get(c, "Al", 0.0),
        pH = pH(rs.state, rs.model), portlandite = get(m, "Portlandite", 0.0), afm,
        ettringite = get(m, "ettringite", 0.0), csh = 1000 * ustrip(gel.mass) * 100 / rs.recipe.binder_mass,
        phases = sort([k for (k, v) in m if v > 1.0e-6]),
        Ca_Si = gel.elements[:Ca] / gel.elements[:Si], certified = rs.certificate.optimal,
    )
end

"""
    be92_events(rows) -> Vector{NamedTuple}

Where each phase of the paste is gone along a leaching sequence of
[`be92_leaching`](@ref): for portlandite, the AFm phases and ettringite, the
first step at which none is left, the water that has passed by then, and the
calcium of the solution, the pH and the Ca/Si of the C-S-H at that step.
"""
function be92_events(rows)
    out = NamedTuple[]
    for (name, field) in (("portlandite", :portlandite), ("AFm", :afm), ("ettringite", :ettringite))
        k = findfirst(r -> getproperty(r, field) < 1.0e-6, rows)
        k === nothing && continue
        r = rows[k]
        push!(out, (; phase = name, water = r.water, Ca = r.Ca, pH = r.pH, Ca_Si = r.Ca_Si))
    end
    return out
end
