# =============================================================================
#  delayed_ettringite.jl — ettringite, monosulfate and the pore solution
#  against temperature
#
#  Lothenbach et al. (2007) hydrated a sulfate-resisting CEM I 52.5 N HTS
#  (SRPC) at w/c 0.4 at 5, 20 and 50 °C and measured its pore solution at 28
#  and 150 days and its phases: at 50 °C monosulfate replaces part of the
#  ettringite, and the sulfate of the solution rises sixfold. Heated above that,
#  as a precast element is in steam curing, a paste loses its ettringite, and
#  forms it again once cold: the delayed ettringite formation. This file
#  computes the paste of the companion paper (Lothenbach et al. 2008), at its
#  degree of hydration after 150 days, at each temperature and back. The page
#  docs/src/examples/delayed_ettringite.md compares it with the measurements.
# =============================================================================

using ChemistryLab
using DynamicQuantities
using OptimaSolver

isdefined(@__MODULE__, :l08t_state) || include(joinpath(@__DIR__, "lothenbach2008_temperature.jl"))

const DE_MODEL = cemdata18_activity_model(:KOH)

"""
    de_state(T_C; start = nothing) -> RecipeState

The SRPC paste at `T_C` °C, its clinker at the degree 150 days give it at
20 °C, certified, from `start` when given.
"""
function de_state(T_C; start = nothing, cs = l08t_system())
    rs, cert = equilibrate_certified(l08t_recipe("SRPC", T_C), cs; model = DE_MODEL, start)
    cert.optimal || error("de_state: the paste at $T_C °C did not certify.")
    return rs
end

"""
    de_observables(rs) -> NamedTuple

The dissolved sulfate, aluminum, iron, calcium, sodium and potassium (mol/L)
and the hydroxide of the paste `rs`, and its ettringite and monosulfate
(g per 100 g of cement).
"""
function de_observables(rs)
    st = rs.state
    cs = st.system
    n = ustrip.(us"mol", st.n)
    V_L = ustrip(uconvert(us"L", volume(st).liquid))
    dissolved(e) = sum(n[k] * Float64(get(atoms(cs.species[k]), e, 0)) for k in cs.idx_aqueous) / V_L
    at(s) = (j = findfirst(x -> symbol(x) == s, cs.species); j === nothing ? 0.0 : n[j])
    mass(s) = (j = findfirst(x -> symbol(x) == s, cs.species); j === nothing ? 0.0 : n[j] * ustrip(us"g/mol", cs.species[j][:M]))
    return (;
        SO4 = dissolved(:S), Al = dissolved(:Al), Fe = dissolved(:Fe), Ca = dissolved(:Ca),
        Na = dissolved(:Na), K = dissolved(:K), OH = at("OH-") / V_L,
        ettringite = mass("ettringite"), monosulfate = mass("monosulphate12"),
    )
end

"""
    de_measured(element, T_C; age = 150) -> Float64

The concentration of `element` (`"SO4"`, `"Al"`, …) Lothenbach et al. (2007)
measured in the pore solution of the SRPC at `T_C` °C and `age` days (mol/L),
`NaN` where their figure does not show it.
"""
function de_measured(element, T_C; age = 150)
    t = literature_table("Lothenbach2007", "pore_solution_srpc")
    i = findfirst(k -> t.element[k] == element && t.temperature_C[k] == T_C && ustrip(us"d", t.age_d[k]) == age, eachindex(t.element))
    return i === nothing ? NaN : ustrip(us"mol/L", t.concentration[i])
end
