# =============================================================================
#  lhopital2016_aluminum.jl — aluminum uptake by C-S-H, L'Hôpital et al. (2016a)
#
#  L'Hôpital, Lothenbach, Kulik and Scrivener (2016a) synthesized C-(A-)S-H at a
#  Ca/Si from 0.6 to 1.6 and an Al/Si from 0 to 0.33, 2 g of lime, silica fume
#  and monocalcium aluminate in 90 mL of water, equilibrated at 20 °C for 182
#  days (some for 364 or 546), and reported the C-S-H and the other solids by
#  mass balance (Appendix B) and the solutions (Appendix D). This file computes
#  the same batches at equilibrium with CNASH_ss (Myers et al. 2014), the one
#  model of the gel the package ships that takes aluminum, and with CSHQ, which
#  takes none, the hydrates of Cemdata18 beside either. The page
#  docs/src/tutorials/validation_aluminum_uptake.md compares them.
#
#  Neither model was fitted on these data: CNASH_ss was derived from earlier
#  data on synthetic C-(N-)A-S-H (Myers et al. 2014), CSHQ on aluminum-free
#  C-S-H.
# =============================================================================

isdefined(@__MODULE__, :lh16_system) || include(joinpath(@__DIR__, "lhopital2016_alkali.jl"))

const LH16A = "LHopital2016a"
lh16a_table(name; where...) = literature_table(LH16A, name; where...)

"""
    lh16a_system(gel) -> ChemicalSystem

The system of an alkali-free batch with the gel `gel`, `"CNASH_ss"` or
`"CSHQ"`, on Cemdata18, and the phases of `LH16_PURE`.
"""
function lh16a_system(gel)
    gel == "CSHQ" && return lh16_system("CSHQ", "none")
    gel == "CNASH_ss" || throw(ArgumentError("lh16a_system: the gel is \"CNASH_ss\" or \"CSHQ\"; got \"$gel\""))
    ss = only(p for p in build_solid_solutions(datapath("solid_solutions.toml"), LH16_DB) if name(p) == gel)
    members = [symbol(m) for m in ss.end_members]
    sp = speciation(
        collect(values(LH16_DB)), vcat(LH16_PURE, members);
        aggregate_state = [AS_AQUEOUS], exclude_species = split("H2@ O2@ CH4@"),
    )
    return ChemicalSystem(sp, CEMDATA_PRIMARIES; solid_solutions = [ss])
end

"""
    lh16a_state(cs, Ca_Si, Al_Si) -> ChemicalState

The batch of Appendix A at the target `Ca_Si` and `Al_Si`, its lime, silica
fume and monocalcium aluminate in 90 mL of water, at 20 °C.

ASSUMED: the 90 mL are 90 g of water; the solids are the oxides they are
written as, the silica fume amorphous silica.
"""
function lh16a_state(cs, Ca_Si, Al_Si)
    mix = lh16a_table("mixing_proportions"; Ca_Si_target = Ca_Si, Al_Si_target = Al_Si)
    g(q) = ustrip(us"g", only(q))
    st = ChemicalState(cs; T = literature_value(LH16A, "temperature"))
    set_quantity!(st, "H2O@", ustrip(u"mL", literature_value(LH16A, "solution_volume")) * u"g")
    set_quantity!(st, "Lim", g(mix.CaO) * u"g")
    set_quantity!(st, "Amor-Sl", g(mix.SiO2) * u"g")
    g(mix.CaO_Al2O3) > 0 && set_quantity!(st, "CA", g(mix.CaO_Al2O3) * u"g")
    return st
end

# The other solids of Appendix B, by the species of the system that each counts.
const LH16A_OTHER = (
    katoite = ("C3AH6", "C3AS0.84H4.32"), straetlingite = ("straetlingite",), portlandite = ("Portlandite",),
)

"""
    lh16a_observables(eq, gel) -> NamedTuple

What Appendices B and D report, on the equilibrium `eq`: the Ca/Si and Al/Si of
the gel, the other solids in wt.% of the solids, and the dissolved Si, Ca and Al
in mmol/L and the pH.

ASSUMED: the other solids are weighed against every solid of the equilibrium,
the water of the gel and of the hydrates included, where the authors weighed
them by TGA on the freeze-dried powder.
"""
function lh16a_observables(eq, gel)
    cs = eq.system
    n = ustrip.(us"mol", eq.n)
    V_L = ustrip(uconvert(us"L", volume(eq).liquid))
    dissolved(e) = 1000 * sum(n[k] * Float64(get(atoms_charge(cs.species[k]), e, 0)) for k in cs.idx_aqueous) / V_L
    e = solid_solution_totals(eq, gel).elements
    g(s) = haskey(cs.dict_species, s) ? ustrip(uconvert(us"g", mass(eq, cs.dict_species[s]))) : 0.0
    solids = ustrip(uconvert(us"g", mass(eq).solid))
    other = map(names -> 100 * sum(g, names) / solids, LH16A_OTHER)
    return (;
        Ca_Si = e[:Ca] / e[:Si], Al_Si = get(e, :Al, 0.0) / e[:Si], other...,
        Si = dissolved(:Si), Ca = dissolved(:Ca), Al = dissolved(:Al), pH = pH(eq, lh16_model("none")),
    )
end

"""
    lh16a_compute(gel) -> Dict

Every batch of Appendix A computed with `gel`, keyed by its target Ca/Si and
Al/Si. Each series of one Ca/Si is walked up from the aluminum-free batch, each
solve starting from the last ([`equilibrate_path`](@ref)).
"""
function lh16a_compute(gel)
    cs = lh16a_system(gel)
    A = Float64.(cs.SM.A)
    model = lh16_model("none")
    mix = lh16a_table("mixing_proportions")
    out = Dict{Tuple{Float64, Float64}, Any}()
    for ca in unique(mix.Ca_Si_target)
        series = sort([al for (c, al) in zip(mix.Ca_Si_target, mix.Al_Si_target) if c == ca])
        states = [lh16a_state(cs, ca, al) for al in series]
        eqs, certs = equilibrate_path(first(states), [A * ustrip.(us"mol", st.n) for st in states]; model)
        for (al, eq, cert) in zip(series, eqs, certs)
            out[(ca, al)] = (; lh16a_observables(eq, gel)..., certified = cert.optimal, state = eq)
        end
    end
    return out
end

"""
    lh16a_measured(Ca_Si, Al_Si, column) -> Union{Float64, Missing}

The value of `column` of Appendix B for the batch at the target `Ca_Si` and
`Al_Si`: on its sample hydrated for 182 days, or for 364 where the batch has no
other; `missing` where it was not measured.
"""
function lh16a_measured(Ca_Si, Al_Si, column)
    t = lh16a_table("gel_composition")
    rows = [k for k in eachindex(t.time_d) if t.Ca_Si_target[k] == Ca_Si && t.Al_Si_target[k] == Al_Si && !ismissing(t.H2O_Si[k])]
    isempty(rows) && return missing
    k = argmin(k -> ustrip(us"d", t.time_d[k]), rows)
    return getproperty(t, column)[k]
end

"""
    lh16a_measured_solution(Ca_Si, Al_Si, element) -> Vector{Float64}

The concentrations of `element` measured in the solutions of the batch at the
target `Ca_Si` and `Al_Si`, one per time, in mmol/L, those below the detection
limit left out; `"pH"` for the pH.
"""
function lh16a_measured_solution(Ca_Si, Al_Si, element)
    if element == "pH"
        t = lh16a_table("pH")
        return Float64[t.pH[i] for i in eachindex(t.pH) if t.Ca_Si_target[i] == Ca_Si && t.Al_Si_target[i] == Al_Si]
    end
    t = lh16a_table("pore_solution")
    return Float64[
        ustrip(us"mmol/L", t.concentration[i]) for i in eachindex(t.element)
            if t.Ca_Si_target[i] == Ca_Si && t.Al_Si_target[i] == Al_Si && t.element[i] == element && t.qualifier[i] == "measured"
    ]
end

"""
    lh16a_aluminum_split(eq, gel) -> NamedTuple

Where the aluminum of the equilibrium `eq` is, in percent of the whole: in the
gel `gel`, in katoite (with its siliceous member), in strätlingite, in
microcrystalline gibbsite, in the other solids and in solution.
"""
function lh16a_aluminum_split(eq, gel)
    cs = eq.system
    n = ustrip.(us"mol", eq.n)
    al(i) = n[i] * Float64(get(atoms(cs.species[i]), :Al, 0))
    members = Set(symbol(m) for p in cs.solid_solutions for m in p.end_members)
    groups = (katoite = LH16A_OTHER.katoite, straetlingite = LH16A_OTHER.straetlingite, gibbsite = ("AlOHmic",))
    held = Dict(k => 0.0 for k in (:gel, :katoite, :straetlingite, :gibbsite, :other, :solution))
    for i in eachindex(n)
        s = symbol(cs.species[i])
        k = i in cs.idx_aqueous ? :solution : s in members ? :gel :
            something(findfirst(names -> s in names, groups), :other)
        held[k] += al(i)
    end
    total = sum(values(held))
    return NamedTuple{(:gel, :katoite, :straetlingite, :gibbsite, :other, :solution)}(
        Tuple(100 * held[k] / total for k in (:gel, :katoite, :straetlingite, :gibbsite, :other, :solution))
    )
end
