# =============================================================================
#  seawater_flushing.jl — a hydrated CEM I flushed with 100 L of seawater
#
#  De Weerdt and Justnes (2015) dripped 100 L of Trondheim fjord seawater, in
#  some two thousand fillings of a 50 mL Soxhlet chamber over six weeks, through
#  50 g of a ground, well-hydrated CEM I 42.5 R paste (w/c 0.4, about 90 %
#  reacted), and measured the elements of the paste dried at 105 °C before and
#  after (their Fig. 4), and its phases. This file computes the exposure as
#  successive renewals of the pore solution by that seawater, with Cemdata18
#  and the phases of the seawater page. The page
#  docs/src/examples/seawater_flushing.md compares the two. Nothing is fitted.
#
#  ASSUMED, each stated on the page:
#  - the 4.1 % of limestone, as calcite, is part of the CaO of the oxide
#    analysis (Table 1), and is taken out of the clinker;
#  - the clinker 90 % reacted (the authors' estimate), at 20 °C, w/c 0.4;
#  - the seawater of Table 3 has its charge closed on sodium (as printed, its
#    anions exceed its cations by 0.071 eq/L), 1 mL weighing 1 g of water;
#  - the 21 % of bound water is per mass of ignited paste, so that the 50 g of
#    moist paste (38 % evaporable water) hold 25.6 g of cement; each filling is
#    then 195 mL per 100 g of cement, 100 L is 3.9e5 mL per 100 g, and the moist
#    paste retains 74 g of pore water per 100 g of cement when it is dried;
#  - each filling replaces the whole pore solution.
# =============================================================================

using ChemistryLab
using DynamicQuantities
using OptimaSolver

isdefined(@__MODULE__, :sw_system) || include(joinpath(@__DIR__, "seawater.jl"))

const DJ_KEY = "DeWeerdtJustnes2015"
dj_value(name) = ustrip(literature_value(DJ_KEY, name))
dj_table(name) = literature_table(DJ_KEY, name)

# The elements Fig. 4 reports.
const DJ_ELEMENTS = (:Mg, :Na, :Cl, :S, :Al, :Fe, :K, :Ca)
_dj_M(e) = ustrip(us"g/mol", Species(String(e))[:M])

"""
    dj_cement() -> NamedTuple

The masses (g per 100 g of cement) of the clinker and of the limestone, and the
clinker's oxide fractions: the oxides of Table 1, the CaO of the limestone taken
out of the CaO of the analysis, the free CaO counted in the CaO.
"""
function dj_cement()
    t = dj_table("cement_oxides")
    ox = Dict(String(o) => ustrip(p) for (o, p) in zip(t.oxide, t.percent))
    limestone = pop!(ox, "limestone")
    pop!(ox, "CaO_free")
    # The CaO the calcite carries, in the CaO of the analysis.
    ox["CaO"] -= limestone * ustrip(us"g/mol", Species("CaO")[:M]) / ustrip(us"g/mol", Species("CaCO3")[:M])
    clinker = sum(values(ox))
    return (; clinker, limestone, oxides = Dict(k => v / clinker for (k, v) in ox))
end

"""
    dj_recipe(; T = 293.15) -> Recipe

The CEM I 42.5 R of De Weerdt and Justnes (2015): its clinker, 90 % reacted, and
its limestone as calcite, at w/c 0.4 and `T` (K).
"""
function dj_recipe(; T = 293.15)
    c = dj_cement()
    clinker = Material(
        "CEM I 42.5 R (De Weerdt and Justnes 2015), clinker", :cement;
        constituents = AbstractConstituent[
            OxideConstituent("clinker", c.oxides; mass_fraction = 1.0, extent = dj_value("degree_of_reaction_percent") / 100),
        ],
    )
    limestone = Material("limestone", :other; constituents = AbstractConstituent[MineralConstituent(SW_DB["Cal"]; mass_fraction = 1.0)])
    tot = c.clinker + c.limestone
    return Recipe(clinker => c.clinker / tot, limestone => c.limestone / tot; w_b = dj_value("water_cement_ratio"), T = T * u"K")
end

"""
    dj_sample() -> NamedTuple

What the paper's numbers give for the sample, per 100 g of cement: the cement in
the 50 g of moist paste (g), one filling and the whole exposure (mL), and the
pore water the moist paste retains (g). The bound water is read per mass of
ignited paste (ASSUMED); `cement_moist` is the cement if it were per mass of
moist paste instead.
"""
function dj_sample()
    moist = ustrip(us"g", literature_value(DJ_KEY, "paste_mass"))
    evap = dj_value("evaporable_water_percent") / 100
    bound = dj_value("bound_water_percent") / 100
    dry = moist * (1 - evap)
    cement = dry / (1 + bound)
    cement_moist = dry - bound * moist
    filling = ustrip(us"mL", literature_value(DJ_KEY, "chamber_volume")) / cement * 100
    total = ustrip(us"mL", literature_value(DJ_KEY, "seawater_total")) / cement * 100
    retained = moist * evap / cement * 100
    return (; cement, cement_moist, filling, total, retained)
end

"""
    dj_seawater(cs) -> Vector{Float64}

The budget of 1 mL of the Trondheim seawater of Table 3 in the primaries of
`cs`: its ions, sodium closing the charge, and 1 g of water.
"""
function dj_seawater(cs)
    t = dj_table("seawater")
    k = findfirst(==("Trondheim fjord"), t.water)
    g = Dict(e => ustrip(us"g/L", getproperty(t, e)[k]) for e in (:Cl, :Na, :Mg, :S, :K, :Ca))
    mol = Dict(e => g[e] / _dj_M(e) / 1000 for e in keys(g))          # mol per mL
    ion = Dict(:Ca => "Ca+2", :Mg => "Mg+2", :K => "K+", :S => "SO4-2", :Cl => "Cl-")
    z = Dict(:Ca => 2, :Mg => 2, :K => 1, :S => -2, :Cl => -1)
    na = -sum(z[e] * mol[e] for e in keys(ion))                       # closes the charge
    col(s) = Float64.(cs.SM.A[:, findfirst(x -> symbol(x) == s, cs.species)])
    Mw = ustrip(us"g/mol", cs.dict_species["H2O@"][:M])
    return sum(mol[e] .* col(ion[e]) for e in keys(ion)) .+ na .* col("Na+") .+ (1 / Mw) .* col("H2O@")
end

"""
    dj_charge_gap() -> NamedTuple

The cations and anions of 1 L of the Trondheim seawater as printed (eq/L), and
the sodium that closes its charge (mol/L) against the sodium printed.
"""
function dj_charge_gap()
    t = dj_table("seawater")
    k = findfirst(==("Trondheim fjord"), t.water)
    m = Dict(e => ustrip(us"g/L", getproperty(t, e)[k]) / _dj_M(e) for e in (:Cl, :Na, :Mg, :S, :K, :Ca))
    cations = m[:Na] + m[:K] + 2m[:Mg] + 2m[:Ca]
    anions = m[:Cl] + 2m[:S]
    return (; cations, anions, Na_printed = m[:Na], Na_closed = m[:Na] + anions - cations)
end

"""
    dj_elements(state; retained = 0.0) -> Dict{Symbol, Float64}

The mass (g) of each element of Fig. 4 in the solids of `state`, plus the
solutes of `retained` g of its pore water, which drying leaves in the sample.
"""
function dj_elements(state; retained = 0.0)
    cs = state.system
    n = ustrip.(us"mol", state.n)
    iw = findfirst(x -> symbol(x) == "H2O@", cs.species)
    out = Dict(e => 0.0 for e in DJ_ELEMENTS)
    aqueous = Set(cs.idx_aqueous)
    for (i, sp) in enumerate(cs.species)
        f = i in aqueous ? (i == iw ? 0.0 : retained / (n[iw] * ustrip(us"g/mol", cs.species[iw][:M]))) : 1.0
        iszero(f) && continue
        for (e, k) in atoms(sp)
            haskey(out, e) && (out[e] += f * k * n[i] * _dj_M(e))
        end
    end
    return out
end

"""
    dj_flush(steps; T = 293.15) -> NamedTuple

The paste of [`dj_recipe`](@ref) equilibrated, then `steps` renewals of its pore
solution by equal fillings of the seawater adding up to the whole exposure of
[`dj_sample`](@ref): the initial state, the states after each renewal, and the
cumulative volume (mL per 100 g of cement).
"""
function dj_flush(steps; T = 293.15)
    cs = sw_system()
    rs, cert = equilibrate_certified(dj_recipe(; T), cs; model = SW_MODEL)
    cert.optimal || error("dj_flush: the paste did not certify.")
    V = dj_sample().total / steps
    p = leach(rs, steps; solution = V .* dj_seawater(cs))
    return (; initial = rs, states = p.states, volume = V .* (1:steps))
end

"""
    dj_measured() -> NamedTuple

The mass ratios to calcium of Fig. 4, original and exposed, and the fraction of
the calcium the exposed sample retained, referred to its aluminum, which the
exposure does not move (the authors read the rise of Al, Fe and S in % of the
dry sample as the loss of mass the calcium carries away).
"""
function dj_measured()
    t = dj_table("bulk_elements")
    v(col, e) = ustrip(getproperty(t, col)[findfirst(==(String(e)), t.element)])
    ratios(col) = Dict(e => v(col, e) / v(col, :Ca) for e in DJ_ELEMENTS)
    retained = (v(:exposed, :Ca) / v(:exposed, :Al)) / (v(:original, :Ca) / v(:original, :Al))
    return (; original = ratios(:original), exposed = ratios(:exposed), retained)
end

"""
    dj_path(steps = 2000; T = 293.15) -> Vector{NamedTuple}

The paste along the exposure, one row per renewal: the cumulative volume (mL per
100 g of cement), whether the equilibrium certified, the fraction of its calcium
retained (referred to its aluminum), the mass (g per 100 g of cement) of each
element of Fig. 4 in its solids, its mass ratios to calcium in the solids, the
same with the pore water the moist sample retains, and the amounts (mol per
100 g of cement) of its solids.
"""
function dj_path(steps = 2000; T = 293.15)
    r = dj_flush(steps; T)
    s = dj_sample()
    el0 = dj_elements(r.initial.state)
    caal0 = el0[:Ca] / el0[:Al]
    rows = NamedTuple[]
    for (k, st) in enumerate(r.states)
        el = dj_elements(st.state)
        elr = dj_elements(st.state; retained = s.retained)
        cs = st.state.system
        n = ustrip.(us"mol", st.state.n)
        solids = Dict(symbol(cs.species[i]) => n[i] for i in eachindex(n) if !(i in cs.idx_aqueous) && n[i] > 1.0e-9)
        push!(
            rows, (;
                V = r.volume[k], certified = st.certificate.optimal, retained = (el[:Ca] / el[:Al]) / caal0,
                elements = el, ratios = Dict(e => el[e] / el[:Ca] for e in DJ_ELEMENTS),
                ratios_pw = Dict(e => elr[e] / elr[:Ca] for e in DJ_ELEMENTS),
                solids,
            )
        )
    end
    return (; initial = el0, rows)
end

"""
    dj_match(rows, retained) -> NamedTuple

The volume at which the calcium retained along `rows` falls to `retained`, and
the mass ratios there, both interpolated linearly between the two renewals that
bracket it.
"""
function dj_match(rows, retained)
    k = findfirst(r -> r.retained <= retained, rows)
    (k === nothing || k == 1) && error("dj_match: the path does not cross $retained.")
    a, b = rows[k - 1], rows[k]
    w = (a.retained - retained) / (a.retained - b.retained)
    lerp(x, y) = x + w * (y - x)
    return (;
        V = lerp(a.V, b.V),
        ratios = Dict(e => lerp(a.ratios[e], b.ratios[e]) for e in DJ_ELEMENTS),
        ratios_pw = Dict(e => lerp(a.ratios_pw[e], b.ratios_pw[e]) for e in DJ_ELEMENTS),
        before = a.solids, after = b.solids,
    )
end
