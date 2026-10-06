# =============================================================================
#  balonis2010_chloride.jl — the AFm phases of a model Portland system as
#  chloride is added
#
#  Balonis et al. (2010) equilibrated at 25 °C a mixture of 0.01 mol of C3A,
#  0.015 mol of portlandite and 0.01 mol of calcium sulfate in 60 mL of water,
#  without and with 0.0075 mol of calcite, with CaCl2 up to 2Cl/Al2O3 = 1, and
#  identified the phases by XRD and analyzed the solutions. The sulfate AFm
#  turns into Kuzel's salt, then into Friedel's salt; with calcite the
#  monocarbonate turns into Friedel's salt directly. This file computes the
#  same titration with Cemdata18, on the system of the companion page on
#  temperature (Balonis 2019, scripts/balonis2019_temperature.jl). The page
#  docs/src/examples/chloride_afm.md compares them.
# =============================================================================

using ChemistryLab
using DynamicQuantities
using OptimaSolver

isdefined(@__MODULE__, :b19_system) || include(joinpath(@__DIR__, "balonis2019_temperature.jl"))

const BA10 = "Balonis2010"
ba10_value(name) = ustrip(literature_value(BA10, name))
const BA10_MODEL = cemdata18_activity_model(:KOH)

"""
    ba10_state(cs; calcite = false, T_C = 25) -> ChemicalState

The mixture of Balonis et al. (2010): C3A, portlandite, the sulfate as
anhydrite and 60 g of water, with their calcite when `calcite`.
"""
function ba10_state(cs; calcite = false, T_C = ba10_value("temperature_C"))
    st = ChemicalState(cs; T = (T_C + 273.15) * u"K")
    set_quantity!(st, "H2O@", ustrip(us"cm^3", literature_value(BA10, "water")) * u"g")
    set_quantity!(st, "C3A", ba10_value("c3a") * u"mol")
    set_quantity!(st, "Portlandite", ba10_value("portlandite") * u"mol")
    set_quantity!(st, "Anh", ba10_value("caso4") * u"mol")
    calcite && set_quantity!(st, "Cal", ba10_value("calcite") * u"mol")
    return st
end

"""
    ba10_titration(ratios; calcite = false) -> Vector{NamedTuple}

The mixture with CaCl2 at each 2Cl/Al2O3 of `ratios` (increasing), each
equilibrium started from the previous one: the amounts of the AFm phases, of
ettringite and portlandite per mole of Al2O3, the calcium and chloride of the
solution (mmol/L), the pH, the volume of the solids (cm³), and whether each
certified.
"""
function ba10_titration(ratios; calcite = false)
    cs = b19_system()
    st = ba10_state(cs; calcite)
    A = Float64.(cs.SM.A)
    b0 = A * ustrip.(us"mol", st.n)
    col(s) = A[:, findfirst(x -> symbol(x) == s, cs.species)]
    al = ba10_value("c3a")
    out = NamedTuple[]
    prev = st
    for r in ratios
        cacl2 = r * al
        b = b0 .+ cacl2 .* col("Ca+2") .+ 2cacl2 .* col("Cl-")
        eq, cert = equilibrate_certified(prev; model = BA10_MODEL, b)
        push!(out, (; ratio = r, certified = cert.optimal, _ba10_observables(eq, al)...))
        prev = eq
    end
    return out
end

function _ba10_observables(eq, al)
    cs = eq.system
    n = ustrip.(us"mol", eq.n)
    # Every instance of a member: a solid solution split in two names its
    # second instance's members `name#2`.
    at(s) = sum((n[j] for j in eachindex(n) if symbol(cs.species[j]) == s || startswith(symbol(cs.species[j]), s * "#")); init = 0.0)
    V_L = ustrip(uconvert(us"L", volume(eq).liquid))
    dissolved(e) = 1000 * sum(n[k] * Float64(get(atoms(cs.species[k]), e, 0)) for k in cs.idx_aqueous) / V_L
    return (;
        monosulfate = at("monosulphate12") / al, kuzel = at("C4AsClH12") / al, friedel = at("C4AClH10") / al,
        monocarbonate = at("monocarbonate") / al, ettringite = (at("ettringite03_ss") + at("tricarboalu03")) / 3 / al,
        portlandite = at("Portlandite"), calcite = at("Cal"), Ca = dissolved(:Ca), Cl = dissolved(:Cl),
        pH = pH(eq, BA10_MODEL), solids = ustrip(uconvert(us"cm^3", volume(eq).solid)),
    )
end

"""
    ba10_phases(row) -> String

The phases present in a row of [`ba10_titration`](@ref), in the abbreviations of
the XRD figures of Balonis et al. (2010): E, Ms, Ks, Fs, Mc, P, Cc.
"""
function ba10_phases(row)
    out = String[]
    for (f, lab) in (
            (:ettringite, "E"), (:monosulfate, "Ms"), (:kuzel, "Ks"), (:friedel, "Fs"),
            (:monocarbonate, "Mc"), (:portlandite, "P"), (:calcite, "Cc"),
        )
        getproperty(row, f) > 1.0e-6 && push!(out, lab)
    end
    return join(out, ", ")
end

"""
    ba10_boundaries(rows) -> NamedTuple

Along a titration of [`ba10_titration`](@ref): the 2Cl/Al2O3 at which the
sulfate AFm is gone, at which Kuzel's salt first forms and at which it is gone
again, `NaN` where it never does.
"""
function ba10_boundaries(rows)
    first_r(f) = (i = findfirst(f, rows); i === nothing ? NaN : rows[i].ratio)
    ks = findfirst(r -> r.kuzel > 1.0e-6, rows)
    gone = ks === nothing ? NaN : first_r(r -> r.ratio > rows[ks].ratio && r.kuzel < 1.0e-6)
    return (; monosulfate_gone = first_r(r -> r.monosulfate < 1.0e-6), kuzel_first = ks === nothing ? NaN : rows[ks].ratio, kuzel_gone = gone)
end
