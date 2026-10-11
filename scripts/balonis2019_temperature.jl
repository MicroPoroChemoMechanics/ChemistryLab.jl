# =============================================================================
#  balonis2019_temperature.jl — the chloride AFm phases from 0 to 99 °C
#
#  Balonis (2019) calculated four model mixtures of 0.01 mol C3A, 0.015 mol
#  portlandite and 60 ml of water with sulfate (SO3/Al2O3 = 1), chloride
#  (2Cl/Al2O3 = 0.5 or 1) and, in two of them, calcite (CO2/Al2O3 = 0.75), from
#  0 to 99 °C, and stated where Kuzel's salt, Friedel's salt, monocarbonate and
#  ettringite give way to monosulfate (Section 3, Figs. 5, 6, 11 and 12). This
#  file computes the same mixtures with the Cemdata18 records of the same
#  hydrates. The page docs/src/examples/chloride_temperature.md compares the
#  two calculations.
# =============================================================================

using ChemistryLab
using DynamicQuantities
using OptimaSolver
using OrderedCollections

const B19 = "Balonis2019"
b19_value(name) = ustrip(literature_value(B19, name))
const B19_DB = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-thermofun.json"); verbose = false))

# The phases the article's figures show, by their Cemdata18 records: Friedel's
# and Kuzel's salts, ettringite with its carbonate analog and monosulfate with
# the hydroxy-AFm as the two non-ideal solid solutions of Cemdata18,
# monocarbonate and hemicarbonate, calcite, portlandite and the calcium
# sulfates. C3A is in the list for the mixture to be given in it, and dissolves.
const B19_PURE = split("C3A Portlandite Gp Anh Cal monocarbonate hemicarbonate C4AClH10 C4AsClH12")
const B19_SS = ("AFm_SO4_OH", "AFt_SO4_CO3")

"""
    b19_system() -> ChemicalSystem

The aqueous species of Cemdata18 for Ca, Al, S, C, Cl, O and H, the phases of
`B19_PURE`, and the two solid solutions of `B19_SS`, with the non-ideal
parameters Cemdata18 gives them.

ASSUMED, where the article differs: Friedel's salt is pure, where the article
mixes it ideally with the hydroxy-AFm and with monocarbonate; hydrogarnet is
left out, as it is absent from the article's figures.
"""
function b19_system()
    ss = [p for p in build_solid_solutions(datapath("solid_solutions.toml"), B19_DB) if name(p) in B19_SS]
    members = [symbol(s) for p in ss for s in p.end_members]
    sp = speciation(
        collect(values(B19_DB)), vcat(B19_PURE, members);
        aggregate_state = [AS_AQUEOUS], exclude_species = split("H2@ O2@ CH4@"),
    )
    return ChemicalSystem(sp, CEMDATA_PRIMARIES; solid_solutions = ss)
end

"""
    b19_state(cs, figure, T_C) -> ChemicalState

The mixture of `figure` ("Fig. 5", "Fig. 6", "Fig. 11" or "Fig. 12") at `T_C`
°C: the C3A, the portlandite and the water of the captions, the sulfate as
anhydrite, the chloride as dissolved CaCl2, the carbonate as calcite.

ASSUMED: the 60 ml of water are 60 g.
"""
function b19_state(cs, figure, T_C)
    t = literature_table(B19, "temperature_scans")
    i = findfirst(==(figure), t.figure)
    c3a = b19_value("c3a")
    st = ChemicalState(cs; T = (T_C + 273.15) * u"K")
    set_quantity!(st, "H2O@", ustrip(us"cm^3", literature_value(B19, "water")) * u"g")
    set_quantity!(st, "C3A", c3a * u"mol")
    set_quantity!(st, "Portlandite", b19_value("portlandite") * u"mol")
    set_quantity!(st, "Anh", b19_value("so3_al2o3") * c3a * u"mol")
    cacl2 = ustrip(t.cl2_al2o3[i]) * c3a
    set_quantity!(st, "Ca+2", cacl2 * u"mol")
    set_quantity!(st, "Cl-", 2cacl2 * u"mol")
    set_quantity!(st, "Cal", ustrip(t.co2_al2o3[i]) * c3a * u"mol")
    return st
end

# The species of each phase as the article names them.
const B19_PHASES = OrderedDict(
    "Kuzel's salt" => ("C4AsClH12",), "Friedel's salt" => ("C4AClH10",),
    "ettringite" => ("ettringite03_ss", "tricarboalu03"), "monosulfate" => ("monosulphate12", "C4AH13"),
    "monocarbonate" => ("monocarbonate",), "hemicarbonate" => ("hemicarbonate",),
    "calcite" => ("Cal",), "portlandite" => ("Portlandite",), "C3A" => ("C3A",),
    "calcium sulfate" => ("Gp", "Anh"),
)

"""
    b19_solids(cs, figure, T_C; model = cemdata18_activity_model(:KOH)) -> OrderedDict

The volume, cm³, of each phase of `B19_PHASES` at the certified equilibrium of
the mixture of `figure` at `T_C` °C, and their total under `"total"`.
"""
function b19_solids(cs, figure, T_C; model = cemdata18_activity_model(:KOH))
    st = b19_state(cs, figure, T_C)
    eq, cert = equilibrate_certified(st; model, b = budget(st))
    cert.optimal || error("b19_solids: the mixture of $figure at $T_C °C is not certified")
    V = ustrip(uconvert(us"cm^3", volume(eq).total))
    solid = Set(symbol(s) for s in cs.species if aggregate_state(s) != AS_AQUEOUS)
    out = OrderedDict{String, Float64}(ph => 0.0 for ph in keys(B19_PHASES))
    for (k, f) in volume_fractions(eq)
        k in solid || continue
        # A second instance of a solid solution names its members `<member>#2`.
        member = first(split(k, '#'))
        ph = findfirst(m -> member in m, B19_PHASES)
        ph === nothing && (ph = "other")
        out[ph] = get(out, ph, 0.0) + V * f
    end
    out["total"] = sum(values(out))
    return out
end

"""
    b19_scan(cs, figure; temperatures = vcat(0:5:95, 99)) -> (; figure, temperatures, volumes)

The phases of `figure` at each temperature, °C.
"""
function b19_scan(cs, figure; temperatures = vcat(0:5:95, 99), model = cemdata18_activity_model(:KOH))
    # A loop, not a comprehension: inferring the element type of a
    # comprehension over the certified search compiled for minutes and took
    # gigabytes.
    volumes = OrderedDict{String, Float64}[]
    for T in temperatures
        push!(volumes, b19_solids(cs, figure, float(T); model))
    end
    return (; figure, temperatures = collect(temperatures), volumes)
end

"""
    b19_volumes(scan, T) -> OrderedDict

The volumes of `scan` at the temperature `T`, °C, one of its temperatures.
"""
b19_volumes(scan, T) = scan.volumes[findfirst(==(T), scan.temperatures)]

"""
    b19_transition(scan, phase, event) -> Union{Tuple, Nothing}

Between which two temperatures of `scan`, °C, `phase` stops being present on
heating (`event = "last"`) or starts (`"first"`), present meaning above 1e-3
cm³: the pair `(lo, hi)`, `hi = nothing` if it is still present at the last
temperature (or `lo = nothing` if already present at the first); `nothing` if
it is never present. Near a change of assemblage the certified search can take
half a minute an equilibrium, so the change is bracketed rather than bisected.
"""
function b19_transition(scan, phase, event)
    present = [v[phase] > 1.0e-3 for v in scan.volumes]
    any(present) || return nothing
    T = scan.temperatures
    if event == "last"
        k = findlast(present)
        return (T[k], k == length(T) ? nothing : T[k + 1])
    end
    k = findfirst(present)
    return (k == 1 ? nothing : T[k - 1], T[k])
end
