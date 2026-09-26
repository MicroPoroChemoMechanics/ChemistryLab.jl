# =============================================================================
#  gruyaert2010.jl — the pastes of Gruyaert, Robeyst & De Belie (2010)
#
#  A CEM I 52.5 N and a blast-furnace slag, pastes at w/b = 0.5 with 0 to 85 %
#  slag, whose degrees of hydration were measured by image analysis (their
#  Table 5), whose heat was measured by isothermal calorimetry (Tables 2 and 6)
#  and whose bound water by thermogravimetry (Fig. 7), all transcribed in
#  data/literature/Gruyaert2010.json. The pages on isothermal calorimetry and
#  on thermogravimetry include this file and compute, at the measured degrees of
#  hydration, what those two instruments would read.
#
#  Written once here so that the two pages cannot drift apart. The assumptions
#  are the file's, and each is stated where it is made.
# =============================================================================

using ChemistryLab
using DynamicQuantities
using OptimaSolver
using OrderedCollections

const G10_DB = Dict(
    symbol(s) => s for s in build_species(datapath("cemdata18-thermofun.json"); verbose = false)
)

gruyaert_oxides(material) = literature_row("Gruyaert2010", "oxides", material)

# ASSUMED: the clinker phases of the cement from its oxides by Bogue's
# calculation, the lime bound to the sulfate (as gypsum) and to the carbon
# dioxide (as calcite) set aside first. The MgO, 0.9 to 1.0 %, is left out:
# CEMDATA18 carries no periclase, and a phase without an enthalpy cannot enter
# a heat balance. The article reports no alkali.
const G10_CLINKER = ("C3S", "C2S", "C3A", "C4AF")
const G10_OXIDES = ("CaO", "SiO2", "Al2O3", "Fe2O3")

_molar(s) = ustrip(us"g/mol", G10_DB[s][:M])
_oxide_molar(ox) = ustrip(us"g/mol", Species(ox)[:M])

# The element each oxide carries, and how many of its atoms per formula.
const G10_OXIDE_ELEMENT = Dict("CaO" => (:Ca, 1), "SiO2" => (:Si, 1), "Al2O3" => (:Al, 2), "Fe2O3" => (:Fe, 2))

# The mass fraction of the oxide `ox` in the phase `ph`, from its formula and the
# molar masses of the library.
function _oxide_fraction(ph, ox)
    el, per = G10_OXIDE_ELEMENT[ox]
    return Float64(get(atoms(G10_DB[ph]), el, 0)) / per * _oxide_molar(ox) / _molar(ph)
end

"""
    bogue(material) -> (; clinker, gypsum, calcite)

Mass fractions of the four clinker phases, of gypsum and of calcite in the cement
`material` of Table 1, from its oxides.
"""
function bogue(material)
    r = gruyaert_oxides(material)
    x(ox) = ustrip(getproperty(r, Symbol(ox))) / 100
    gypsum = x("SO3") * _molar("Gp") / _oxide_molar("SO3")
    calcite = x("CO2") * _molar("Cal") / _oxide_molar("CO2")
    cao = x("CaO") - x("SO3") * _oxide_molar("CaO") / _oxide_molar("SO3") -
        x("CO2") * _oxide_molar("CaO") / _oxide_molar("CO2")
    # The phases that reproduce the available oxides.
    M = [_oxide_fraction(ph, ox) for ox in G10_OXIDES, ph in G10_CLINKER]
    y = M \ [cao, x("SiO2"), x("Al2O3"), x("Fe2O3")]
    return (; clinker = OrderedDict(zip(G10_CLINKER, y)), gypsum, calcite)
end

const G10_PURE = split(
    "C3S C2S C3A C4AF Gp Anh Cal Portlandite ettringite monosulphate12 " *
        "monocarbonate hemicarbonate C4AH13 C3AH6 C3FH6 straetlingite " *
        "hydrotalcite Mg2AlC0.5OH Brc FeOOHmic AlOHmic Amor-Sl"
)
const G10_SS = [
    "CSHQ" => ["CSHQ-JenD", "CSHQ-JenH", "CSHQ-TobD", "CSHQ-TobH"],
    "C3(AF)S0.84H" => ["C3AFS0.84H4.32", "C3FS0.84H4.32"],
]

"""
    gruyaert_system() -> ChemicalSystem

The pure phases of the cement pages, CSHQ without its alkali end members (the
article reports no alkali), and the siliceous hydrogarnet. The slag's sulfur
is reported as SO3 and is taken as sulfate.
"""
function gruyaert_system()
    members = reduce(vcat, last.(G10_SS))
    sp = speciation(
        collect(values(G10_DB)), vcat(G10_PURE, members);
        aggregate_state = [AS_AQUEOUS], exclude_species = split("H2@ O2@ CH4@"),
    )
    # `enthalpy` skips a species without an enthalpy of formation, and the heat
    # would leave it out without a word.
    none = [symbol(s) for s in sp if !haskey(s, :ΔₐH⁰)]
    isempty(none) || error("no enthalpy of formation for " * join(none, ", "))
    ss = [SolidSolutionPhase(n, [G10_DB[m] for m in ms]) for (n, ms) in G10_SS]
    return ChemicalSystem(sp, CEMDATA_PRIMARIES; solid_solutions = ss)
end

"""
    gruyaert_paste(cs; slag, alpha_cement, alpha_slag, batch = "CAL", binder = 100.0)
        -> (; initial, b, eq, certificate)

100 g of binder at w/b = 0.5 with the fraction `slag` of slag, the cement reacted
to `alpha_cement` and the slag to `alpha_slag`, at 20 °C. `batch` is the material
of Table 1: "CAL" for the calorimetry, "TG/b" for the thermogravimetry after 2
days. The unreacted cement and slag are not part of the equilibrium.

ASSUMED: the four clinker phases at the one degree of hydration the image
analysis gives the cement; the gypsum and the calcite all in the equilibrium,
which dissolves what is not stable.

`initial` holds the reacted cement and the water, as species: its enthalpy is
that of the anhydrous reactants. The slag enters `b` through its oxides and has
no species, since a glass has no formula and no enthalpy in any database.
"""
function gruyaert_paste(cs; slag, alpha_cement, alpha_slag, batch = "CAL", binder = 100.0)
    wb = ustrip(literature_value("Gruyaert2010", "water_binder_ratio"))
    c = bogue("OPC-" * batch)
    m_cement = binder * (1 - slag)
    st = ChemicalState(cs; T = 293.15u"K")
    for (ph, f) in c.clinker
        set_quantity!(st, ph, alpha_cement * m_cement * f / _molar(ph) * u"mol")
    end
    set_quantity!(st, "Gp", m_cement * c.gypsum / _molar("Gp") * u"mol")
    set_quantity!(st, "Cal", m_cement * c.calcite / _molar("Cal") * u"mol")
    set_quantity!(st, "H2O@", binder * wb / _molar("H2O@") * u"mol")
    b = Float64.(cs.SM.A) * ustrip.(us"mol", st.n)
    if slag > 0 && alpha_slag > 0
        r = gruyaert_oxides("BFS-" * batch)
        ox = Dict(o => ustrip(getproperty(r, Symbol(o))) / 100 for o in ("CaO", "SiO2", "Al2O3", "Fe2O3", "MgO", "SO3"))
        b .+= oxide_budget(ox, cs.SM.primaries; mass = alpha_slag * slag * binder * u"g")
    end
    model = HKFActivityModel(å = 0.0, Ḃ = G10_BDOT, Kₙ = 0.0)
    eq, cert = equilibrate_certified(st; model, b)
    return (; initial = st, b, eq, certificate = cert)
end

# The B-dot of the cement pages, identified from the activity coefficients a
# GEM-Selektor run printed for |z| = 1 and 2.
const G10_BDOT = let g = ChemistryLab.JSON.parsefile(
        joinpath(pkgdir(ChemistryLab), "test", "reference", "gems_cemdata18_portland.json")
    )
    lg1, lg2 = log10(g["gamma"]["z1"]), log10(g["gamma"]["z2"])
    (lg1 + (lg1 - lg2) / 3) / g["ionic_strength_mol_per_kg"]
end
