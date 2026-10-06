# =============================================================================
#  seawater.jl — a hydrated CEM I in seawater, in zero dimensions
#
#  De Weerdt, Justnes and Geiker (2014) cored a concrete of CEM I 42.5 R with
#  limestone filler exposed for 16 years to the sea, measured the phases and
#  the elements against the depth, and computed with GEMS the paste in contact
#  with increasing volumes of seawater: 100 g of cement (89 g of CEM I, 11 g of
#  limestone), 42 g of water, 70 % of the CEM I reacted (their Section 3.3 and
#  Fig. 20). This file computes the same titration with Cemdata18, the
#  chloride AFm phases, thaumasite, brucite, hydrotalcite and M-S-H allowed.
#  The page docs/src/examples/seawater.md compares the two calculations and
#  the zones of the core. Nothing is fitted.
#
#  ASSUMED: the seawater of their Table 1, its charge closed on sodium (as
#  printed, its anions exceed its cations by 0.04 eq/L), 1 mL weighing 1 g of
#  water; 20 °C, the temperature of the calculation being unstated.
# =============================================================================

using ChemistryLab
using DynamicQuantities
using OptimaSolver

const SW_SUBSTANCES = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
const SW_DB = Dict(symbol(s) => s for s in SW_SUBSTANCES)
const SW_KEY = "DeWeerdt2014"
sw_value(name) = ustrip(literature_value(SW_KEY, name))
sw_table(name) = literature_table(SW_KEY, name)

# The phases of a Portland paste, and those seawater can form beside them.
const SW_ADDED = ["C4AClH10", "C4AsClH12", "thaumasite", "M075SH", "M15SH"]
sw_system() = phase_list_system("Portland paste (Lothenbach and Winnefeld 2006)", SW_SUBSTANCES; add = SW_ADDED)

"""
    sw_recipe(; T = 293.15) -> Recipe

89 g of the CEM I 42.5 R of Table 3, 70 % of it reacted, 11 g of limestone
taken as calcite, and 42 g of water, at `T` (K).
"""
function sw_recipe(; T = 293.15)
    t = sw_table("cement_oxides")
    ox = Dict(String(o) => ustrip(p) / 100 for (o, p) in zip(t.oxide, t.percent) if o != "LOI")
    total = sum(values(ox))
    cement = Material(
        "CEM I 42.5 R (De Weerdt et al. 2014)", :cement;
        constituents = AbstractConstituent[
            OxideConstituent("clinker", Dict(k => v / total for (k, v) in ox); mass_fraction = 1.0, extent = sw_value("model_opc_reacted") / 100),
        ],
    )
    limestone = Material("limestone", :other; constituents = AbstractConstituent[MineralConstituent(SW_DB["Cal"]; mass_fraction = 1.0)])
    c, l, w = sw_value("model_cement"), sw_value("model_limestone"), sw_value("model_water")
    return Recipe(cement => c / (c + l), limestone => l / (c + l); w_b = w / (c + l), T = T * u"K")
end

"""
    sw_seawater(cs) -> Vector{Float64}

The budget of 1 mL of the seawater of Table 1 in the primaries of `cs`: its
ions, sodium closing the charge, and 1 g of water.
"""
function sw_seawater(cs)
    t = sw_table("seawater")
    g = Dict(String(e) => ustrip(us"g/L", c) for (e, c) in zip(t.element, t.concentration))
    M(e) = ustrip(us"g/mol", Species(e)[:M])
    mol = Dict(e => g[e] / M(e) / 1000 for e in keys(g))          # mol per mL
    ion = Dict("Ca" => "Ca+2", "Mg" => "Mg+2", "K" => "K+", "S" => "SO4-2", "Cl" => "Cl-", "C" => "HCO3-")
    z = Dict("Ca" => 2, "Mg" => 2, "K" => 1, "S" => -2, "Cl" => -1, "C" => -1)
    na = -sum(z[e] * mol[e] for e in keys(ion))                  # closes the charge
    col(s) = Float64.(cs.SM.A[:, findfirst(x -> symbol(x) == s, cs.species)])
    Mw = ustrip(us"g/mol", cs.dict_species["H2O@"][:M])
    return sum(mol[e] .* col(ion[e]) for e in keys(ion)) .+ na .* col("Na+") .+ (1 / Mw) .* col("H2O@")
end

"""
    sw_charge_gap() -> NamedTuple

The cations and anions of 1 L of the seawater as printed (eq/L), and the
sodium that closes its charge (mol/L) against the sodium printed.
"""
function sw_charge_gap()
    t = sw_table("seawater")
    g = Dict(String(e) => ustrip(us"g/L", c) for (e, c) in zip(t.element, t.concentration))
    M(e) = ustrip(us"g/mol", Species(e)[:M])
    m = Dict(e => g[e] / M(e) for e in keys(g))
    cations = m["Na"] + m["K"] + 2m["Mg"] + 2m["Ca"]
    anions = m["Cl"] + 2m["S"] + m["C"]
    return (; cations, anions, Na_printed = m["Na"], Na_closed = m["Na"] + anions - cations)
end

const SW_MODEL = cemdata18_activity_model(:NaOH)

"""
    sw_titration(volumes; T = 293.15) -> Vector{NamedTuple}

The paste in contact with each of `volumes` (mL of seawater per 100 g of
cement, increasing), each equilibrium started from the previous one: the volume,
whether it certified, and the mass (g per 100 g of cement) of the phases that
mark the sequence.
"""
function sw_titration(volumes; T = 293.15)
    cs = sw_system()
    rs, cert = equilibrate_certified(sw_recipe(; T), cs; model = SW_MODEL)
    cert.optimal || error("sw_titration: the paste did not certify.")
    sea = sw_seawater(cs)
    prev = rs.state
    out = NamedTuple[]
    for V in volumes
        eq, c = equilibrate_certified(prev; model = SW_MODEL, b = rs.b .+ V .* sea)
        push!(out, (; V, certified = c.optimal, _sw_masses(eq)...))
        prev = eq
    end
    return out
end

const SW_PHASES = (
    friedel = "C4AClH10", kuzel = "C4AsClH12", brucite = "Brc", monocarbonate = "monocarbonate",
    thaumasite = "thaumasite", portlandite = "Portlandite", calcite = "Cal", ettringite = "ettringite",
    gypsum = "Gp", hydrotalcite = "hydrotalcite", msh_low = "M075SH", msh_high = "M15SH",
)

function _sw_masses(eq)
    cs = eq.system
    n = ustrip.(us"mol", eq.n)
    mass(s) = (j = findfirst(x -> symbol(x) == s, cs.species); j === nothing ? 0.0 : n[j] * ustrip(us"g/mol", cs.species[j][:M]))
    return NamedTuple{keys(SW_PHASES)}(Tuple(mass(s) for s in values(SW_PHASES)))
end

"""
    sw_transitions(rows) -> Dict

For each phase of the sequence, the volume (mL per 100 g of cement) at which it
first appears and the last at which it is present, `NaN` where it never does.
"""
function sw_transitions(rows)
    out = Dict{Symbol, NamedTuple}()
    for p in keys(SW_PHASES)
        present = [getproperty(r, p) > 1.0e-6 for r in rows]
        i, j = findfirst(present), findlast(present)
        out[p] = (; first = i === nothing ? NaN : rows[i].V, last = j === nothing ? NaN : rows[j].V)
    end
    return out
end
