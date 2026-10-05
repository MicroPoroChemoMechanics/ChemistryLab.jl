# =============================================================================
#  wiebe_gaddy1940_co2.jl — carbon dioxide in water under pressure
#
#  Wiebe and Gaddy (1940) measured the solubility of carbon dioxide in pure
#  water from 12 to 40 °C and from 25 to 500 atm (their Table I). This file
#  computes it at equilibrium, the gas ideal or following the equation of state
#  of Peng and Robinson (1976), the dissolved carbon dioxide and its ions from
#  Cemdata18 (the HKF records of slop98), and splits the logarithm of the
#  dissolved amount into the terms each piece of physics contributes. The page
#  docs/src/tutorials/validation_co2_solubility.md compares them.
#
#  Nothing is fitted: the critical constants of the gas are those of
#  phreeqc.dat, the standard energies those of the database.
# =============================================================================

using ChemistryLab
using DynamicQuantities

const WG40 = "WiebeGaddy1940"
const WG40_DB = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-thermofun.json"); verbose = false))
const WG40_AQUEOUS = ["H2O@", "H+", "OH-", "CO2@", "HCO3-", "CO3-2"]

"""
    wg40_measured() -> NamedTuple

Table I of Wiebe and Gaddy (1940): the temperature in °C (`T_C`) and in K
(`T`), the total pressure in Pa (`P`), and the dissolved carbon dioxide in mol
per kg of water (`m`), from the volume at 0 °C and 1 atm they report through
the molar volume of the ideal gas there.
"""
function wg40_measured()
    t = literature_table(WG40, "co2_solubility")
    # 0 °C and the standard atmosphere, 101325 Pa by definition.
    V_stp = R_GAS * 273.15 / 101325
    T_C = Float64.(t.temperature_C)
    return (; T_C, T = T_C .+ 273.15, P = ustrip.(us"Pa", t.pressure), m = ustrip.(us"m^3/kg", t.solubility) ./ V_stp)
end

"""
    wg40_system(gas) -> ChemicalSystem

Water, its ions and the dissolved carbon dioxide over a gas of carbon dioxide,
ideal (`gas = :ideal`) or following the equation of state of Peng and Robinson
with the constants of phreeqc.dat (`gas = :real`).
"""
function wg40_system(gas::Symbol)
    co2 = gas === :real ? peng_robinson(WG40_DB["CO2"]) :
        gas === :ideal ? WG40_DB["CO2"] :
        throw(ArgumentError("wg40_system: the gas is :ideal or :real; got :$gas"))
    return ChemicalSystem(vcat([WG40_DB[k] for k in WG40_AQUEOUS], [co2]), ["H2O@", "H+", "CO3-2", "Zz"])
end

"""
    wg40_dissolved(cs, T, P; model = HKFActivityModel()) -> Float64

The carbon dissolved in one kilogram of water, in mol/kg, at equilibrium with
twenty moles of carbon dioxide at `T` (K) and `P` (Pa): enough that the gas
remains at every condition of the table, even where an ideal gas would dissolve
some nine moles per kilogram.
"""
function wg40_dissolved(cs, T, P; model = HKFActivityModel())
    st = ChemicalState(cs; T = T * u"K", P = P * u"Pa")
    set_quantity!(st, "H2O@", 1.0u"kg")
    set_quantity!(st, "CO2", 20.0u"mol")
    eq, cert = equilibrate_certified(st; model)
    cert.optimal || error("wg40_dissolved: no certified equilibrium at $T K and $P Pa")
    n = Dict(symbol(s) => ustrip(us"mol", q) for (s, q) in zip(cs.species, eq.n))
    n["CO2"] > 1 || error("wg40_dissolved: the gas is exhausted at $T K and $P Pa")
    kg_water = n["H2O@"] * ustrip(us"kg/mol", WG40_DB["H2O@"][:M])
    return (n["CO2@"] + n["HCO3-"] + n["CO3-2"]) / kg_water
end

"""
    wg40_terms(T, P) -> NamedTuple

The natural logarithm of the molality of the dissolved carbon dioxide in an
ideal solution under the pure gas, split into its terms,
``\\ln m = \\ln K_H + \\ln(P/P^\\circ) + \\ln\\varphi - \\Delta g_{aq}``:
`henry`, the constant of CO2(g) = CO2(aq) at `T` and ``P^\\circ``; `pressure`,
``\\ln(P/P^\\circ)``, all an ideal gas has; `fugacity`, the logarithm of the
fugacity coefficient of the gas; and `volume`, minus the rise of the standard
Gibbs energy of the dissolved carbon dioxide from ``P^\\circ`` to `P` over ``RT``,
the work of its partial molar volume.
"""
function wg40_terms(T, P)
    g(k, Pv) = WG40_DB[k][:ΔₐG⁰](T = T, P = Pv) / (R_GAS * T)
    st = ChemicalState(ChemicalSystem([peng_robinson(WG40_DB["CO2"])], ["CO2"]); T = T * u"K", P = P * u"Pa", n = [1.0u"mol"])
    return (;
        henry = g("CO2", P_STANDARD) - g("CO2@", P_STANDARD),
        pressure = log(P / P_STANDARD),
        fugacity = log(fugacity_coefficients(st)["CO2"]),
        volume = -(g("CO2@", P) - g("CO2@", P_STANDARD)),
    )
end

"""
    wg40_compute() -> NamedTuple

Every point of Table I: the measured molality and those computed with the gas
ideal and real, and the terms of `wg40_terms`.
"""
function wg40_compute()
    d = wg40_measured()
    n = length(d.T)
    ideal, real_gas = zeros(n), zeros(n)
    terms = Vector{NamedTuple{(:henry, :pressure, :fugacity, :volume), NTuple{4, Float64}}}(undef, n)
    # Behind an inference barrier: inferred from here, the whole certified
    # equilibrium is inferred again inside this function, which took the
    # compiler more than ten minutes and 4 GB; called through the barrier, it is
    # compiled once on its own, in under a minute.
    cs_ideal, cs_real = Base.inferencebarrier(wg40_system(:ideal)), Base.inferencebarrier(wg40_system(:real))
    for i in 1:n
        ideal[i] = wg40_dissolved(cs_ideal, d.T[i], d.P[i])
        real_gas[i] = wg40_dissolved(cs_real, d.T[i], d.P[i])
        terms[i] = wg40_terms(d.T[i], d.P[i])
    end
    return (; d..., ideal, real = real_gas, terms)
end
