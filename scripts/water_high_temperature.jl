# =============================================================================
#  water_high_temperature.jl — water and quartz from 0 to 1000 °C
#
#  Two equilibria of the solvent with nothing but itself and one mineral, from
#  0 to 1000 °C and up to 5000 bar, computed from the standard states of the
#  slop98 database (the equation of state of water for the solvent, HKF for the
#  aqueous species, heat-capacity functions for quartz): the ionization constant
#  of water, against the formulation Bandura and Lvov (2006) fitted to the
#  measurements (their Table 4), and the solubility of quartz, against the
#  equation Manning (1994) fitted to the measurements from 25 °C and 1 bar to
#  20 kbar. The page docs/src/tutorials/validation_high_temperature.md compares
#  them.
#
#  Nothing is fitted: the constants are those of the database.
# =============================================================================

using ChemistryLab
using DynamicQuantities

const HW_DB = Dict(symbol(s) => s for s in build_species(datapath("slop98-inorganic-thermofun.json"); verbose = false))

# The standard Gibbs energy of a species at T (K) and P (Pa), over RT ln 10.
hw_g10(name, T, P) = HW_DB[name][:ΔₐG⁰](T = T, P = P) / (R_GAS * T * log(10))

"""
    hw_pkw(T, P) -> Float64

The negative decimal logarithm of the ionization constant of water,
H₂O = H⁺ + OH⁻, on the molality scale, at `T` (K) and `P` (Pa).
"""
hw_pkw(T, P) = hw_g10("H+", T, P) + hw_g10("OH-", T, P) - hw_g10("H2O@", T, P)

"""
    hw_quartz_log_m(T, P) -> Float64

The decimal logarithm of the molality of dissolved silica in equilibrium with
quartz, quartz = SiO₂(aq), the activity coefficient of the neutral solute taken
as one: the decimal logarithm of the equilibrium constant.
"""
hw_quartz_log_m(T, P) = hw_g10("Qtz", T, P) - hw_g10("SiO2@", T, P)

"""
    hw_density(T, P) -> Float64

The density of water at `T` (K) and `P` (Pa), in g/cm³, by the equation of state
of the solvent.
"""
hw_density(T, P) = ustrip(us"g/cm^3", HW_DB["H2O@"][:M] / (HW_DB["H2O@"][:V⁰](T = T, P = P) * us"m^3/mol"))

"""
    hw_manning(T, P) -> Float64

The decimal logarithm of the molality of dissolved silica by the equation of
Manning (1994), at `T` (K), with the density of water of the solvent at `P` (Pa).
"""
function hw_manning(T, P)
    c = Dict(
        zip(
            literature_table("Manning1994", "log_k_coefficients").coefficient,
            ustrip.(literature_table("Manning1994", "log_k_coefficients").value)
        )
    )
    ρ = hw_density(T, P)
    return c["a"] + c["b"] / T + c["c"] / T^2 + c["d"] / T^3 + (c["e"] + c["f"] / T + c["g"] / T^2) * log10(ρ)
end

"""
    hw_bandura_lvov() -> NamedTuple

The ionization constants of Table 4 of Bandura and Lvov (2006) at the states
where the HKF equations hold, water denser than 350 kg/m³ and at most 5000 bar:
the temperature (°C and K), the pressure (Pa), their pKw and this package's.
"""
function hw_bandura_lvov()
    t = literature_table("BanduraLvov2006", "pKw")
    T_C, P, ref = Float64[], Float64[], Float64[]
    for (Tc, p, pk) in zip(t.temperature_C, ustrip.(us"Pa", t.pressure), t.pKw)
        T = Tc + 273.15
        p <= 5.0e8 || continue
        hw_density(T, p) >= 0.35 || continue
        push!(T_C, Tc); push!(P, p); push!(ref, pk)
    end
    ours = [hw_pkw(Tc + 273.15, p) for (Tc, p) in zip(T_C, P)]
    return (; T_C, T = T_C .+ 273.15, P, ref, ours)
end
