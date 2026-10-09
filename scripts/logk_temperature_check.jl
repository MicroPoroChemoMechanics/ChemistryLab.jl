# =============================================================================
#  logk_temperature_check.jl — Cemdata18's equilibrium constants from 5 to 90 °C
#
#  The temperature dependence of twelve equilibrium constants as Cemdata18
#  computes them (the HKF equations of state of its aqueous species, the heat
#  capacities of its solids, water from its own equation of state), against the
#  analytical expressions of the PHREEQC database, fits to measured constants
#  whose sources the file names in its comments where it names them. The file is
#  PHREEQC's phreeqc.dat of release 3.7.3, which datapath obtains. The page
#  docs/src/tutorials/validation_logk_temperature.md runs it.
# =============================================================================

using ChemistryLab
using DynamicQuantities

const LKT_DB = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-thermofun.json"); verbose = false))
const LKT_PHREEQC = datapath("phreeqc.dat")

# Each reaction: the line or the phase name that heads its record in
# phreeqc.dat, and the same reaction in Cemdata18 symbols, products positive.
# The aqueous silica of Cemdata18 is SiO2@, which PHREEQC writes H4SiO4: the
# same constant, water at unit activity.
const LKT_REACTIONS = [
    "water" => ("H2O = OH- + H+", Dict("OH-" => 1, "H+" => 1, "H2O@" => -1)),
    "bicarbonate" => ("CO3-2 + H+ = HCO3-", Dict("HCO3-" => 1, "CO3-2" => -1, "H+" => -1)),
    "carbon dioxide (aq)" => ("CO3-2 + 2 H+ = CO2 + H2O", Dict("CO2@" => 1, "H2O@" => 1, "CO3-2" => -1, "H+" => -2)),
    "bisulfate" => ("SO4-2 + H+ = HSO4-", Dict("HSO4-" => 1, "SO4-2" => -1, "H+" => -1)),
    "CO2(g)" => ("CO2(g)", Dict("CO2@" => 1, "CO2" => -1)),
    "calcite" => ("Calcite", Dict("Ca+2" => 1, "CO3-2" => 1, "Cal" => -1)),
    "aragonite" => ("Aragonite", Dict("Ca+2" => 1, "CO3-2" => 1, "Arg" => -1)),
    "dolomite" => ("Dolomite", Dict("Ca+2" => 1, "Mg+2" => 1, "CO3-2" => 2, "Ord-Dol" => -1)),
    "gypsum" => ("Gypsum", Dict("Ca+2" => 1, "SO4-2" => 1, "H2O@" => 2, "Gp" => -1)),
    "anhydrite" => ("Anhydrite", Dict("Ca+2" => 1, "SO4-2" => 1, "Anh" => -1)),
    "quartz" => ("Quartz", Dict("SiO2@" => 1, "Qtz" => -1)),
    "amorphous silica" => ("SiO2(a)", Dict("SiO2@" => 1, "Amor-Sl" => -1)),
]

"""
    lkt_analytic(header) -> Vector{Float64}

The six coefficients of the analytical expression of the record headed by
`header` in phreeqc.dat, `log K = A₁ + A₂T + A₃/T + A₄ log₁₀T + A₅/T² + A₆T²`,
T in kelvin, those the record does not give at zero. Where a record gives two
expressions (gypsum), the second, which its comment calls the better fit.
"""
function lkt_analytic(header)
    lines = readlines(LKT_PHREEQC)
    k = findfirst(l -> strip(first(split(l, '#'))) == header, lines)
    k === nothing && error("lkt_analytic: no record \"$header\" in phreeqc.dat")
    coefficients = nothing
    # A phase is headed by its name, its reaction on the next line.
    reaction_next = !occursin('=', header)
    for l in lines[(k + 1):end]
        s = strip(first(split(l, '#')))
        isempty(s) && continue
        if !startswith(s, "-")
            reaction_next && (reaction_next = false; continue)
            # The record ends where the next one starts: a line that is not an option.
            break
        end
        option, rest... = split(s)
        if lowercase(option) in ("-analytic", "-analytical", "-analytical_expression", "-a_e")
            a = zeros(6)
            for (i, x) in enumerate(rest)
                a[i] = parse(Float64, x)
            end
            coefficients = a
        end
    end
    coefficients === nothing && error("lkt_analytic: \"$header\" has no analytical expression")
    return coefficients
end

lkt_phreeqc(a, T) = a[1] + a[2] * T + a[3] / T + a[4] * log10(T) + a[5] / T^2 + a[6] * T^2

"""
    lkt_cemdata(ν, T) -> Float64

log₁₀ K of the reaction `ν` (symbol => coefficient, products positive) at `T`
kelvin and 1 bar, from the standard Gibbs energies of Cemdata18.
"""
function lkt_cemdata(ν, T)
    ΔG = sum(c * ustrip(us"J/mol", LKT_DB[s][:ΔₐG⁰](T = T * u"K", P = 1.0e5u"Pa"; unit = true)) for (s, c) in ν)
    return -ΔG / (R_GAS * T * log(10))
end

"""
    lkt_compare(; temperatures = 5:5:90) -> Vector{NamedTuple}

For each reaction: log K of both at 25 °C, and the difference of the two once
each is taken from its own value at 25 °C (what the temperature dependence
alone contributes), at 50 °C and at its largest over `temperatures`, °C.
"""
function lkt_compare(; temperatures = 5:5:90)
    return map(LKT_REACTIONS) do (name, (header, ν))
        a = lkt_analytic(header)
        p(T) = lkt_phreeqc(a, T + 273.15)
        c(T) = lkt_cemdata(ν, T + 273.15)
        shift(T) = (c(T) - c(25)) - (p(T) - p(25))
        worst = argmax(T -> abs(shift(T)), temperatures)
        (; name, phreeqc25 = p(25), cemdata25 = c(25), at50 = shift(50), worst_T = worst, worst = shift(worst))
    end
end
