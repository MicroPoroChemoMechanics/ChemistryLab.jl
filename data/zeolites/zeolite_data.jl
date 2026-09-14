# Standard thermodynamic data of Na- and K-based zeolites at 25 °C.
#
# TRANSCRIBED, NOT DERIVED. Every number below is copied from a published
# table; none is estimated, interpolated or adjusted here.
#
#   Na series (14 phases) — Table 6 of
#     Ma, B. & Lothenbach, B. (2020). Synthesis, characterization, and
#     thermodynamic study of selected Na-based zeolites.
#     Cement and Concrete Research 135, 106111.
#     doi:10.1016/j.cemconres.2020.106111        [CrossRef-verified]
#
#   K series (14 phases) — Table 5 of
#     Ma, B. & Lothenbach, B. (2021). Synthesis, characterization, and
#     thermodynamic study of selected K-based zeolites.
#     Cement and Concrete Research 148, 106537.
#     doi:10.1016/j.cemconres.2021.106537        [CrossRef-verified]
#
# Both papers state that log Ksp is referred to AlO2-, Ca+2, Na+, K+, SiO2@,
# Cl-, OH-, NO3- and H2O — the CEMDATA18 primary species. That is what makes
# the merge legitimate, and `build_database.jl` checks it rather than assuming
# it: it recomputes each log Ksp from the CEMDATA18 aqueous Gibbs energies and
# refuses to emit a database if any disagrees.
#
# `origin` records the paper's own footnote for S⁰ and Cp⁰: `:measured` where
# the authors determined it, `:additivity` where they built it from the
# elementary (hydro)oxide components, `:literature` where they carried over a
# published value. It is metadata, never used in a calculation.

"One zeolite as published: energies in kJ/mol, S⁰/Cp⁰ in J/(mol·K), V⁰ in cm³/mol."
struct ZeoliteRecord
    symbol::String
    name::String
    formula::String          # ThermoFun formula syntax
    logKsp::Float64
    logKsp_err::Float64
    ΔfG⁰::Float64
    ΔfH⁰::Float64
    S⁰::Float64
    Cp⁰::Float64
    V⁰::Float64
    products::Dict{String, Float64}   # aqueous products of congruent dissolution
    doi::String
    origin::Symbol
end

const DOI_NA = "10.1016/j.cemconres.2020.106111"
const DOI_K = "10.1016/j.cemconres.2021.106537"

# ── Na-based zeolites (Ma & Lothenbach 2020, Table 6) ────────────────────────
const ZEOLITES_NA = [
    ZeoliteRecord(
        "ANA-Na", "analcime (Na)", "Na2(Al2Si4)O12(H2O)2",
        -26.8, 0.8, -6139.7, -6575.84, 469.0, 425.0, 194.84,
        Dict("Na+" => 2.0, "AlO2-" => 2.0, "SiO2@" => 4.0, "H2O@" => 2.0),
        DOI_NA, :literature,
    ),
    ZeoliteRecord(
        "GIS-LSP-Na", "gismondine, low-silica P (Na)", "Na2(Al2Si2)O8(H2O)3.8",
        -19.6, 0.6, -4858.72, -5314.82, 374.0, 384.0, 153.49,
        Dict("Na+" => 2.0, "AlO2-" => 2.0, "SiO2@" => 2.0, "H2O@" => 3.8),
        DOI_NA, :additivity,
    ),
    ZeoliteRecord(
        "PHI-Na", "phillipsite (Na)", "Na2.5(Al2.5Si5.5)O16(H2O)5",
        -39.4, 1.2, -8717.83, -9438.72, 692.0, 620.0, 304.74,
        Dict("Na+" => 2.5, "AlO2-" => 2.5, "SiO2@" => 5.5, "H2O@" => 5.0),
        DOI_NA, :additivity,
    ),
    ZeoliteRecord(
        "PHI-NaK", "phillipsite (Na,K)", "Na1.5K(Al2.5Si5.5)O16(H2O)5",
        -39.9, 1.2, -8741.26, -9461.67, 707.0, 626.0, 304.74,
        Dict("Na+" => 1.5, "K+" => 1.0, "AlO2-" => 2.5, "SiO2@" => 5.5, "H2O@" => 5.0),
        DOI_NA, :additivity,
    ),
    ZeoliteRecord(
        "LTA-Na", "zeolite LTA (Na)", "Na1.98(Al1.98Si2.02)O8(H2O)5.31",
        -18.2, 0.6, -5203.75, -5701.89, 584.0, 513.0, 186.95,
        Dict("Na+" => 1.98, "AlO2-" => 1.98, "SiO2@" => 2.02, "H2O@" => 5.31),
        DOI_NA, :additivity,
    ),
    ZeoliteRecord(
        "LTA-4A-Na", "zeolite 4A (Na)", "Na2(Al2Si2)O8(H2O)4.5",
        -20.5, 0.6, -5029.88, -5486.36, 536.0, 475.0, 187.0,
        Dict("Na+" => 2.0, "AlO2-" => 2.0, "SiO2@" => 2.0, "H2O@" => 4.5),
        DOI_NA, :additivity,
    ),
    ZeoliteRecord(
        "SOD-OH-Na", "sodalite, hydroxy (Na)", "Na8(Al6Si6)O24(OH)2(H2O)2",
        -65.2, 2.0, -13221.4, -14120.1, 943.0, 895.0, 424.74,
        Dict("Na+" => 8.0, "AlO2-" => 6.0, "SiO2@" => 6.0, "OH-" => 2.0, "H2O@" => 2.0),
        DOI_NA, :measured,
    ),
    ZeoliteRecord(
        "SOD-Cl-Na", "sodalite, chloride (Na)", "Na8(Al6Si6)O24Cl2",
        -69.4, 2.1, -12719.1, -13473.4, 848.0, 812.0, 421.53,
        Dict("Na+" => 8.0, "AlO2-" => 6.0, "SiO2@" => 6.0, "Cl-" => 2.0),
        DOI_NA, :measured,
    ),
    ZeoliteRecord(
        "CAN-NO3-Na", "cancrinite, nitrate (Na)", "Na8(Al6Si6)O24(NO3)2(H2O)4",
        -64.8, 1.9, -13600.8, -14717.6, 1149.0, 1119.0, 435.96,
        Dict("Na+" => 8.0, "AlO2-" => 6.0, "SiO2@" => 6.0, "NO3-" => 2.0, "H2O@" => 4.0),
        DOI_NA, :additivity,
    ),
    ZeoliteRecord(
        "CHA-Na", "chabazite (Na)", "Na2(Al2Si4)O12(H2O)6",
        -31.9, 1.0, -7117.55, -7808.31, 548.0, 578.0, 249.95,
        Dict("Na+" => 2.0, "AlO2-" => 2.0, "SiO2@" => 4.0, "H2O@" => 6.0),
        DOI_NA, :additivity,
    ),
    ZeoliteRecord(
        "FAU-X-Na", "zeolite X, faujasite (Na)", "Na2(Al2Si2.5)O9(H2O)6.2",
        -21.9, 0.7, -5857.79, -6456.94, 566.0, 586.0, 195.8,
        Dict("Na+" => 2.0, "AlO2-" => 2.0, "SiO2@" => 2.5, "H2O@" => 6.2),
        DOI_NA, :literature,
    ),
    ZeoliteRecord(
        "FAU-Y-Na", "zeolite Y, faujasite (Na)", "Na2(Al2Si4)O12(H2O)8",
        -29.5, 0.9, -7578.22, -8352.62, 734.0, 739.0, 282.94,
        Dict("Na+" => 2.0, "AlO2-" => 2.0, "SiO2@" => 4.0, "H2O@" => 8.0),
        DOI_NA, :literature,
    ),
    ZeoliteRecord(
        "NAT-Na", "natrolite (Na)", "Na2(Al2Si3)O10(H2O)2",
        -26.6, 0.8, -5305.15, -5707.02, 360.0, 359.0, 169.36,
        Dict("Na+" => 2.0, "AlO2-" => 2.0, "SiO2@" => 3.0, "H2O@" => 2.0),
        DOI_NA, :measured,
    ),
    ZeoliteRecord(
        "MOR-Na", "mordenite (Na)", "Na0.72(Al0.72Si5.28)O12(H2O)2.71",
        -22.5, 0.7, -5955.95, -6442.4, 388.0, 405.0, 210.59,
        Dict("Na+" => 0.72, "AlO2-" => 0.72, "SiO2@" => 5.28, "H2O@" => 2.71),
        DOI_NA, :additivity,
    ),
]

# ── K-based zeolites (Ma & Lothenbach 2021, Table 5) ─────────────────────────
const ZEOLITES_K = [
    ZeoliteRecord(
        "LEU-K", "leucite (K)", "K2(Al2Si4)O12",
        -27.6, 0.83, -5711.07, -6048.59, 360.0, 328.0, 177.35,
        Dict("K+" => 2.0, "AlO2-" => 2.0, "SiO2@" => 4.0),
        DOI_K, :literature,
    ),
    ZeoliteRecord(
        "GIS-LSP-K", "gismondine, low-silica P (K)", "K2(Al2Si2)O8(H2O)2",
        -19.6, 0.59, -4472.95, -4814.83, 364.0, 309.0, 140.26,
        Dict("K+" => 2.0, "AlO2-" => 2.0, "SiO2@" => 2.0, "H2O@" => 2.0),
        DOI_K, :measured,
    ),
    ZeoliteRecord(
        "GIS-P1-K", "gismondine P1 (K)", "K1.67(Al1.67Si2.33)O8(H2O)1.9",
        -21.2, 0.64, -4367.11, -4699.8, 347.0, 299.0, 140.34,
        Dict("K+" => 1.67, "AlO2-" => 1.67, "SiO2@" => 2.33, "H2O@" => 1.9),
        DOI_K, :measured,
    ),
    ZeoliteRecord(
        "PHI-K", "phillipsite (K)", "K2.5(Al2.5Si5.5)O16(H2O)5",
        -42.6, 1.28, -8787.69, -9546.58, 598.0, 639.0, 312.19,
        Dict("K+" => 2.5, "AlO2-" => 2.5, "SiO2@" => 5.5, "H2O@" => 5.0),
        DOI_K, :measured,
    ),
    ZeoliteRecord(
        "LTA-K", "zeolite LTA (K)", "K2(Al2Si2)O8(H2O)3.3",
        -20.5, 0.62, -4786.42, -5218.42, 365.0, 371.0, 186.82,
        Dict("K+" => 2.0, "AlO2-" => 2.0, "SiO2@" => 2.0, "H2O@" => 3.3),
        DOI_K, :measured,
    ),
    ZeoliteRecord(
        "CHA-K", "chabazite (K)", "K2(Al2Si4)O12(H2O)4",
        -32.3, 0.97, -6686.63, -7228.69, 607.0, 564.0, 252.91,
        Dict("K+" => 2.0, "AlO2-" => 2.0, "SiO2@" => 4.0, "H2O@" => 4.0),
        DOI_K, :measured,
    ),
    ZeoliteRecord(
        "FAU-X-K", "zeolite X, faujasite (K)", "K2.03(Al2.03Si2.47)O9(H2O)6.04",
        -22.5, 0.68, -5872.72, -6453.89, 618.0, 577.0, 223.48,
        Dict("K+" => 2.03, "AlO2-" => 2.03, "SiO2@" => 2.47, "H2O@" => 6.04),
        DOI_K, :measured,
    ),
    ZeoliteRecord(
        "FAU-Y-K", "zeolite Y, faujasite (K)", "K2.18(Al2.18Si3.82)O12(H2O)7.72",
        -32.35, 0.97, -7619.01, -8374.57, 772.0, 745.0, 291.27,
        Dict("K+" => 2.18, "AlO2-" => 2.18, "SiO2@" => 3.82, "H2O@" => 7.72),
        DOI_K, :measured,
    ),
    ZeoliteRecord(
        "NAT-tetra-K", "tetranatrolite (K)", "K2(Al2Si3)O10(H2O)2",
        -25.27, 0.76, -5338.72, -5731.87, 416.0, 370.0, 186.51,
        Dict("K+" => 2.0, "AlO2-" => 2.0, "SiO2@" => 3.0, "H2O@" => 2.0),
        DOI_K, :measured,
    ),
    ZeoliteRecord(
        "NAT-K", "natrolite (K)", "K2(Al2Si3)O10(H2O)2",
        -26.35, 0.79, -5344.89, -5738.04, 416.0, 370.0, 186.51,
        Dict("K+" => 2.0, "AlO2-" => 2.0, "SiO2@" => 3.0, "H2O@" => 2.0),
        DOI_K, :measured,
    ),
    ZeoliteRecord(
        "MOR-K", "mordenite (K)", "K0.65(Al0.65Si5.35)O12(H2O)2.3",
        -22.0, 0.66, -5851.31, -6323.1, 346.0, 388.0, 190.87,
        Dict("K+" => 0.65, "AlO2-" => 0.65, "SiO2@" => 5.35, "H2O@" => 2.3),
        DOI_K, :measured,
    ),
    ZeoliteRecord(
        "STI-K", "stilbite (K)", "K2.2(Al2.2Si6.8)O18(H2O)4.8",
        -45.2, 1.36, -9505.55, -10301.12, 630.0, 669.0, 316.67,
        Dict("K+" => 2.2, "AlO2-" => 2.2, "SiO2@" => 6.8, "H2O@" => 4.8),
        DOI_K, :additivity,
    ),
    ZeoliteRecord(
        "HEU-K", "heulandite (K)", "K2.22(Al2.22Si6.78)O18(H2O)4.7",
        -45.15, 1.35, -9487.07, -10289.26, 586.0, 665.0, 324.8,
        Dict("K+" => 2.22, "AlO2-" => 2.22, "SiO2@" => 6.78, "H2O@" => 4.7),
        DOI_K, :measured,
    ),
    ZeoliteRecord(
        "CLI-K", "clinoptilolite (K)", "K1.01(Al1.01Si4.99)O12(H2O)2.3",
        -26.8, 0.8, -5978.26, -6448.47, 378.0, 395.0, 191.26,
        Dict("K+" => 1.01, "AlO2-" => 1.01, "SiO2@" => 4.99, "H2O@" => 2.3),
        DOI_K, :measured,
    ),
]

const ZEOLITES = vcat(ZEOLITES_NA, ZEOLITES_K)
