# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2026 Jean-François Barthélémy (Cerema, UMR MCD)

# How the dissolved species interact: the activity model of the pore solution.

"""
    ActivityChoice

The activity model of the aqueous phase.
"""
abstract type ActivityChoice <: Choice end

"Extended Debye–Hückel with the parameters Cemdata18 recommends, for a KOH or NaOH solution."
struct Cemdata18Activity <: ActivityChoice
    electrolyte::String
    temperature_dependent::Bool
end
kind(::Type{Cemdata18Activity}) = "cemdata18"
fields(::Type{Cemdata18Activity}) = [
    Field("electrolyte", :choice; default = "KOH", options = ["KOH", "NaOH"]),
    Field("temperature_dependent", :bool; default = true, level = :expert),
]
activity_lines(x::Cemdata18Activity) = [
    "# How the dissolved species interact: extended Debye–Hückel, with the",
    "# parameters Cemdata18 recommends for a pore solution dominated by $(x.electrolyte).",
    "model = cemdata18_activity_model(:$(x.electrolyte)$(x.temperature_dependent ? "" : "; temperature_dependent = false"))",
]

"Every activity coefficient equal to one: a dilute solution."
struct DiluteActivity <: ActivityChoice end
kind(::Type{DiluteActivity}) = "dilute"
fields(::Type{DiluteActivity}) = Field[]
activity_lines(::DiluteActivity) = ["# A dilute solution: every activity coefficient is one.", "model = DiluteSolutionModel()"]

"The equation of Davies."
struct DaviesActivity <: ActivityChoice
    b::Float64
end
kind(::Type{DaviesActivity}) = "davies"
fields(::Type{DaviesActivity}) = [Field("b", :number; default = 0.3, min = 0, max = 1)]
activity_lines(x::DaviesActivity) = ["# The equation of Davies.", "model = DaviesActivityModel(; b = $(julia_number(x.b)))"]

"Truesdell–Jones, with the parameters of PHREEQC's database."
struct TruesdellJonesActivity <: ActivityChoice end
kind(::Type{TruesdellJonesActivity}) = "truesdell_jones"
fields(::Type{TruesdellJonesActivity}) = Field[]
activity_lines(::TruesdellJonesActivity) = [
    "# Truesdell–Jones, with the ion-size parameters of PHREEQC's database.",
    "model = TruesdellJonesActivityModel(; parameters = phreeqc_gamma_parameters(datapath(\"phreeqc.dat\")))",
]

"The B-dot model of the LLNL database."
struct LLNLActivity <: ActivityChoice end
kind(::Type{LLNLActivity}) = "llnl"
fields(::Type{LLNLActivity}) = Field[]
activity_lines(::LLNLActivity) = [
    "# The B-dot model of the LLNL database.",
    "model = database_activity_model(import_database(datapath(\"llnl.dat\"))[2])",
]

"Specific ion interaction theory, with the coefficients of PHREEQC's sit.dat."
struct SITActivity <: ActivityChoice end
kind(::Type{SITActivity}) = "sit"
fields(::Type{SITActivity}) = Field[]
activity_lines(::SITActivity) = [
    "# Specific ion interaction theory, with the coefficients of PHREEQC's sit.dat.",
    "model = SITActivityModel(; parameters = build_sit_parameters(datapath(\"sit.dat\")))",
]

"Pitzer's model, with the binary parameters of Reardon (1990)."
struct PitzerActivity <: ActivityChoice
    missing_pairs::String
end
kind(::Type{PitzerActivity}) = "pitzer"
fields(::Type{PitzerActivity}) = [Field("missing_pairs", :choice; default = "refuse", options = ["refuse", "zero"])]
activity_lines(x::PitzerActivity) = [
    "# Pitzer's model, with the binary parameters of Reardon (1990); a pair of ions",
    "# without parameters is $(x.missing_pairs == "zero" ? "taken as not interacting" : "refused").",
    "model = PitzerActivityModel(; parameters = build_pitzer_parameters(datapath(\"pitzer-reardon1990.toml\")), missing_pairs = :$(x.missing_pairs))",
]

members(::Type{ActivityChoice}) =
    (Cemdata18Activity, DiluteActivity, DaviesActivity, TruesdellJonesActivity, LLNLActivity, SITActivity, PitzerActivity)
