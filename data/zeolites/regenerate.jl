# Build `data/cemdata18-zeolites.json`: CEMDATA18 plus the Empa zeolite series.
#
#   julia --project=docs data/zeolites/regenerate.jl
#
# The file this writes is VENDORED — it ships with the package, and this script
# exists so that it can be reproduced and audited rather than trusted. Run it to
# check that the shipped database is what the two published tables say it is.
#
# What this does, and — as importantly — what it refuses to do.
#
# CEMDATA18 is copied through **verbatim**: all 228 substances, the reactions,
# the elements, byte-for-byte the values ChemistryLab already ships. Nothing is
# recalibrated, nothing is overwritten. The 28 zeolites of `zeolite_data.jl` are
# appended under new symbols, so `natrolite` and `NAT-Na` coexist and a caller
# chooses; `README.md` in this directory records why their `ΔfG⁰` differ.
#
# Merging two thermodynamic datasets is only defensible if they share a
# reference state, and the usual failure is silent: an offset on Na+ of a few
# kJ/mol moves every dissolution equilibrium by an order of magnitude without
# any solver complaining. So this script *proves* the compatibility instead of
# asserting it. Each paper publishes both ΔfG⁰ and log Ksp for the same phase,
# referred to the CEMDATA18 primary species. Recomputing one from the other
# through CEMDATA18's own aqueous Gibbs energies closes the loop: agreement can
# only happen if the reference states coincide and the transcription is exact.
# A disagreement anywhere aborts the build.

import ChemistryLab
using ChemistryLab: parse_formula

# `JSON` rather than a dependency of its own: it is already what the package
# reads ThermoFun files with, so this script needs no environment beyond the one
# ChemistryLab itself resolves.
const JSON = ChemistryLab.JSON

const HERE = @__DIR__
const SRC = joinpath(pkgdir(ChemistryLab), "data", "cemdata18-thermofun.json")
const OUT = joinpath(pkgdir(ChemistryLab), "data", "cemdata18-zeolites.json")

include(joinpath(HERE, "zeolite_data.jl"))

const R = 8.31446261815324          # J/(mol·K), CODATA
const T₀ = 298.15                   # K
const RTln10 = R * T₀ * log(10)     # J/mol
const LOGK_TOL = 0.05               # log units
const CM3_PER_MOL_TO_J_PER_BAR = 0.1

# ── the CEMDATA18 side ───────────────────────────────────────────────────────

"Gibbs energies of formation (J/mol) of every CEMDATA18 substance, by symbol."
function cemdata_gibbs(db)
    g = Dict{String, Float64}()
    for s in db["substances"]
        v = get(get(s, "sm_gibbs_energy", Dict()), "values", nothing)
        v === nothing || isempty(v) || (g[String(s["symbol"])] = Float64(v[1]))
    end
    return g
end

# ── check 1: the dissolution reaction balances ───────────────────────────────
#
# Reading the stoichiometry back out of the formula also proves the formula
# string parses under the same parser ChemistryLab will use to load it — a
# formula this script accepted but the package could not read would fail much
# later and much less legibly.

const AQUEOUS_FORMULA = Dict(
    "Na+" => "Na|+|", "K+" => "K|+|", "AlO2-" => "AlO2|-|",
    "SiO2@" => "SiO2", "H2O@" => "H2O", "OH-" => "OH|-|",
    "Cl-" => "Cl|-|", "NO3-" => "NO3|-|",
)

const AQUEOUS_CHARGE = Dict(
    "Na+" => 1.0, "K+" => 1.0, "AlO2-" => -1.0, "SiO2@" => 0.0,
    "H2O@" => 0.0, "OH-" => -1.0, "Cl-" => -1.0, "NO3-" => -1.0,
)

function check_balance(z::ZeoliteRecord)
    lhs = Dict{Symbol, Float64}(k => Float64(v) for (k, v) in parse_formula(z.formula))
    rhs = Dict{Symbol, Float64}()
    charge = 0.0
    for (sp, ν) in z.products
        for (el, n) in parse_formula(AQUEOUS_FORMULA[sp])
            el === :Zz && continue
            rhs[el] = get(rhs, el, 0.0) + ν * n
        end
        charge += ν * AQUEOUS_CHARGE[sp]
    end
    # the solid dissolves congruently, so H2O@ on the right also feeds back the
    # structural water on the left: compare totals element by element
    for el in union(keys(lhs), keys(rhs))
        el === :Zz && continue
        Δ = get(lhs, el, 0.0) - get(rhs, el, 0.0)
        abs(Δ) > 1.0e-9 && error(
            "$(z.symbol): dissolution is not balanced in $el " *
                "(solid $(get(lhs, el, 0.0)), products $(get(rhs, el, 0.0)))",
        )
    end
    abs(charge) > 1.0e-9 &&
        error("$(z.symbol): dissolution products carry a net charge of $charge")
    return nothing
end

# ── check 2: log Ksp recomputed through CEMDATA18 ────────────────────────────

function recomputed_logK(z::ZeoliteRecord, G::Dict{String, Float64})
    ΔrG = sum(ν * G[sp] for (sp, ν) in z.products) - z.ΔfG⁰ * 1000
    return -ΔrG / RTln10
end

# ── emission ─────────────────────────────────────────────────────────────────

function substance_entry(z::ZeoliteRecord)
    return Dict{String, Any}(
        "name" => z.name,
        "symbol" => z.symbol,
        "formula" => z.formula,
        "formula_charge" => 0,
        "aggregate_state" => Dict("3" => "AS_CRYSTAL"),
        "class_" => Dict("0" => "SC_COMPONENT"),
        "Tst" => T₀,
        "Pst" => 100000,
        "sm_gibbs_energy" => Dict("values" => [z.ΔfG⁰ * 1000]),
        "sm_enthalpy" => Dict("values" => [z.ΔfH⁰ * 1000]),
        "sm_entropy_abs" => Dict("values" => [z.S⁰]),
        "sm_heat_capacity_p" => Dict("values" => [z.Cp⁰]),
        "sm_volume" => Dict("values" => [z.V⁰ * CM3_PER_MOL_TO_J_PER_BAR]),
        "datasources" => [z.doi],
        # provenance kept on the entry itself, so a phase can never be separated
        # from the paper it came from
        "zeolite_provenance" => Dict(
            "doi" => z.doi,
            "log_Ksp_298K" => z.logKsp,
            "log_Ksp_error" => z.logKsp_err,
            "S_Cp_origin" => String(z.origin),
            "transcribed_from" => "published table, verbatim",
        ),
    )
end

function main()
    db = JSON.parsefile(SRC; dicttype = Dict{String, Any})
    G = cemdata_gibbs(db)
    n_cem = length(db["substances"])
    existing = Set(String(s["symbol"]) for s in db["substances"])

    println("CEMDATA18 source : $(n_cem) substances")
    println("zeolites to add  : $(length(ZEOLITES))\n")

    worst = 0.0
    println(rpad("zeolite", 14), lpad("log Ksp", 9), lpad("recomputed", 12), lpad("Δ", 8))
    for z in ZEOLITES
        z.symbol in existing &&
            error("$(z.symbol) would overwrite a CEMDATA18 substance — refusing")
        check_balance(z)
        lk = recomputed_logK(z, G)
        Δ = lk - z.logKsp
        worst = max(worst, abs(Δ))
        flag = abs(Δ) > LOGK_TOL ? "  ✗" : ""
        println(
            rpad(z.symbol, 14), lpad(round(z.logKsp; digits = 2), 9),
            lpad(round(lk; digits = 2), 12), lpad(round(Δ; digits = 3), 8), flag,
        )
    end
    println()
    if worst > LOGK_TOL
        error(
            "log Ksp mismatch of $(round(worst, digits = 3)) log units exceeds the " *
                "tolerance of $LOGK_TOL: the two datasets do not share a reference " *
                "state, or a value was mistranscribed. No database written.",
        )
    end
    println("largest discrepancy: $(round(worst, digits = 4)) log units — within $LOGK_TOL\n")

    append!(db["substances"], substance_entry.(ZEOLITES))
    db["thermodataset"] = "cemdata18-zeolites"
    db["zeolite_extension"] = Dict(
        "base" => "cemdata18 (Lothenbach et al. 2019, doi:10.1016/j.cemconres.2018.04.018)",
        "base_substances" => n_cem,
        "added_substances" => length(ZEOLITES),
        "sources" => [DOI_NA, DOI_K],
        "verification" => Dict(
            "method" => "log Ksp recomputed from CEMDATA18 aqueous Gibbs energies",
            "tolerance_log_units" => LOGK_TOL,
            "worst_discrepancy_log_units" => worst,
        ),
        "note" => "CEMDATA18 entries are copied unchanged; nothing is overwritten.",
    )
    open(OUT, "w") do io
        JSON.print(io, db, 2)
    end
    println("wrote $OUT  ($(length(db["substances"])) substances)")
    return nothing
end

abspath(PROGRAM_FILE) == (@__FILE__) && main()
