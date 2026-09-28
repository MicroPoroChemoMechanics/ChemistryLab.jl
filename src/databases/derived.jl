# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using SHA: sha256

# ── Databases built by ChemistryLab from a published one ─────────────────────
#
# Two databases extend Cemdata18 with phases it does not carry: 28 zeolites
# (Ma & Lothenbach 2020, 2021) and a chloride end member of CSHQ fitted on the
# sorption tests of Hirao et al. (2005). Only what ChemistryLab adds lives in the
# package — the published zeolite data in `data/literature/`, the fitted
# parameter in `data/chloride/cshq_cl.json`. The database itself is assembled on
# first use from the downloaded Cemdata18 file, whose entries are copied through
# unchanged, and cached under a key made of the SHA-256 of everything it was
# built from, so a new release of either side rebuilds it.

"""
    DerivedDatabase

A database ChemistryLab builds from a published one (`base`, a name of
[`THIRD_PARTY_DATABASES`](@ref)) and data of its own (`inputs`, paths under the
package's `data/`). `build(base_path, out_path)` writes it.
"""
struct DerivedDatabase{B}
    name::String
    base::String
    description::String
    inputs::Vector{String}
    build::B
end

# Bumped whenever a builder changes what it writes, so that cached builds are
# redone.
const _DERIVED_VERSION = "1"

function _derived_key(d::DerivedDatabase, base_path)
    parts = [_DERIVED_VERSION, d.name, _digest(base_path)]
    append!(parts, (_digest(joinpath(pkgdir(@__MODULE__), "data", i)) for i in d.inputs))
    return bytes2hex(sha256(join(parts, "\n")))
end

function _derived_path(d::DerivedDatabase; download::Bool = true)
    base = database_path(d.base; download)
    key = _derived_key(d, base)
    dir = joinpath(database_cache(), "derived")
    out = joinpath(dir, d.name)
    keyfile = out * ".key"
    isfile(out) && isfile(keyfile) && strip(read(keyfile, String)) == key && return out
    mkpath(dir)
    tmp = out * ".part"
    try
        d.build(base, tmp)
        mv(tmp, out; force = true)
    finally
        rm(tmp; force = true)
    end
    write(keyfile, key)
    return out
end

function _write_json_atomically(path, db)
    open(path, "w") do io
        JSON.print(io, db, 2)
    end
    return path
end

# ── the zeolites of Ma & Lothenbach ──────────────────────────────────────────

const _ZEOLITE_SOURCES = ("MaLothenbach2020", "MaLothenbach2021")

# The published log Ksp is recomputed from the published ΔfG⁰ through the base's
# aqueous Gibbs energies; both come from the same laboratory and must share a
# reference state, which this checks for every phase before anything is added.
const _ZEOLITE_LOGK_TOLERANCE = 0.05
const _CM3_PER_MOL_TO_J_PER_BAR = 0.1

const _AQUEOUS_FORMULA = Dict(
    "Na+" => "Na|+|", "K+" => "K|+|", "AlO2-" => "AlO2|-|", "SiO2@" => "SiO2",
    "H2O@" => "H2O", "OH-" => "OH|-|", "Cl-" => "Cl|-|", "NO3-" => "NO3|-|",
)
const _AQUEOUS_CHARGE = Dict(
    "Na+" => 1.0, "K+" => 1.0, "AlO2-" => -1.0, "SiO2@" => 0.0,
    "H2O@" => 0.0, "OH-" => -1.0, "Cl-" => -1.0, "NO3-" => -1.0,
)

"""
    zeolite_records() -> Vector{NamedTuple}

The 28 zeolites of Ma & Lothenbach (2020, Na series; 2021, K series) as the two
papers publish them, read from `data/literature/MaLothenbach2020.json` and
`MaLothenbach2021.json`: energies in kJ/mol, S⁰ and Cp⁰ in J/(mol K), V⁰ in
cm³/mol, and the products of their congruent dissolution.
"""
function zeolite_records()
    return reduce(vcat, _zeolite_records(key) for key in _ZEOLITE_SOURCES)
end

function _zeolite_records(key)
    d = JSON.parsefile(datapath("literature", key * ".json"); dicttype = Dict{String, Any})
    t = d["tables"]["zeolites"]
    c = Dict(name => i for (i, name) in enumerate(t["columns"]))
    products = Dict{String, Dict{String, Float64}}()
    for (sym, sp, ν) in d["tables"]["dissolution_products"]["rows"]
        get!(products, sym, Dict{String, Float64}())[sp] = Float64(ν)
    end
    num(r, name) = Float64(r[c[name]])
    doi = String(d["source"]["doi"])
    return [
        (;
            symbol = String(r[c["symbol"]]), name = String(r[c["name"]]),
            formula = String(r[c["formula"]]),
            logKsp = num(r, "log_Ksp"), logKsp_err = num(r, "log_Ksp_error"),
            ΔfG⁰ = num(r, "dfG"), ΔfH⁰ = num(r, "dfH"), S⁰ = num(r, "S"),
            Cp⁰ = num(r, "Cp"), V⁰ = num(r, "V"),
            products = products[r[c["symbol"]]], doi, origin = String(r[c["S_Cp_origin"]]),
        ) for r in t["rows"]
    ]
end

function _check_zeolite_balance(z)
    lhs = Dict{Symbol, Float64}(k => Float64(v) for (k, v) in parse_formula(z.formula))
    rhs = Dict{Symbol, Float64}()
    charge = 0.0
    for (sp, ν) in z.products
        for (el, n) in parse_formula(_AQUEOUS_FORMULA[sp])
            el === :Zz && continue
            rhs[el] = get(rhs, el, 0.0) + ν * n
        end
        charge += ν * _AQUEOUS_CHARGE[sp]
    end
    for el in union(keys(lhs), keys(rhs))
        el === :Zz && continue
        abs(get(lhs, el, 0.0) - get(rhs, el, 0.0)) > 1.0e-9 && error(
            "$(z.symbol): its dissolution is not balanced in $el " *
                "(solid $(get(lhs, el, 0.0)), products $(get(rhs, el, 0.0)))",
        )
    end
    abs(charge) > 1.0e-9 && error("$(z.symbol): its dissolution products carry a net charge of $charge")
    return nothing
end

"""
    zeolite_logK_check(base_path) -> Vector{NamedTuple}

For every zeolite of [`zeolite_records`](@ref), the published log Ksp at 25 °C
and the one recomputed from its published ΔfG⁰ through the aqueous Gibbs
energies of the base database. Their agreement, within 0.05 log units, is what
shows the two datasets share a reference state; the build refuses otherwise.
"""
function zeolite_logK_check(base_path)
    db = JSON.parsefile(base_path; dicttype = Dict{String, Any})
    G = Dict{String, Float64}()
    for s in db["substances"]
        v = get(get(s, "sm_gibbs_energy", Dict()), "values", nothing)
        v === nothing || isempty(v) || (G[String(s["symbol"])] = Float64(v[1]))
    end
    RTln10 = R_GAS * 298.15 * log(10)
    return map(zeolite_records()) do z
        ΔrG = sum(ν * G[sp] for (sp, ν) in z.products) - z.ΔfG⁰ * 1000
        recomputed = -ΔrG / RTln10
        (; symbol = z.symbol, published = z.logKsp, recomputed, difference = recomputed - z.logKsp)
    end
end

function _zeolite_entry(z)
    return Dict{String, Any}(
        "name" => z.name,
        "symbol" => z.symbol,
        "formula" => z.formula,
        "formula_charge" => 0,
        "aggregate_state" => Dict("3" => "AS_CRYSTAL"),
        "class_" => Dict("0" => "SC_COMPONENT"),
        "Tst" => 298.15,
        "Pst" => 100000,
        "sm_gibbs_energy" => Dict("values" => [z.ΔfG⁰ * 1000]),
        "sm_enthalpy" => Dict("values" => [z.ΔfH⁰ * 1000]),
        "sm_entropy_abs" => Dict("values" => [z.S⁰]),
        "sm_heat_capacity_p" => Dict("values" => [z.Cp⁰]),
        "sm_volume" => Dict("values" => [z.V⁰ * _CM3_PER_MOL_TO_J_PER_BAR]),
        "datasources" => [z.doi],
        # The provenance stays on the entry, so that a phase is never separated
        # from the paper it came from; the dissolution products are the ones the
        # check above used, transcribed from the paper rather than inferred.
        "zeolite_provenance" => Dict(
            "doi" => z.doi,
            "log_Ksp_298K" => z.logKsp,
            "log_Ksp_error" => z.logKsp_err,
            "S_Cp_origin" => z.origin,
            "transcribed_from" => "published table, verbatim",
            "dissolution_products" => z.products,
        ),
    )
end

function _build_zeolites(base_path, out_path)
    db = JSON.parsefile(base_path; dicttype = Dict{String, Any})
    existing = Set(String(s["symbol"]) for s in db["substances"])
    n_base = length(db["substances"])
    records = zeolite_records()
    for z in records
        z.symbol in existing && error("$(z.symbol) would overwrite a Cemdata18 substance; nothing built.")
        _check_zeolite_balance(z)
    end
    check = zeolite_logK_check(base_path)
    worst = maximum(abs(c.difference) for c in check)
    worst > _ZEOLITE_LOGK_TOLERANCE && error(
        "a zeolite's log Ksp recomputed through the base database differs from the published one " *
            "by $(round(worst; digits = 3)) log units (tolerance $(_ZEOLITE_LOGK_TOLERANCE)): the two " *
            "datasets do not share a reference state, or a value was mistranscribed. Nothing built.",
    )
    append!(db["substances"], _zeolite_entry.(records))
    db["thermodataset"] = "cemdata18-zeolites"
    db["zeolite_extension"] = Dict(
        "base" => "cemdata18 (Lothenbach et al. 2019, doi:10.1016/j.cemconres.2018.04.018)",
        "base_substances" => n_base,
        "added_substances" => length(records),
        "sources" => unique(z.doi for z in records),
        "verification" => Dict(
            "method" => "log Ksp recomputed from the aqueous Gibbs energies of the base",
            "tolerance_log_units" => _ZEOLITE_LOGK_TOLERANCE,
            "worst_discrepancy_log_units" => worst,
        ),
        "note" => "The base entries are copied unchanged; nothing is overwritten.",
    )
    return _write_json_atomically(out_path, db)
end

# ── the chloride end member of CSHQ ──────────────────────────────────────────

"""
    cshq_chloride_entry(db, δ; provenance = nothing) -> Dict

The ThermoFun record of `CSHQ-Cl` = (CaCl2)0.5, the chloride end member of CSHQ,
built on the pattern of the `NaSiOH` record of `db` (a parsed Cemdata18
database). `δ` (J/mol) is the Gibbs energy of forming it from ½ Ca²⁺ and Cl⁻ at
293.15 K and 1 bar, the one parameter fitted by `data/chloride/regenerate.jl`.
Its entropy and volume are those of the ions it is made of, so that its
dissolution has neither; its heat capacity is zero.
"""
function cshq_chloride_entry(db, δ; provenance = nothing)
    ext = _chloride_extension()
    T_ref = Float64(ext["reference_temperature_K"])
    template = only(s for s in db["substances"] if s["symbol"] == ext["template"])
    template["Tst"] == T_ref ||
        error("$(ext["template"]) is no longer referred to $T_ref K; the chloride end member must be redone.")
    composition = [(String(k), Float64(v)) for (k, v) in ext["composition"]]
    # Through the default JSON objects, which the species reader expects, whatever
    # dictionary type `db` was parsed with.
    records = JSON.parse(JSON.json([s for s in db["substances"] if s["symbol"] in first.(composition)]))
    byname = Dict(symbol(s) => s for s in build_species(DataFrame(Tables.dictrowtable(records))))
    at_ref(s, key, unit) = ustrip(unit, s[key](T = T_ref * u"K", P = 1.0e5u"Pa"; unit = true))
    along(key, unit) = sum(ν * at_ref(byname[s], key, unit) for (s, ν) in composition)
    G = along(:ΔₐG⁰, us"J/mol") + δ
    S = along(:S⁰, us"J/(mol*K)")
    H = along(:ΔₐH⁰, us"J/mol") + δ
    e = JSON.parse(JSON.json(template))
    # The template's molar mass is NaSiOH's, not this member's; the package
    # computes molar masses from formulas, so the field is dropped, not rewritten.
    haskey(e, "mass_per_mole") && delete!(e, "mass_per_mole")
    e["name"] = ext["name"]
    e["symbol"] = ext["symbol"]
    e["formula"] = ext["formula"]
    e["sm_gibbs_energy"]["values"] = [G]
    e["sm_enthalpy"]["values"] = [H]
    e["sm_entropy_abs"]["values"] = [S]
    e["sm_heat_capacity_p"]["values"] = [0.0]
    e["sm_volume"]["values"] = [along(:V⁰, us"J/(bar*mol)")]
    for m in e["TPMethods"]
        haskey(m, "m_heat_capacity_ft_coeffs") || continue
        m["m_heat_capacity_ft_coeffs"]["values"] = zero.(m["m_heat_capacity_ft_coeffs"]["values"])
    end
    e["datasources"] = ["Hirao2005: fitted by ChemistryLab, data/chloride/regenerate.jl"]
    provenance === nothing || (e["chloride_provenance"] = provenance)
    return e
end

_chloride_extension() = JSON.parsefile(datapath("chloride", "cshq_cl.json"); dicttype = Dict{String, Any})

function _build_chloride(base_path, out_path)
    db = JSON.parsefile(base_path; dicttype = Dict{String, Any})
    ext = _chloride_extension()
    sym = ext["symbol"]
    any(s -> s["symbol"] == sym, db["substances"]) && error("$sym would overwrite a Cemdata18 substance; nothing built.")
    push!(db["substances"], cshq_chloride_entry(db, Float64(ext["delta_J_per_mol"]); provenance = ext["provenance"]))
    db["thermodataset"] = "cemdata18-chloride"
    db["chloride_extension"] = Dict(
        "base" => "cemdata18 (Lothenbach et al. 2019, doi:10.1016/j.cemconres.2018.04.018)",
        "added_substances" => [sym],
        "solid_solution" => "CSHQ_Cl in data/solid_solutions.toml",
        "note" => "The base entries are copied unchanged; nothing is overwritten. The added end member is fitted, not measured: see its chloride_provenance.",
    )
    return _write_json_atomically(out_path, db)
end

"""
    DERIVED_DATABASES

The databases ChemistryLab builds from a published one, by file name: the
Cemdata18 zeolite extension and the Cemdata18 chloride extension.
"""
const DERIVED_DATABASES = Dict(
    d.name => d for d in (
            DerivedDatabase(
                "cemdata18-zeolites.json", "cemdata18-thermofun.json",
                "the 28 zeolites of Ma & Lothenbach (2020, 2021)",
                ["literature/MaLothenbach2020.json", "literature/MaLothenbach2021.json"],
                _build_zeolites,
            ),
            DerivedDatabase(
                "cemdata18-chloride.json", "cemdata18-thermofun.json",
                "the chloride end member of CSHQ fitted on Hirao et al. (2005)",
                ["chloride/cshq_cl.json"],
                _build_chloride,
            ),
        )
)
