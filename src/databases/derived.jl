# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using SHA: sha256

# ── Databases built by ChemistryLab from a published one ─────────────────────
#
# Three databases extend Cemdata18 with phases it does not carry: 28 zeolites
# (Ma & Lothenbach 2020, 2021); a chloride end member of CSHQ fitted on the
# sorption tests of Hirao et al. (2005), with the Fe-Friedel's salt that the
# Cemdata18 paper tabulates and its ThermoFun export lacks; and the CASH+ model of
# C-S-H. Only what ChemistryLab adds lives in the package — the published data in
# `data/literature/`, the fitted parameter in `data/chloride/cshq_cl.json`. The
# database itself is assembled on first use from the downloaded Cemdata18 file,
# whose entries are copied through unchanged, and cached under a key made of the
# SHA-256 of everything it was built from, so a new release of either side
# rebuilds it.

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
const _DERIVED_VERSION = "3"

function _derived_key(d::DerivedDatabase, base_path)
    parts = [_DERIVED_VERSION, d.name, _digest(base_path)]
    append!(parts, (_digest(joinpath(pkgdir(@__MODULE__), "data", i)) for i in d.inputs))
    return bytes2hex(sha256(join(parts, "\n")))
end

# Under the lock of `remote.jl`: two threads asking for the same derived
# database build it once, not twice into the same file.
function _derived_path(d::DerivedDatabase; download::Bool = true)
    return lock(_DATABASE_LOCK) do
        base = database_path(d.base; download)
        key = _derived_key(d, base)
        dir = joinpath(database_cache(), "derived")
        out = joinpath(dir, d.name)
        keyfile = out * ".key"
        isfile(out) && isfile(keyfile) && strip(read(keyfile, String)) == key && return out
        mkpath(dir)
        tmp = _partial_path(out)
        try
            d.build(base, tmp)
            mv(tmp, out; force = true)
        finally
            rm(tmp; force = true)
        end
        write(keyfile, key)
        out
    end
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

The 28 zeolites of [MaLothenbach2020, MaLothenbach2021](@citet) (Na and K series) as the two
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

# ── Fe-Friedel's salt ────────────────────────────────────────────────────────
#
# Table 1 of Cemdata18 lists Fe-Friedel's salt and Table 2 its solubility, but
# the ThermoFun export of the database does not carry it. Its record is written
# from the table, on the pattern of Friedel's salt, and the log Ks0 of Table 2 is
# recomputed through the base's aqueous Gibbs energies before it is added: the
# same check as for the zeolites, with the same tolerance.

const _FE_FRIEDEL = "C4FCl2H10"
const _FE_FRIEDEL_FORMULA = "Ca4Fe|3|2Cl2(OH)12(H2O)4"

"""
    fe_friedel_entry(db) -> Dict

The ThermoFun record of Fe-Friedel's salt, Ca₄Fe₂Cl₂(OH)₁₂·4H₂O, from Table 1 of
[Lothenbach2019](@cite) as transcribed in `data/literature/Lothenbach2019.json`
(table `solid_standard_properties`), on the pattern of the `C4AClH10` record of
`db`, a parsed Cemdata18 database: ΔfG⁰, ΔfH⁰, S⁰ and V⁰ at 298.15 K and 1 bar,
and the heat capacity `Cp = a₀ + a₁T + a₂T⁻² + a₃T^(-1/2)`.
"""
function fe_friedel_entry(db)
    r = literature_row("Lothenbach2019", "solid_standard_properties", _FE_FRIEDEL)
    template = only(s for s in db["substances"] if s["symbol"] == "C4AClH10")
    template["Tst"] == 298.15 || error("C4AClH10 is no longer referred to 298.15 K; Fe-Friedel's salt must be redone.")
    # The four coefficients are in SI units already: their values are the table's.
    a = Float64[ustrip(r.a0), ustrip(r.a1), ustrip(r.a2), ustrip(r.a3)]
    T = 298.15
    e = JSON.parse(JSON.json(template))
    # As for the chloride end member, the package computes molar masses from
    # formulas, and the template's is Friedel's salt's.
    haskey(e, "mass_per_mole") && delete!(e, "mass_per_mole")
    e["name"] = "C4FCl2H10 Fe-Friedel's salt"
    e["symbol"] = _FE_FRIEDEL
    e["formula"] = _FE_FRIEDEL_FORMULA
    e["sm_gibbs_energy"]["values"] = [ustrip(us"J/mol", r.dfG)]
    e["sm_enthalpy"]["values"] = [ustrip(us"J/mol", r.dfH)]
    e["sm_entropy_abs"]["values"] = [ustrip(us"J/(mol*K)", r.S)]
    e["sm_heat_capacity_p"]["values"] = [a[1] + a[2] * T + a[3] / T^2 + a[4] / sqrt(T)]
    e["sm_volume"]["values"] = [ustrip(us"cm^3/mol", r.V) * _CM3_PER_MOL_TO_J_PER_BAR]
    for m in e["TPMethods"]
        haskey(m, "m_heat_capacity_ft_coeffs") || continue
        c = m["m_heat_capacity_ft_coeffs"]["values"]
        c .= 0
        c[1:4] .= a
    end
    e["datasources"] = ["Lothenbach2019: Table 1, after Dilnesa (2012)"]
    e["literature_provenance"] = Dict(
        "doi" => literature("Lothenbach2019").source["doi"],
        "table" => "Table 1, row C4FCl2H10; log Ks0 of Table 2",
        "transcribed_from" => "data/literature/Lothenbach2019.json, table solid_standard_properties",
    )
    return e
end

"""
    fe_friedel_logK_check(base_path) -> NamedTuple

The log Ks0 of Fe-Friedel's salt at 25 °C printed in Table 2 of
[Lothenbach2019](@cite) and the one recomputed from its ΔfG⁰ of Table 1 through
the aqueous Gibbs energies of the base database, over the products of
`data/literature/Lothenbach2019.json` (table `dissolution_products`), the water
recovered from the hydrogen balance.
"""
function fe_friedel_logK_check(base_path)
    db = JSON.parsefile(base_path; dicttype = Dict{String, Any})
    G = Dict(String(s["symbol"]) => Float64(s["sm_gibbs_energy"]["values"][1]) for s in db["substances"])
    products = literature_table("Lothenbach2019", "dissolution_products"; phase = _FE_FRIEDEL)
    ν = Dict(zip(products.species, products.coefficient))
    # Of the products, only OH⁻ carries hydrogen.
    ν["H2O@"] = (parse_formula(_FE_FRIEDEL_FORMULA)[:H] - ν["OH-"]) / 2
    t2 = literature_table("Lothenbach2019", "solubility_products")
    published = t2.log_Ks0[only(findall(==("Fe-Friedel's salt"), t2.printed))]
    ΔfG = ustrip(us"J/mol", literature_row("Lothenbach2019", "solid_standard_properties", _FE_FRIEDEL).dfG)
    recomputed = -(sum(n * G[sp] for (sp, n) in ν) - ΔfG) / (R_GAS * 298.15 * log(10))
    return (; published, recomputed, difference = recomputed - published)
end

function _build_chloride(base_path, out_path)
    db = JSON.parsefile(base_path; dicttype = Dict{String, Any})
    ext = _chloride_extension()
    sym = ext["symbol"]
    for s in (sym, _FE_FRIEDEL)
        any(x -> x["symbol"] == s, db["substances"]) && error("$s would overwrite a Cemdata18 substance; nothing built.")
    end
    check = fe_friedel_logK_check(base_path)
    abs(check.difference) > _ZEOLITE_LOGK_TOLERANCE && error(
        "the log Ks0 of Fe-Friedel's salt recomputed through the base database differs from the " *
            "published one by $(round(check.difference; digits = 3)) log units (tolerance " *
            "$(_ZEOLITE_LOGK_TOLERANCE)): a value was mistranscribed. Nothing built.",
    )
    push!(db["substances"], cshq_chloride_entry(db, Float64(ext["delta_J_per_mol"]); provenance = ext["provenance"]))
    push!(db["substances"], fe_friedel_entry(db))
    db["thermodataset"] = "cemdata18-chloride"
    db["chloride_extension"] = Dict(
        "base" => "cemdata18 (Lothenbach et al. 2019, doi:10.1016/j.cemconres.2018.04.018)",
        "added_substances" => [sym, _FE_FRIEDEL],
        "solid_solution" => "CSHQ_Cl and Friedel_AlFe in data/solid_solutions.toml",
        "verification" => Dict(
            "method" => "log Ks0 of Fe-Friedel's salt recomputed from the aqueous Gibbs energies of the base",
            "tolerance_log_units" => _ZEOLITE_LOGK_TOLERANCE,
            "discrepancy_log_units" => check.difference,
        ),
        "note" => "The base entries are copied unchanged; nothing is overwritten. The CSHQ " *
            "end member is fitted, not measured: see its chloride_provenance. Fe-Friedel's " *
            "salt is the one of Table 1 of the Cemdata18 paper, which its ThermoFun export lacks.",
    )
    return _write_json_atomically(out_path, db)
end

# ── the CASH+ core model of C-S-H ────────────────────────────────────────────

"""
    cashplus_entries() -> Vector{Dict}

The ThermoFun records of the thirty-three end-members of the CASH+ model of C-S-H
with its extension to the alkali and alkaline-earth metals, the heat capacity
taken as constant:

  - the six of the core model, from `data/literature/Kulik2022.json`: G° and H°
    of its Table 8, S°, Cp° and V° of its Table 4, except the H° of TSvh, which
    is the one [Miron2022a](@citet) reprint (Table 8 transposes two digits);
  - the six with sodium or potassium, from `data/literature/Miron2022a.json`
    (Table A1), with the G° and H° of TCNh and TCKh fine-tuned by
    [Miron2022b; Table 5](@citet) for cement pore solutions;
  - the twenty-one with Li, Rb, Cs, Mg, Sr, Ba or Ra, from the same Table A1.
"""
function cashplus_entries()
    core = literature_table("Kulik2022", "cashplus_standard_properties")
    alkali = literature_table("Miron2022a", "alkali_standard_properties")
    tuned = literature_table("Miron2022b", "fine_tuned_end_members")
    formulas = literature_table("Miron2022a", "cashplus_full_end_members")
    doi(key) = literature(key).source["doi"]
    rows = Dict{String, Any}[]
    function record(name, G, H, S, Cp, V, sources)
        k = findfirst(==(name), formulas.end_member)
        return Dict{String, Any}(
            "name" => "CASH+ end-member $name",
            "symbol" => name,
            "formula" => formulas.formula[k],
            "formula_charge" => 0,
            "aggregate_state" => Dict("3" => "AS_CRYSTAL"),
            "class_" => Dict("0" => "SC_COMPONENT"),
            "Tst" => 298.15,
            "Pst" => 100000,
            "sm_gibbs_energy" => Dict("values" => [ustrip(us"J/mol", G)]),
            "sm_enthalpy" => Dict("values" => [ustrip(us"J/mol", H)]),
            "sm_entropy_abs" => Dict("values" => [ustrip(us"J/(mol*K)", S)]),
            "sm_heat_capacity_p" => Dict("values" => [ustrip(us"J/(mol*K)", Cp)]),
            "sm_volume" => Dict("values" => [ustrip(us"cm^3/mol", V) * _CM3_PER_MOL_TO_J_PER_BAR]),
            "datasources" => doi.(sources),
        )
    end
    for r in eachindex(core.end_member)
        name = core.end_member[r]
        H, src = name == "TSvh" ? (literature_value("Miron2022a", "TSvh_H"), ["Kulik2022", "Miron2022a"]) :
            (core.H[r], ["Kulik2022"])
        push!(rows, record(name, core.G[r], H, core.S[r], core.Cp[r], core.V[r], src))
    end
    for r in eachindex(alkali.end_member)
        name = alkali.end_member[r]
        t = findfirst(==(name), tuned.end_member)
        G, H, src = t === nothing ? (alkali.G[r], alkali.H[r], ["Miron2022a"]) :
            (tuned.G[t], tuned.H[t], ["Miron2022a", "Miron2022b"])
        push!(rows, record(name, G, H, alkali.S[r], alkali.Cp[r], alkali.V[r], src))
    end
    ext = literature_table("Miron2022a", "extension_standard_properties")
    for r in eachindex(ext.end_member)
        push!(rows, record(ext.end_member[r], ext.G[r], ext.H[r], ext.S[r], ext.Cp[r], ext.V[r], ["Miron2022a"]))
    end
    return rows
end

# The aqueous species the CASH+ extension needs and Cemdata18 lacks, with the
# standard properties Miron et al. (2022a) give them (their Table A1): the
# cations of lithium, rubidium, cesium, barium and radium, and the Ca(OH)2@
# complex they derived and kept when fitting the extension. Built on the pattern
# of CaSiO3@, a solute of constant heat capacity.
const _CASHPLUS_IONS = ("Li+" => "Li+", "Rb+" => "Rb+", "Cs+" => "Cs+", "Ba2+" => "Ba+2", "Ra2+" => "Ra+2")

function _cashplus_aqueous_entries(db)
    template = only(s for s in db["substances"] if s["symbol"] == "CaSiO3@")
    doi = literature("Miron2022a").source["doi"]
    function entry(symbol, formula, charge, G, H, S, Cp, V)
        e = JSON.parse(JSON.json(template))
        haskey(e, "mass_per_mole") && delete!(e, "mass_per_mole")
        e["name"], e["symbol"], e["formula"], e["formula_charge"] = symbol, symbol, formula, charge
        e["sm_gibbs_energy"]["values"] = [ustrip(us"J/mol", G)]
        e["sm_enthalpy"]["values"] = [ustrip(us"J/mol", H)]
        e["sm_entropy_abs"]["values"] = [ustrip(us"J/(mol*K)", S)]
        e["sm_heat_capacity_p"]["values"] = [ustrip(us"J/(mol*K)", Cp)]
        e["sm_volume"]["values"] = [ustrip(us"cm^3/mol", V) * _CM3_PER_MOL_TO_J_PER_BAR]
        for m in e["TPMethods"]
            haskey(m, "m_heat_capacity_ft_coeffs") || continue
            m["m_heat_capacity_ft_coeffs"]["values"] = vcat(ustrip(us"J/(mol*K)", Cp), zeros(length(m["m_heat_capacity_ft_coeffs"]["values"]) - 1))
        end
        e["datasources"] = [doi]
        return e
    end
    ms = literature_table("Miron2022a", "master_species")
    aq = literature_table("Miron2022a", "aqueous_species")
    out = Dict{String, Any}[]
    for (printed, symbol) in _CASHPLUS_IONS
        i = findfirst(==(printed), ms.species)
        push!(out, entry(symbol, symbol, occursin("2", printed) ? 2 : 1, ms.G[i], ms.H[i], ms.S[i], ms.Cp[i], ms.V[i]))
    end
    for i in eachindex(aq.species)
        push!(out, entry(aq.species[i], aq.formula[i] * "@", 0, aq.G[i], aq.H[i], aq.S[i], aq.Cp[i], aq.V[i]))
    end
    return out
end

# The aqueous complex CaSiO3@ refitted together with the core model (Kulik et al.
# 2022, Table 9, accepted variant): the model and this complex set the Si of a
# C-S-H solution together, and neither holds with the other's Cemdata18 value.
function _cashplus_casio3!(entry)
    q(name, unit) = ustrip(unit, literature_value("Kulik2022", name))
    entry["sm_gibbs_energy"]["values"] = [q("CaSiO3_aq_G", us"J/mol")]
    entry["sm_enthalpy"]["values"] = [q("CaSiO3_aq_H", us"J/mol")]
    entry["sm_entropy_abs"]["values"] = [q("CaSiO3_aq_S", us"J/(mol*K)")]
    Cp = q("CaSiO3_aq_Cp", us"J/(mol*K)")
    entry["sm_heat_capacity_p"]["values"] = [Cp]
    entry["sm_volume"]["values"] = [q("CaSiO3_aq_V", us"cm^3/mol") * _CM3_PER_MOL_TO_J_PER_BAR]
    for m in entry["TPMethods"]
        haskey(m, "m_heat_capacity_ft_coeffs") || continue
        c = m["m_heat_capacity_ft_coeffs"]["values"]
        all(iszero, c[2:end]) || error("CaSiO3@ no longer has a constant heat capacity in the base; nothing built.")
        c[1] = Cp
    end
    entry["datasources"] = [literature("Kulik2022").source["doi"]]
    return entry
end

function _build_cashplus(base_path, out_path)
    db = JSON.parsefile(base_path; dicttype = Dict{String, Any})
    added = cashplus_entries()
    for e in added
        any(s -> s["symbol"] == e["symbol"], db["substances"]) &&
            error("$(e["symbol"]) would overwrite a Cemdata18 substance; nothing built.")
    end
    _cashplus_casio3!(only(s for s in db["substances"] if s["symbol"] == "CaSiO3@"))
    aqueous = _cashplus_aqueous_entries(db)
    for e in aqueous
        any(s -> s["symbol"] == e["symbol"], db["substances"]) &&
            error("$(e["symbol"]) would overwrite a Cemdata18 substance; nothing built.")
    end
    append!(added, aqueous)
    append!(db["substances"], added)
    db["thermodataset"] = "cemdata18-cashplus"
    db["cashplus_extension"] = Dict(
        "base" => "cemdata18 (Lothenbach et al. 2019, doi:10.1016/j.cemconres.2018.04.018)",
        "added_substances" => [e["symbol"] for e in added],
        "replaced_substances" => ["CaSiO3@"],
        "sources" => [
            "Kulik, Miron & Lothenbach (2022), doi:$(literature("Kulik2022").source["doi"])",
            "Miron, Kulik, Yan, Tits & Lothenbach (2022), doi:$(literature("Miron2022a").source["doi"])",
            "Miron, Kulik & Lothenbach (2022), doi:$(literature("Miron2022b").source["doi"])",
        ],
        "solid_solutions" => "CASH+, CASH+NK and CASH+ext in data/solid_solutions.toml",
        "note" => "The base entries are copied unchanged except CaSiO3@, whose standard " *
            "properties are those the core model was fitted with (Table 9, accepted variant). " *
            "The aqueous species added (Li+, Rb+, Cs+, Ba+2, Ra+2, Ca(OH)2@) carry the " *
            "properties of Miron et al. (2022a, Table A1), with a constant heat capacity.",
    )
    return _write_json_atomically(out_path, db)
end

"""
    DERIVED_DATABASES

The databases ChemistryLab builds from a published one, by file name: the
Cemdata18 zeolite extension, the Cemdata18 chloride extension, and Cemdata18 with
the CASH+ model of C-S-H and its extension to the alkali and alkaline-earth
metals.
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
                "the chloride end member of CSHQ fitted on Hirao et al. (2005), and Fe-Friedel's salt of the Cemdata18 paper",
                ["chloride/cshq_cl.json", "literature/Lothenbach2019.json"],
                _build_chloride,
            ),
            DerivedDatabase(
                "cemdata18-cashplus.json", "cemdata18-thermofun.json",
                "the CASH+ model of C-S-H of Kulik et al. (2022), with the alkali and alkaline-earth end-members of Miron et al. (2022a, b)",
                ["literature/Kulik2022.json", "literature/Miron2022a.json", "literature/Miron2022b.json"],
                _build_cashplus,
            ),
        )
)
