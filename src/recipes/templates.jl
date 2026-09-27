# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using TOML

# ── Materials of published mixes ──────────────────────────────────────────────

# What a phase named by an XRD-Rietveld table is, in the database or by its
# oxides. The database symbol when Cemdata18 has the phase; the formula,
# entered through its oxides, when it has none (periclase, mullite, dolomite).
# "CaO + Ca(OH)2", free lime and portlandite reported together, is taken as
# lime: in a clinker it is mostly free lime, and both end as portlandite.
const RIETVELD_PHASES = Dict{String, Union{String, Tuple{Symbol, String}}}(
    "C3S" => "C3S", "C2S" => "C2S", "C3A" => "C3A", "C4AF" => "C4AF",
    "Calcite" => "Cal", "Anhydrite" => "Anh", "Gypsum" => "Gp",
    "Arcanite" => "K2SO4", "Quartz" => "Qtz", "Portlandite" => "Portlandite",
    "Periclase" => (:oxides, "MgO"),
    "Mullite" => (:oxides, "Al6Si2O13"),
    "Dolomite" => (:oxides, "CaMgC2O6"),
    "CaO + Ca(OH)2" => (:oxides, "CaO"),
)

# The oxide fractions of a formula, per the cation of each oxide.
function _formula_oxides(f::AbstractString)
    comp = composition(Formula(f))
    M = ustrip(us"g/mol", calculate_molar_mass(comp))
    names = Dict(
        :Ca => "CaO", :Mg => "MgO", :Al => "Al2O3", :Si => "SiO2", :Fe => "Fe2O3",
        :C => "CO2", :S => "SO3", :K => "K2O", :Na => "Na2O", :Ti => "TiO2", :H => "H2O"
    )
    out = OrderedDict{String, Float64}()
    for (el, k) in comp
        el === :O && continue
        haskey(names, el) || throw(ArgumentError("no oxide is defined here for $el in $f"))
        ox = names[el]
        _, per = _oxide_cation(ox)
        out[ox] = Float64(k) / per * _oxide_molar_mass(ox) / M
    end
    return out
end

"""
    material_template(name, species) -> Material
    material_templates() -> Vector{String}

A material of a published mix, built from `data/recipe_templates.toml` and the
literature files it names: every number is read from `data/literature`, and the
file holds only references and choices. `species` maps database symbols to
species (a `Dict`, or a [`ChemicalSystem`](@ref) holding them all).
`material_templates()` lists the names.

An entry is built one of four ways:

  - `phases = "<key>:<table>:<material>"`: from an XRD-Rietveld table, each
    phase a [`MineralConstituent`](@ref) or, when the database has no record of
    it, an [`OxideConstituent`](@ref) (see `ChemistryLab.RIETVELD_PHASES`); the
    `Amorphous` entry, when given, is the glass by difference from `analysis`
    ([`reactive_part`](@ref));
  - `analysis = …` and `route = "bogue"`: Bogue's calculation on the oxide
    analysis ([`bogue`](@ref)), the sulfate carried by `sulfate` (gypsum by
    default), with the oxides no phase takes as one oxide constituent unless
    `minor_oxides = false`;
  - `analysis = …` and `route = "oxides"`: the analysis alone, one oxide
    constituent;
  - `analysis = …` and `crystalline = "<key>:<table>:<material>"`: the glass by
    difference, reacting, beside its crystals, inert.

Every constituent of a template reacts completely unless the entry says
otherwise; the extents of a real mix are the caller's, set with
`with_extents`.
"""
function material_template(name::AbstractString, species)
    e = _template_entry(name)
    kind = Symbol(get(e, "kind", "other"))
    src = get(e, "source", nothing)
    if haskey(e, "phases")
        ph = _literature_phases(e["phases"])
        crystals = OrderedDict(k => v for (k, v) in ph if k != "Amorphous")
        cons = AbstractConstituent[_rietveld_constituent(p, w, species) for (p, w) in crystals]
        if haskey(ph, "Amorphous")
            haskey(e, "analysis") || error("template \"$name\": an amorphous part needs `analysis` to be found by difference")
            glass = _glass(e["analysis"], crystals, species)
            push!(cons, OxideConstituent("glass", glass.oxides; mass_fraction = glass.mass_fraction))
        end
        return Material(name, kind; constituents = cons, source = src)
    elseif get(e, "route", "") == "bogue"
        ox = literature_oxides(_split3(e["analysis"])...)
        bg = bogue(ox, species; sulfate = get(e, "sulfate", "Gp"))
        cons = AbstractConstituent[MineralConstituent(species[p]; mass_fraction = f) for (p, f) in bg.phases if f > 0]
        if get(e, "minor_oxides", true) && !isempty(bg.remainder)
            tot = sum(values(bg.remainder))
            push!(cons, OxideConstituent("minor oxides", OrderedDict(k => v / tot for (k, v) in bg.remainder); mass_fraction = tot))
        end
        return Material(name, kind; constituents = cons, source = src)
    elseif get(e, "route", "") == "oxides"
        return oxide_material(name, _split3(e["analysis"])...; kind)
    elseif haskey(e, "crystalline")
        crystals = _literature_phases(e["crystalline"])
        pop!(crystals, "Amorphous", nothing)
        cons = AbstractConstituent[_rietveld_constituent(p, w, species; extent = 0.0) for (p, w) in crystals]
        glass = _glass(e["analysis"], crystals, species)
        push!(cons, OxideConstituent("glass", glass.oxides; mass_fraction = glass.mass_fraction))
        return Material(name, kind; constituents = cons, source = src)
    end
    error("template \"$name\": give `phases`, or `analysis` with `route` or `crystalline`")
end

material_templates() = [e["name"] for e in _templates()]

"""
    with_extents(material, extents::AbstractDict; material_extent = nothing) -> Material

`material` with the extents of its constituents set by name (`"C3S" => 0.7`,
`"glass" => TabulatedExtent(…)`), the others unchanged, and the extent of the
material as a whole replaced when `material_extent` is given.
"""
function with_extents(m::Material, extents::AbstractDict; material_extent = nothing)
    cons = AbstractConstituent[]
    for c in m.constituents
        e = get(extents, c.name, nothing)
        push!(cons, e === nothing ? c : _with_extent(c, _as_extent(e)))
    end
    unknown = setdiff(String.(keys(extents)), [c.name for c in m.constituents])
    isempty(unknown) || throw(ArgumentError("with_extents: $(m.name) has no constituent $(join(unknown, ", "))."))
    return Material(m.name, m.kind, cons, material_extent === nothing ? m.extent : _as_extent(material_extent), m.source)
end
_with_extent(c::MineralConstituent, e) = MineralConstituent{typeof(c.species), typeof(e)}(c.name, c.species, c.mass_fraction, e)
_with_extent(c::OxideConstituent, e) = OxideConstituent{typeof(e)}(c.name, c.oxides, c.mass_fraction, e, c.density, c.enthalpy, c.source)

function _rietveld_constituent(phase, w, species; extent = 1.0)
    target = get(RIETVELD_PHASES, phase, nothing)
    target === nothing && error("no database phase or formula is known here for the Rietveld phase \"$phase\"")
    if target isa String
        return MineralConstituent(species[target]; mass_fraction = w, extent, name = phase)
    end
    return OxideConstituent(phase, _formula_oxides(last(target)); mass_fraction = w, extent)
end

function _glass(analysis, crystals, species)
    ox = literature_oxides(_split3(analysis)...)
    # Mass balance on the oxides the analysis gives, the crystals entered through
    # the same map as their constituents.
    rest = copy(ox)
    for (ph, w) in crystals
        t = RIETVELD_PHASES[ph]
        c = t isa String ? oxide_content(species[t], keys(rest)) :
            OrderedDict(k => get(_formula_oxides(last(t)), k, 0.0) for k in keys(rest))
        for k in keys(rest)
            rest[k] -= w * c[k]
        end
    end
    frac = 1 - sum(values(crystals); init = 0.0)
    return (; oxides = OrderedDict(k => max(v, 0.0) / frac for (k, v) in rest), mass_fraction = frac)
end

function _literature_phases(ref)
    key, table, material = _split3(ref)
    t = literature_table(key, table)
    out = OrderedDict{String, Float64}()
    for (m, p, v) in zip(t.material, t.phase, t.percent)
        m == material && (out[String(p)] = ustrip(v) / 100)
    end
    isempty(out) && error("$ref names no rows")
    return out
end

function _split3(ref::AbstractString)
    p = split(ref, ":")
    length(p) == 3 || error("expected \"<key>:<table>:<material>\"; got \"$ref\"")
    return String.(p)
end

_templates() = get(TOML.parsefile(resolve_data_path("recipe_templates.toml")), "material", Any[])
function _template_entry(name)
    for e in _templates()
        e["name"] == name && return e
    end
    throw(KeyError("no material template \"$name\"; the templates are: " * join(material_templates(), "; ")))
end
