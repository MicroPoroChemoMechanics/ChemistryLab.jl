# =============================================================================
#  shi_2016.jl — the carbonated mortars of Shi et al. (2016)
#
#  A white Portland cement (CEM I 52.5 N), alone (P) and with 31.9 % of its mass
#  replaced by limestone (L), metakaolin (M) or both (ML), hydrated 91 days at
#  w/b 0.5 and 20 °C, then carbonated. Their materials (Table 1, Section 2.1),
#  blends (Table 2), degrees of hydration (Table 3) and the CaO of their hydrates,
#  computed and measured (Table 5), are transcribed in
#  data/literature/Shi2016.json. The carbonation page, its test and the generator
#  of the GEMS3K replay include this file.
#
#  The assumptions are the file's, and each is stated where it is made.
# =============================================================================

using ChemistryLab
using DynamicQuantities
using OptimaSolver
using OrderedCollections

include(joinpath(@__DIR__, "validation_common.jl"))

const SHI16 = "Shi2016"
const SHI16_SUBSTANCES = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
const SHI16_DB = Dict(symbol(s) => s for s in SHI16_SUBSTANCES)

shi16_value(name) = ustrip(literature_value(SHI16, name))
shi16_table(name) = literature_table(SHI16, name)

"""The age at which the mortars were carbonated, in days (91; Section 2.3)."""
const SHI16_AGE = ustrip(us"d", literature_value(SHI16, "hydration_time"))

"""The four mortars of Table 2."""
const SHI16_MIXES = ["P", "L", "ML", "M"]

# The phases the pastes may form: the phase list of a Portland paste, whose
# products hold what the metakaolin adds (straetlingite) and what carbonation
# forms (calcite, gypsum, amorphous silica, aluminum and iron hydroxides).
const SHI16_PHASES = "Portland paste (Lothenbach and Winnefeld 2006)"

"""The chemical system of the phase list."""
shi16_system() = phase_list_system(SHI16_PHASES, SHI16_SUBSTANCES)

"""
    shi16_extent(mix, phase) -> ConstantExtent

The degree of reaction of `phase` (`"alite"`, `"belite"`, `"MK"`) in `mix` when
the mortars were carbonated, at 91 days. Table 3 gives it at 28 and 180 days
only, and Section 2.2 reads from the table that little hydration takes place
after 91 days: the value at 180 days is taken. Interpolated linearly between
the two instead, the portlandite of the cement alone comes out 9 % below the
one measured.
"""
function shi16_extent(mix, phase)
    t = shi16_table("degree_of_hydration")
    rows = [k for k in eachindex(t.mix) if t.mix[k] == mix && t.phase[k] == phase]
    isempty(rows) && throw(ArgumentError("Table 3 gives no degree of hydration of $phase in $mix."))
    k = rows[argmax(ustrip.(t.time_d[rows]))]
    return ConstantExtent(ustrip(t.percent[k]) / 100)
end

"""
    shi16_cement(mix) -> Material

The white Portland cement, its alite and belite reacting as Table 3 measured in
`mix`. The paper gives no degree of reaction for the aluminate: it is taken as
reacted by 91 days, as the gypsum and the free lime are. The minor oxides
(magnesia, alkalis, the sulfate beyond the gypsum) are released with the clinker,
at the mean extent of its three phases weighted by their contents.
"""
function shi16_cement(mix)
    wpc = material_template("white Portland cement (Shi 2016)", SHI16_DB)
    alite, belite = shi16_extent(mix, "alite"), shi16_extent(mix, "belite")
    w = Dict(c.name => c.mass_fraction for c in wpc.constituents)
    clinker = w["C3S"] + w["C2S"] + w["C3A"]
    minors(t) = (w["C3S"] * extent(alite, t) + w["C2S"] * extent(belite, t) + w["C3A"]) / clinker
    return with_extents(wpc, Dict("C3S" => alite, "C2S" => belite, "minor oxides" => minors))
end

"""The metakaolin, reacting as Table 3 measured it in `mix`."""
shi16_metakaolin(mix) = with_extents(material_template("metakaolin (Shi 2016)", SHI16_DB), Dict{String, Any}();
                                     material_extent = shi16_extent(mix, "MK"))

"""The limestone: its calcite available, the rest of its analysis inert."""
shi16_limestone() = with_extents(material_template("limestone (Shi 2016)", SHI16_DB), Dict("minor oxides" => 0.0))

"""
    shi16_recipe(mix) -> Recipe

The paste of the mortar `mix`: its binder in the proportions of Table 2, at
w/b 0.5 and 20 °C. The sand is left out, as inert quartz.
"""
function shi16_recipe(mix)
    m = shi16_table("mixes")
    i = findfirst(==(mix), m.mix)
    wpc, mk, ls = ustrip(m.wpc[i]) / 100, ustrip(m.mk[i]) / 100, ustrip(m.ls[i]) / 100
    parts = Pair{Material, Float64}[shi16_cement(mix) => wpc]
    mk > 0 && push!(parts, shi16_metakaolin(mix) => mk)
    ls > 0 && push!(parts, shi16_limestone() => ls)
    return Recipe(parts...; w_b = shi16_value("water_binder_ratio"), T = shi16_value("curing_temperature") * u"K")
end

"""
    shi16_ignited_mortar(mix) -> Float64

The mass (g) of ignited mortar per 100 g of binder: the sand, three times the
binder (Section 2.2), and the binder less its loss on ignition (Table 1). Table 5
counts the CaO per 100 g of it.
"""
function shi16_ignited_mortar(mix)
    m = shi16_table("mixes")
    i = findfirst(==(mix), m.mix)
    c = shi16_table("chemical_composition")
    loi(mat) = ustrip(only(c.percent[k] for k in eachindex(c.material) if c.material[k] == mat && c.oxide[k] == "LOI"))
    binder_loi = (ustrip(m.wpc[i]) * loi("wPc") + ustrip(m.mk[i]) * loi("MK") + ustrip(m.ls[i]) * loi("LS")) / 100
    return 100 * shi16_value("sand_binder_ratio") + 100 - binder_loi
end

"""What Table 5 gives for `phase` of `mix` (mol CaO per 100 g of ignited mortar), computed by the authors."""
function shi16_computed(mix, phase)
    t = shi16_table("cao_computed")
    return ustrip(only(t.cao[k] for k in eachindex(t.mix) if t.mix[k] == mix && t.phase[k] == phase))
end

"""What Table 5 gives for `quantity` of `mix` (mol CaO per 100 g of ignited mortar), measured by TGA."""
function shi16_measured(mix, quantity)
    t = shi16_table("cao_measured")
    return ustrip(only(t.cao[k] for k in eachindex(t.mix) if t.mix[k] == mix && t.quantity[k] == quantity))
end

"""
    shi16_hydrate_cao(rs) -> OrderedDict

The calcium of each hydrate of the paste `rs`, in mol per 100 g of binder: what
Table 5 calls the CO2 binding capacity, the calcium a full carbonation can turn
into calcium carbonate. `"C-S-H"` is the gel's; `"total"` their sum. The calcite
and the unreacted binder are left out, as the table leaves them out.
"""
function shi16_hydrate_cao(rs)
    cs = rs.state.system
    n = ustrip.(us"mol", rs.state.n)
    scale = 100 / rs.recipe.binder_mass
    out = OrderedDict{String, Float64}()
    gel = Set(symbol.(end_members(only(p for p in cs.solid_solutions if name(p) == "CSHQ"))))
    for i in cs.idx_crystal
        sym = symbol(cs.species[i])
        sym in ("Cal", "Gp", "Anh", "hemihydrate", "C3S", "C2S", "C3A", "C4AF") && continue
        ca = get(atoms(cs.species[i]), :Ca, 0) * n[i] * scale
        ca > 0 || continue
        key = sym in gel ? "C-S-H" : sym
        out[key] = get(out, key, 0.0) + ca
    end
    out["total"] = sum(values(out); init = 0.0)
    return out
end

"""
The amounts of carbon dioxide added to each paste, in g per 100 g of binder: the
range of the authors' calculations (Section 2.4.5), in steps of 5 g.
"""
shi16_co2_grams() = collect(0.0:5.0:shi16_value("modeled_co2_max"))

"""
    shi16_carbonated_budget(cs, b, grams) -> Vector{Float64}

The budget `b` of a paste of 100 g of binder with `grams` of carbon dioxide added,
as dissolved `CO2@`: what [`carbonate`](@ref) adds at each step.
"""
function shi16_carbonated_budget(cs, b, grams)
    j = findfirst(s -> symbol(s) == "CO2@", cs.species)
    M = ustrip(us"g/mol", cs.species[j][:M])
    return b .+ (grams / M) .* Float64.(cs.SM.A[:, j])
end
