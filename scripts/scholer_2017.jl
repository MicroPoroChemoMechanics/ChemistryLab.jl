# =============================================================================
#  scholer_2017.jl — the early pore solutions of Schöler et al. (2017)
#
#  Forty-eight pore solutions of a CEM I 52.5 R alone and blended with slag, fly
#  ash, limestone or quartz, analyzed during the first six hours of hydration
#  (their Table 6), with the saturation indices the authors computed from them
#  (their Table 7), transcribed in data/literature/Scholer2017.json. This file
#  speciates each solution at its measured pH and computes the same indices; the
#  validation page and its test include it.
# =============================================================================

using ChemistryLab
using DynamicQuantities
using OptimaSolver

const S17 = "Scholer2017"
s17_table(name) = literature_table(S17, name)

# The solids whose indices Table 7 gives, by their Cemdata18 records, whose
# solubility products are those the authors used (Table 4 and Cemdata07):
# monosulfoaluminate, ettringite, portlandite and gypsum. Monocarboaluminate
# needs a carbonate activity and the table gives no carbon.
const S17_SOLIDS = ["monosulphate12" => "Ms", "ettringite" => "E", "Portlandite" => "CH", "Gp" => "Gp"]

# No redox reaction is followed: the sulfur stays sulfate, and the dissolved
# gases are left out.
const S17_EXCLUDE = ["HS-", "S-2", "H2S@", "SO3-2", "HSO3-", "S2O3-2", "H2@", "O2@"]

"""
    s17_systems() -> (; aqueous, full, substances)

The aqueous species of Cemdata18 for the elements of Table 6 (Al, Ca, S, K, Na,
Si), and the same with the solids of Table 7 beside them, which gives the
composition of each solid in the primaries of the solution.
"""
function s17_systems()
    substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
    # The solids, and the alkalis and silicon, which no solid holds but every
    # solution does.
    sp = speciation(substances, vcat(first.(S17_SOLIDS), ["K+", "Na+", "SiO2@"]); aggregate_state = [AS_AQUEOUS], exclude_species = S17_EXCLUDE)
    aq = [s for s in sp if aggregate_state(s) == AS_AQUEOUS]
    return (;
        aqueous = ChemicalSystem(aq, CEMDATA_PRIMARIES),
        full = ChemicalSystem(sp, CEMDATA_PRIMARIES),
        substances,
    )
end

"""
    s17_solution(cs, system, age; model) -> (; state, certificate)

The pore solution of `system` at `age` (hours) as Table 6 gives it, in a
kilogram of water (the table's mmol per liter taken as mmol per kilogram), at the
measured pH: the concentrations fix the elements, and hydroxide is added or taken
until the pH is the one measured, so that the charge the analysis leaves
unbalanced is carried by it.
"""
function s17_solution(cs, system::AbstractString, age::Real; model)
    ps = s17_table("pore_solution")
    at(k) = ps.system[k] == system && ustrip(us"hr", ps.age[k]) ≈ age
    rows = [k for k in eachindex(ps.system) if at(k)]
    isempty(rows) && throw(ArgumentError("s17_solution: Table 6 gives no solution $system at $age h"))
    # mol per kg of water, the table's mmol per liter taken as mmol per kilogram;
    # an analyte not detected is absent from the table, and held at a trace.
    mol(el) = (k = findfirst(k -> ps.analyte[k] == el, rows); k === nothing ? 1.0e-12 : ustrip(us"mol/m^3", ps.concentration[rows[k]]) / 1000)
    st = ChemicalState(cs)
    Mw = ustrip(us"kg/mol", only(s for s in cs.species if symbol(s) == "H2O@")[:M])
    set_quantity!(st, "H2O@", (1 / Mw)u"mol")
    for (el, sp) in ("Ca" => "Ca+2", "Al" => "AlO2-", "SO4" => "SO4-2", "K" => "K+", "Na" => "Na+", "Si" => "SiO2@")
        set_quantity!(st, sp, mol(el)u"mol")
    end
    b = Float64.(conservation_matrix(cs)) * ustrip.(us"mol", st.n)
    eq, cert = equilibrate_certified(st; model, b, constraint = FixedpH(s17_pH(system, age); titrant = "OH-"))
    return (; state = eq, certificate = cert)
end

"""The pH Table 6 gives for `system` at `age` hours, measured."""
function s17_pH(system, age)
    t = s17_table("pore_solution_pH")
    return Float64(t.pH[findfirst(k -> t.system[k] == system && ustrip(us"hr", t.age[k]) ≈ age, eachindex(t.system))])
end

"""
    s17_indices(eq, systems; model) -> Dict

The saturation indices log10(IAP/K) of the solids of Table 7 in the solution `eq`,
and of the Ca-rich C-S-H of Table 4 (its reaction and log K = -17.2, since Cemdata18
holds no such end-member): the potentials of the primaries of the solution
against the energy of each solid, as `saturation_indices` computes them.
"""
function s17_indices(eq, systems; model)
    cs, full = systems.aqueous, systems.full
    T = 293.15
    RT = ChemistryLab.R_GAS * T
    g(sp) = ustrip(us"J/mol", sp[:ΔₐG⁰](T = T * u"K", P = 1.0e5u"Pa"; unit = true)) / RT
    lna = log_activities(eq, model)
    # The potential of each primary, from the solution.
    y = map(full.SM.primaries) do p
        s = symbol(p)
        k = findfirst(x -> symbol(x) == s, cs.species)
        k === nothing ? 0.0 : g(cs.species[k]) + lna[s]
    end
    A = Float64.(conservation_matrix(full))
    out = Dict{String, Float64}()
    for (record, label) in S17_SOLIDS
        j = findfirst(s -> symbol(s) == record, full.species)
        out[label] = (sum(A[c, j] * y[c] for c in eachindex(y)) - g(full.species[j])) / log(10)
    end
    # The C-S-H of Table 4: (CaO)1.5(SiO2)(H2O)2.5 + H2O = 1.5 Ca+2 + Si(OH)4 + 3 OH-,
    # with Si(OH)4 = SiO2@ + 2 H2O in the database's species.
    la(s) = lna[s] / log(10)
    out["CSH"] = 1.5 * la("Ca+2") + la("SiO2@") + 2 * la("H2O@") + 3 * la("OH-") - la("H2O@") + 17.2
    return out
end

"""The index of `label` in Table 7 for `system` at `age` hours (`nothing` when not given)."""
function s17_published(system, age, label)
    t = s17_table("saturation_indices")
    i = findfirst(k -> t.system[k] == system && ustrip(us"hr", t.age[k]) ≈ age && t.index[k] == label, eachindex(t.system))
    return i === nothing ? nothing : Float64(t.value[i])
end

s17_systems_listed() = unique(s17_table("pore_solution").system)
s17_ages() = sort(unique(ustrip.(us"hr", s17_table("pore_solution").age)))
