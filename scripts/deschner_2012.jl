# =============================================================================
#  deschner_2012.jl — the pore solutions of Deschner et al. (2012)
#
#  The pore solutions of a CEM I 42.5 N alone, and blended with 50 % of one of
#  two siliceous fly ashes, of quartz, or of fly ash and limestone, analyzed from
#  one hour to 550 days (their Figs. 13-16), with the effective saturation
#  indices the authors computed from them (their Table 3), in
#  data/literature/Deschner2012.json. This file speciates each solution at its
#  measured hydroxide and computes the same indices; the validation page and its
#  test include it.
# =============================================================================

using ChemistryLab
using DynamicQuantities
using OptimaSolver

const D12 = "Deschner2012"
d12_table(name) = literature_table(D12, name)

# The solids of Table 3 whose indices the analyses determine, by their Cemdata18
# records, with the number of ions each dissolves into (Ca+2, Al(OH)4-,
# SiO(OH)3-, OH-, SO4-2, water left out), by which the authors divide the index.
# Hemicarbonate and monocarbonate need a carbonate activity, and the pore
# solutions were not analyzed for carbon.
const D12_SOLIDS = [
    "Portlandite" => ("CH", 3),
    "Gp" => ("Gypsum", 2),
    "ettringite" => ("Ettr", 15),
    "straetlingite" => ("Strätl", 6),
    "monosulphate12" => ("MS", 11),
]

# No redox reaction is followed: the sulfur stays sulfate, and the dissolved
# gases are left out.
const D12_EXCLUDE = ["HS-", "S-2", "H2S@", "SO3-2", "HSO3-", "S2O3-2", "H2@", "O2@"]

# The analytes of the figures and the species that carry them in the solution.
const D12_SPECIES = ["Ca" => "Ca+2", "Al" => "AlO2-", "S" => "SO4-2", "K" => "K+", "Na" => "Na+", "Si" => "SiO2@"]

"""
    d12_systems() -> (; aqueous, full, substances)

The aqueous species of Cemdata18 for the elements of the figures (Al, Ca, S, K,
Na, Si), and the same with the solids of Table 3 beside them, which gives the
composition of each solid in the primaries of the solution.
"""
function d12_systems()
    substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
    sp = speciation(
        substances, vcat(first.(D12_SOLIDS), ["K+", "Na+", "SiO2@"]);
        aggregate_state = [AS_AQUEOUS], exclude_species = D12_EXCLUDE,
    )
    aq = [s for s in sp if aggregate_state(s) == AS_AQUEOUS]
    return (;
        aqueous = ChemicalSystem(aq, CEMDATA_PRIMARIES),
        full = ChemicalSystem(sp, CEMDATA_PRIMARIES),
        substances,
    )
end

"""
    d12_concentrations(system, age) -> Dict{String, Float64}

What the figures give for `system` at `age` (hours), in mol per liter, by analyte
(`"Ca"`, `"Al"`, `"S"`, `"K"`, `"Na"`, `"Si"`, `"OH-"`); an analyte not plotted
at that age is absent.
"""
function d12_concentrations(system::AbstractString, age::Real)
    ps = d12_table("pore_solution")
    return Dict(
        ps.analyte[k] => ustrip(us"mol/m^3", ps.concentration[k]) / 1000
            for k in eachindex(ps.system) if ps.system[k] == system && ustrip(us"hr", ps.age[k]) ≈ age
    )
end

"""
    d12_solution(cs, system, age; model) -> (; state, certificate)

The pore solution of `system` at `age` (hours) as the figures give it, in a
kilogram of water (their mmol per liter taken as mmol per kilogram), at the
measured hydroxide. The concentrations fix the elements, and hydroxide is added
or taken until its molality is the one measured, so that the charge the analysis
leaves unbalanced is carried by it.

The measured quantity is a concentration, and a constraint holds an activity, so
the activity held is the measured molality times the activity coefficient of the
last solve, until the two agree to 1e-8.
"""
function d12_solution(cs, system::AbstractString, age::Real; model, T = 296.15)
    c = d12_concentrations(system, age)
    haskey(c, "OH-") || throw(ArgumentError("d12_solution: no hydroxide is plotted for $system at $age h"))
    st = ChemicalState(cs; T = T * u"K")
    Mw = ustrip(us"kg/mol", only(s for s in cs.species if symbol(s) == "H2O@")[:M])
    set_quantity!(st, "H2O@", (1 / Mw)u"mol")
    # an analyte not plotted is held at a trace
    for (el, sp) in D12_SPECIES
        set_quantity!(st, sp, get(c, el, 1.0e-12)u"mol")
    end
    b = budget(st)
    γ = 1.0
    local eq, cert
    from = st
    for _ in 1:20
        eq, cert = equilibrate_certified(from; model, b, constraint = FixedActivity("OH-", γ * c["OH-"]; titrant = "OH-"))
        γ_new = activity_coefficients(eq, model)["OH-"]
        converged = abs(log(γ_new / γ)) < 1.0e-8
        γ = γ_new
        converged && break
        from = eq          # the next activity is close: start from this answer
    end
    return (; state = eq, certificate = cert)
end

"""
    d12_indices(eq, systems, present; model, T = 296.15) -> Dict

The effective saturation indices of the solids of Table 3 in the solution `eq`:
log10(IAP/K) from the potentials of the primaries of the solution against the
energy of each solid, divided by the number of ions of its dissolution. A solid
one of whose analytes is not in `present` (aluminum before it is plotted,
calcium at 550 days) is left out.
"""
function d12_indices(eq, systems, present; model, T = 296.15)
    cs, full = systems.aqueous, systems.full
    RT = ChemistryLab.R_GAS * T
    g(sp) = ustrip(us"J/mol", sp[:ΔₐG⁰](T = T * u"K", P = 1.0e5u"Pa"; unit = true)) / RT
    lna = log_activities(eq, model)
    y = map(full.SM.primaries) do p
        s = symbol(p)
        k = findfirst(x -> symbol(x) == s, cs.species)
        k === nothing ? 0.0 : g(cs.species[k]) + lna[s]
    end
    A = Float64.(conservation_matrix(full))
    needs = Dict(
        "Portlandite" => ("Ca",), "Gp" => ("Ca", "S"), "ettringite" => ("Ca", "Al", "S"),
        "straetlingite" => ("Ca", "Al", "Si"), "monosulphate12" => ("Ca", "Al", "S")
    )
    out = Dict{String, Float64}()
    for (record, (label, ions)) in D12_SOLIDS
        all(el -> el in present, needs[record]) || continue
        j = findfirst(s -> symbol(s) == record, full.species)
        out[label] = (sum(A[c, j] * y[c] for c in eachindex(y)) - g(full.species[j])) / log(10) / ions
    end
    return out
end

"""The effective index of `label` in Table 3 for `system` at `age` hours (`nothing` when not given)."""
function d12_published(system, age, label)
    t = d12_table("effective_saturation_indices")
    i = findfirst(k -> t.system[k] == system && ustrip(us"hr", t.age[k]) ≈ age && t.phase[k] == label, eachindex(t.system))
    return i === nothing ? nothing : Float64(t.value[i])
end

d12_systems_listed() = unique(d12_table("pore_solution").system)
d12_ages() = sort(unique(ustrip.(us"hr", d12_table("pore_solution").age)))
