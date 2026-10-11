# =============================================================================
#  deschner_2013.jl — the pore solutions of Deschner et al. (2013), 7 to 80 °C
#
#  The pore solutions of a CEM I blended with 50 % of quartz powder or of a
#  siliceous fly ash, cured sealed at 7, 23, 40, 50 and 80 °C and analyzed from
#  one day to 250 days (their Tables A.1 and A.2), with the effective saturation
#  indices the authors computed from them (their Table B.1), in
#  data/literature/Deschner2013.json. Each solution is speciated at its measured
#  hydroxide and at its curing temperature, and the same indices are computed:
#  the method of scripts/deschner_2012.jl, whose functions are used here, carried
#  over the temperatures of the run.
# =============================================================================

isdefined(@__MODULE__, :d12_indices) || include(joinpath(@__DIR__, "deschner_2012.jl"))

const D13 = "Deschner2013"
d13_table(name) = literature_table(D13, name)

# The analytes of Tables A.1 and A.2 and the species that carry them. Iron and
# chloride, analyzed here and not in 2012, are kept: they add to the ionic
# strength and complex nothing the indices read.
const D13_SPECIES = [
    "Ca" => "Ca+2", "Al" => "AlO2-", "S" => "SO4-2", "K" => "K+", "Na" => "Na+",
    "Si" => "SiO2@", "Fe" => "FeO2-", "Cl" => "Cl-",
]

# The four solids of Table B.1, with the labels of the table.
const D13_LABELS = Dict("Portlandite" => "CH", "ettringite" => "Ettr", "monosulphate12" => "Ms", "straetlingite" => "Strätl")

"""
    d13_systems() -> (; aqueous, full, substances)

The aqueous species of Cemdata18 for the analytes of Tables A.1 and A.2, and the
same with the solids of `d12_systems` beside them, with iron and chlorine
added.
"""
function d13_systems()
    substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
    # Every solid of `D12_SOLIDS`, gypsum included, so that `d12_indices` finds
    # each of them; Table B.1 gives no index for gypsum, and none is reported.
    sp = speciation(
        substances, vcat(first.(D12_SOLIDS), ["K+", "Na+", "SiO2@", "FeO2-", "Cl-"]);
        aggregate_state = [AS_AQUEOUS], exclude_species = D12_EXCLUDE,
    )
    aq = [s for s in sp if aggregate_state(s) == AS_AQUEOUS]
    return (; aqueous = ChemicalSystem(aq, CEMDATA_PRIMARIES), full = ChemicalSystem(sp, CEMDATA_PRIMARIES), substances)
end

"""
    d13_concentrations(system, T, age) -> Dict{String, Float64}

What Tables A.1 and A.2 give for `system` cured at `T` (K) after `age` days, in
mol per liter, by analyte; an analyte not available is absent.
"""
function d13_concentrations(system::AbstractString, T::Real, age::Real)
    ps = d13_table("pore_solution")
    return Dict(
        ps.analyte[k] => ustrip(us"mol/m^3", ps.concentration[k]) / 1000
            for k in eachindex(ps.system)
            if ps.system[k] == system && ustrip(us"K", ps.temperature[k]) ≈ T && ustrip(us"d", ps.age[k]) ≈ age
    )
end

"""
    d13_solution(cs, system, T, age; model) -> (; state, certificate)

The pore solution of `system` cured at `T` after `age` days, in a kilogram of
water, at the measured hydroxide and at `T`, as `d12_solution` sets it up: the
activity of hydroxide held is the measured molality times the activity
coefficient of the last solve, until the two agree to 1e-8.
"""
function d13_solution(cs, system::AbstractString, T::Real, age::Real; model)
    c = d13_concentrations(system, T, age)
    haskey(c, "OH-") || throw(ArgumentError("d13_solution: no hydroxide for $system at $T K after $age d"))
    st = ChemicalState(cs; T = T * u"K")
    Mw = ustrip(us"kg/mol", only(s for s in cs.species if symbol(s) == "H2O@")[:M])
    set_quantity!(st, "H2O@", (1 / Mw)u"mol")
    for (el, sp) in D13_SPECIES
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
        from = eq
    end
    return (; state = eq, certificate = cert)
end

"""
    d13_indices(eq, systems, T; model) -> Dict

The effective saturation indices of the solids of Table B.1 in the solution `eq`
at `T`, by `d12_indices`, under the labels of Table B.1.
"""
function d13_indices(eq, systems, T; model)
    present = ("Ca", "Al", "S", "Si")
    ix = d12_indices(eq, systems, present; model, T)
    label12 = Dict(label => record for (record, (label, _)) in D12_SOLIDS)
    return Dict(D13_LABELS[label12[l]] => v for (l, v) in ix if haskey(D13_LABELS, label12[l]))
end

"""The effective index of `label` in Table B.1 for `system` at `T` after `age` days (`nothing` when not given)."""
function d13_published(system, T, age, label)
    t = d13_table("effective_saturation_indices")
    i = findfirst(
        k -> t.system[k] == system && ustrip(us"K", t.temperature[k]) ≈ T && ustrip(us"d", t.age[k]) ≈ age && t.phase[k] == label,
        eachindex(t.system),
    )
    return i === nothing ? nothing : Float64(t.value[i])
end

"""The (system, temperature, age) of every solution Table B.1 gives indices for."""
function d13_solutions()
    t = d13_table("effective_saturation_indices")
    return unique((t.system[k], ustrip(us"K", t.temperature[k]), ustrip(us"d", t.age[k])) for k in eachindex(t.system))
end
