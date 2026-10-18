# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2026 Jean-François Barthélémy (Cerema, UMR MCD)

# From the `study` a script ends its calculation with to the plain numbers the
# page plots and tabulates. Dispatch on the type of `study.output` decides what
# is read; a new kind of calculation adds a method.

# JSON has no NaN or Inf: a value that is not finite (an enthalpy left unknown,
# say) is sent as null and shown as missing. These are display values of a
# finished Float64 run, hence the conversion.
finite(x::Real) = isfinite(x) ? Float64(x) : nothing

"A labeled series: the keys of an ordered dictionary and their values."
series(d::AbstractDict) = Dict{String, Any}("labels" => string.(collect(keys(d))), "values" => [finite(v) for v in values(d)])

"A quantity as a number in unit `u`: converted first, whatever unit it carries."
in_unit(u, x) = ustrip(u, uconvert(u, x))

"A parameter of a path as a number: days, a fraction, an amount, a step."
parameter_number(x::Real) = Float64(x)

"""
    extract(study) -> Dict{String, Any}

The results of a study, by the type of its `output`.
"""
extract(study) = extract_output(study.output, study)

"The name under which the aqueous species of a state are counted together."
const PORE_SOLUTION = "pore solution"

"""
The volume fractions of a state, its aqueous species counted together as the
pore solution (where the first of them comes), or why they are unknown.
"""
function fractions_entry(rs::RecipeState)
    v = volume(rs)
    if isempty(v.missing)
        aqueous = Set(symbol(s) for s in rs.state.system.species if aggregate_state(s) == AS_AQUEOUS)
        labels, values = String[], Float64[]
        for (k, f) in volume_fractions(rs)
            key = k in aqueous ? PORE_SOLUTION : k
            i = findfirst(==(key), labels)
            i === nothing ? (push!(labels, key); push!(values, f)) : (values[i] += f)
        end
        return Dict{String, Any}("labels" => labels, "values" => [finite(v) for v in values])
    end
    return Dict{String, Any}(
        "reason" => "No sourced density for the unreacted part of $(join(v.missing, ", ")): " *
            "its volume, and every volume fraction, is unknown.",
    )
end

function extract_output(rs::RecipeState, study)
    out = Dict{String, Any}("kind" => "equilibrium")
    out["binder_mass_g"] = finite(rs.recipe.binder_mass)
    out["phase_masses"] = series(phase_masses(rs))
    out["volume_fractions"] = fractions_entry(rs)
    p = porosity(rs)
    out["porosity"] = Dict{String, Any}("liquid" => finite(p.liquid), "void" => finite(p.void), "total" => finite(p.total))
    ps = pore_solution(rs)
    out["pore_solution"] = Dict{String, Any}("pH" => finite(ps.pH), "molalities" => series(ps.elements))
    out["ionic_strength"] = finite(ionic_strength(rs.state))
    out["bound_water"] = finite(bound_water(rs))
    si = saturation_indices(rs.state, rs.model)
    out["saturation_indices"] = series(sort(si; byvalue = true, rev = true))
    out["certificate"] = certificate_entry(study.certificate)
    return out
end

"""
The results of a path: for each phase and each dissolved element, a series over
the parameter; the pH, porosity, bound water and certification at each point.
"""
function extract_output(p::ProcessResult, study)
    states = collect(p)
    masses = [phase_masses(rs) for rs in states]
    phases = unique(k for d in masses for k in keys(d))
    solutions = [pore_solution(rs) for rs in states]
    elements = unique(k for s in solutions for k in keys(s.elements))
    fractions = [fractions_entry(rs) for rs in states]
    known = all(f -> haskey(f, "labels"), fractions)
    return Dict{String, Any}(
        "kind" => "path", "label" => string(p.label),
        "parameter" => [parameter_number(x) for x in p.parameter],
        "binder_mass_g" => finite(first(states).recipe.binder_mass),
        "phase_masses" => Dict{String, Any}(
            "labels" => phases, "series" => [[finite(get(d, k, 0.0)) for d in masses] for k in phases],
        ),
        "volume_fractions" => known ? volume_series(fractions) : Dict{String, Any}("reason" => first(f["reason"] for f in fractions if haskey(f, "reason"))),
        "pH" => [finite(s.pH) for s in solutions],
        "molalities" => Dict{String, Any}(
            "labels" => string.(elements), "series" => [[finite(get(s.elements, e, 0.0)) for s in solutions] for e in elements],
        ),
        "porosity" => [finite(porosity(rs).total) for rs in states],
        "bound_water" => [finite(bound_water(rs)) for rs in states],
        "certified" => [rs.certificate !== nothing && hasproperty(rs.certificate, :optimal) && rs.certificate.optimal === true for rs in states],
    )
end

"""
The results of a kinetic run: the degrees of reaction in time, from the steps of
the integration; at each age to report, the phases, the pore solution, the bound
water and the volumes of its certified equilibrium; the heat released, and the
temperature of a semi-adiabatic cell.
"""
function extract_output(o::NamedTuple{(:problem, :solution, :ages, :states)}, study)
    kp, sol, states = o.problem, o.solution, o.states
    day = 86400.0
    binder = Float64(study.recipe.binder_mass)
    times = sampled(sol.t)
    degrees = degrees_of_hydration(sol, kp; times)
    masses = [crystal_masses(st) for st in states]
    phases = unique(first(p) for v in masses for p in v)
    lookup = [Dict(v) for v in masses]
    solutions = [pore_solution(st, study.model) for st in states]
    elements = unique(k for s in solutions for k in keys(s.elements))
    fractions = [kinetic_fractions(st, kp.initial_state) for st in states]
    out = Dict{String, Any}(
        "kind" => "kinetics", "label" => "time", "parameter" => Float64.(o.ages), "binder_mass_g" => binder,
        "retcode" => string(sol.retcode),
        "time" => times ./ day,
        "degrees" => Dict{String, Any}("labels" => collect(keys(degrees)), "series" => [[finite(x) for x in v] for v in values(degrees)]),
        "phase_masses" => Dict{String, Any}(
            "labels" => phases, "series" => [[finite(get(d, k, 0.0)) for d in lookup] for k in phases],
        ),
        "volume_fractions" => volume_series(fractions),
        "pH" => [finite(s.pH) for s in solutions],
        "molalities" => Dict{String, Any}(
            "labels" => string.(elements), "series" => [[finite(get(s.elements, e, 0.0)) for s in solutions] for e in elements],
        ),
        # the pores: what the pore solution fills, and what the reactions emptied
        "porosity" => [finite(fraction_of(f, PORE_SOLUTION) + fraction_of(f, "void")) for f in fractions],
        "bound_water" => [finite(in_unit(us"g", bound_water(st)) / binder) for st in states],
    )
    out["heat"] = heat_entry(kp.calorimeter, sol, kp, o, binder)
    return out
end

"At most `n` of the instants of a solution, the last one kept, the start left out (a log axis)."
function sampled(t::AbstractVector; n::Integer = 300)
    ts = filter(>(0), Float64.(t))
    length(ts) <= n && return ts
    k = ceil(Int, length(ts) / n)
    return unique!(vcat(ts[1:k:end], [last(ts)]))
end

"The mass (g) of each solid of a state above a microgram, largest first."
function crystal_masses(st)
    cs = st.system
    out = [symbol(cs.species[i]) => in_unit(us"g", st.n[i] * cs.species[i][:M]) for i in cs.idx_crystal]
    return sort!(filter!(p -> last(p) > 1.0e-6, out); by = last, rev = true)
end

"The fraction a series of volume fractions gives `key`, zero when it has none."
function fraction_of(f, key)
    i = findfirst(==(key), f["labels"])
    return i === nothing ? 0.0 : f["values"][i]
end

"The volume fractions of a state of a run, relative to the volume at its start, the aqueous species together."
function kinetic_fractions(st, start)
    aqueous = Set(symbol(st.system.species[i]) for i in st.system.idx_aqueous)
    labels, values = String[], Float64[]
    for (k, f) in volume_fractions(st; reference = start)
        key = k in aqueous ? PORE_SOLUTION : k
        i = findfirst(==(key), labels)
        i === nothing ? (push!(labels, key); push!(values, f)) : (values[i] += f)
    end
    return Dict{String, Any}("labels" => labels, "values" => values)
end

# The heat released (J per g of binder) at the ages, from the certified states
# for an isothermal run; the cell's own balance in time, and its temperature,
# for a semi-adiabatic one. A species without an enthalpy of formation leaves
# the heat unknown, and the reason is said.
function heat_entry(::Nothing, sol, kp, o, binder)
    isempty(o.ages) && return Dict{String, Any}("time" => Float64[], "heat" => Float64[])
    missing_h = missing_enthalpy(kp.initial_state)
    isempty(missing_h) || return Dict{String, Any}("reason" => "No enthalpy of formation for $(join(missing_h, ", ")): the heat is unknown.")
    _, Q, _ = heat_release(sol, kp; times = o.ages .* 86400.0, states = o.states, reference = kp.initial_state)
    return Dict{String, Any}("time" => Float64.(o.ages), "heat" => [finite(q / binder) for q in Q])
end
function heat_entry(cal::SemiAdiabaticCalorimeter, sol, kp, o, binder)
    t, Q = cumulative_heat(sol, cal)
    keep = findall(>(0), t)
    _, T = temperature_profile(sol, cal; times = t[keep])
    return Dict{String, Any}(
        "time" => t[keep] ./ 86400.0, "heat" => [finite(q / binder) for q in Q[keep]],
        "temperature" => [finite(x - ustrip(T_ZERO_CELSIUS_Q)) for x in T],
    )
end

"Volume fractions of the points of a path, as one series per phase."
function volume_series(fractions)
    labels = unique(l for f in fractions for l in f["labels"])
    lookup = [Dict(zip(f["labels"], f["values"])) for f in fractions]
    return Dict{String, Any}("labels" => labels, "series" => [[get(d, l, 0.0) for d in lookup] for l in labels])
end

"The fields of a certificate the page shows, those it has."
function certificate_entry(cert)
    out = Dict{String, Any}()
    for key in (:optimal, :scope, :balance, :worst_supersaturation, :within_activity_range)
        hasproperty(cert, key) || continue
        v = getproperty(cert, key)
        out[String(key)] = v isa Bool ? v : v isa Real ? finite(v) : string(v)
    end
    return out
end
