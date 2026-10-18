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
