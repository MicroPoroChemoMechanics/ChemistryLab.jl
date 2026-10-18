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

"""
    extract(study) -> Dict{String, Any}

The results of a study, by the type of its `output`.
"""
extract(study) = extract_output(study.output, study)

function extract_output(rs::RecipeState, study)
    out = Dict{String, Any}("kind" => "equilibrium")
    out["binder_mass_g"] = finite(rs.recipe.binder_mass)
    out["phase_masses"] = series(phase_masses(rs))
    v = volume(rs)
    if isempty(v.missing)
        out["volume_fractions"] = series(volume_fractions(rs))
    else
        out["volume_fractions"] = Dict{String, Any}(
            "reason" => "No sourced density for the unreacted part of $(join(v.missing, ", ")): " *
                "its volume, and every volume fraction, is unknown.",
        )
    end
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
