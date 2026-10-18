# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2026 Jean-François Barthélémy (Cerema, UMR MCD)

# The conditions of the mix: water, binder mass, temperature and pressure.

"""
    MixChoice

The conditions of the mix.
"""
abstract type MixChoice <: Choice end

"The water-to-binder ratio, the binder mass, the temperature and the pressure."
struct Mix <: MixChoice
    w_b::Float64
    temperature::Float64
    binder_mass::Float64
    pressure::Float64
end
kind(::Type{Mix}) = "mix"
fields(::Type{Mix}) = [
    Field("w_b", :number; default = 0.5, min = 0.01, max = 100),
    Field("temperature", :number; default = 20.0, unit = "°C", min = 0, max = 99),
    Field("binder_mass", :number; default = 100.0, unit = "g", min = 1.0e-9, level = :expert),
    Field("pressure", :number; default = 1.0, unit = "bar", min = 0.5, max = 1000, level = :expert),
]
members(::Type{MixChoice}) = (Mix,)

"The lines that set the conditions, before the materials (their laws may read them)."
mix_lines(x::Mix) = [
    "# The conditions of the mix, set first: the laws of the degrees of reaction",
    "# may depend on them.",
    "w_b = $(julia_number(x.w_b))                        # water-to-binder ratio",
    "T = T_ZERO_CELSIUS_Q + $(julia_number(x.temperature))u\"K\"   # $(julia_number(x.temperature)) °C, in kelvins",
]

"The keywords of `Recipe` for the conditions."
function recipe_keywords(x::Mix)
    kws = ["w_b", "T"]
    x.pressure == 1 || push!(kws, "P = $(julia_number(x.pressure))u\"bar\"")
    x.binder_mass == 100 || push!(kws, "binder_mass = $(julia_number(x.binder_mass))u\"g\"")
    return join(kws, ", ")
end
