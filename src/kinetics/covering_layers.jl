# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# ── covering_layers.jl ────────────────────────────────────────────────────────
#
# A thin layer that forms on a mineral while a phase is stable and blocks part of
# its sites, and opens again, from the sites already free, once that phase is not
# stable: a passivation counted in sites. The layer is a species of its own,
# whose amount is the moles of sites it blocks.

"""
    CoveringLayer(layer, substrate, area, Γ; blockable = 1, floor = 1e-6) -> CoveringLayer

A layer, the species named `layer`, covering the surface of the solid
`substrate`, whose area model is `area`, at `Γ` blockable sites per m² (mol/m²,
a `Quantity` or a plain `Real`). The amount of `layer` counts the sites it
blocks, so that the sites of the substrate and the fraction left free are

```math
N = \\Gamma\\, \\mathcal{A}(n_s, n_{s,0}),
\\qquad
\\varphi = 1 - \\frac{n_\\ell}{N} ,
```

`φ` floored at `floor` (continuously with its slope). At most a fraction
`blockable` of the sites can be blocked, the rest — the defects — staying
open. [`covered`](@ref) multiplies a rate law of the substrate by `φ`, and
[`covering_rate`](@ref) is the law of the layer itself.
"""
struct CoveringLayer{A <: AbstractSurfaceModel, T <: Real}
    layer::String
    substrate::String
    area::A
    Γ::T
    blockable::T
    floor::T
end

function CoveringLayer(layer::AbstractString, substrate::AbstractString, area::AbstractSurfaceModel, Γ; blockable::Real = 1, floor::Real = 1.0e-6)
    g = _area_si(us"mol/m^2", Γ, "CoveringLayer Γ")
    g > 0 || throw(ArgumentError("CoveringLayer: Γ must be positive; got $g mol/m²."))
    0 < blockable <= 1 || throw(ArgumentError("CoveringLayer: `blockable` is a fraction in (0, 1]; got $blockable."))
    0 < floor < 1 || throw(ArgumentError("CoveringLayer: `floor` is a fraction in (0, 1); got $floor."))
    v = promote(g, float(blockable), float(floor))
    return CoveringLayer{typeof(area), eltype(v)}(String(layer), String(substrate), area, v...)
end

# The sites of the substrate at the amounts of the run, mol.
function _layer_capacity(c::CoveringLayer, n, n₀, M)
    ns = max(n[c.substrate], oneunit(n[c.substrate]) * 1.0e-30)
    return c.Γ * total_area(c.area, ns, max(n₀[c.substrate], ns), M)
end

# The free fraction, floored continuously with its slope: φ_min + c(1 − φ_min − θ).
function _free_fraction(c::CoveringLayer, n, n₀, M)
    θ = max(n[c.layer], zero(n[c.layer])) / _layer_capacity(c, n, n₀, M)
    return c.floor + _c1_cut(1 - c.floor - θ, c.floor)
end

function _substrate_molar_mass(cs::ChemicalSystem, c::CoveringLayer)
    for s in (c.layer, c.substrate)
        haskey(cs.dict_species, s) || throw(ArgumentError("CoveringLayer: \"$s\" is not a species of the system."))
    end
    return _molar_mass_si(cs.dict_species[c.substrate])
end

"""
    covered(rate, cs, layer::CoveringLayer) -> KineticFunc

The rate law `rate` of the substrate of `layer`, multiplied by the fraction of
its sites the layer leaves free, `φ` ([`CoveringLayer`](@ref)): a dissolution
that only the open sites carry.
"""
function covered(rate, cs::ChemicalSystem, layer::CoveringLayer)
    M = _substrate_molar_mass(cs, layer)
    f = (T, P, t, n, lna, n₀) -> rate(T, P, t, n, lna, n₀) * _free_fraction(layer, n, n₀, M)
    return KineticFunc(f, _rate_refs(rate), u"mol/s")
end

"""
    covering_rate(layer::CoveringLayer, cs, rxn; k_form, τ_open, opening = nothing) -> KineticFunc

The rate law of the layer, `rxn` being the dissolution of the species `layer`
counts, `Ω` its saturation ratio. While the phase of the layer is stable it
forms on the blockable sites still open; once it is not, it opens from the sites
already free, a logistic recession:

```math
r = \\begin{cases}
  -\\,k_\\text{form}\\,\\bigl(b N - n_\\ell\\bigr)^+ \\dfrac{\\Omega - 1}{\\Omega + 1} & \\Omega \\ge 1 \\\\[1ex]
  \\dfrac{n_\\ell}{\\tau_\\text{open}} \\Bigl(1 - \\dfrac{n_\\ell}{N}\\Bigr)^+ (1 - \\Omega)\\; g & \\Omega < 1
\\end{cases}
\\qquad [\\text{mol/s}],
```

negative while it forms, the reaction being written as a dissolution; `b` is
`blockable`, `N` the sites of the substrate, `k_form` in s⁻¹, `τ_open` in s.
At a fixed `Ω < 1` the blocked fraction `θ = n_ℓ/N` follows
`dθ/dt = −θ(1 − θ)(1 − Ω)/τ`. The two branches meet at zero at `Ω = 1`,
continuously but not with their slopes. `opening`, a function
`(T, P, t, n, lna, n₀) -> g` of the run, slows the opening — a species that
holds the layer in place, say — and is one when not given.
"""
function covering_rate(
        layer::CoveringLayer, cs::ChemicalSystem, rxn::AbstractReaction;
        k_form::Real, τ_open::Real, opening = nothing,
    )
    (k_form >= 0 && τ_open > 0) || throw(ArgumentError("covering_rate: k_form ≥ 0 and τ_open > 0 are required."))
    M = _substrate_molar_mass(cs, layer)
    terms = _saturation_terms(cs, rxn, "covering_rate")
    g = opening === nothing ? _AlwaysOpen() : opening
    return _covering_law(layer, terms, M, k_form, τ_open, g)
end

struct _AlwaysOpen end
(::_AlwaysOpen)(T, P, t, n, lna, n₀) = true

function _covering_law(layer, terms, M, k_form, τ_open, g)
    rate = (T, P, t, n, lna, n₀) -> begin
        Ω = exp(_log_saturation(terms, T, P, lna))
        N = _layer_capacity(layer, n, n₀, M)
        nℓ = max(n[layer.layer], zero(n[layer.layer]))
        # Not symbolic: the branch is on the stability of the layer's phase.
        _primal(Ω) >= 1 ? -k_form * max(layer.blockable * N - nℓ, zero(nℓ)) * (Ω - 1) / (Ω + 1) :
            nℓ / τ_open * max(1 - nℓ / N, zero(nℓ)) * (1 - Ω) * g(T, P, t, n, lna, n₀)
    end
    return KineticFunc(rate, (T = T_STANDARD_Q, P = P_STANDARD_Q), u"mol/s")
end
