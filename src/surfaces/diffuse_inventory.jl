# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# ── What a diffuse layer holds: the Donnan approach ──────────────────────────
#
# A `DiffuseLayer` family raises a potential from its charge, and the ions that
# screen that charge are left implicit: the solution carries the counter-charge
# and the surface is not neutral on its own. That is PHREEQC's default too. Its
# `SURFACE -Donnan` makes the layer explicit, and this file reproduces it: the
# surface speciation and its Gouy-Chapman potential are unchanged, and a layer of
# water of fixed thickness is added on the surface, holding each solute at an
# average Boltzmann enrichment chosen so that the layer's charge balances the
# surface's. The layer's solutes are withdrawn from the solution; its water is
# added, as PHREEQC adds it when a solution and a surface are reacted together.

"""
    DonnanLayer(; thickness = 1.0e-8u"m")

The ions of the diffuse layer that balances a charged surface, averaged over a
layer of water of fixed `thickness`: the Donnan approach that
[AppeloWersin2007](@cite) added to PHREEQC as `SURFACE -Donnan`, which
[`equilibrate_donnan`](@ref) solves with and [`diffuse_layer_contents`](@ref)
reads.

On each support carrying a [`DiffuseLayer`](@ref) family, of area `𝒜`, the layer
holds `W_D = ρ 𝒜 t` of water (`ρ = 1000 kg/m³`, `t` the thickness) and each
solute `i` of the solution at

```math
n_i^{D} = m_i \\, W_D \\, e^{-z_i \\tilde\\psi_D},
\\qquad
\\frac{\\mathcal{A}\\,\\sigma}{F} + W_D \\sum_i z_i\\, m_i\\, e^{-z_i \\tilde\\psi_D} = 0 ,
```

with `m_i` the molality of the free solution and `ψ̃_D = FΨ_D/RT` the average
potential of the layer, the one whose charge balances the surface charge `σ`.
It is not the surface potential: the surface complexes still see the potential
Gouy-Chapman gives, and `ψ̃_D` is the smaller average over the layer. The
left-hand side decreases in `ψ̃_D`, so the root is unique.

The solutes of the layer come out of the system's totals, and the free solution
holds the rest. Where its water comes from is a convention, and
[`equilibrate_donnan`](@ref) offers PHREEQC's two: added to the solution's,
which PHREEQC does when a solution and a surface are reacted together, or taken
from it, which PHREEQC does for a surface equilibrated beforehand and which is
the one a closed pore solution obeys. The **excess** of a solute,
`n_i^D − m_i W_D`, is what the layer holds beyond its water at the molality of
the free solution: positive for a counter-ion, negative for a co-ion.

`thickness` accepts a length or a number of meters. PHREEQC's default is 10 nm;
a layer thicker than the pores it lines is not a physical statement, and the
fixed point of [`equilibrate_donnan`](@ref) stops converging well before it.
"""
struct DonnanLayer{T <: Real}
    thickness::T
    function DonnanLayer(thickness::Real)
        thickness > 0 || throw(ArgumentError("thickness must be positive; got $thickness m."))
        t = float(thickness)
        return new{typeof(t)}(t)
    end
end
DonnanLayer(; thickness = 1.0e-8u"m") =
    DonnanLayer(thickness isa Real ? thickness : ustrip(us"m", thickness))

"""
    diffuse_layer_contents(state, layer::DonnanLayer) -> NamedTuple

What the Donnan layers of `state` hold, one per support carrying a
[`DiffuseLayer`](@ref) family (see [`DonnanLayer`](@ref)):

  - `amounts`: moles of each species of the system in the layers, zero but for
    the aqueous solutes and the solvent, whose entry is the layers' water;
  - `excess`: the same, less what the layers' water holds at the molality of the
    free solution — the ions the layers add to the solution's totals;
  - `water`: the layers' water, kg;
  - `potential`: `ψ̃_D = FΨ_D/RT` of each layer, by support name;
  - `charge`: the surface charge each layer balances, moles of charge.

`state` is a solution with its surfaces: the molalities are those of its free
water, so the layers' water is not in the solvent of `state`.
"""
function diffuse_layer_contents(state::ChemicalState, layer::DonnanLayer)
    cs = state.system
    n = [ustrip(us"mol", x) for x in state.n]
    iw = only(cs.idx_solvent)
    W = n[iw] * ustrip(us"kg/mol", cs.species[iw][:M])
    W > 0 || throw(ArgumentError("diffuse_layer_contents: the state has no solvent."))
    solutes = cs.idx_solutes
    z = Float64[charge(cs.species[i]) for i in solutes]
    m = n[solutes] ./ W
    idx = Dict(symbol(sp) => k for (k, sp) in enumerate(cs.species))

    # In the number type of the state, of the layer and of the areas: a
    # thickness, an area or an amount being differentiated carries its
    # derivative into what the layers hold.
    groups = _diffuse_layer_groups(cs)
    R = promote_type(
        Float64, eltype(n), typeof(layer.thickness),
        (typeof(float(f.model.area)) for fs in values(groups) for f in fs)...,
    )
    amounts = zeros(R, length(n))
    excess = zeros(R, length(n))
    potential = Dict{String, R}()
    charges = Dict{String, R}()
    water = zero(R)
    for (support, families) in groups
        area = only(unique(f.model.area for f in families))
        σ = sum(
            charge(sp) * n[idx[symbol(sp)]] for f in families for sp in vcat([f.free_site], f.complexes)
        )
        W_D = 1000 * area * layer.thickness
        ψ = _donnan_potential(σ, z, m, W_D)
        for (k, i) in enumerate(solutes)
            held = m[k] * W_D * exp(-z[k] * ψ)
            amounts[i] += held
            excess[i] += held - m[k] * W_D
        end
        amounts[iw] += W_D / ustrip(us"kg/mol", cs.species[iw][:M])
        potential[support] = ψ
        charges[support] = σ
        water += W_D
    end
    return (; amounts, excess, water, potential, charge = charges)
end

# The DiffuseLayer families of a system, by the support they share: one layer,
# and one area, per support.
function _diffuse_layer_groups(cs::ChemicalSystem)
    groups = OrderedDict{String, Vector{Any}}()
    for f in cs.site_families
        f.model isa DiffuseLayer || continue
        push!(get!(groups, f.support.name, Any[]), f)
    end
    isempty(groups) && throw(
        ArgumentError(
            "a Donnan layer balances the charge of a DiffuseLayer family, and the " *
                "system has none."
        )
    )
    return groups
end

# The average potential ψ̃ whose layer balances σ moles of surface charge:
# σ + W Σ z_i m_i exp(-z_i ψ̃) = 0, decreasing in ψ̃, by Newton from the
# symmetric-electrolyte root, steps capped at one unit as PHREEQC caps them.
#
# On dual numbers the derivative is exact: once the residual vanishes, the next
# Newton step carries the implicit-function derivative `−f_θ/f_ψ`, and the loop
# stops on a step below the tolerance, one step after the residual vanished.
function _donnan_potential(σ::Real, z::AbstractVector, m::AbstractVector, W::Real)
    iszero(σ) && return zero(promote_type(typeof(σ), eltype(m), typeof(W)))
    I = 0.5 * sum(z[k]^2 * m[k] for k in eachindex(z))
    ψ = asinh(σ / (2 * W * max(I, eps())))
    for _ in 1:100
        # The two sums in one loop: a generator here captured `ψ`, which the
        # iteration reassigns, and Julia boxed it.
        sf = zero(promote_type(typeof(ψ), eltype(m), eltype(z)))
        sdf = sf
        for k in eachindex(z)
            e = m[k] * exp(-z[k] * ψ)
            sf += z[k] * e
            sdf += z[k]^2 * e
        end
        f = σ + W * sf
        df = -W * sdf
        step = clamp(-f / df, -1.0, 1.0)
        ψ += step
        abs(step) < 1.0e-14 * max(1.0, abs(ψ)) && return ψ
    end
    throw(
        ErrorException(
            "the Donnan potential did not converge for a surface charge of $(_plain(σ)) mol; " *
                "the layer's solution cannot balance it."
        )
    )
end

"""
    equilibrate_donnan(state, layer::DonnanLayer; model = DiluteSolutionModel(),
                       water = :added, b = nothing, maxiter = 200, rtol = 1e-8,
                       kwargs...)
        -> NamedTuple

The equilibrium of `state` with the ions of its diffuse layers counted, as
PHREEQC's `SURFACE -Donnan` counts them (see [`DonnanLayer`](@ref)).

The layers take their solutes out of the solution's totals, and what they take
depends on the solution they leave. Each step solves the solution and its
surfaces with [`equilibrate_certified`](@ref) on the totals less what the
layers held, from the previous answer, until the layers' contents stop moving
by more than `rtol` relative. A layer that holds `r` times what the free
solution holds of a solute would make a plain substitution oscillate with a
factor `r`, and diverge beyond one; each step therefore moves that solute only
the fraction `1/(1 + r)` of the way, which cancels the factor. `b` is the system's total
budget, the layers included; by default the one `state` holds.

`water` says where the layers' water comes from. `:added`, the default, adds it
to the solution's: PHREEQC's convention for a solution and a surface reacted
together, and the one its `-Donnan` results are compared against. `:taken`
takes it out of the solution's water, so that the free water is the total less
the layers': the convention of a closed system, and the one to use when the
layers hold a noticeable share of the water, which `:added` would otherwise add
to the system and dilute it by.
Keywords other than these are passed to [`equilibrate_certified`](@ref).

Returns `(state, certificate, layer, iterations)`: the free solution with its
surfaces and solids, the certificate of its last solve, and
[`diffuse_layer_contents`](@ref) of that state. An iteration that does not
settle within `maxiter` steps is an error that gives the last relative change.
The more of the water the layers hold, the slower they settle: a layer holding
70 % of it takes about a hundred steps.
"""
function equilibrate_donnan(
        state::ChemicalState, layer::DonnanLayer;
        model::AbstractActivityModel = DiluteSolutionModel(), water::Symbol = :added,
        b = nothing, maxiter::Integer = 200, rtol::Real = 1.0e-8, kwargs...,
    )
    _refuse_state_keywords(kwargs, "equilibrate_donnan")
    water in (:added, :taken) || throw(
        ArgumentError("water must be :added or :taken; got :$water.")
    )
    cs = state.system
    _diffuse_layer_groups(cs)
    # In the number type of the data: a budget, a state or a layer being
    # differentiated is carried through the fixed point, whose derivative is
    # that of the iterations reaching it.
    A = conservation_matrix(cs)
    b_total = b === nothing ? A * [ustrip(us"mol", x) for x in state.n] : collect(b)
    iw = only(cs.idx_solvent)
    held = zeros(promote_type(Float64, eltype(b_total), typeof(layer.thickness)), length(cs.species))
    current = state
    last = Inf
    for it in 1:maxiter
        eq, cert = equilibrate_certified(current; model, b = b_total - A * held, kwargs...)
        contents = diffuse_layer_contents(eq, layer)
        target = copy(contents.amounts)
        water === :added && (target[iw] = 0.0)
        # The model or the temperature may carry duals the budget does not.
        held = convert(Vector{promote_type(eltype(held), eltype(target))}, held)
        change = maximum(abs, target - held)
        last = change / max(maximum(abs, target), eps())
        last <= rtol &&
            return (; state = eq, certificate = cert, layer = contents, iterations = it)
        # What the layers hold of each species against what the free solution
        # holds: the factor a plain substitution would oscillate with.
        free = [ustrip(us"mol", x) for x in eq.n]
        for i in eachindex(held)
            r = free[i] > 0 ? target[i] / free[i] : 0.0
            held[i] += (target[i] - held[i]) / (1 + r)
        end
        current = eq
    end
    throw(
        ErrorException(
            "equilibrate_donnan: the diffuse layer did not settle in $maxiter steps " *
                "(last relative change $(round(_plain(last); sigdigits = 3)), rtol $rtol); a layer " *
                "holding a large share of the water may have no fixed point."
        )
    )
end
