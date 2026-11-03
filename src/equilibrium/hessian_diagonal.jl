# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# ── The diagonal of the Gibbs Hessian, without the Jacobian ───────────────────
#
# The interior point needs `∂μᵢ/∂nᵢ` at every iteration, and only that: the
# diagonal of `∇²G = ∂μ/∂n` (Gibbs–Duhem). Read off a forward-mode Jacobian of
# `μ` it costs `ns / chunk` evaluations of the whole activity model on dual
# numbers — six for 66 species — for `ns` numbers out of `ns²`. On a cement
# paste with the C-S-H surface of Guo et al. that Jacobian was 30 % of the time
# of an equilibrium and 40 % of what it allocated.
#
# The structure of the models makes the diagonal cheap, exactly:
#
#  - a solute's `ln aᵢ = ln10 log₁₀γᵢ(I) + ln mᵢ` depends on its own amount through
#    `mᵢ = nᵢ/(n_w M_w)` and the ionic strength `I`, so
#    `∂ln aᵢ/∂nᵢ = 1/nᵢ + (∂ln aᵢ/∂I) zᵢ²/(2 n_w M_w)`, with `∂ln aᵢ/∂I` for every
#    solute in one pass of the model's own solute terms on a dual `I`;
#  - the solvent's row is one pass of the model's own solvent term on a dual `n_w`;
#  - the gas phase, the solid solutions and the site families depend on the
#    amounts of their own members (and on `I`, which none of them changes): their
#    diagonal is that of a forward-mode Jacobian over those members only, with
#    whatever model they follow, their electrostatics included.
#
# Each piece is the model's own code differentiated, not a formula written a
# second time, so the diagonal cannot drift from the activities. An amount at or
# below the activity floor is evaluated at the floor, and its row is then flat:
# its entry is zero, as the Jacobian of `max.(n, floor)` makes it.

"""
    _hessian_diagonal_hook(μ, q) -> Union{Nothing, Function}

The exact diagonal of the Gibbs Hessian `∂μᵢ/∂nᵢ`, floored at zero, as a function
`(hf, n) -> hf` filling `hf` at the amounts `n` and the parameters `q`, for an
interior point to read in place of a forward-mode Jacobian of `μ`; `nothing` when
the potentials `μ` do not provide one (see [`_lna_diagonal!`](@ref)).

Floored at zero as the interior point floors the diagonal it reads off the
Jacobian: a pure phase has none, and the barrier supplies the curvature there.
"""
_hessian_diagonal_hook(μ, q) = nothing
function _hessian_diagonal_hook(μ::_Potentials, q)
    _has_lna_diagonal(μ.lna) || return nothing
    lna = μ.lna
    return function (hf, n)
        _lna_diagonal!(hf, lna, n, q)
        return hf .= max.(hf, zero(eltype(hf)))
    end
end

"""
    _has_lna_diagonal(lna) -> Bool

Whether [`_lna_diagonal!`](@ref) has a method for the activity closure `lna`:
the Debye–Hückel family ([`HKFActivityModel`](@ref), [`DaviesActivityModel`](@ref),
[`TruesdellJonesActivityModel`](@ref), [`LLNLActivityModel`](@ref)), whatever
gas, solid solutions and site families the system holds.
"""
_has_lna_diagonal(lna) = false
_has_lna_diagonal(::_DebyeHuckelLna) = true

"""
    _lna_diagonal!(d, lna, n, p) -> d

Write into `d` the diagonal `∂ln aᵢ/∂nᵢ` of the Jacobian of the log activities
`lna(n, p)`, without forming the Jacobian; equal to it to rounding.
"""
function _lna_diagonal!(d::AbstractVector, f::_DebyeHuckelLna, n::AbstractVector, p)
    ϵ = p.ϵ
    fl = _activity_floor(p)
    _n = _Floored(n, fl)
    fill!(d, zero(eltype(d)))
    _aqueous_diagonal!(d, f, n, _n, fl, p, ϵ)
    _mixing_diagonal!(d, f, n, _n, fl, p, ϵ)
    return d
end

# Whether the activities read the amount `x` itself, rather than the floor.
_above_floor(x, fl) = isless(fl, x)

function _aqueous_diagonal!(d, f::_DebyeHuckelLna, n, _n, fl, p, ϵ)
    aq = f.aq
    A, B = _coefficients(f.kernel, p)
    n_w = _n[aq.idx_solvent]
    denom_mol = n_w * aq.M_w
    I = _ionic_strength(aq, _n, denom_mol)

    # ∂ln aᵢ/∂I of every solute: the model's solute terms on a dual `I`, each
    # derivative written into `d` as it is computed.
    TI = typeof(ForwardDiff.Tag(_aqueous_diagonal!, typeof(I)))
    ID = ForwardDiff.Dual{TI}(I, one(I))
    _solute_lna!(_DerivativeSink{TI}(d), f.kernel, aq, _n, denom_mol, ID, sqrt(ID + ϵ), A, B)
    @inbounds for i in aq.idx_solutes
        d[i] = _above_floor(n[i], fl) ? inv(_n[i]) + d[i] * aq.zv[i]^2 / denom_mol / 2 : zero(eltype(d))
    end

    # The solvent row on a dual `n_w`, its kilograms and the ionic strength
    # following.
    w = aq.idx_solvent
    if _above_floor(n[w], fl)
        Tw = typeof(ForwardDiff.Tag(_solvent_row, typeof(n_w)))
        n_wD = ForwardDiff.Dual{Tw}(n_w, one(n_w))
        denomD = n_wD * aq.M_w
        IwD = _ionic_strength(aq, _n, denomD)
        rw = _solvent_row(f.kernel, aq, _n, n_wD, denomD, IwD, sqrt(IwD + ϵ), A, B, ϵ)
        d[w] = ForwardDiff.extract_derivative(Tw, rw)
    end
    return d
end

# The gas phase, the solid solutions and the site families: the diagonal of the
# Jacobian of their activities with respect to their own members, the other
# amounts held at their (floored) values.
function _mixing_diagonal!(d, f::_DebyeHuckelLna, n, _n, fl, p, ϵ)
    idx = f.idx_mixing
    isempty(idx) && return d
    block = let f = f, _n = _n, fl = fl, p = p, ϵ = ϵ
        v -> _mixing_block(f, _n, fl, p, ϵ, v)
    end
    J = ForwardDiff.jacobian(block, n[idx])
    @inbounds for (k, i) in enumerate(idx)
        d[i] = J[k, k]
    end
    return d
end

function _mixing_block(f::_DebyeHuckelLna{K, G, X, M}, _n, fl, p, ϵ, v::AbstractVector{T}) where {K, G, X, M, T}
    idx = f.idx_mixing
    nD = Vector{T}(undef, length(_n))
    @inbounds for i in eachindex(nD)
        nD[i] = _n[i]
    end
    @inbounds for (k, i) in enumerate(idx)
        nD[i] = max(v[k], fl)
    end
    # in the number type of the amounts, the parameters and the data, as the
    # activities themselves are
    out = zeros(promote_type(T, _number_type_of(p), M), length(_n))
    isempty(f.idx_gas) || _gas_lna!(out, nD, f.idx_gas, p, f.gas_mix)
    _mixing_lna!(out, nD, f.mix, p, ϵ)
    return out[idx]
end

# The amounts `max.(n, floor)` the activities are evaluated at, read lazily.
struct _Floored{T, V <: AbstractVector, F} <: AbstractVector{T}
    n::V
    floor::F
end
_Floored(n::AbstractVector, floor) = _Floored{promote_type(eltype(n), typeof(floor)), typeof(n), typeof(floor)}(n, floor)
Base.size(x::_Floored) = size(x.n)
Base.@propagate_inbounds Base.getindex(x::_Floored, i::Int) = max(x.n[i], x.floor)

# A vector the solute terms write their dual values into, keeping only the
# derivative along the tag `T`: the derivatives of all the solutes in one pass,
# without a vector of dual numbers. Its element type is that of the derivatives
# it keeps; it is written to, never read.
struct _DerivativeSink{T, E, V <: AbstractVector{E}} <: AbstractVector{E}
    d::V
end
_DerivativeSink{T}(d::AbstractVector) where {T} = _DerivativeSink{T, eltype(d), typeof(d)}(d)
Base.size(s::_DerivativeSink) = size(s.d)
Base.@propagate_inbounds function Base.setindex!(s::_DerivativeSink{T}, v, i::Int) where {T}
    s.d[i] = ForwardDiff.extract_derivative(T, v)
    return s
end
