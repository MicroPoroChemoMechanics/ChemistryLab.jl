# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using LinearAlgebra

# ── Which parameters a measurement can actually determine ────────────────────
#
# Fitting six numbers to a curve and reporting six numbers are different acts.
# `scripts/hydration_calibration.jl` worked this out for calorimetry and stated
# the result plainly: over six candidates the singular values of `∂Q/∂log θ` came
# out at `[420, 100, 60, 6.3, 1.4, 0.20]`, a factor of nine between the third and
# the fourth, so the measurement determined three COMBINATIONS and not six
# numbers. What is here is that reasoning taken out of one script and made to
# work on any forward model.

"""
    log_sensitivity(forward, θ; relstep = 0.05) -> Matrix

`∂y/∂log θⱼ` by central differences, with `y = forward(θ)`.

Differentiating with respect to the **logarithm** is what makes the columns
comparable: a rate constant and a dimensionless exponent have no common unit,
and a matrix mixing `∂y/∂k` with `∂y/∂n` has a singular spectrum that says more
about the units than about the data.

Central differences rather than automatic differentiation, because a forward
model is often a solver whose parameters do not carry duals — `hydration_calibration.jl`
names exactly why for its own: the integrator casts to `Float64` on the way in.
`forward` is free to be AD-clean; this does not require it, and costs `2n`
evaluations.

`relstep` is a **relative** step, so a parameter near zero needs rescaling before
it is passed here — which is the same condition as the logarithm being defined.
"""
function log_sensitivity(forward, θ; relstep::Real = 0.05)
    p = collect(float.(θ))
    any(iszero, p) && throw(
        ArgumentError(
            "a relative step is undefined at zero, and so is the logarithm this " *
                "differentiates against. Rescale or shift the parameter first."
        ),
    )
    y0 = forward(p)
    J = Matrix{Float64}(undef, length(y0), length(p))
    for j in eachindex(p)
        up = copy(p); up[j] *= (1 + relstep)
        dn = copy(p); dn[j] *= (1 - relstep)
        J[:, j] = (forward(up) .- forward(dn)) ./ (2 * relstep)
    end
    return J
end

"""
    struct Identifiability

What a measurement can and cannot determine about a parameter vector, at one
point.

# Fields

  - `J`: the log-sensitivity matrix `∂y/∂log θ`.
  - `S`, `U`, `V`: its singular value decomposition. `V[:, k]` is the parameter
    combination the k-th singular value belongs to.
  - `condition`: `S[1]/S[end]`.
  - `rank`: how many directions the data constrain, by the largest gap in the
    spectrum — see [`identifiable_rank`](@ref).
  - `correlation`: the parameter correlation matrix from `(JᵀJ)⁻¹`.
  - `stderr`: approximate **relative** standard errors, `σ√diag((JᵀJ)⁻¹)`, or
    `nothing` when no observation was given to get `σ` from.
  - `rmse`: the residual root-mean-square, or `nothing`.
  - `names`: parameter names, for reading the output.

These are **linearized** errors at one point. They say which numbers in a fit
deserve to be quoted; they are not confidence intervals, and a model that is
nonlinear in its parameters — which is most of them — will have a likelihood
that is not the ellipse this describes.
"""
struct Identifiability{T}
    J::Matrix{Float64}
    U::Matrix{Float64}
    S::Vector{Float64}
    V::Matrix{Float64}
    condition::Float64
    rank::Int
    correlation::Matrix{Float64}
    stderr::Union{Nothing, Vector{Float64}}
    rmse::Union{Nothing, Float64}
    names::T
end

"""
    identifiability(forward, θ; observed = nothing, relstep = 0.05,
                    names = nothing, gap = 10.0) -> Identifiability

How much of `θ` the output of `forward` determines.

`observed` turns the linearized errors into numbers: without it there is no
residual to scale them by, and `stderr` comes back `nothing` rather than
pretending.

# Reading it

A `condition` of a few means every direction is constrained. Several orders of
magnitude means the data determine a **combination** of parameters, and
`V[:, end]` names which one — the entries of that column say which parameters
trade off against each other.

`correlation` is the sharper instrument when two parameters are collinear:
`hydration_calibration.jl` found `k₁` and `n₁` correlated at −0.985, which says
the data see a product and not its factors, and therefore that only one of the
two can be fitted. **Which one to keep is a modeling judgement, not a
statistical one** — there the rate constant was kept and the shape exponent
fixed, because a rate constant is the quantity a different clinker plausibly
changes.

See also: [`identifiable_rank`](@ref), [`as_traced`](@ref).
"""
function identifiability(
        forward, θ; observed = nothing, relstep::Real = 0.05,
        names = nothing, gap::Real = 10.0,
    )
    J = log_sensitivity(forward, θ; relstep)
    F = svd(J)
    JtJ = J' * J
    C = try
        inv(JtJ)
    catch
        pinv(JtJ)
    end
    d = sqrt.(abs.(diag(C)))
    correlation = C ./ (d * d')
    rmse = if observed === nothing
        nothing
    else
        # Written out rather than reaching for `Statistics` — one mean is not a
        # dependency's worth of surface.
        r = forward(collect(float.(θ))) .- observed
        sqrt(sum(abs2, r) / length(r))
    end
    stderr = rmse === nothing ? nothing : rmse .* d
    return Identifiability(
        J, Matrix(F.U), Vector(F.S), Matrix(F.V),
        F.S[1] / max(F.S[end], eps()),
        identifiable_rank(F.S; gap),
        correlation, stderr, rmse,
        names === nothing ? ["θ$i" for i in eachindex(θ)] : collect(names),
    )
end

"""
    identifiable_rank(S; gap = 10.0) -> Int

How many directions a singular spectrum constrains: the number of singular
values before the first **ratio** larger than `gap`.

A rank is read off a gap rather than a threshold because a threshold has units
and a gap does not. `[420, 100, 60, 6.3, 1.4, 0.20]` has its largest ratio
between the third and the fourth — a factor of nine at `gap = 5`, and the reason
that calibration fits three parameters and not six.

Returns `length(S)` when no gap exceeds `gap`, which is the honest answer when
the spectrum is flat: everything is constrained, or nothing distinguishes what
is not.
"""
function identifiable_rank(S::AbstractVector; gap::Real = 10.0)
    isempty(S) && return 0
    for k in 1:(length(S) - 1)
        S[k + 1] <= 0 && return k
        S[k] / S[k + 1] > gap && return k
    end
    return length(S)
end

identifiable_rank(id::Identifiability; gap::Real = 10.0) =
    identifiable_rank(id.S; gap)

"""
    as_traced(id::Identifiability, θ; source, rank_limit = true) -> Vector{Traced}

The fitted parameters, each carrying `PROV_FITTED`, the dataset it was fitted
to, and its linearized standard error as [`uncertainty`](@ref).

This is where an identification meets [`Traced`](@ref): a number that came out
of an optimization is not the same kind of claim as one somebody measured, and
`is_evidence` returns `false` for it on purpose.

With `rank_limit`, parameters beyond `id.rank` are marked `PROV_PLACEHOLDER`
instead — because a number the data did not constrain is a value the fit had to
put somewhere, not a value the fit determined. Pass `false` to report them all
as fitted and take the responsibility.
"""
function as_traced(
        id::Identifiability, θ; source::AbstractString, rank_limit::Bool = true,
    )
    σ = id.stderr
    return [
        Traced(
            float(θ[i]),
            (rank_limit && i > id.rank) ? PROV_PLACEHOLDER : PROV_FITTED,
            source;
            uncertainty = σ === nothing ? nothing : σ[i],
        ) for i in eachindex(θ)
    ]
end

function Base.show(io::IO, ::MIME"text/plain", id::Identifiability)
    println(io, "Identifiability of ", length(id.S), " parameters")
    println(io, "  singular values  ", round.(id.S; sigdigits = 3))
    println(io, "  condition        ", round(id.condition; sigdigits = 4))
    println(io, "  constrained      ", id.rank, " of ", length(id.S), " directions")
    id.rmse === nothing || println(io, "  residual RMSE    ", round(id.rmse; sigdigits = 4))
    worst, pair = 0.0, (0, 0)
    for i in eachindex(id.S), j in (i + 1):length(id.S)
        abs(id.correlation[i, j]) > worst &&
            ((worst, pair) = (abs(id.correlation[i, j]), (i, j)))
    end
    return iszero(pair[1]) ? nothing : print(
            io, "  most collinear   ", id.names[pair[1]], " / ", id.names[pair[2]],
            "  r = ", round(id.correlation[pair[1], pair[2]]; digits = 3),
        )
end

Base.show(io::IO, id::Identifiability) =
    print(io, "Identifiability(", length(id.S), " parameters, rank ", id.rank, ")")
