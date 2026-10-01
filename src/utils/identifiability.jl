# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using LinearAlgebra
using ForwardDiff

# ── Which parameters a measurement can actually determine ────────────────────
#
# Fitting six numbers to a curve and reporting six numbers are different acts.
# `scripts/hydration_calibration.jl` worked this out for calorimetry and stated
# the result plainly: over six candidates, the standard errors of the singular
# directions of `∂Q/∂log θ` say that the measurement determines three
# COMBINATIONS and not six numbers. What is here is that reasoning taken out of
# one script and made to work on any forward model.

"""
    log_sensitivity(forward, θ) -> Matrix

`∂y/∂log θⱼ = θⱼ ∂y/∂θⱼ`, with `y = forward(θ)`, exact, by forward-mode
differentiation.

Differentiating with respect to the **logarithm** is what makes the columns
comparable: a rate constant and a dimensionless exponent have no common unit,
and a matrix mixing `∂y/∂k` with `∂y/∂n` has a singular spectrum that says more
about the units than about the data.

`forward` must accept `ForwardDiff` dual numbers and return a vector. Every
forward model of this package does, a kinetic run included: an equilibrium
inside it is differentiated by the implicit-function theorem at its certified
answer. The cost is one evaluation of `forward` on dual numbers per chunk of
parameters (up to twelve at a time), and the derivatives carry no step: an
exact degeneracy between two parameters shows as a singular value at the
rounding of the computation, not as a condition number of a few hundred.

The logarithm is undefined at zero, and so is this: rescale or shift a
parameter that can vanish first.

`relstep` set the difference step until 0.28.2. It is accepted and ignored.
"""
function log_sensitivity(forward, θ; relstep = nothing)
    relstep === nothing || Base.depwarn(
        "`relstep` is ignored: `log_sensitivity` differentiates exactly, by " *
            "forward mode, since ChemistryLab 0.28.3.",
        :log_sensitivity,
    )
    p = collect(float.(θ))
    any(iszero, p) && throw(
        ArgumentError(
            "the logarithm this differentiates against is undefined at zero. " *
                "Rescale or shift the parameter first."
        ),
    )
    # `ForwardDiff.jacobian` tags its duals with the function it differentiates,
    # so a forward model that differentiates internally, or a caller who
    # differentiates this, keeps its own perturbations apart.
    J = ForwardDiff.jacobian(q -> collect(forward(q)), p)
    return J .* transpose(p)
end

"""
    struct Identifiability

What a measurement can and cannot determine about a parameter vector, at one
point.

# Fields

  - `J`: the log-sensitivity matrix `∂y/∂log θ`.
  - `S`, `U`, `V`: its singular value decomposition. `V[:, k]` is the parameter
    combination the k-th singular value belongs to.
  - `condition`: `S[1]/S[end]`, or `Inf` when the problem is exactly
    rank-deficient — see below.
  - `rank`: how many **directions** the data constrain — see
    [`identifiable_rank`](@ref): with an observation, those whose standard error
    in `log θ` is below `tol`; without one, those before the largest gap in the
    spectrum that exceeds `gap`. It is a count of directions and not of
    parameters; [`null_participation`](@ref) is what says which parameters those
    directions leave undetermined.
  - `correlation`: the parameter correlation matrix from `(JᵀJ)⁻¹`.
  - `stderr`: approximate **relative** standard errors, `σ√diag((JᵀJ)⁻¹)`, or
    `nothing` when no observation was given to get `σ` from.
  - `rmse`: the residual root-mean-square, or `nothing`.
  - `names`: parameter names, for reading the output.

These are **linearized** errors at one point. They say which numbers in a fit
deserve to be quoted; they are not confidence intervals, and a model that is
nonlinear in its parameters — which is most of them — will have a likelihood
that is not the ellipse this describes.

# Fewer observations than parameters

The problem is then exactly rank-deficient, and that case is reported rather
than smoothed over: `condition` is `Inf`, and `null_participation` accounts for
the directions the data cannot see at all — not only the ones they see badly.
`correlation` and `stderr`, on the other hand, come from a pseudo-inverse of a
singular `JᵀJ`, so they are finite where the truth is not, and **understate**
the error on any parameter with weight in that null space. The participation is
what to read there; it is also what marks those parameters `PROV_PLACEHOLDER` in
[`as_traced`](@ref).
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
    nparameters(id::Identifiability) -> Int

How many parameters were differentiated against.

Not `length(id.S)`: the singular spectrum is as long as the SMALLER of the two
dimensions, so with fewer observations than parameters it is shorter than the
parameter vector. Every count of parameters reads this, and every count of
directions reads `length(id.S)` — conflating the two is what made an
under-determined fit report a shorter answer than it was asked about.
"""
nparameters(id::Identifiability) = size(id.V, 1)

"""
    identifiability(forward, θ; observed = nothing,
                    names = nothing, gap = 5.0, tol = 1.0) -> Identifiability

How much of `θ` the output of `forward` determines.

`observed` turns the linearized errors into numbers: without it there is no
residual to scale them by, and `stderr` comes back `nothing` rather than
pretending. It also decides how the rank is read: by the standard error of each
direction (`tol`) when there is a residual, by the gap in the spectrum (`gap`)
when there is not ([`identifiable_rank`](@ref)). The sensitivity is
[`log_sensitivity`](@ref)'s, exact; `relstep` is accepted and ignored, as there.

# Reading it

A `condition` of a few means every direction is constrained. Several orders of
magnitude means the data determine a **combination** of parameters, and
`V[:, end]` names which one — the entries of that column say which parameters
trade off against each other.

`correlation` is the sharper instrument when two parameters are collinear:
`hydration_calibration.jl` found `k₁` and `n₁` correlated at −0.96, which says
the data see a product and not its factors, and therefore that only one of the
two can be fitted. **Which one to keep is a modeling judgement, not a
statistical one** — there the rate constant was kept and the shape exponent
fixed, because a rate constant is the quantity a different clinker plausibly
changes.

See also: [`identifiable_rank`](@ref), [`as_traced`](@ref).
"""
function identifiability(
        forward, θ; observed = nothing, relstep = nothing,
        names = nothing, gap::Real = 5.0, tol::Real = 1.0,
    )
    J = Matrix{Float64}(log_sensitivity(forward, θ; relstep))
    # `full` only when the thin `V` would be SHORT OF COLUMNS. With fewer
    # observations than parameters the thin factorization returns a `V` of size
    # `n_par × n_obs`, which silently omits the `n_par - n_obs` directions the
    # data cannot see at all — precisely the ones worth naming. Asking for the
    # complete factorization then costs nothing extra, because the `U` that grows
    # with it is `n_obs × n_obs` and `n_obs` is the small dimension in that case.
    F = svd(J; full = size(J, 1) < size(J, 2))
    # The covariance `(JᵀJ)⁻¹` from the factorization, not by inverting `JᵀJ`,
    # which squares the condition number: an exact degeneracy, which exact
    # derivatives give as one, has `cond(J)` near 1e10, `JᵀJ` beyond the
    # arithmetic, and an inverse whose correlations then have no sign. A
    # direction the spectrum does not carry (fewer observations than parameters)
    # or carries at an exact zero is left out, as the pseudo-inverse would.
    keep = [k for k in eachindex(F.S) if F.S[k] > 0]
    Vk = F.V[:, keep]
    C = Vk * Diagonal(inv.(F.S[keep] .^ 2)) * transpose(Vk)
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
        # An exactly rank-deficient problem has an infinite condition number, and
        # `S[1] / S[end]` over the singular values that EXIST would report a
        # finite one — the structurally missing directions are the worst ones.
        size(F.V, 2) > length(F.S) ? Inf : F.S[1] / max(F.S[end], eps()),
        identifiable_rank(F.S; gap, rmse, tol, floor = _rank_floor(J, F.S)),
        correlation, stderr, rmse,
        names === nothing ? ["θ$i" for i in eachindex(θ)] : collect(names),
    )
end

"""
    identifiable_rank(S; gap = 5.0, rmse = nothing, tol = 1.0, floor = 0.0) -> Int

How many directions a singular spectrum `S` of `∂y/∂log θ` constrains.

With the residual `rmse` of a fit, a direction is constrained when the data fix
it to better than a factor `exp(tol)`: its standard error in `log θ`,
`rmse / S[k]`, is below `tol`. That is a ratio of two quantities in the units of
`y`, so it has none. A singular value at or below `floor`, the rounding of the
factorization, is no direction at all, whatever the residual.

Without a residual, the noise is not known, and the rank is read off the
spectrum alone: the number of singular values before its largest **ratio**
exceeding `gap`. A threshold would have units; a gap does not.

# Why the residual decides when there is one

A gap separates singular values, not what the data see from what they do not:
whether a direction is seen depends on the noise, which the spectrum does not
carry. Two cases of this repository part the two rules:

  - `scripts/hydration_calibration.jl`, six rate parameters against one heat
    curve. The spectrum is `[424, 108, 63.6, 6.46, 1.98, 0.152]` and the residual
    26 J/g. The standard errors of the six directions are then 0.06, 0.24, 0.41,
    4.0, 13 and 170 in `log θ`: three combinations are determined, and the other
    three are known to within factors of 55 and more. The largest ratio, 13, is
    between the fifth singular value and the sixth, so the gap answers five.
  - Three separated peaks of a synthetic thermogram (`test/thermogram.jl`), a
    midpoint of 400 to 950 K and a width of 12 to 20 K each, fitted exactly.
    Every direction is determined. A log-sensitivity carries the size of its
    parameter, so the three midpoints come out 15 to 33 times above the three
    widths, and the largest ratio, 5.8, falls between the two groups: the gap
    answers three.

The default `gap` of 5 was set on the first spectrum when it was measured by
central differences, `[420, 100, 60, 6.3, 1.4, 0.20]`, whose largest ratio, 9.5,
was between the third and the fourth. Exact derivatives moved the smallest
singular value by a quarter, enough to move the largest ratio, and the rule with
it.

# Which gap, when there are several

Without a residual, the cut is at the **largest** qualifying ratio, not the
first one. A rate law measured against its own shrinking-core exponent has a
ratio of 6 between its first two singular values and a third at the rounding of
the arithmetic: only the sum of two exponents is visible. Cutting at the first
ratio would say that the amplitude alone is determined, when the amplitude and
the sum are, which is **two**. A factor of six is ordinary conditioning; a
singular value at the rounding is a structure.

Returns `length(S)` when no ratio exceeds `gap`, which is the honest answer for
a flat spectrum: everything is constrained, or nothing distinguishes what is
not.
"""
function identifiable_rank(
        S::AbstractVector; gap::Real = 5.0, rmse = nothing, tol::Real = 1.0, floor::Real = 0.0,
    )
    isempty(S) && return 0
    rmse === nothing || return count(s -> s > floor && rmse < tol * s, S)
    best, cut = zero(float(gap)), 0
    for k in 1:(length(S) - 1)
        # A singular value at or below zero is a direction that does not exist,
        # and no finite ratio describes it: cut there and stop.
        S[k + 1] <= 0 && return k
        r = S[k] / S[k + 1]
        if r > gap && r > best
            best, cut = r, k
        end
    end
    return cut == 0 ? length(S) : cut
end

identifiable_rank(id::Identifiability; gap::Real = 5.0, tol::Real = 1.0) =
    identifiable_rank(id.S; gap, rmse = id.rmse, tol, floor = _rank_floor(id.J, id.S))

# The rounding of a singular value decomposition of `J`: the tolerance of
# `LinearAlgebra.rank`.
_rank_floor(J, S) = maximum(size(J)) * eps(isempty(S) ? 0.0 : float(S[1]))

"""
    null_participation(id::Identifiability) -> Vector{Float64}

For each parameter, how much of it lies in the directions the data do **not**
constrain: `Σ_{k > rank} V[i,k]²`, which is between 0 and 1.

# Why this and not "beyond the rank"

A rank of 2 out of 3 says the data constrain two *directions in parameter
space*. It says nothing about which *parameters* are determined, because the
directions are the singular vectors and the parameter order is whatever the
caller packed. Reading the rank as "the first two parameters are fine" confuses
an index with a subspace, and for a model where two parameters trade off it
picks the wrong one — the one that happened to be listed second.

The participation is the honest version. In a model where `a` and `c` enter only
as their product, it comes out `[0.5, 0, 0.5]`: neither `a` nor `c` is
determined on its own, `b` is, and no ordering was consulted to say so.
"""
function null_participation(id::Identifiability)
    n = nparameters(id)
    id.rank >= size(id.V, 2) && return zeros(Float64, n)
    return [sum(abs2, @view id.V[i, (id.rank + 1):size(id.V, 2)]) for i in 1:n]
end

"""
    as_traced(id::Identifiability, θ; source, null_threshold = 0.1)
        -> Vector{Traced}

The fitted parameters, each carrying `PROV_FITTED`, the dataset it was fitted
to, and its linearized standard error as [`uncertainty`](@ref).

This is where an identification meets [`Traced`](@ref): a number that came out
of an optimization is not the same kind of claim as one somebody measured, and
`is_evidence` returns `false` for it on purpose.

A parameter with more than `null_threshold` of its weight in the directions the
data do not constrain is marked `PROV_PLACEHOLDER` instead — a value the fit had
to leave somewhere is not a value the fit determined. That test is
[`null_participation`](@ref) and **not** the parameter's position relative to
the rank: the rank counts directions, the position is the caller's packing
order, and confusing the two flags whichever of a trading-off pair was listed
second. Pass `null_threshold = Inf` to report them all as fitted and take the
responsibility.
"""
function as_traced(
        id::Identifiability, θ; source::AbstractString, null_threshold::Real = 0.1,
    )
    σ = id.stderr
    p = null_participation(id)
    return [
        Traced(
            float(θ[i]),
            p[i] > null_threshold ? PROV_PLACEHOLDER : PROV_FITTED,
            source;
            uncertainty = σ === nothing ? nothing : σ[i],
        ) for i in eachindex(θ)
    ]
end

function Base.show(io::IO, ::MIME"text/plain", id::Identifiability)
    np = nparameters(id)
    println(io, "Identifiability of ", np, " parameters")
    println(io, "  singular values  ", round.(id.S; sigdigits = 3))
    println(io, "  condition        ", round(id.condition; sigdigits = 4))
    println(io, "  constrained      ", id.rank, " of ", size(id.V, 2), " directions")
    id.rmse === nothing || println(io, "  residual RMSE    ", round(id.rmse; sigdigits = 4))
    worst, pair = 0.0, (0, 0)
    for i in 1:np, j in (i + 1):np
        abs(id.correlation[i, j]) > worst &&
            ((worst, pair) = (abs(id.correlation[i, j]), (i, j)))
    end
    return iszero(pair[1]) ? nothing : print(
            io, "  most collinear   ", id.names[pair[1]], " / ", id.names[pair[2]],
            "  r = ", round(id.correlation[pair[1], pair[2]]; digits = 3),
        )
end

Base.show(io::IO, id::Identifiability) =
    print(io, "Identifiability(", nparameters(id), " parameters, rank ", id.rank, ")")
