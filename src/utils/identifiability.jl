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
`forward` is free to be AD-clean; this does not require it, and costs exactly
`2n` evaluations: the output is sized from the first perturbed call rather than
from an extra unperturbed one. On a forward model that is a solver, that saved
call is a whole solve.

`relstep` is a **relative** step, so a parameter near zero needs rescaling before
it is passed here — which is the same condition as the logarithm being defined.

!!! warning "A coarse step can hide an exact degeneracy"
    The default 5 % is chosen so a noisy forward model still gives a usable
    derivative, and it is coarse. On a rate law that goes as `(1-ξ)^{n₃}`, a 5 %
    step moves `n₃ = 3.3` by 0.165 in the exponent, which is far enough that the
    second-order error differs between two parameters that are **exactly**
    collinear — and the collinearity then shows up as a condition number of 80
    rather than of 3 × 10⁵. Measured on that case, the third singular value goes
    from 2.0e-7 at `relstep = 0.05` to 5.1e-11 at 0.01 and is unchanged below.

    So a condition number of a few hundred is not evidence that a model is well
    posed. **Refine `relstep` and see whether the answer moves**; if it does, the
    coarse one was measuring the differencing and not the model.
"""
function log_sensitivity(forward, θ; relstep::Real = 0.05)
    p = collect(float.(θ))
    any(iszero, p) && throw(
        ArgumentError(
            "a relative step is undefined at zero, and so is the logarithm this " *
                "differentiates against. Rescale or shift the parameter first."
        ),
    )
    # Sized from the first column rather than from an extra unperturbed call. A
    # forward model here is often a solver, so that call is a whole solve, and
    # the first perturbation already says how long the output is.
    J = nothing
    for j in eachindex(p)
        up = copy(p); up[j] *= (1 + relstep)
        dn = copy(p); dn[j] *= (1 - relstep)
        col = (forward(up) .- forward(dn)) ./ (2 * relstep)
        J === nothing && (J = Matrix{Float64}(undef, length(col), length(p)))
        J[:, j] = col
    end
    return J === nothing ? Matrix{Float64}(undef, 0, 0) : J
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
  - `rank`: how many **directions** the data constrain, by the largest gap in
    the spectrum that exceeds `gap` — see [`identifiable_rank`](@ref). It is a
    count of directions
    and not of parameters; [`null_participation`](@ref) is what says which
    parameters those directions leave undetermined.
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
    identifiability(forward, θ; observed = nothing, relstep = 0.05,
                    names = nothing, gap = 5.0) -> Identifiability

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
        names = nothing, gap::Real = 5.0,
    )
    J = log_sensitivity(forward, θ; relstep)
    # `full` only when the thin `V` would be SHORT OF COLUMNS. With fewer
    # observations than parameters the thin factorization returns a `V` of size
    # `n_par × n_obs`, which silently omits the `n_par - n_obs` directions the
    # data cannot see at all — precisely the ones worth naming. Asking for the
    # complete factorization then costs nothing extra, because the `U` that grows
    # with it is `n_obs × n_obs` and `n_obs` is the small dimension in that case.
    F = svd(J; full = size(J, 1) < size(J, 2))
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
        # An exactly rank-deficient problem has an infinite condition number, and
        # `S[1] / S[end]` over the singular values that EXIST would report a
        # finite one — the structurally missing directions are the worst ones.
        size(F.V, 2) > length(F.S) ? Inf : F.S[1] / max(F.S[end], eps()),
        identifiable_rank(F.S; gap),
        correlation, stderr, rmse,
        names === nothing ? ["θ$i" for i in eachindex(θ)] : collect(names),
    )
end

"""
    identifiable_rank(S; gap = 5.0) -> Int

How many directions a singular spectrum constrains: the number of singular
values before the first **ratio** larger than `gap`.

A rank is read off a gap rather than a threshold because a threshold has units
and a gap does not.

# Where the default comes from

From the one case in this repository where the answer is known independently.
`scripts/hydration_calibration.jl` measured `[420, 100, 60, 6.3, 1.4, 0.20]` and
concluded, on the correlation structure, that the data determine three
combinations — and the largest ratio in that spectrum is **9.5**, between the
third singular value and the fourth. A default of 10 would have returned 6 on
the very case the rule exists for, which is how this default came to be 5 rather
than a round number chosen for looking careful.

# Which gap, when there are several

The cut is at the **largest** qualifying ratio, not the first one. That
distinction has a case behind it. Measuring a rate law against its own
shrinking-core exponent gives `[1.62e-5, 2.71e-6, 5.1e-11]`: a ratio of 6.0 and
then one of fifty thousand. The first rule cut at 6.0 and answered **one**
determined direction, which says the amplitude alone is visible; the truth is
that the amplitude and one exponent combination are both determined and only
their split is not, which is **two**. A factor of six is ordinary conditioning.
A factor of fifty thousand is a structure.

Returns `length(S)` when no ratio exceeds `gap`, which is the honest answer for
a flat spectrum: everything is constrained, or nothing distinguishes what is
not.
"""
function identifiable_rank(S::AbstractVector; gap::Real = 5.0)
    isempty(S) && return 0
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

identifiable_rank(id::Identifiability; gap::Real = 5.0) =
    identifiable_rank(id.S; gap)

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
