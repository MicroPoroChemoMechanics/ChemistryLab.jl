# ── The linear-programming initial approximation ─────────────────────────────
#
# A cold cement is the hard case for any Gibbs minimization: from the cast state
# — all the mass in the reactants, every product at the floor — the answer is
# tens of solids away, and an interior-point method spends its whole budget
# discovering WHICH solids before it can start refining how much of each.
#
# The observation that fixes it is not ours. GEM-Selektor's automatic initial
# approximation solves a **simplified problem first**: drop the mixing entropy,
# keep only the linear part of the Gibbs energy, and minimize
#
#     min  gᵀn    subject to    A n = b,   n ≥ 0
#
# by the simplex method [Kulik2013]. That is a linear program, its optimum is a
# **vertex**, and a vertex of this polyhedron carries at most one nonzero per
# component — which is to say it IS a phase assemblage. Handed to the nonlinear
# solver, it replaces "discover the assemblage" with "refine this one".
#
# The approximation is crude by construction: without the mixing entropy no
# solution phase can exist, so the vertex puts pure end-members where the answer
# has a C-S-H. That does not matter. It is a starting point, and the certificate
# still judges the answer.

"""
    lp_vertex(A, b, g; maxit, tol) -> Union{Vector{Float64}, Nothing}

A vertex of `{n ≥ 0 : A n = b}` minimizing `gᵀn`, by a two-phase revised
simplex, or `nothing` if the polyhedron is empty or the iteration does not
finish.

Bland's rule is used throughout rather than the steepest-descent rule: it is
slower per solve and it **cannot cycle**, and a starting point that occasionally
fails to terminate is worse than one that is occasionally suboptimal.

Rows of `A` that are linearly dependent are tolerated: phase one drives the
artificial variables of redundant rows to zero and leaves them in the basis at
zero, which is harmless here because the result is only a starting point.
"""
function lp_vertex(
        A::AbstractMatrix{<:Real}, b::AbstractVector{<:Real},
        g::AbstractVector{<:Real}; maxit::Int = 10_000, tol::Float64 = 1.0e-9,
    )
    m, n = size(A)
    length(b) == m || throw(DimensionMismatch("b has length $(length(b)), expected $m"))
    length(g) == n || throw(DimensionMismatch("g has length $(length(g)), expected $n"))

    # Phase one needs `b ≥ 0`, so flip the rows that are not. A cement's H⁺ row
    # is routinely negative, which is why this is not a formality.
    #
    # `_primal` strips dual numbers here rather than propagating them, and that
    # is deliberate: this returns a STARTING POINT. A vertex is piecewise
    # constant in the data, so its derivative is zero almost everywhere and
    # undefined on the pivots — differentiating it would be meaningless, and
    # carrying `Dual`s through a simplex tableau would be expensive for nothing.
    # The sensitivity of the ANSWER is obtained where it belongs, by
    # differentiating the converged equilibrium.
    Aw = Matrix{Float64}(_primal.(A))
    bw = Vector{Float64}(_primal.(b))
    @inbounds for i in 1:m
        if bw[i] < 0
            bw[i] = -bw[i]
            @views Aw[i, :] .= .-Aw[i, :]
        end
    end

    # Tableau: [Aw I | bw], the identity being the artificial variables.
    tab = [Aw Matrix{Float64}(LinearAlgebra.I, m, m) bw]
    basis = collect((n + 1):(n + m))
    ncol = n + m

    function pivot!(r, c)
        p = tab[r, c]
        @views tab[r, :] ./= p
        @inbounds for i in 1:m
            i == r && continue
            f = tab[i, c]
            f == 0 && continue
            @views tab[i, :] .-= f .* tab[r, :]
        end
        basis[r] = c
        return nothing
    end

    # Cost row of a phase, expressed in the current basis.
    function reduced_costs(c)
        z = zeros(Float64, ncol + 1)
        @inbounds for i in 1:m
            ci = basis[i] <= length(c) ? c[basis[i]] : 0.0
            ci == 0 && continue
            @views z .+= ci .* tab[i, :]
        end
        return z
    end

    function run_phase!(c, allowed)
        for _ in 1:maxit
            z = reduced_costs(c)
            # Bland: the LOWEST-INDEX improving column, which is what forbids
            # cycling.
            enter = 0
            @inbounds for j in 1:ncol
                allowed(j) || continue
                cj = j <= length(c) ? c[j] : 0.0
                if cj - z[j] < -tol
                    enter = j
                    break
                end
            end
            enter == 0 && return true              # optimal
            # Ratio test, ties broken by the lowest basis index (Bland again).
            leave, best = 0, Inf
            @inbounds for i in 1:m
                aij = tab[i, enter]
                aij > tol || continue
                r = tab[i, end] / aij
                if r < best - tol || (r < best + tol && leave != 0 && basis[i] < basis[leave])
                    best, leave = r, i
                end
            end
            leave == 0 && return false             # unbounded
            pivot!(leave, enter)
        end
        return false                               # ran out of iterations
    end

    # Phase one: minimize the sum of the artificials.
    c1 = [zeros(Float64, n); ones(Float64, m)]
    run_phase!(c1, _ -> true) || return nothing
    infeas = sum(i -> basis[i] > n ? tab[i, end] : 0.0, 1:m)
    infeas > 1.0e-7 * max(1.0, maximum(abs, bw; init = 1.0)) && return nothing

    # Drive whatever artificials remain out of the basis where a real column can
    # replace them; a row that cannot be pivoted is redundant and is left alone.
    @inbounds for i in 1:m
        basis[i] > n || continue
        j = findfirst(k -> abs(tab[i, k]) > tol, 1:n)
        j === nothing || pivot!(i, j)
    end

    # Phase two: the real objective, artificials barred from re-entering.
    c2 = [Vector{Float64}(_primal.(g)); zeros(Float64, m)]
    run_phase!(c2, j -> j <= n) || return nothing

    x = zeros(Float64, n)
    @inbounds for i in 1:m
        basis[i] <= n && (x[basis[i]] = max(tab[i, end], 0.0))
    end
    return x
end

"""
    _lp_start(des, state, b, ϵ, verbose) -> Union{ChemicalState, Nothing}

The vertex of [`lp_vertex`](@ref) as a `ChemicalState`, or `nothing` if the
linear program is infeasible or does not terminate.

Amounts are floored at `ϵ` rather than left at zero: a species at exactly zero
is outside the domain of the log-activities the nonlinear solver evaluates, so a
vertex handed over verbatim would be rejected at its first residual.
"""
function _lp_start(des, state::ChemicalState, b, ϵ::Float64, verbose::Bool)
    return try
        p = _build_params(state; ϵ = ϵ)
        g = Float64.(_primal.(p.ΔₐG⁰overRT))
        x = lp_vertex(des.A, b, g)
        x === nothing && return nothing
        n = max.(x, ϵ)
        ChemicalState(
            state.system, n .* u"mol";
            T = temperature(state), P = pressure(state),
        )
    catch err
        verbose && @info "the LP initial approximation did not run" err
        nothing
    end
end
