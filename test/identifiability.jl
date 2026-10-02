# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using ChemistryLab
using ChemistryLab: value
using ForwardDiff
using LinearAlgebra
using Test

@testsection "which parameters a measurement can determine" begin

    # A model with a COLLINEARITY PUT IN ON PURPOSE, so the method has a known
    # right answer rather than a plausible one: `a` and `c` enter only as their
    # product, so no measurement of `y` can separate them, whatever the data.
    _t = range(0, 5; length = 40)
    _decay(p) = @. p[1] * p[3] * exp(-p[2] * _t)
    _θ = [2.0, 0.7, 1.0]

    @testset "an exact collinearity is found, and named" begin
        id = identifiability(_decay, _θ; names = ["a", "b", "c"], observed = _decay(_θ))

        # Two directions carried, one empty — the spectrum falls off a cliff.
        @test id.rank == 2
        @test identifiable_rank(id) == id.rank
        @test id.S[1] / id.S[2] < 10
        @test id.S[3] / id.S[2] < 1.0e-10
        @test id.condition > 1.0e10
        # The fit is exact and leaves a residual of zero, which is no noise
        # level: σ is the floor, √eps times the root-mean-square of the curve,
        # and the direction at the rounding stays out.
        y = _decay(_θ)
        @test id.rmse == 0
        @test id.noise ≈ sqrt(eps()) * sqrt(sum(abs2, y) / length(y))

        # And the empty direction is the a-versus-c trade-off, which is the
        # answer put in: equal and opposite in `a` and `c`, nothing in `b`.
        v = id.V[:, end]
        @test abs(v[2]) < 1.0e-8
        @test abs(abs(v[1]) - abs(v[3])) < 1.0e-8
        @test v[1] * v[3] < 0

        # The correlation says the same thing more directly, which is why it is
        # the sharper instrument when two parameters trade off.
        # A trade-off: the two move against each other.
        @test id.correlation[1, 3] < -0.999
        @test abs(id.correlation[1, 2]) < 0.999
        @test occursin("a / c", sprint(show, MIME"text/plain"(), id))
    end

    @testset "a well-conditioned case says so" begin
        # Drop the redundant parameter and everything is determined.
        f(p) = @. p[1] * exp(-p[2] * _t)
        id = identifiability(f, [2.0, 0.7]; names = ["a", "b"])
        @test id.rank == 2
        @test id.condition < 100
        # No observation was given, so there is no residual to scale errors by
        # and it says `nothing` rather than inventing one.
        @test id.stderr === nothing
        @test id.rmse === nothing
    end

    @testset "the noise level is the instrument's, or the fit's on n − p" begin
        f(q) = @. q[1] * exp(-q[2] * _t)
        θ = [2.0, 0.7]
        y = f(θ)
        # A residual of known size, on 40 − 2 degrees of freedom.
        r = 1.0e-3 .* cos.(3 .* _t)
        id = identifiability(f, θ; observed = y .- r)
        @test id.rmse ≈ sqrt(sum(abs2, r) / 40) rtol = 1.0e-12
        @test id.noise ≈ sqrt(sum(abs2, r) / 38) rtol = 1.0e-12
        @test id.stderr ≈ id.noise .* sqrt.(diag(inv(transpose(id.J) * id.J))) rtol = 1.0e-8
        # Given, it takes precedence over the residual and decides the rank: a
        # noise level above every singular value leaves nothing determined.
        loud = identifiability(f, θ; observed = y .- r, noise = 2 * id.S[1])
        @test loud.noise == 2 * id.S[1]
        @test loud.rank == 0
        @test loud.stderr ≈ (loud.noise / id.noise) .* id.stderr
        @test occursin("noise level σ", sprint(show, MIME"text/plain"(), loud))
        # Without an observation or a noise level, the gap, and the report says so.
        bare = identifiability(f, θ)
        @test bare.noise === nothing
        @test occursin("by the gap", sprint(show, MIME"text/plain"(), bare))
        # A noise level without an observation still gives the errors.
        quiet = identifiability(f, θ; noise = 1.0e-3)
        @test quiet.rmse === nothing
        @test quiet.stderr ≈ (1.0e-3 / id.noise) .* id.stderr
    end

    @testset "a rank is read off the standard errors, or else a gap" begin
        # The spectrum `hydration_calibration.jl` measures over six candidates,
        # with the noise level of that fit (J/g): the residual of 25.94 J/g over
        # thirty instants, on the 24 degrees of freedom six parameters leave. The
        # standard errors of the six directions are 0.07, 0.27, 0.46, 4.5, 15 and
        # 190 in log θ: three combinations are determined, the others to within
        # factors of 89 and more.
        S = [424.41, 107.88, 63.619, 6.4596, 1.9803, 0.15242]
        σ = 25.9354 * sqrt(30 / 24)
        @test identifiable_rank(S; noise = σ) == 3
        # The spectrum alone does not say where the noise is, and its largest
        # ratio is between the fifth and the sixth.
        @test identifiable_rank(S) == 5
        @test identifiable_rank(S; gap = 20.0) == 6     # no gap that large
        # The tolerance is a factor in θ: at e⁵ the fourth direction is in.
        @test identifiable_rank(S; noise = σ, tol = 5.0) == 4
        # A noise level is a standard deviation, and a zero one is no
        # measurement: an exact fit gets its floor from `identifiability`.
        @test_throws ArgumentError identifiable_rank(S; noise = 0.0)
        @test_throws ArgumentError identifiability(x -> x, [1.0]; noise = -1.0)
        # A flat spectrum constrains everything, and saying `length` is the
        # honest answer rather than a smaller number chosen to look careful.
        @test identifiable_rank([1.0, 0.9, 0.8]; gap = 10.0) == 3
        @test identifiable_rank([1.0, 0.0]; gap = 10.0) == 1
        @test identifiable_rank(Float64[]) == 0
        @test identifiable_rank(Float64[]; noise = 1.0) == 0

        # THE LARGEST QUALIFYING GAP, NOT THE FIRST. Measured spectrum of a rate
        # law against its own shrinking-core exponent (the manual page "Where
        # the numbers come from"): a ratio of 6.0 and then one of a billion.
        # Cutting at the first answers one determined direction — the amplitude
        # alone — when the amplitude AND the sum of the two exponents are
        # determined and only their split is not.
        real = [1.6187e-5, 2.7103e-6, 1.5616e-15]
        @test real[1] / real[2] > 5                    # the 6.0 does qualify
        @test identifiable_rank(real) == 2             # and is not where the cut is
        # Order matters only through the size of the ratio, so a spectrum whose
        # cliff comes first still cuts there.
        @test identifiable_rank([1.0e-5, 1.0e-10, 2.0e-11]) == 1
    end

    @testset "fewer observations than parameters is reported, not smoothed" begin
        # Two numbers measured, three asked about: one direction is invisible by
        # COUNTING, before any conditioning argument. The thin factorization
        # `svd` returns by default has a `V` of size 3×2 and does not carry that
        # direction at all, so reading the parameter count off the singular
        # spectrum answered about two parameters when asked about three.
        g(p) = [p[1] + p[2] + p[3], p[1] - p[2]]
        id = identifiability(g, [1.0, 1.0, 1.0]; names = ["a", "b", "c"])

        @test length(id.S) == 2              # directions the data can rank
        @test size(id.V) == (3, 3)           # but three parameters to place
        @test id.rank == 2
        @test id.condition == Inf            # exactly rank-deficient, and says so

        # The invisible direction is `(1, 1, -2)/√6`, so `c` carries four sixths
        # of it and `a` and `b` one sixth each. One null direction means the
        # participations sum to exactly one.
        q = null_participation(id)
        @test length(q) == 3
        @test sum(q) ≈ 1.0 atol = 1.0e-10
        @test q ≈ [1 / 6, 1 / 6, 2 / 3] atol = 1.0e-8

        # All three are above the default threshold, so none is reported as
        # determined — which is the honest answer to two equations in three
        # unknowns.
        out = as_traced(id, [1.0, 1.0, 1.0]; source = "two measurements")
        @test length(out) == 3
        @test all(t -> provenance(t) === PROV_PLACEHOLDER, out)

        # And the printed form counts parameters, not singular values.
        @test occursin("Identifiability of 3 parameters", sprint(show, MIME"text/plain"(), id))
        @test occursin("3 parameters", sprint(show, id))
    end

    @testset "the sensitivity is exact, and needs a nonzero parameter" begin
        # Refused rather than silently differentiating against nothing: the
        # logarithm is undefined at zero.
        @test_throws ArgumentError log_sensitivity(_decay, [2.0, 0.0, 1.0])

        # Exact, by forward mode: `θⱼ ∂y/∂θⱼ` in closed form.
        J = log_sensitivity(_decay, _θ)
        a, b, c = _θ
        e = exp.(-b .* _t)
        @test J[:, 1] ≈ a * c .* e rtol = 1.0e-14
        @test J[:, 2] ≈ -a * b * c .* _t .* e rtol = 1.0e-14
        @test J[:, 3] ≈ a * c .* e rtol = 1.0e-14

        # So the exact collinearity of `a` and `c` gives two IDENTICAL columns,
        # and a singular value at the rounding of the computation.
        @test J[:, 1] == J[:, 3]

        # One evaluation of the forward model, on dual numbers, for up to twelve
        # parameters; and `identifiability` one more only when it has an
        # observation to form a residual against.
        calls = Ref(0)
        counted(q) = (calls[] += 1; _decay(q))
        log_sensitivity(counted, _θ)
        @test calls[] == 1
        calls[] = 0
        identifiability(counted, _θ)
        @test calls[] == 1
        calls[] = 0
        identifiability(counted, _θ; observed = _decay(_θ))
        @test calls[] == 2

        # The step of the differences it replaced is accepted and ignored.
        @test (@test_deprecated log_sensitivity(_decay, _θ; relstep = 0.05)) == J

        # And differentiated again, it gives the exact second derivative: the
        # tags keep the two levels apart.
        d2 = ForwardDiff.derivative(β -> log_sensitivity(_decay, [a, β, c])[end, 2], b)
        t = _t[end]
        @test d2 ≈ -a * c * exp(-b * t) * t * (1 - b * t) rtol = 1.0e-12
    end

    @testset "a fitted number is not a measured one" begin
        id = identifiability(_decay, _θ; names = ["a", "b", "c"], observed = _decay(_θ))
        out = as_traced(id, _θ; source = "a synthetic fit")
        @test length(out) == 3
        @test all(!is_evidence, out)             # fitted is never evidence
        @test value(out[1]) == _θ[1]
        @test uncertainty(out[1]) !== nothing

        # WHICH PARAMETERS ARE UNDETERMINED IS A SUBSPACE QUESTION, not a
        # question of position. `a` and `c` enter only as their product, so
        # NEITHER is determined alone and `b` is — and reading "beyond the
        # rank 2" as "the third one" would have flagged `c` and cleared `a`,
        # which is the wrong answer arrived at by a plausible route.
        p = null_participation(id)
        @test p[1] ≈ 0.5 atol = 1.0e-6
        @test p[2] < 1.0e-8
        @test p[3] ≈ 0.5 atol = 1.0e-6
        @test provenance(out[1]) === PROV_PLACEHOLDER
        @test provenance(out[2]) === PROV_FITTED
        @test provenance(out[3]) === PROV_PLACEHOLDER

        # Unless the caller takes the responsibility explicitly.
        all_fitted = as_traced(id, _θ; source = "a synthetic fit", null_threshold = Inf)
        @test all(t -> provenance(t) === PROV_FITTED, all_fitted)

        # And a well-conditioned problem has nothing in the null space at all.
        f(q) = @. q[1] * exp(-q[2] * _t)
        clean = identifiability(f, [2.0, 0.7]; names = ["a", "b"])
        @test all(iszero, null_participation(clean))
        @test all(
            t -> provenance(t) === PROV_FITTED,
            as_traced(clean, [2.0, 0.7]; source = "a clean fit"),
        )

        # And the report over the lot says what it is made of.
        r = provenance_report(out)
        @test r.weakest === PROV_PLACEHOLDER
        @test !r.all_evidence
        @test r.counts[PROV_FITTED] == 1
        @test r.counts[PROV_PLACEHOLDER] == 2
    end
end
