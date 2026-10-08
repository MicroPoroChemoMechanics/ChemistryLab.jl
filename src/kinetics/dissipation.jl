# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# ── what a rate law did with the second law ─────────────────────────────────

"""
    reaction_affinities(sol, kp::KineticsProblem; times = sol.t, states = nothing)
        -> Matrix{Float64}

The affinity of each kinetic reaction of `kp` (columns) at each instant of
`times` (rows), in J/mol:

```math
\\mathcal{A}_j = -\\sum_i \\nu_{ij}\\,\\mu_i = -\\sum_i \\nu_{ij}
\\left(\\Delta_a G^\\circ_i + RT \\ln a_i\\right) = -RT\\ln\\Omega_j
```

with `ν_ij` the coefficients of the reaction as the problem holds it, the
species it controls at `−1`: a positive affinity drives the reaction the way a
positive rate runs it.

The composition is the certified one, [`speciated_states`](@ref), passed as
`states` when already computed; the activities are those the rate laws of the
run read, at the temperature of each state. A reaction with a participant that
carries no standard Gibbs energy, a glass given a formula for its rate law, has
no affinity: its column is `NaN`.

See also [`dissipation`](@ref).
"""
function reaction_affinities(sol, kp::KineticsProblem; times = sol.t, states = nothing)
    return _affinities_and_rates(sol, kp, times, states; rates = false).A
end

"""
    dissipation(sol, kp::KineticsProblem; times = sol.t, states = nothing,
                lnΩ_tol = 1e-6, rate_rtol = 1e-6) -> NamedTuple

How each kinetic reaction of a run stands with the second law: at each instant,
its affinity `𝒜_j` ([`reaction_affinities`](@ref)), its rate `r_j` evaluated by
its own law on the certified composition, and the product `𝒜_j r_j`, the power
it dissipates, which must not be negative.

The fields are `times`, `affinity` (J/mol), `rate` (mol/s), `power` (`𝒜 r`, W),
`entropy_production` (`Σ_j 𝒜_j r_j / T`, W/K, one value per instant) and
`violations`, one named tuple `(time, reaction, affinity, rate)` per instant and
reaction where `𝒜_j r_j < 0` while `|ln Ω_j| > lnΩ_tol` and
`|r_j| > rate_rtol · maxₜ |r_j|`.

# What it tells

A law written as `k F (1 − Ω)` with `F ≥ 0`, such as [`transition_state`](@ref)
or [`sorption_rate`](@ref), cannot appear in `violations`: its rate has the sign
of the affinity by construction. A law that does not read the composition, a
Parrott–Killoh or Waller degree of reaction, can: it goes on dissolving a phase
the solution has become supersaturated with, and the run then dissipates a
negative power. So can a law whose reverse rate is not tied to the equilibrium
constant of the database: it stops where `Ω ≠ 1`, and between `Ω = 1` and its
own stopping point it runs against its affinity.

`lnΩ_tol` keeps out of `violations` the instants where the reaction sits at
equilibrium to the accuracy of the speciation, whose sign of `1 − Ω` is then
noise; `rate_rtol` keeps out those where the law has stopped, its rate a rounding
of zero, at an `Ω` that is not one. Stopping there is allowed: a rate may vanish
while the affinity does not, which is a metastable state. Such a stop, at an
`Ω` other than one, is visible in `affinity` at the last instants.
"""
function dissipation(
        sol, kp::KineticsProblem; times = sol.t, states = nothing,
        lnΩ_tol::Real = 1.0e-6, rate_rtol::Real = 1.0e-6,
    )
    res = _affinities_and_rates(sol, kp, times, states; rates = true)
    A, r, Ts = res.A, res.r, res.T
    power = A .* r
    production = [sum((power[i, j] for j in axes(power, 2) if isfinite(power[i, j])); init = 0.0) / Ts[i] for i in axes(power, 1)]
    names = [String(kr.reaction.symbol) for kr in kp.kinetic_reactions]
    violations = NamedTuple{(:time, :reaction, :affinity, :rate), Tuple{Float64, String, Float64, Float64}}[]
    r_max = [maximum((abs(x) for x in r[:, j] if isfinite(x)); init = 0.0) for j in axes(r, 2)]
    for i in axes(power, 1), j in axes(power, 2)
        isfinite(power[i, j]) || continue
        lnΩ = -A[i, j] / (R_GAS * Ts[i])
        power[i, j] < 0 && abs(lnΩ) > lnΩ_tol && abs(r[i, j]) > rate_rtol * r_max[j] &&
            push!(violations, (time = Float64(times[i]), reaction = names[j], affinity = A[i, j], rate = r[i, j]))
    end
    return (
        times = collect(Float64, times), affinity = A, rate = r, power = power,
        entropy_production = production, violations = violations,
    )
end

# The affinities and, when asked, the rates of every reaction at every instant,
# on the certified compositions, with the temperature of each.
function _affinities_and_rates(sol, kp, times, states; rates::Bool)
    sts = states === nothing ? speciated_states(sol, kp; times) : states
    length(sts) == length(times) || throw(
        DimensionMismatch("$(length(sts)) states for $(length(times)) instants.")
    )
    p = sol.prob.p
    sys = kp.system
    g_fns = [haskey(sp, :ΔₐG⁰) ? sp[:ΔₐG⁰] : nothing for sp in sys.species]
    M = length(kp.kinetic_reactions)
    A = fill(NaN, length(times), M)
    r = fill(NaN, length(times), M)
    Ts = zeros(length(times))
    n0 = StateView(Float64[_plain(x) for x in p.n_initial_full], p.species_index)
    for (i, st) in enumerate(sts)
        T = Float64(_plain(ustrip(us"K", temperature(st))))
        P = Float64(_plain(ustrip(us"Pa", pressure(st))))
        Ts[i] = T
        n = Float64[_plain(ustrip(us"mol", x)) for x in st.n]
        lna = Float64[_plain(x) for x in p.lna_fn(n, _lna_params(p, T))]
        μ = [isnothing(g) ? NaN : g(; T = T, P = P, unit = false) + R_GAS * T * lna[k] for (k, g) in enumerate(g_fns)]
        for (j, kr) in enumerate(kp.kinetic_reactions)
            acc = 0.0
            for (k, ν) in enumerate(kr.stoich)
                iszero(ν) && continue
                acc -= ν * μ[k]
            end
            A[i, j] = acc
            if rates
                r[i, j] = Float64(
                    _plain(
                        kr.rate_fn(
                            T, P, Float64(times[i]), StateView(n, p.species_index),
                            StateView(lna, p.species_index), n0,
                        )
                    )
                )
            end
        end
    end
    return (A = A, r = r, T = Ts)
end
