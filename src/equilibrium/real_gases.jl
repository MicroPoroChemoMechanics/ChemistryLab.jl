# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# ── Real gases: the equation of state of Peng and Robinson (1976) ────────────
#
# A gas phase is ideal unless its species carry critical constants, which
# `peng_robinson` attaches to them. The constants travel with the species, so
# that every system rebuilt from them (an equilibrium partition, a second
# instance, a phase list) keeps the gas model without being told.

# Ω_a and Ω_b, the values the equation takes at the critical point, where the
# cubic in Z has a triple root Z_c: matching (Z − Z_c)³ coefficient by
# coefficient gives Z_c = (1 − B)/3, A = 3Z_c² + 3B² + 2B, and for B the cubic
# 64B³ + 6B² + 12B − 1 = 0. Peng and Robinson print the two roots rounded,
# 0.45724 and 0.07780 (their Eqs. 9 and 10); rounded, they move the triple root
# by about 0.01, the cube root of the rounding error, so the critical point of
# the equation would no longer be the one it is given. The coefficients of κ(ω)
# are those of their Eq. 18.
const _PR_ΩB = let B = 0.0778
    for _ in 1:8
        B -= (64B^3 + 6B^2 + 12B - 1) / (192B^2 + 12B + 12)
    end
    B
end
const _PR_ZC = (1 - _PR_ΩB) / 3
const _PR_ΩA = 3_PR_ZC^2 + 3_PR_ΩB^2 + 2_PR_ΩB
const _PR_KAPPA = (0.37464, 1.54226, -0.26992)

"""
    peng_robinson(s::Species; T_c, P_c, ω, kij = Pair{Symbol, Float64}[]) -> Species
    peng_robinson(s::Species; kij = Pair{Symbol, Float64}[]) -> Species

The gas `s` described by the equation of state of [PengRobinson1976](@citet),
through its critical temperature `T_c`, critical pressure `P_c` and acentric
factor `ω`, with the binary interaction parameters `kij`, a list of
`Symbol(partner) => k`, zero for a partner it does not name. Without the three
constants, they are those PHREEQC ships in `phreeqc.dat` for the gas of the same
symbol ([ParkhurstAppelo2013](@citet)), which `literature_table("ParkhurstAppelo2013",
"gas_critical_constants")` lists.

A gas phase is described as a whole: either every gas of a system carries
critical constants and the phase follows the equation of state, or none does
and the phase is ideal. The activity of a gas is then
``\\ln a_i = \\ln y_i + \\ln\\varphi_i(T, P, \\mathbf{y}) + \\ln(P/P^\\circ)``, its
fugacity coefficient ``\\varphi_i`` given by the equation of state, and the
standard state stays the ideal gas at ``P^\\circ``.

`T_c` and `P_c` are quantities, or plain numbers in K and Pa. The copy shares
the thermodynamic functions of `s`; only its properties gain the constants, in
SI units, under `:T_c` (K), `:P_c` (Pa), `:ω` and `:kij`.

# Examples

```julia
co2 = peng_robinson(db["CO2"])                       # PHREEQC's constants
co2 = peng_robinson(db["CO2"]; T_c = 304.128u"K", P_c = 73.773u"bar", ω = 0.2239)
n2 = peng_robinson(db["N2"]; kij = [:CO2 => -0.02])
```
"""
function peng_robinson(
        s::Species{T}; T_c = nothing, P_c = nothing, ω = nothing, kij = Pair{Symbol, Float64}[],
    ) where {T}
    aggregate_state(s) == AS_GAS || throw(
        ArgumentError("peng_robinson: $(symbol(s)) is not a gas; the equation of state describes a gas phase.")
    )
    given = (T_c !== nothing, P_c !== nothing, ω !== nothing)
    if !any(given)
        t = literature_table("ParkhurstAppelo2013", "gas_critical_constants")
        i = findfirst(==(symbol(s)), t.gas)
        i === nothing && throw(
            ArgumentError(
                "peng_robinson: phreeqc.dat gives no critical constants for $(symbol(s)); it gives " *
                    "them for $(join(t.gas, ", ")). Pass `T_c`, `P_c` and `ω`.",
            )
        )
        T_c, P_c, ω = t.T_c[i], t.P_c[i], t.omega[i]
    elseif !all(given)
        throw(ArgumentError("peng_robinson: give `T_c`, `P_c` and `ω` together, or none of them."))
    end
    Tc = _pr_si(T_c, u"K")
    Pc = _pr_si(P_c, u"Pa")
    (Tc > 0 && Pc > 0) || throw(ArgumentError("peng_robinson: T_c and P_c must be positive; got $T_c and $P_c."))
    props = copy(s.properties)
    props[:T_c] = Tc
    props[:P_c] = Pc
    props[:ω] = Float64(ω)
    props[:kij] = Pair{Symbol, Float64}[Symbol(first(p)) => Float64(last(p)) for p in kij]
    return Species{T}(s.name, s.symbol, s.formula, s.aggregate_state, s.class, props)
end

# A value in SI units: a plain number is taken as given in them (K, Pa), a
# quantity is converted, and refused when its dimension is not the expected one.
_pr_si(x::Real, _) = Float64(x)
function _pr_si(x::DynamicQuantities.AbstractQuantity, unit)
    q = uexpand(x)
    dimension(q) == dimension(unit) || throw(DimensionError(q, unit))
    return Float64(ustrip(q))
end

"""
    _PengRobinsonMixing

The Peng-Robinson description of the gas phase of a system: per gas, in the
order of `idx_gas`, the critical temperature (K) and pressure (Pa) and the
coefficient ``\\kappa(\\omega)``, and the matrix of binary interaction
parameters.
"""
struct _PengRobinsonMixing
    T_c::Vector{Float64}
    P_c::Vector{Float64}
    κ::Vector{Float64}
    kij::Matrix{Float64}
end

"""
    _gas_mixing(cs) -> Union{Nothing, _PengRobinsonMixing}

The mixing model of the gas phase of `cs`: `nothing` for an ideal mixture, a
Peng-Robinson description when its gases carry critical constants. A phase in
which some gases carry them and others do not is refused: there is no mixing
rule between an equation of state and its absence.
"""
function _gas_mixing(cs)
    idx = cs.idx_gas
    isempty(idx) && return nothing
    has = [haskey(properties(cs.species[i]), :T_c) for i in idx]
    any(has) || return nothing
    all(has) || throw(
        ArgumentError(
            "the gas phase is ideal or follows the equation of state of Peng and Robinson as a " *
                "whole: $(join([symbol(cs.species[i]) for (i, h) in zip(idx, has) if h], ", ")) carry " *
                "critical constants and $(join([symbol(cs.species[i]) for (i, h) in zip(idx, has) if !h], ", ")) " *
                "do not. Give them with `peng_robinson`.",
        )
    )
    # `properties`, not `s[:key]`, which looks among the atoms first.
    pr = [properties(cs.species[i]) for i in idx]
    Tc = [Float64(q[:T_c]) for q in pr]
    Pc = [Float64(q[:P_c]) for q in pr]
    κ = [_PR_KAPPA[1] + _PR_KAPPA[2] * q[:ω] + _PR_KAPPA[3] * q[:ω]^2 for q in pr]
    sym = [Symbol(symbol(cs.species[i])) for i in idx]
    m = length(idx)
    kij = zeros(m, m)
    for (a, q) in enumerate(pr), (partner, k) in q[:kij]
        b = findfirst(==(partner), sym)
        b === nothing && continue
        kij[a, b] != 0 && kij[a, b] != k && throw(
            ArgumentError("the interaction parameter between $(sym[a]) and $(sym[b]) is given twice, as $(kij[a, b]) and $k.")
        )
        kij[a, b] = kij[b, a] = k
    end
    return _PengRobinsonMixing(Tc, Pc, κ, kij)
end

# The parameters a and b of each gas at T, in SI units (Pa m⁶/mol², m³/mol).
function _pr_ab(mix::_PengRobinsonMixing, T)
    R = R_GAS
    a = [_PR_ΩA * (R * Tc)^2 / Pc * (1 + κ * (1 - sqrt(T / Tc)))^2 for (Tc, Pc, κ) in zip(mix.T_c, mix.P_c, mix.κ)]
    b = [_PR_ΩB * R * Tc / Pc for (Tc, Pc) in zip(mix.T_c, mix.P_c)]
    return a, b
end

# The dimensionless residual Gibbs energy of the phase per mole, G_res/(nRT),
# at the compressibility factor Z, for the reduced parameters A and B.
_pr_g_res(Z, A, B) = Z - 1 - log(Z - B) -
    A / (2 * sqrt(2) * B) * log((Z + (1 + sqrt(2)) * B) / (Z + (1 - sqrt(2)) * B))

# The cubic in Z and its derivative.
_pr_cubic(Z, A, B) = Z^3 - (1 - B) * Z^2 + (A - 3B^2 - 2B) * Z - (A * B - B^2 - B^3)
_pr_cubic′(Z, A, B) = 3Z^2 - 2 * (1 - B) * Z + (A - 3B^2 - 2B)

"""
    _pr_roots(a, b) -> Vector{Float64}

The real roots above `b` of the cubic in `Z` at the reduced parameters `a`, `b`
(plain numbers), in closed form, each sharpened by one Newton step: one root in
a single-phase region, three where the equation has a vapor and a liquid root.
"""
function _pr_roots(a, b)
    # Z³ + c₂Z² + c₁Z + c₀, depressed by Z = t − c₂/3 and solved in closed form.
    c2, c1, c0 = -(1 - b), a - 3b^2 - 2b, -(a * b - b^2 - b^3)
    p = c1 - c2^2 / 3
    q = 2c2^3 / 27 - c2 * c1 / 3 + c0
    Δ = (q / 2)^2 + (p / 3)^3
    roots = if Δ > 0
        s = sqrt(Δ)
        [cbrt(-q / 2 + s) + cbrt(-q / 2 - s) - c2 / 3]
    elseif p == 0
        # The triple root of the critical point.
        [-c2 / 3]
    else
        r = 2 * sqrt(-p / 3)
        θ = acos(clamp(3q / (p * r), -1.0, 1.0)) / 3
        [r * cos(θ - 2π * k / 3) - c2 / 3 for k in 0:2]
    end
    return [_pr_newton(z, a, b) for z in roots if z > b]
end

# One Newton step on the cubic, kept only where it is defined: at a double or
# triple root the derivative vanishes with the cubic, and the step is 0/0.
function _pr_newton(Z, A, B)
    step = _pr_cubic(Z, A, B) / _pr_cubic′(Z, A, B)
    return isfinite(_plain(step)) ? Z - step : Z
end

"""
    _pr_z(A, B) -> Z

The compressibility factor: of the real roots of the cubic above `B`, the one
of least Gibbs energy, which is the stable phase where the equation has a vapor
and a liquid root. The roots are found on the values, then lifted into the dual
numbers of `A` and `B` by two Newton steps on the cubic: the value does not move
and the derivatives of the root are exact to second order, as the
implicit-function theorem gives them.
"""
function _pr_z(A, B)
    a, b = _plain(A), _plain(B)
    roots = _pr_roots(a, b)
    isempty(roots) && throw(DomainError((a, b), "the cubic of Peng and Robinson has no root above B."))
    z0 = roots[argmin([_pr_g_res(z, a, b) for z in roots])]
    Z = z0 + zero(A) + zero(B)
    for _ in 1:2
        Z = _pr_newton(Z, A, B)
    end
    return Z
end

"""
    _pr_ln_phi(mix, y, T, P) -> (lnφ, Z)

The fugacity coefficients of the gases at the mole fractions `y`, the
temperature `T` (K) and the pressure `P` (Pa), and the compressibility factor of
the phase, by the van der Waals mixing rules,
``a = \\sum_{ij} y_i y_j (1 - k_{ij})\\sqrt{a_i a_j}`` and ``b = \\sum_i y_i b_i``
([PengRobinson1976](@citet), their Eq. 19 and the mixture form of their
fugacity coefficient).
"""
function _pr_ln_phi(mix::_PengRobinsonMixing, y, T, P)
    R = R_GAS
    ai, bi = _pr_ab(mix, T)
    m = length(y)
    aij(i, j) = (1 - mix.kij[i, j]) * sqrt(ai[i] * ai[j])
    s = [sum(y[i] * aij(i, k) for i in 1:m) for k in 1:m]
    a = sum(y[k] * s[k] for k in 1:m)
    b = sum(y[k] * bi[k] for k in 1:m)
    A = a * P / (R * T)^2
    B = b * P / (R * T)
    Z = _pr_z(A, B)
    L = log((Z + (1 + sqrt(2)) * B) / (Z + (1 - sqrt(2)) * B))
    lnφ = [
        bi[k] / b * (Z - 1) - log(Z - B) - A / (2 * sqrt(2) * B) * (2 * s[k] / a - bi[k] / b) * L
            for k in 1:m
    ]
    return lnφ, Z
end

"""
    _pr_residual_gibbs(mix, n, T, P) -> G_res/(RT)

The residual Gibbs energy of the gas phase over ``RT``, for the amounts `n` (in
the order of the mixture), at `T` (K) and `P` (Pa). Its derivatives with
respect to the amounts are the logarithms of the fugacity coefficients, which
the tests check `_pr_ln_phi` against.
"""
function _pr_residual_gibbs(mix::_PengRobinsonMixing, n, T, P)
    R = R_GAS
    ai, bi = _pr_ab(mix, T)
    N = sum(n)
    y = n ./ N
    m = length(y)
    a = sum(y[i] * y[j] * (1 - mix.kij[i, j]) * sqrt(ai[i] * ai[j]) for i in 1:m, j in 1:m)
    b = sum(y[k] * bi[k] for k in 1:m)
    A = a * P / (R * T)^2
    B = b * P / (R * T)
    return N * _pr_g_res(_pr_z(A, B), A, B)
end

"""
    _pr_phase_volume(mix, n, T, P) -> V

The volume of the gas phase, ``V = Z N R T / P`` in m³, for the amounts `n` of
its gases (mol, in the order of the mixture), at `T` (K) and `P` (Pa): the
volume of the ideal gas scaled by the compressibility factor of the mixture.
"""
function _pr_phase_volume(mix::_PengRobinsonMixing, n, T, P)
    N = sum(n; init = zero(eltype(n)))
    iszero(N) && return zero(promote_type(typeof(N), typeof(T), typeof(P)))
    _, Z = _pr_ln_phi(mix, n ./ N, T, P)
    return Z * N * R_GAS * T / P
end

# The gases of `state`, their fugacity coefficients and the compressibility
# factor of their phase, at the composition the activity models use: a gas
# alone is the whole phase, the members of a mixture floored as there.
function _gas_phase_state(state::ChemicalState)
    cs = state.system
    idx = cs.idx_gas
    isempty(idx) && throw(ArgumentError("the system has no gas phase."))
    mix = _gas_mixing(cs)
    mix === nothing && return (idx, ones(length(idx)), 1.0)
    T = _plain(ustrip(us"K", temperature(state)))
    P = _plain(ustrip(us"Pa", pressure(state)))
    n = [_plain(ustrip(us"mol", state.n[i])) + _AMOUNT_FLOOR for i in idx]
    y = length(idx) == 1 ? [1.0] : n ./ sum(n)
    lnφ, Z = _pr_ln_phi(mix, y, T, P)
    return (idx, exp.(lnφ), Z)
end

"""
    fugacity_coefficients(state::ChemicalState) -> OrderedDict{String, Float64}

The fugacity coefficient ``\\varphi_i`` of every gas of `state`, at its
temperature, pressure and gas composition: 1 in an ideal gas phase, the value
of the equation of state of Peng and Robinson when the gases carry critical
constants ([`peng_robinson`](@ref)). The fugacity of a gas is
``\\varphi_i y_i P``.

```julia
co2 = peng_robinson(db["CO2"])
st = ChemicalState(ChemicalSystem([co2], ["CO2"]); T = 298.15u"K", P = 50u"bar", n = [1.0u"mol"])
fugacity_coefficients(st)["CO2"]             # 0.741
```

See also: [`compressibility_factor`](@ref).
"""
function fugacity_coefficients(state::ChemicalState)
    idx, φ, _ = _gas_phase_state(state)
    return OrderedDict{String, Float64}(symbol(state.system.species[i]) => φ[k] for (k, i) in enumerate(idx))
end

"""
    compressibility_factor(state::ChemicalState) -> Float64

The compressibility factor ``Z = PV/(NRT)`` of the gas phase of `state`: 1 for
an ideal gas phase, the root of the cubic of Peng and Robinson of least Gibbs
energy when the gases carry critical constants ([`peng_robinson`](@ref)). The
gas phase occupies ``Z N R T/P``.

See also: [`fugacity_coefficients`](@ref).
"""
compressibility_factor(state::ChemicalState) = _gas_phase_state(state)[3]
