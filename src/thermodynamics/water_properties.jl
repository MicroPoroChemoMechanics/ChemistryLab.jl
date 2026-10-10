# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)
# Portions of this file are Julia ports adapted from the Reaktoro C++ library
# (https://github.com/reaktoro/reaktoro), Copyright © 2014-2024 Allan Leal,
# distributed under the LGPL-2.1-or-later.

using ForwardDiff

# Water thermodynamic and electrostatic properties using:
# - the HGK equation of state of Haar et al. (1984) for density and derivatives;
# - the dielectric constant model of Johnson and Norton (1991) and Born functions;
# - the g-function of Shock et al. (1992) for electrostatic corrections.
#
# The functions accept dual numbers, nested ones included. The density is solved
# on plain numbers and lifted into the duals (`water_density_hgk`), so symbolic
# numbers do not pass through it.

# ============================================================
#  Public structs
# ============================================================

"""
    WaterThermoProps{T<:Real}

Density `D` (kg/m³) of liquid water and its partial derivatives w.r.t. T (K) and P (Pa),
computed from the equation of state of [Haar1984](@citet).
"""
struct WaterThermoProps{T <: Real}
    D::T
    DT::T
    DP::T
    DTT::T
    DTP::T
    DPP::T
end

"""
    WaterElectroProps{T<:Real}

Dielectric constant `ε` of water and Born functions `Z, Y, Q, X, U, N`,
from [JohnsonNorton1991](@citet).
"""
struct WaterElectroProps{T <: Real}
    epsilon::T
    bornZ::T
    bornY::T
    bornQ::T
    bornX::T
    bornU::T
    bornN::T
end

"""
    HKFGState{T<:Real}

The g-function of [Shock1992](@citet) and its first- and second-order partial derivatives
w.r.t. T and P.
"""
struct HKFGState{T <: Real}
    g::T
    gT::T
    gP::T
    gTT::T
    gTP::T
    gPP::T
end

"""
    SpeciesElectroPropsHKF{T<:Real}

Born coefficient `ω` for an aqueous species and its partial derivatives,
from the HKF model of [Helgeson1981](@citet) extended by [Shock1992](@citet).
"""
struct SpeciesElectroPropsHKF{T <: Real}
    w::T
    wT::T
    wP::T
    wTT::T
    wTP::T
    wPP::T
end

# ============================================================
#  Internal Helmholtz accumulator
# ============================================================

struct _HGKHelm{T <: Real}
    A::T
    AD::T
    AT::T
    ADD::T
    ATD::T
    ATT::T
    ADDD::T
    ATDD::T
    ATTD::T
    ATTT::T
end

function _HGKHelm(::Type{T}) where {T <: Real}
    return _HGKHelm{T}(ntuple(_ -> zero(T), 10)...)
end

function Base.:+(a::_HGKHelm, b::_HGKHelm)
    return _HGKHelm(
        a.A + b.A,
        a.AD + b.AD,
        a.AT + b.AT,
        a.ADD + b.ADD,
        a.ATD + b.ATD,
        a.ATT + b.ATT,
        a.ADDD + b.ADDD,
        a.ATDD + b.ATDD,
        a.ATTD + b.ATTD,
        a.ATTT + b.ATTT,
    )
end

# ============================================================
#  HGK reference constants
# ============================================================

const _HGK_Tk = 647.27   # reference temperature  (K)
const _HGK_Dk = 317.763  # reference density      (kg/m³)
const _HGK_Hk = 69595.89 # reference Helmholtz    (J/kg)

# ============================================================
#  HGK auxiliary function 0 — ideal-gas contribution (no density dependence)
# ============================================================

function _hgk_helm0(t::T, ::T) where {T <: Real}
    A0 = (
        -0.130840393653e+2,
        -0.85702042094e+2,
        0.765192919131e-2,
        -0.620600116069e+0,
        -0.106924329402e+2,
        -0.280671377296e+1,
        0.119843634845e+3,
        -0.823907389256e+2,
        0.555864146443e+2,
        -0.31069812298e+2,
        0.136200239305e+2,
        -0.457116129409e+1,
        0.115382128188e+1,
        -0.214242224683e+0,
        0.282800597384e-1,
        -0.250384152737e-2,
        0.132952679669e-3,
        -0.319277411208e-5,
    )

    ln_t = log(t)

    # First two terms (A0[1]+A0[2]*t)*log(t)   (C++ indices 0,1)
    A = (A0[1] + A0[2] * t) * ln_t
    AT = A0[1] / t + A0[2] * (ln_t + 1)
    ATT = -A0[1] / (t * t) + A0[2] / t
    ATTT = 2 * A0[1] / (t^3) - A0[2] / (t * t)

    # Remaining terms A0[j]*t^(j-5), C++ i=2..17 → Julia j=3..18, exponent = i-4 = j-5
    for j in 3:18
        p = j - 5
        aux = A0[j] * t^p
        A += aux
        AT += aux * p / t
        ATT += aux * p * (p - 1) / (t * t)
        ATTT += aux * p * (p - 1) * (p - 2) / (t^3)
    end

    z = zero(T)
    return _HGKHelm(A, z, AT, z, z, ATT, z, z, z, ATTT)
end

# ============================================================
#  HGK auxiliary function 1 — first-order density correction
# ============================================================

function _hgk_helm1(t::T, d::T) where {T <: Real}
    A1 = (0.15383053e+1, -0.81048367e+0, -0.68305748e+1, 0.0e+0, 0.86756271e+0)

    A = AT = ATT = ATTT = zero(T)

    # C++ loop i=0..4: aux = d*A1[i]*t^(1-i)
    # Julia j=1..5: i = j-1, exponent = 1-i = 2-j
    for j in 1:5
        i = j - 1
        p = 1 - i
        aux = d * A1[j] * t^p
        A += aux
        AT -= aux * (i - 1) / t
        ATT += aux * (i - 1) * i / (t * t)
        ATTT -= aux * (i - 1) * i * (i + 1) / (t^3)
    end

    AD = A / d
    ATD = AT / d
    ATTD = ATT / d

    z = zero(T)
    return _HGKHelm(A, AD, AT, z, ATD, ATT, z, z, ATTD, ATTT)
end

# ============================================================
#  HGK auxiliary function 2 — second-order density correction
# ============================================================

function _hgk_helm2(t::T, d::T) where {T <: Real}
    A20 = 0.42923415e+1
    yc = (0.59402227e-1, -0.28128238e-1, 0.56826674e-3, -0.27987451e-3)

    t3 = t^(-3)
    t5 = t^(-5)
    ln_t = log(t)

    # y and derivatives (∂y/∂d = y_r, ∂y/∂t = y_t, etc.)
    y = d * (yc[1] + yc[2] * ln_t + yc[3] * t3 + yc[4] * t5)
    y_r = y / d
    y_t = d * (yc[2] - 3 * yc[3] * t3 - 5 * yc[4] * t5) / t
    y_tt = d * (-yc[2] + 12 * yc[3] * t3 + 30 * yc[4] * t5) / (t * t)
    y_rt = y_t / d
    # y_rrr = 0, y_rrt = 0
    y_rtt = y_tt / d
    y_ttt = d * (2 * yc[2] - 60 * yc[3] * t3 - 210 * yc[4] * t5) / (t^3)

    x = 1 / (1 - y)
    x2 = x * x
    x_r = y_r * x2
    x_t = y_t * x2
    x_rr = 2 * y_r * x_r * x          # y_rr=0
    x_tt = y_tt * x2 + 2 * y_t * x_t * x
    x_rt = y_rt * x2 + 2 * y_r * x_t * x
    x_rrr = 2 * y_r * (x_rr * x + x_r * x_r)  # y_rrr=0
    x_rrt = y_rt * x_r + x * (y_rt * x_r + y_r * x_rt)   # simplified from C++
    # direct translation: y_rrt*x2 + 2*(y_rt*x_r + y_rr*x_t) + 2*y_r*(x_rt*x + x_r*x_t)
    # with y_rr=0, y_rrt=0:
    x_rrt = 2 * y_rt * x_r + 2 * y_r * (x_rt * x + x_r * x_t)
    x_rtt = y_rtt * x2 + 4 * y_rt * x_t * x + 2 * y_r * (x_tt * x + x_t * x_t)
    x_ttt = y_ttt * x2 + 4 * y_tt * x_t * x + 2 * y_t * (x_tt * x + x_t * x_t)

    u = log(d * x)
    u_r = x_r / x + 1 / d
    u_t = x_t / x
    u_rr = x_rr / x - (x_r / x)^2 - 1 / (d * d)
    u_rt = x_rt / x - x_r * x_t / (x * x)
    u_tt = x_tt / x - (x_t / x)^2
    u_rrr = x_rrr / x - 3 * x_rr * x_r / (x * x) + 2 * (x_r / x)^3 + 2 / (d^3)
    u_rrt = x_rrt / x - (2 * x_rt * x_r + x_rr * x_t) / (x * x) + 2 * x_r * x_r * x_t / (x^3)
    u_rtt = x_rtt / x - (2 * x_rt * x_t + x_tt * x_r) / (x * x) + 2 * x_t * x_t * x_r / (x^3)
    u_ttt = x_ttt / x - 3 * x_tt * x_t / (x * x) + 2 * (x_t / x)^3

    c1 = -130.0 / 3
    c2 = 169.0 / 6
    c3 = -14.0

    A = A20 * t * (u + c1 * x + c2 * x * x + c3 * y)
    AD = A20 * t * (u_r + c1 * x_r + 2 * c2 * x * x_r + c3 * y_r)
    AT = A20 * t * (u_t + c1 * x_t + 2 * c2 * x * x_t + c3 * y_t) + A / t
    ADD = A20 * t * (u_rr + c1 * x_rr + 2 * c2 * (x * x_rr + x_r * x_r))  # y_rr=0
    ATD = A20 * t * (u_rt + c1 * x_rt + 2 * c2 * (x * x_rt + x_r * x_t) + c3 * y_rt) + AD / t
    ATT = A20 * t * (u_tt + c1 * x_tt + 2 * c2 * (x * x_tt + x_t * x_t) + c3 * y_tt) +
        2 * (AT / t - A / (t * t))
    ADDD = A20 * t * (u_rrr + c1 * x_rrr + 2 * c2 * (3 * x_r * x_rr + x * x_rrr))  # y_rrr=0
    ATDD = A20 * t * (u_rrt + c1 * x_rrt + 2 * c2 * (x_t * x_rr + 2 * x_r * x_rt + x * x_rrt)) +
        ADD / t  # y_rrt=0
    ATTD = A20 * t * (
        u_rtt + c1 * x_rtt + 2 * c2 * (x_r * x_tt + 2 * x_t * x_rt + x * x_rtt) +
            c3 * y_rtt
    ) +
        2 * (ATD - AD / t) / t
    ATTT = A20 * t * (u_ttt + c1 * x_ttt + 2 * c2 * (3 * x_t * x_tt + x * x_ttt) + c3 * y_ttt) +
        3 * (ATT - 2 * AT / t + A / (t * t)) / t

    return _HGKHelm(A, AD, AT, ADD, ATD, ATT, ADDD, ATDD, ATTD, ATTT)
end

# ============================================================
#  HGK auxiliary function 3 — main polynomial series
# ============================================================

function _hgk_helm3(t::T, d::T) where {T <: Real}
    ki = (
        1, 1, 1, 1, 2, 2, 2, 2,
        3, 3, 3, 3, 4, 4, 4, 4,
        5, 5, 5, 5, 6, 6, 6, 6,
        7, 7, 7, 7, 9, 9, 9, 9,
        3, 3, 1, 5,
    )
    li = (
        1, 2, 4, 6, 1, 2, 4, 6,
        1, 2, 4, 6, 1, 2, 4, 6,
        1, 2, 4, 6, 1, 2, 4, 6,
        1, 2, 4, 6, 1, 2, 4, 6,
        0, 3, 3, 3,
    )
    A3 = (
        -0.76221190138079e+1,
        0.32661493707555e+2,
        0.11305763156821e+2,
        -0.10015404767712e+1,
        0.12830064355028e+3,
        -0.28371416789846e+3,
        0.24256279839182e+3,
        -0.99357645626725e+2,
        -0.12275453013171e+4,
        0.23077622506234e+4,
        -0.16352219929859e+4,
        0.58436648297764e+3,
        0.42365441415641e+4,
        -0.78027526961828e+4,
        0.38855645739589e+4,
        -0.91225112529381e+3,
        -0.90143895703666e+4,
        0.15196214817734e+5,
        -0.39616651358508e+4,
        -0.72027511617558e+3,
        0.1114712670599e+5,
        -0.1741206525221e+5,
        0.99918281207782e+3,
        0.33504807153854e+4,
        -0.64752644922631e+4,
        0.98323730907847e+4,
        0.83877854108422e+3,
        -0.27919349903103e+4,
        0.11112410081192e+4,
        -0.17287587261807e+4,
        -0.36233262795423e+3,
        0.61139429010144e+3,
        0.32968064728562e+2,
        0.10411239605066e+3,
        -0.3822587471259e+2,
        -0.20307478607599e+3,
    )

    z0_hgk = 0.317763

    z = 1 - exp(-z0_hgk * d)
    z_r = z0_hgk * (1 - z)
    z_rr = -z0_hgk * z_r
    z_rrr = -z0_hgk * z_rr

    A = AD = AT = ADD = ATD = ATT = ADDD = ATDD = ATTD = ATTT = zero(T)

    for i in 1:36
        k = ki[i]
        l = li[i]

        λ = A3[i] * t^(-l) * z^k
        λ_r = k * z_r * λ / z
        λ_t = -l * λ / t
        λ_rr = λ_r * (z_rr / z_r + λ_r / λ - z_r / z)
        λ_rt = λ_r * λ_t / λ
        λ_tt = λ_t * (λ_t / λ - 1 / t)
        λ_rrr = λ_rr * (z_rr / z_r + λ_r / λ - z_r / z) +
            λ_r * (
            z_rrr / z_r - (z_rr / z_r)^2 + λ_rr / λ - (λ_r / λ)^2 -
                z_rr / z + (z_r / z)^2
        )
        λ_rrt = -(λ_r / λ)^2 * λ_t + (λ_rr * λ_t + λ_rt * λ_r) / λ
        λ_rtt = -(λ_t / λ)^2 * λ_r + (λ_tt * λ_r + λ_rt * λ_t) / λ
        λ_ttt = λ_tt * (λ_t / λ - 1 / t) +
            λ_t * (λ_tt / λ - (λ_t / λ)^2 + 1 / (t * t))

        A += λ
        AD += λ_r
        AT += λ_t
        ADD += λ_rr
        ATD += λ_rt
        ATT += λ_tt
        ADDD += λ_rrr
        ATDD += λ_rrt
        ATTD += λ_rtt
        ATTT += λ_ttt
    end

    return _HGKHelm(A, AD, AT, ADD, ATD, ATT, ADDD, ATDD, ATTD, ATTT)
end

# ============================================================
#  HGK auxiliary function 4 — critical region correction
# ============================================================

function _hgk_helm4(t::T, d::T) where {T <: Real}
    mi_arr = (2, 2, 2, 4)
    ni_arr = (0, 2, 0, 0)
    alpha_arr = (34.0, 40.0, 30.0, 1050.0)
    beta_arr = (20000.0, 20000.0, 40000.0, 25.0)
    ri_arr = (0.10038928e+1, 0.10038928e+1, 0.10038928e+1, 0.48778492e+1)
    ti_arr = (0.98876821e+0, 0.98876821e+0, 0.99124013e+0, 0.41713659e+0)
    A4_arr = (-0.32329494e-2, -0.24139355e-1, 0.79027651e-3, -0.13362857e+1)

    A = AD = AT = ADD = ATD = ATT = ADDD = ATDD = ATTD = ATTT = zero(T)

    for i in 1:4
        ri = ri_arr[i]
        ti_val = ti_arr[i]
        m = mi_arr[i]
        n = ni_arr[i]
        α = alpha_arr[i]
        β = beta_arr[i]
        A4i = A4_arr[i]

        delta = (d - ri) / ri
        tau = (t - ti_val) / ti_val
        delta_r = 1 / ri
        tau_t = 1 / ti_val

        delta_m = delta^m
        delta_n = (n == 0) ? one(T) : delta^n

        ψ = (n - α * m * delta_m) * delta_r / delta
        ψ_r = -(n + α * m * (m - 1) * delta_m) * (delta_r / delta)^2
        ψ_rr = (2 * n - α * m * (m - 1) * (m - 2) * delta_m) * (delta_r / delta)^3

        θ = A4i * delta_n * exp(-α * delta_m - β * tau * tau)
        θ_r = ψ * θ
        θ_t = -2 * β * tau * tau_t * θ
        θ_rr = ψ_r * θ + ψ * θ_r
        θ_tt = 2 * β * (2 * β * tau * tau - 1) * tau_t * tau_t * θ
        θ_rt = -2 * β * tau * tau_t * θ_r
        θ_rrr = ψ_rr * θ + 2 * ψ_r * θ_r + ψ * θ_rr
        θ_rrt = ψ_r * θ_r + ψ * θ_rt
        θ_rtt = ψ * θ_tt
        θ_ttt = -2 * β * (2 * tau_t * tau_t * θ_t + tau * tau_t * θ_tt)

        A += θ
        AD += θ_r
        AT += θ_t
        ADD += θ_rr
        ATD += θ_rt
        ATT += θ_tt
        ADDD += θ_rrr
        ATDD += θ_rrt
        ATTD += θ_rtt
        ATTT += θ_ttt
    end

    return _HGKHelm(A, AD, AT, ADD, ATD, ATT, ADDD, ATDD, ATTD, ATTT)
end

# ============================================================
#  Full HGK Helmholtz free energy (dimensionless → dimensional)
# ============================================================

"""
    water_helmholtz_hgk(T_K, D_kgm3) -> _HGKHelm

Compute the specific Helmholtz free energy of water and its derivatives at temperature
`T_K` (K) and density `D_kgm3` (kg/m³) using the equation of state of
[Haar1984](@citet).

All outputs are in SI units (J/kg per derivative w.r.t. appropriate variables).
AD-compatible (ForwardDiff-safe).
"""
function water_helmholtz_hgk(T_K::T, D_kgm3::T) where {T <: Real}
    t = T_K / _HGK_Tk
    d = D_kgm3 / _HGK_Dk

    res = _hgk_helm0(t, d) + _hgk_helm1(t, d) + _hgk_helm2(t, d) +
        _hgk_helm3(t, d) + _hgk_helm4(t, d)

    # Convert dimensionless derivatives to dimensional (SI) units
    Hk = _HGK_Hk
    Tk = _HGK_Tk
    Dk = _HGK_Dk

    return _HGKHelm(
        res.A * Hk,
        res.AD * Hk / Dk,
        res.AT * Hk / Tk,
        res.ADD * Hk / (Dk * Dk),
        res.ATD * Hk / (Dk * Tk),
        res.ATT * Hk / (Tk * Tk),
        res.ADDD * Hk / (Dk * Dk * Dk),
        res.ATDD * Hk / (Dk * Dk * Tk),
        res.ATTD * Hk / (Dk * Tk * Tk),
        res.ATTT * Hk / (Tk * Tk * Tk),
    )
end

water_helmholtz_hgk(T_K::Real, D_kgm3::Real) =
    water_helmholtz_hgk(promote(T_K, D_kgm3)...)

# ============================================================
#  Newton solver for water density given T and P
# ============================================================

"""
    water_density_hgk(T_K, P_Pa; D0 = 1000.0) -> density (kg/m³)

The density of water at the temperature `T_K` (K) and the pressure `P_Pa` (Pa)
by the HGK equation of state: the liquid where a liquid root exists, otherwise
the only root, that of the vapor or of the supercritical fluid.

The root of ``P_{HGK}(T, \\rho) = P`` is found on plain numbers, to ``10^{-12}`` in
the density: by Newton's method from `D0`, and, when that does not land on the
liquid, from a density no liquid reaches, where the pressure of the equation is
increasing and convex in the density so that the iteration descends onto the
largest root without overshooting it; failing a liquid root, from the side of the
dilute gas. It is then lifted into the dual numbers of `T_K` and `P_Pa` by two
Newton steps: the value does not move, and the derivatives are exact to second
order, as the implicit-function theorem gives them.
"""
function water_density_hgk(T_K::Real, P_Pa::Real; D0::Real = 1000.0)
    R = promote_type(typeof(T_K), typeof(P_Pa))
    R <: ForwardDiff.Dual || return _water_density_hgk(T_K, P_Pa, Float64(D0))
    return _WATER_DENSITY_MEMO(_water_density_hgk, (T_K, P_Pa, Float64(D0)), R)
end
water_density_hgk(T_K::Float64, P_Pa::Float64; D0::Real = 1000.0) =
    _WATER_DENSITY_MEMO(_water_density_hgk, (T_K, P_Pa, Float64(D0)))

const _WATER_DENSITY_MEMO = _TPMemo{NTuple{3, Float64}, Float64}()

function _water_density_hgk(T_K, P_Pa, D0::Float64)
    ρ = _hgk_density_root(_hgk_value(T_K), _hgk_value(P_Pa), D0) + zero(T_K) + zero(P_Pa)
    for _ in 1:2
        h = water_helmholtz_hgk(promote(T_K, ρ)...)
        ρ -= (ρ^2 * h.AD - P_Pa) / (2ρ * h.AD + ρ^2 * h.ADD)
    end
    return ρ
end

# The value of a number, through every level of dual numbers.
_hgk_value(x::ForwardDiff.Dual) = _hgk_value(ForwardDiff.value(x))
_hgk_value(x::Real) = Float64(x)

# The pressure of the equation of state at (T, ρ) and its derivative in ρ, or
# NaN beyond the density at which its repulsive term diverges (a logarithm of a
# negative number there).
function _hgk_pressure(T::Float64, ρ::Float64)
    h = try
        water_helmholtz_hgk(T, ρ)
    catch e
        e isa DomainError || rethrow()
        return NaN, NaN
    end
    return ρ^2 * h.AD, 2ρ * h.AD + ρ^2 * h.ADD
end

# Newton's method on P_HGK(T, ρ) = P from ρ₀, while the pressure increases with
# the density: the root and the slope there, or NaN where the iteration leaves
# the branch it started on or does not converge in `maxiter` steps. A step beyond
# the domain of the equation is halved until it falls inside.
function _hgk_newton(T::Float64, P::Float64, ρ₀::Float64; maxiter::Int = 200)
    ρ = ρ₀
    p, dp = _hgk_pressure(T, ρ)
    isfinite(p) || return NaN, 0.0
    for _ in 1:maxiter
        # Converged when the density is good to 1e-12: the residual against the
        # bulk modulus ρ ∂P/∂ρ, which is what the rounding of a liquid's pressure,
        # a small difference of large terms, scales with.
        abs(p - P) <= 1.0e-12 * (abs(P) + ρ * abs(dp)) && return ρ, dp
        dp > 0 || return NaN, dp
        ρ_new = ρ - (p - P) / dp
        ρ_new > 0 || (ρ_new = ρ / 2)
        p_new, dp_new = _hgk_pressure(T, ρ_new)
        for _ in 1:60
            isfinite(p_new) && break
            ρ_new = (ρ + ρ_new) / 2
            p_new, dp_new = _hgk_pressure(T, ρ_new)
        end
        isfinite(p_new) || return NaN, 0.0
        ρ, p, dp = ρ_new, p_new, dp_new
    end
    return NaN, 0.0
end

# The density no liquid water reaches on the range of the equation, from which
# the descent onto the liquid root starts, below the divergence of its repulsive
# term from 0 to 1000 °C; and the density above which a root is that of the
# liquid or of the supercritical fluid (twice the critical density, a vapor root
# never exceeding the critical density).
const _HGK_DENSITY_CEILING = 1500.0
const _HGK_LIQUID_DENSITY = 600.0

# The root is a Newton solve of the equation of state, asked for at the same
# values of (T, P) by every species and every state an evaluation builds, on dual
# numbers as well (their root is solved on the values): see `_TPMemo`.
_hgk_density_root(T::Float64, P::Float64, D0::Float64) = _HGK_ROOT_MEMO(_hgk_density_root_uncached, (T, P, D0))
const _HGK_ROOT_MEMO = _TPMemo{NTuple{3, Float64}, Float64}()

function _hgk_density_root_uncached(T::Float64, P::Float64, D0::Float64)
    ρ, dp = _hgk_newton(T, P, D0)
    isfinite(ρ) && dp > 0 && ρ >= _HGK_LIQUID_DENSITY && return ρ
    # From above, where the pressure grows without bound towards the divergence.
    ρ_liquid, dp_liquid = _hgk_newton(T, P, _HGK_DENSITY_CEILING)
    isfinite(ρ_liquid) && dp_liquid > 0 && return ρ_liquid
    isfinite(ρ) && dp > 0 && return ρ
    # The dilute side: the gas constant of water read off the equation itself,
    # P → ρ R T as ρ → 0, gives the ideal-gas density to start from.
    ρ_dilute = 1.0e-8
    R = water_helmholtz_hgk(T, ρ_dilute).AD * ρ_dilute / T
    ρ_gas, dp_gas = _hgk_newton(T, P, P / (R * T))
    isfinite(ρ_gas) && dp_gas > 0 && return ρ_gas
    throw(DomainError((T, P), "water_density_hgk: no density of the equation of state of water at $T K and $P Pa."))
end

# ============================================================
#  WaterThermoProps from T and P
# ============================================================

"""
    water_thermo_props(T_K, P_Pa) -> WaterThermoProps

Compute density and its partial derivatives for liquid water at `T_K` (K), `P_Pa` (Pa),
using the equation of state of [Haar1984](@citet).

AD-compatible (ForwardDiff-safe).
"""
function water_thermo_props(T_K::T, P_Pa::T) where {T <: Real}
    D = water_density_hgk(T_K, P_Pa)
    h = water_helmholtz_hgk(T_K, D)

    AD = h.AD
    ADD = h.ADD
    ADDD = h.ADDD
    ATD = h.ATD
    ATDD = h.ATDD  # ∂³f/∂T∂D²  (helmholtzTDD)
    ATTD = h.ATTD  # ∂³f/∂T²∂D  (helmholtzTTD)

    PD = 2 * D * AD + D * D * ADD
    PT = D * D * ATD
    PDD = 2 * AD + 4 * D * ADD + D * D * ADDD
    PTD = 2 * D * ATD + D * D * ATDD  # uses helmholtzTDD, not helmholtzTTD
    PTT = D * D * ATTD                 # uses helmholtzTTD

    DT = -PT / PD
    DP = 1 / PD
    DTT = -DT * DP * (DT * PDD + 2 * PTD + PTT / DT)
    DTP = -DP * DP * (DT * PDD + PTD)
    DPP = -DP * DP * DP * PDD

    return WaterThermoProps(D, DT, DP, DTT, DTP, DPP)
end

water_thermo_props(T_K::Real, P_Pa::Real) =
    water_thermo_props(promote(T_K, P_Pa)...)

# ============================================================
#  Surface tension, IAPWS R1-76(2014)
# ============================================================

"""
    water_surface_tension(T_K) -> surface tension (N/m)

Surface tension of liquid water against its vapor at the absolute temperature
`T_K` (K), from the interpolating equation of the International Association for
the Properties of Water and Steam [IAPWS2014](@cite):

```math
\\sigma = B \\, \\tau^{\\mu} (1 + b \\, \\tau), \\qquad \\tau = 1 - T / T_c,
```

with ``T_c = 647.096`` K, ``B = 235.8`` mN/m, ``b = -0.625`` and ``\\mu = 1.256``.
The equation holds from the triple point to ``T_c``, where the surface tension
vanishes, and the release extrapolates it into the supercooled liquid down to
−25 °C. [IAPWS2014; Table 1](@citet) lists the recommended values beside the
equation's, and from the triple point to 200 °C the two differ by 0.01 mN/m at
most. Above ``T_c`` there is no liquid, and asking is an error.

AD-compatible in `T_K`.

# Examples

```jldoctest
julia> round(1000 * water_surface_tension(298.15); digits = 2)   # mN/m at 25 °C
71.97
```
"""
function water_surface_tension(T_K::Real)
    # IAPWS R1-76(2014), p. 3: these constants define the equation.
    Tc, B, b, μ = 647.096, 235.8e-3, -0.625, 1.256
    T_K <= Tc || throw(DomainError(T_K, "water has no liquid surface above T_c = $Tc K"))
    τ = 1 - T_K / Tc
    return B * τ^μ * (1 + b * τ)
end

# ============================================================
#  Johnson-Norton (1991) dielectric constant and Born functions
# ============================================================

"""
    water_electro_props_jn(T_K, P_Pa, wtp::WaterThermoProps) -> WaterElectroProps

Compute dielectric constant and Born functions Z, Y, Q, X, U, N for water using the
model of [JohnsonNorton1991](@citet).

AD-compatible (ForwardDiff-safe).
"""
function water_electro_props_jn(T_K::T, ::T, wtp::WaterThermoProps) where {T <: Real}
    # Johnson-Norton (1991) coefficients
    a = (
        0.0,
        0.1470333593e+2,
        0.2128462733e+3,
        -0.1154445173e+3,
        0.1955210915e+2,
        -0.833034798e+2,
        0.3213240048e+2,
        -0.6694098645e+1,
        -0.3786202045e+2,
        0.6887359646e+2,
        -0.2729401652e+2,
    )

    kReferenceTemperature = T_STANDARD
    kReferenceDensity = 1000.0

    alpha = -wtp.DT / wtp.D
    beta = wtp.DP / wtp.D
    alphaT = -wtp.DTT / wtp.D + alpha * alpha
    betaT = wtp.DTP / wtp.D + alpha * beta
    betaP = wtp.DPP / wtp.D - beta * beta

    t = T_K / kReferenceTemperature
    r = wtp.D / kReferenceDensity

    # k0..k4 and their t-derivatives (scaled by 1/Tr, 1/Tr²)
    Tr = kReferenceTemperature
    k_vals = (one(T), a[2] / t, a[3] / t + a[4] + a[5] * t, a[6] / t + a[7] * t + a[8] * t * t, a[9] / t / t + a[10] / t + a[11])
    k_t = (zero(T), -a[2] / (t * t) / Tr, (-a[3] / (t * t) + a[5]) / Tr, (-a[6] / (t * t) + a[7] + 2 * a[8] * t) / Tr, (-2 * a[9] / (t * t * t) - a[10] / (t * t)) / Tr)
    k_tt = (zero(T), 2 * a[2] / (t * t * t) / Tr^2, 2 * a[3] / (t * t * t) / Tr^2, (2 * a[6] / (t * t * t) + 2 * a[8]) / Tr^2, (6 * a[9] / (t * t * t * t) + 2 * a[10] / (t * t * t)) / Tr^2)

    epsilon = zero(T)
    epsilonT = zero(T)
    epsilonP = zero(T)
    epsilonTT = zero(T)
    epsilonTP = zero(T)
    epsilonPP = zero(T)

    for i in 1:5
        ri = r^(i - 1)
        ki = k_vals[i]
        ki_t = k_t[i]
        ki_tt = k_tt[i]

        i0 = i - 1  # power of r, same as C++ i (0-based)

        epsilon += ki * ri
        epsilonT += ri * (ki_t - i0 * alpha * ki)
        epsilonP += ri * ki * i0 * beta
        epsilonTT += ri * (
            ki_tt - i0 * (alpha * ki_t + ki * alphaT) -
                i0 * alpha * (ki_t - i0 * alpha * ki)
        )
        epsilonTP += ri * ki * i0 * beta * (ki_t / ki - i0 * alpha + betaT / beta)
        epsilonPP += ri * ki * i0 * beta * (i0 * beta + betaP / beta)
    end

    e2 = epsilon * epsilon

    bornZ = -1 / epsilon
    bornY = epsilonT / e2
    bornQ = epsilonP / e2
    bornU = epsilonTP / e2 - 2 * bornY * bornQ * epsilon
    bornN = epsilonPP / e2 - 2 * bornQ * bornQ * epsilon
    bornX = epsilonTT / e2 - 2 * bornY * bornY * epsilon

    return WaterElectroProps(epsilon, bornZ, bornY, bornQ, bornX, bornU, bornN)
end

water_electro_props_jn(T_K::Real, P_Pa::Real, wtp::WaterThermoProps) =
    water_electro_props_jn(promote(T_K, P_Pa)..., wtp)

# ============================================================
#  Shock et al. (1992) g-function
# ============================================================

"""
    hkf_g_function(T_K, P_Pa, wtp::WaterThermoProps) -> HKFGState

Compute the g-function of [Shock1992](@citet) and its derivatives for the HKF model.
Returns zero state if water density is outside [350, 1000] kg/m³.

AD-compatible (ForwardDiff-safe).
"""
function hkf_g_function(T_K::T, P_Pa::T, wtp::WaterThermoProps) where {T <: Real}
    z = zero(T)
    zero_state = HKFGState(z, z, z, z, z, z)

    # Outside region of validity (Fig. 6 of Shock et al. 1992).
    # Use _primal so that the guard compares physical values only; AD derivatives
    # flow correctly through the returned zero_state when outside the valid range.
    if _primal(wtp.D) > 1000 || _primal(wtp.D) < 350
        return zero_state
    end

    TdegC = _celsius(T_K)
    Pbar = P_Pa * 1.0e-5

    # Region I coefficients (Shock et al. 1992, eqs. 24-31)
    ag1 = -2.037662
    ag2 = 5.747e-3
    ag3 = -6.557892e-6

    bg1 = 6.107361
    bg2 = -1.074377e-2
    bg3 = 1.268348e-5

    ag = ag1 + ag2 * TdegC + ag3 * TdegC * TdegC
    bg = bg1 + bg2 * TdegC + bg3 * TdegC * TdegC

    agT = ag2 + 2 * ag3 * TdegC
    bgT = bg2 + 2 * bg3 * TdegC
    agTT = 2 * ag3
    bgTT = 2 * bg3

    r = wtp.D / 1000

    alpha = -wtp.DT / wtp.D
    beta = wtp.DP / wtp.D
    alphaT = -wtp.DTT / wtp.D + alpha * alpha
    alphaP = -wtp.DTP / wtp.D - alpha * beta
    betaP = wtp.DPP / wtp.D - beta * beta

    g = ag * (1 - r)^bg
    gT = g * (agT / ag + bgT * log(1 - r) + r * alpha * bg / (1 - r))
    gP = -g * r * beta * bg / (1 - r)
    gTT = g * (
        agTT / ag - (agT / ag)^2 + bgTT * log(1 - r) +
            r * alpha * bg / (1 - r) * (2 * bgT / bg + alphaT / alpha - alpha - r * alpha / (1 - r))
    ) +
        gT * gT / g
    gTP = gP * (bgT / bg - alpha - alphaP / beta - r * alpha / (1 - r)) + gP * gT / g
    gPP = gP * (gP / g + beta + betaP / beta + r * beta / (1 - r))

    # Region II correction (Shock et al. 1992, eqs. 32-44).
    # Use _primal for the boundary conditions so AD derivatives flow through
    # the region-II terms when the primal T/P is inside the correction zone.
    if _primal(TdegC) > 155 && _primal(TdegC) < 355 && _primal(Pbar) < 1000
        af1 = 3.66666e+1   # K
        af2 = -1.504956e-10  # Å/bar³
        af3 = 5.01799e-14  # Å/bar⁴

        auxT = (TdegC - 155) / 300
        auxT1 = auxT^4.8
        auxT2 = auxT^16

        auxP = 1000 - Pbar
        auxP1 = auxP^3
        auxP2 = auxP^4

        ft = auxT1 + af1 * auxT2
        ftT = (4.8 * auxT1 + 16.0 * af1 * auxT2) / (300 * auxT)
        ftTT = (18.24 * auxT1 + 240.0 * af1 * auxT2) / (300 * auxT)^2

        fp = af2 * auxP1 + af3 * auxP2
        fpP = -(3 * af2 * auxP1 + 4 * af3 * auxP2) / (auxP * 1.0e+5)
        fpPP = (6 * af2 * auxP1 + 12 * af3 * auxP2) / (auxP * 1.0e+5)^2

        g -= ft * fp
        gT -= fp * ftT
        gP -= ft * fpP
        gTT -= fp * ftTT
        gTP -= ftT * fpP
        gPP -= ft * fpPP
    end

    return HKFGState(g, gT, gP, gTT, gTP, gPP)
end

hkf_g_function(T_K::Real, P_Pa::Real, wtp::WaterThermoProps) =
    hkf_g_function(promote(T_K, P_Pa)..., wtp)

# ============================================================
#  Species electrostatic properties (Born coefficient ω)
# ============================================================

"""
    species_electro_props_hkf(gstate::HKFGState, z, wref) -> SpeciesElectroPropsHKF

Compute the Born coefficient `ω` and its derivatives for an aqueous species with
charge `z` and reference Born coefficient `wref` (J/mol).

Uses the model of [Helgeson1981](@citet), extended by [Shock1992](@citet).
AD-compatible (ForwardDiff-safe).
"""
function species_electro_props_hkf(gstate::HKFGState{T}, z::Real, wref::Real) where {T <: Real}
    # η = 6.94656968×10⁵ J·Å/mol  (= 1.66027×10⁵ cal·Å/mol)
    η = 6.94656968e+5

    if iszero(z)
        # Neutral species: w = wref, all derivatives zero
        return SpeciesElectroPropsHKF(T(wref), zero(T), zero(T), zero(T), zero(T), zero(T))
    end

    g = gstate.g
    gT = gstate.gT
    gP = gstate.gP
    gTT = gstate.gTT
    gTP = gstate.gTP
    gPP = gstate.gPP

    reref = z * z / (wref / η + z / 3.082)
    re = reref + abs(z) * g

    X1 = -η * (abs(z * z * z) / (re * re) - z / (3.082 + g)^2)
    X2 = 2 * η * (z^4 / (re^3) - z / (3.082 + g)^3)

    w = η * (z * z / re - z / (3.082 + g))
    wT = X1 * gT
    wP = X1 * gP
    wTT = X1 * gTT + X2 * gT * gT
    wTP = X1 * gTP + X2 * gT * gP
    wPP = X1 * gPP + X2 * gP * gP

    return SpeciesElectroPropsHKF(w, wT, wP, wTT, wTP, wPP)
end
