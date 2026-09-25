using ForwardDiff
using JSON

# The HKF parameters of Al(OH)2+, read from the shipped aq17 database and brought
# to SI with the reader's own factors: this file tests the thermo functions, and
# a copy of the database typed here would test a transcription instead.
_hkf_aloh2_entry() = only(
    s for s in JSON.parsefile(datapath("aq17-thermofun.json"); dicttype = Dict{String, Any})["substances"]
        if s["symbol"] == "Al(OH)2+"
)
_hkf_aloh2_value(key) = float(only(_hkf_aloh2_entry()[key]["values"]))
function _hkf_aloh2_params()
    s = _hkf_aloh2_entry()
    hkf = only(m for m in s["TPMethods"] if haskey(m, "eos_hkf_coeffs"))["eos_hkf_coeffs"]["values"]
    names = [:a1, :a2, :a3, :a4, :c1, :c2, :wref]
    return vcat(
        [n => float(hkf[i]) * ChemistryLab.HKF_SI_CONVERSIONS[n] for (i, n) in enumerate(names)],
        [
            :z => float(s["formula_charge"]),
            :ΔₐG⁰ => _hkf_aloh2_value("sm_gibbs_energy") * u"J/mol",
            :ΔₐH⁰ => _hkf_aloh2_value("sm_enthalpy") * u"J/mol",
            :S⁰ => _hkf_aloh2_value("sm_entropy_abs") * u"J/(mol*K)",
            :T => float(s["Tst"]) * u"K",
            :P => float(s["Pst"]) * u"Pa",
        ],
    )
end


@testsection "surface tension of water (IAPWS R1-76(2014))" begin
    # The equation against the release's own table: its calculated column to the
    # rounding of the table, and the recommended values within their uncertainty.
    t = literature_table("IAPWS2014", "surface_tension")
    for (θ, σ, Δσ, σ_calc) in zip(t.t_C, t.sigma, t.uncertainty, t.sigma_calc)
        γ = water_surface_tension(273.15 + θ) * u"N/m"
        @test abs(ustrip(u"mN/m", γ - σ_calc)) <= 0.005 + 1.0e-9
        @test abs(ustrip(u"mN/m", γ - σ)) <= ustrip(u"mN/m", Δσ)
    end
    # It vanishes at the critical point, falls with temperature, differentiates,
    # and refuses a temperature at which there is no liquid.
    @test water_surface_tension(647.096) == 0
    @test ForwardDiff.derivative(water_surface_tension, 298.15) < 0
    @test_throws DomainError water_surface_tension(700.0)
end

@testsection "HKF water properties (HGK EOS)" begin
    wtp = water_thermo_props(298.15, 1.0e5)

    # Liquid water density at 25°C, 1 bar ≈ 997 kg/m³
    @test isapprox(wtp.D, 997.0; rtol = 1.0e-2)

    # DT should be negative (water expands when heated)
    @test wtp.DT < 0

    # DP should be positive (density increases under pressure)
    @test wtp.DP > 0
end

@testsection "HKF water electrostatics (Johnson-Norton)" begin
    wtp = water_thermo_props(298.15, 1.0e5)
    wep = water_electro_props_jn(298.15, 1.0e5, wtp)

    # Dielectric constant of water at 25°C ≈ 78.4
    @test isapprox(wep.epsilon, 78.4; rtol = 1.0e-2)

    # Born function Z = -1/ε
    @test isapprox(wep.bornZ, -1 / wep.epsilon; rtol = 1.0e-6)

    # Y = εT/ε² should be negative (ε decreases with T)
    @test wep.bornY < 0
end

@testsection "HKF g-function (Shock 1992)" begin
    wtp = water_thermo_props(298.15, 1.0e5)
    gs = hkf_g_function(298.15, 1.0e5, wtp)

    # g ≈ 0 at reference conditions (density ~997, well inside validity range)
    @test isapprox(gs.g, 0.0; atol = 1.0e-2)

    # Outside validity range → zero state
    wtp_bad = WaterThermoProps(1100.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    gs_bad = hkf_g_function(298.15, 1.0e5, wtp_bad)
    @test iszero(gs_bad.g)
end

@testsection "HKF thermo functions — Al(OH)2+" begin
    params = _hkf_aloh2_params()

    thermo = build_thermo_functions(:solute_hkf88_reaktoro, params)

    @test haskey(thermo, :Cp⁰)
    @test haskey(thermo, :ΔₐH⁰)
    @test haskey(thermo, :S⁰)
    @test haskey(thermo, :ΔₐG⁰)
    @test haskey(thermo, :V⁰)

    @test thermo[:Cp⁰] isa NumericFunc
    @test thermo[:ΔₐG⁰] isa NumericFunc

    # vars must be (:T, :P)
    @test thermo[:Cp⁰].vars == (:T, :P)

    # refs must store T and P as Quantities in SI
    @test haskey(thermo[:Cp⁰].refs, :T)
    @test haskey(thermo[:Cp⁰].refs, :P)
    @test thermo[:Cp⁰].refs.T isa AbstractQuantity
    @test thermo[:Cp⁰].refs.P isa AbstractQuantity
    @test ustrip(thermo[:Cp⁰].refs.T) ≈ 298.15
    @test ustrip(thermo[:Cp⁰].refs.P) ≈ 1.0e5

    # Cp at (Tr, Pr) must match the value the database tabulates beside the coefficients
    @test isapprox(thermo[:Cp⁰](; T = 298.15, P = 1.0e5), _hkf_aloh2_value("sm_heat_capacity_p"); rtol = 1.0e-2)

    # Calling without kwargs must use refs and give the same result
    @test isapprox(thermo[:Cp⁰](), thermo[:Cp⁰](; T = 298.15, P = 1.0e5); rtol = 1.0e-10)

    # G at (Tr, Pr) must match Gf
    G_ref = _hkf_aloh2_value("sm_gibbs_energy")
    @test isapprox(thermo[:ΔₐG⁰](; T = 298.15, P = 1.0e5), G_ref; rtol = 1.0e-3)
    @test isapprox(thermo[:ΔₐG⁰](), G_ref; rtol = 1.0e-3)

    # unit=true returns a Quantity
    @test thermo[:Cp⁰](; T = 298.15, P = 1.0e5, unit = true) isa AbstractQuantity
end

@testsection "HKF ForwardDiff compatibility" begin
    params = _hkf_aloh2_params()
    thermo = build_thermo_functions(:solute_hkf88_reaktoro, params)

    # ∂G/∂T = -S (thermodynamic identity)
    dG_dT = ForwardDiff.derivative(T -> thermo[:ΔₐG⁰](; T = T, P = 1.0e5), 298.15)
    S_val = thermo[:S⁰](; T = 298.15, P = 1.0e5)
    @test isapprox(dG_dT, -S_val; rtol = 1.0e-3)
end

@testsection "NumericFunc arithmetic — scalar" begin
    params = _hkf_aloh2_params()
    thermo = build_thermo_functions(:solute_hkf88_reaktoro, params)
    G = thermo[:ΔₐG⁰]
    T0, P0 = 298.15, 1.0e5
    G0 = G(; T = T0, P = P0)

    # Unary negation
    neg = -G
    @test neg isa NumericFunc
    @test neg(; T = T0, P = P0) ≈ -G0

    # Scalar multiplication (stoichiometric coefficient)
    scaled = 2 * G
    @test scaled isa NumericFunc
    @test scaled(; T = T0, P = P0) ≈ 2 * G0

    scaled2 = G * 3.0
    @test scaled2 isa NumericFunc
    @test scaled2(; T = T0, P = P0) ≈ 3.0 * G0

    # Scalar division
    div_hkf = G / 2.0
    @test div_hkf isa NumericFunc
    @test div_hkf(; T = T0, P = P0) ≈ G0 / 2.0

    # Scalar addition/subtraction (offset in same unit)
    offset = 1000.0
    add_hkf = G + offset
    @test add_hkf(; T = T0, P = P0) ≈ G0 + offset
    sub_hkf = G - offset
    @test sub_hkf(; T = T0, P = P0) ≈ G0 - offset
end

@testsection "NumericFunc arithmetic — HKF op HKF" begin
    params = _hkf_aloh2_params()
    thermo = build_thermo_functions(:solute_hkf88_reaktoro, params)
    G = thermo[:ΔₐG⁰]
    H = thermo[:ΔₐH⁰]
    T0, P0 = 298.15, 1.0e5
    G0 = G(; T = T0, P = P0)
    H0 = H(; T = T0, P = P0)

    sum_hkf = G + H
    @test sum_hkf isa NumericFunc
    @test sum_hkf(; T = T0, P = P0) ≈ G0 + H0

    diff_hkf = G - H
    @test diff_hkf isa NumericFunc
    @test diff_hkf(; T = T0, P = P0) ≈ G0 - H0

    # Linear combination (reaction-like)
    lc = 2 * G + (-1) * H
    @test lc isa NumericFunc
    @test lc(; T = T0, P = P0) ≈ 2 * G0 - H0
end

@testsection "NumericFunc arithmetic — cross-type with SymbolicFunc" begin
    params = _hkf_aloh2_params()
    thermo = build_thermo_functions(:solute_hkf88_reaktoro, params)
    G = thermo[:ΔₐG⁰]
    S = thermo[:S⁰]
    T0, P0 = 298.15, 1.0e5
    G0 = G(; T = T0, P = P0)
    S0 = S(; T = T0, P = P0)

    # SymbolicFunc representing a constant (must carry same unit as G for +/-)
    c = _hkf_aloh2_value("sm_gibbs_energy")
    tf_const = SymbolicFunc(c * u"J/mol")

    add_cross = G + tf_const
    @test add_cross isa NumericFunc
    @test add_cross(; T = T0, P = P0) ≈ G0 + c

    diff_cross = tf_const - G
    @test diff_cross isa NumericFunc
    @test diff_cross(; T = T0, P = P0) ≈ c - G0

    # Division of HKF by a SymbolicFunc of T (as in logK computation)
    R_log10 = ustrip(DynamicQuantities.Constants.R) * log(10)
    tf_T = R_log10 * SymbolicFunc(:T; units = [:T => "K"])
    logK = -G / tf_T
    @test logK isa NumericFunc
    # logK = -G / (R*ln10*T), numerical check
    expected_logK = -G0 / (R_log10 * T0)
    @test isapprox(logK(; T = T0, P = P0), expected_logK; rtol = 1.0e-6)
end

# ── AD smoke tests for water properties (direct) ────────────────────────────

@testsection "ForwardDiff — water_thermo_props" begin
    # ∂(density)/∂T at 298.15 K, 1 bar — must be finite (negative for liquid water)
    dD_dT = ForwardDiff.derivative(T -> water_thermo_props(T, 1.0e5).D, 298.15)
    @test isfinite(dD_dT)
    @test dD_dT < 0   # water expands when heated

    # ∂(density)/∂P at 298.15 K, 1 bar — must be positive (compression)
    dD_dP = ForwardDiff.derivative(P -> water_thermo_props(298.15, P).D, 1.0e5)
    @test isfinite(dD_dP)
    @test dD_dP > 0
end

@testsection "ForwardDiff — water_electro_props_jn" begin
    # ∂(epsilon)/∂T — dielectric constant decreases with temperature
    function eps_of_T(T)
        wtp = water_thermo_props(T, 1.0e5)
        wep = water_electro_props_jn(T, 1.0e5, wtp)
        return wep.epsilon
    end
    deps_dT = ForwardDiff.derivative(eps_of_T, 298.15)
    @test isfinite(deps_dT)
    @test deps_dT < 0   # ε decreases with T
end

@testsection "ForwardDiff — hkf_g_function" begin
    # ∂g/∂T at conditions where g ≈ 0 (reference state)
    function g_of_T(T)
        wtp = water_thermo_props(T, 1.0e5)
        gs = hkf_g_function(T, 1.0e5, wtp)
        return gs.g
    end
    dg_dT = ForwardDiff.derivative(g_of_T, 298.15)
    @test isfinite(dg_dT)
end
