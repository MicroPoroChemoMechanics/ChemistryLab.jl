using Symbolics
using ForwardDiff
using LinearAlgebra: norm

# ── Shared fixtures ───────────────────────────────────────────────────────────

function _nacl_system()
    h2o = Species("H2O"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
    na = Species("Na+"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    cl = Species("Cl-"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    cs = ChemicalSystem([h2o, na, cl])
    return cs, h2o, na, cl
end

# kg/mol, from the library's atomic masses: the same computation that gives the
# solvent of these systems its molar mass, so the two cannot disagree.
const M_W = ustrip(us"kg/mol", calculate_molar_mass(Dict(:H => 2, :O => 1)))

# The ion size of NaCl (Helgeson et al. 1981, Table 2), the common radius the
# structural tests below give every ion.
const å_NaCl = ChemistryLab._NACL_ION_SIZE

function _moles_from_molality(m_NaCl, n_w)
    return [n_w, m_NaCl * n_w * M_W, m_NaCl * n_w * M_W]
end

# ── _hkf_sigma ────────────────────────────────────────────────────────────────

@testsection "_hkf_sigma" begin
    σ = ChemistryLab._hkf_sigma

    # x = 0: Taylor branch → σ(0) = 1
    @test isapprox(σ(0.0), 1.0; atol = 1.0e-12)

    # Small x: Taylor approximation
    x_small = 1.0e-4
    σ_taylor = 1.0 - 1.5 * x_small + 1.8 * x_small^2
    @test isapprox(σ(x_small), σ_taylor; rtol = 1.0e-8)

    # Large x: exact formula at x = 1
    # σ(1) = 3(1 - 2ln2 - 1/2 + 1) = 3(1.5 - 2ln2)
    σ_exact_1 = 3 * (1 - 2 * log(2) - 0.5 + 1)
    @test isapprox(σ(1.0), σ_exact_1; rtol = 1.0e-8)

    # Monotone decreasing
    @test σ(0.1) > σ(0.5) > σ(1.0) > σ(2.0)

    # Continuity across branch threshold (1e-3): both branches agree to O(x³)
    δ = 1.0e-3
    @test isapprox(σ(0.999 * δ), σ(1.001 * δ); rtol = 1.0e-3)
end

# ── hkf_debye_huckel_params ───────────────────────────────────────────────────

@testsection "hkf_debye_huckel_params" begin
    p = hkf_debye_huckel_params(298.15, 1.0e5)

    # The tabulated values at 25 °C / 1 bar of the LLNL model (Parkhurst & Appelo
    # 2013, p. 118), against the package's own water model
    llnl = literature_table("ParkhurstAppelo2013", "llnl_debye_huckel")
    i25 = findfirst(==(25.0), llnl.temperature_C)
    @test isapprox(p.A, llnl.A[i25]; rtol = 1.0e-3)
    @test isapprox(p.B, llnl.B[i25]; rtol = 1.0e-3)

    # A increases with temperature (water structure breaks down)
    p100 = hkf_debye_huckel_params(373.15, 1.0e5)
    @test p100.A > p.A

    # AD smoke-test: gradient of A with respect to T must be finite and positive
    dA_dT = ForwardDiff.derivative(T -> hkf_debye_huckel_params(T, 1.0e5).A, 298.15)
    @test isfinite(dA_dT)
    @test dA_dT > 0

    # AD smoke-test: gradient of B with respect to T
    dB_dT = ForwardDiff.derivative(T -> hkf_debye_huckel_params(T, 1.0e5).B, 298.15)
    @test isfinite(dB_dT)
end

# ── REJ_HKF table ─────────────────────────────────────────────────────────────

@testsection "REJ_HKF" begin
    # Every row of the table, and nothing else, keyed by its PHREEQC formula.
    t = literature_table("Helgeson1981", "effective_electrostatic_radii")
    @test length(REJ_HKF) == length(t.species)
    for (sp, r) in zip(t.species, t.radius_angstrom)
        @test REJ_HKF[sp] === r
    end
    @test all(>(0), values(REJ_HKF))
end

# ── REJ_CHARGE_DEFAULT table ──────────────────────────────────────────────────

@testsection "REJ_CHARGE_DEFAULT" begin
    t = literature_table("Xu2012", "radius_by_charge")
    @test sort!(collect(keys(REJ_CHARGE_DEFAULT))) == sort!(Int.(t.charge))
    for (z, r) in zip(t.charge, t.radius_angstrom)
        @test REJ_CHARGE_DEFAULT[Int(z)] === r
        @test r > 0
    end
    # Roughly monotone: larger |z| → larger radius
    @test REJ_CHARGE_DEFAULT[2] < REJ_CHARGE_DEFAULT[3] < REJ_CHARGE_DEFAULT[4]
    @test REJ_CHARGE_DEFAULT[-1] < REJ_CHARGE_DEFAULT[-2]
end

# ── HKFActivityModel constructors ─────────────────────────────────────────────

@testsection "HKFActivityModel constructors" begin
    m = HKFActivityModel()
    # Every default is the value its source gives.
    llnl = literature_table("ParkhurstAppelo2013", "llnl_debye_huckel")
    i25 = findfirst(==(25.0), llnl.temperature_C)
    @test m.A === llnl.A[i25]
    @test m.B === llnl.B[i25]
    @test m.Ḃ === llnl.B_dot[i25]
    @test m.Kₙ === literature_value("ParkhurstAppelo2013", "uncharged_log_gamma_coefficient")
    @test m.å_default === literature_value("Helgeson1981", "nacl_distance_of_closest_approach")
    @test !m.temperature_dependent

    m2 = HKFActivityModel(A = 0.52, B = 0.33, temperature_dependent = true)
    @test m2.A ≈ 0.52
    @test m2.B ≈ 0.33
    @test m2.temperature_dependent

    # Numeric type promotion: all explicit args must be Float32 for the struct to be Float32
    m3 = HKFActivityModel(A = 0.51f0, B = 0.33f0, Ḃ = 0.041f0, Kₙ = 0.1f0, å_default = 3.72f0)
    @test m3 isa HKFActivityModel{Float32}
end

# ── DaviesActivityModel constructors ──────────────────────────────────────────

@testsection "DaviesActivityModel constructors" begin
    m = DaviesActivityModel()
    @test m.A === ChemistryLab._DH_A_25C
    @test m.b ≈ 0.3
    @test m.bₙ === literature_value("ParkhurstAppelo2013", "uncharged_log_gamma_coefficient")
    @test !m.temperature_dependent

    m2 = DaviesActivityModel(b = 0.2, temperature_dependent = true)
    @test m2.b ≈ 0.2
    @test m2.temperature_dependent
end

# ── HKF: dilution limit (I → 0) ──────────────────────────────────────────────

@testsection "HKF: dilution limit" begin
    cs, h2o, na, cl = _nacl_system()
    model = HKFActivityModel()
    lna = activity_model(cs, model)

    # Very dilute: 1 kg water (≈ 55.5 mol), 1e-6 mol NaCl
    n_w = 1.0 / M_W
    n = [n_w, 1.0e-6, 1.0e-6]
    p = (ΔₐG⁰overRT = zeros(3), T = 298.15, P = 1.0e5, ϵ = 1.0e-30)
    out = lna(n, p)

    # ln aᵢ ≈ ln(mᵢ) as I → 0 (activity coefficient → 1)
    m_ion = 1.0e-6 / (n_w * M_W)
    @test isapprox(out[2], log(m_ion); rtol = 1.0e-3)
    @test isapprox(out[3], log(m_ion); rtol = 1.0e-3)

    # Water activity → 1 (ln a_w → 0) at infinite dilution
    @test isapprox(out[1], 0.0; atol = 1.0e-4)
end

# ── HKF: NaCl 0.1 mol/kg (analytical verification) ───────────────────────────

@testsection "HKF: NaCl 0.1 mol/kg" begin
    cs, h2o, na, cl = _nacl_system()
    model = HKFActivityModel()
    lna = activity_model(cs, model)

    m = 0.1
    n_w = 1.0 / M_W
    n = _moles_from_molality(m, n_w)
    p = (ΔₐG⁰overRT = zeros(3), T = 298.15, P = 1.0e5, ϵ = 1.0e-30)
    out = lna(n, p)

    A = ChemistryLab._DH_A_25C
    B = ChemistryLab._DH_B_25C
    Ḃ = ChemistryLab._BDOT_25C
    I = m   # NaCl 1:1 electrolyte: I = ½ × 2m = m
    sqI = sqrt(I)
    ln10 = log(10.0)

    for (i, (åᵢ, key)) in enumerate([(REJ_HKF["Na+"], 2), (REJ_HKF["Cl-"], 3)])
        log10γ = -A * sqI / (1 + B * åᵢ * sqI) + Ḃ * I   # z=1
        @test isapprox(out[key], ln10 * log10γ + log(m); rtol = 1.0e-6)
    end

    # Activity coefficient < 1 at I = 0.1 mol/kg (electrostatic attraction dominates)
    γ_Na = exp(out[2] - log(m))
    @test γ_Na < 1.0
end

# ── HKF: neutral species salting-out ─────────────────────────────────────────

@testsection "HKF: neutral CO2 salting-out" begin
    h2o = Species("H2O"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
    na = Species("Na+"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    cl = Species("Cl-"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    co2 = Species("CO2"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)   # neutral

    cs = ChemicalSystem([h2o, na, cl, co2])
    model = HKFActivityModel(Kₙ = 0.1)
    lna = activity_model(cs, model)

    m_NaCl = 0.5
    m_co2 = 0.01
    n_w = 1.0 / M_W
    n = [n_w, m_NaCl * n_w * M_W, m_NaCl * n_w * M_W, m_co2 * n_w * M_W]
    p = (ΔₐG⁰overRT = zeros(4), T = 298.15, P = 1.0e5, ϵ = 1.0e-30)
    out = lna(n, p)

    I = m_NaCl   # 1:1 electrolyte
    ln10 = log(10.0)
    @test isapprox(out[4], ln10 * 0.1 * I + log(m_co2); rtol = 1.0e-6)
end

# ── HKF: ionic radius priority lookup ────────────────────────────────────────

@testsection "HKF: ionic radius priority" begin
    model = HKFActivityModel()

    # Priority 1: sp[:å] — explicit value overrides REJ_HKF
    sp_custom = Species("Na+"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    sp_custom[:å] = 5.0
    @test ChemistryLab._hkf_lookup_å(sp_custom, model) ≈ 5.0

    # Priority 2: REJ_HKF — for well-known species
    sp_normal = Species("Na+"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    @test ChemistryLab._hkf_lookup_å(sp_normal, model) ≈ 1.91

    sp_so4 = Species("SO4-2"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    @test ChemistryLab._hkf_lookup_å(sp_so4, model) ≈ 3.15

    # Priority 4: å_default — species with unlisted charge (charge 0, which is neutral)
    # For neutral species zv[i]=0, _hkf_lookup_å is not called for neutrals in activity_model.
    # We verify the model fallback for a species with e.g. charge +5 not in REJ_CHARGE_DEFAULT.
    # We construct such a species via Dict-based Formula constructor.
    sp_exotic = Species(
        Dict(:Zz => 1);
        name = "ExoticIon",
        symbol = "Zz+5",
        aggregate_state = AS_AQUEOUS,
        class = SC_AQSOLUTE,
    )
    # charge(sp_exotic) must be +5 for fallback to å_default
    # If Formula parsing gives a different result, we check REJ_CHARGE_DEFAULT[Int(charge(sp_exotic))]
    z_exotic = Int(charge(sp_exotic))
    expected = haskey(REJ_CHARGE_DEFAULT, z_exotic) ?
        REJ_CHARGE_DEFAULT[z_exotic] : model.å_default
    @test ChemistryLab._hkf_lookup_å(sp_exotic, model) ≈ expected
end

# ── HKF: fallback by charge ±2 ───────────────────────────────────────────────

@testsection "HKF: fallback charge +2" begin
    # Use Ca+2 then modify formula to a fictitious species (not in REJ_HKF)
    # Use Mg+2 already in REJ_HKF for comparison
    h2o = Species("H2O"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
    ca = Species("Ca+2"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    # verify Ca+2 IS in REJ_HKF (priority 2)
    model = HKFActivityModel()
    @test ChemistryLab._hkf_lookup_å(ca, model) ≈ REJ_HKF["Ca+2"]

    # Build a system with a +2 species not in REJ_HKF.
    # We use formula "Sr+2" which IS in REJ_HKF — let's verify Sr+2 lookup
    sr = Species("Sr+2"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    @test ChemistryLab._hkf_lookup_å(sr, model) ≈ REJ_HKF["Sr+2"]
    @test REJ_HKF["Sr+2"] === literature_row("Helgeson1981", "effective_electrostatic_radii", "Sr+2").radius_angstrom
end

# ── HKF: temperature-dependent mode ──────────────────────────────────────────

@testsection "HKF: temperature-dependent mode" begin
    cs, h2o, na, cl = _nacl_system()
    model_fixed = HKFActivityModel(temperature_dependent = false)
    model_tdep = HKFActivityModel(temperature_dependent = true)
    lna_fixed = activity_model(cs, model_fixed)
    lna_tdep = activity_model(cs, model_tdep)

    m = 0.5
    n_w = 1.0 / M_W
    n = _moles_from_molality(m, n_w)

    p25 = (ΔₐG⁰overRT = zeros(3), T = 298.15, P = 1.0e5, ϵ = 1.0e-30)
    p100 = (ΔₐG⁰overRT = zeros(3), T = 373.15, P = 1.0e5, ϵ = 1.0e-30)

    # At 25 °C, fixed A/B ≈ T-dependent (same values)
    @test isapprox(lna_fixed(n, p25)[2], lna_tdep(n, p25)[2]; rtol = 1.0e-3)

    # At 100 °C, A is larger → T-dependent gives different result from fixed
    @test !isapprox(lna_fixed(n, p100)[2], lna_tdep(n, p100)[2]; rtol = 1.0e-2)

    # T-dependent mode produces higher (more negative) activity coefficients at 100 °C
    # (A increases, so |log10γ| is larger)
    @test abs(lna_tdep(n, p100)[2] - log(m)) > abs(lna_fixed(n, p100)[2] - log(m))
end

# ── Gibbs-Duhem consistency ───────────────────────────────────────────────────

@testsection "Gibbs-Duhem consistency (HKF)" begin
    cs, h2o, na, cl = _nacl_system()
    model = HKFActivityModel()
    μ = build_potentials(cs, model)

    m = 0.3
    n_w = 1.0 / M_W
    n0 = _moles_from_molality(m, n_w)
    p = (ΔₐG⁰overRT = zeros(3), T = 298.15, P = 1.0e5, ϵ = 1.0e-30)

    # Finite-difference perturbation along NaCl dissolution direction
    δ = 1.0e-6
    dn = [0.0, δ, -δ]
    μ0 = μ(n0, p)
    μ1 = μ(n0 + dn, p)
    dμ = (μ1 - μ0) / δ

    # Gibbs-Duhem: Σᵢ nᵢ ∂μᵢ/∂ξ ≈ 0
    # Note: the B-dot model uses a charge-weighted average å for the osmotic coefficient,
    # which is an approximation when ion radii differ (Na⁺ å=1.91 Å vs Cl⁻ å=1.81 Å).
    # This introduces a small physical inconsistency (~0.2% for m=0.3 mol/kg NaCl),
    # so the tolerance is relaxed to 5e-3.
    residual = abs(sum(n0 .* dμ)) / max(norm(n0 .* abs.(dμ)), 1.0)
    @test residual < 5.0e-3
end

# ── Gas phase: ideal mixture ──────────────────────────────────────────────────

@testsection "Gas phase: ideal mixture (HKF)" begin
    h2o = Species("H2O"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
    co2 = Species("CO2g"; aggregate_state = AS_GAS, class = SC_GASFLUID)
    n2 = Species("N2g"; aggregate_state = AS_GAS, class = SC_GASFLUID)
    cs = ChemicalSystem([h2o, co2, n2])
    model = HKFActivityModel()
    lna = activity_model(cs, model)

    n_w = 1.0 / M_W
    n = [n_w, 0.3, 0.7]
    p = (ΔₐG⁰overRT = zeros(3), T = 298.15, P = 1.0e5, ϵ = 1.0e-30)
    out = lna(n, p)

    @test isapprox(out[2], log(0.3 / 1.0); rtol = 1.0e-8)   # ln(x_CO2)
    @test isapprox(out[3], log(0.7 / 1.0); rtol = 1.0e-8)   # ln(x_N2)
end

# ── build_potentials with HKF ─────────────────────────────────────────────────

@testsection "build_potentials with HKF" begin
    cs, h2o, na, cl = _nacl_system()
    model = HKFActivityModel()
    μ = build_potentials(cs, model)
    lna = activity_model(cs, model)

    m = 0.1
    n_w = 1.0 / M_W
    n = _moles_from_molality(m, n_w)
    ΔG = [-95.0, -4.0, -2.0]
    p = (ΔₐG⁰overRT = ΔG, T = 298.15, P = 1.0e5, ϵ = 1.0e-30)

    expected = ΔG .+ lna(n, p)
    result = μ(n, p)
    @test isapprox(result, expected; rtol = 1.0e-10)
end

# ── ForwardDiff AD compatibility ──────────────────────────────────────────────

@testsection "ForwardDiff AD compatibility" begin
    # _hkf_sigma: gradient finite at x=0 (Taylor branch) and x=1 (exact branch)
    σ = ChemistryLab._hkf_sigma
    gσ_0 = ForwardDiff.derivative(σ, 0.0)
    gσ_1 = ForwardDiff.derivative(σ, 1.0)
    @test isfinite(gσ_0)
    @test isfinite(gσ_1)
    # σ'(0) = -3/2 from Taylor: dσ/dx|₀ = -(3/2)
    @test isapprox(gσ_0, -1.5; atol = 1.0e-8)

    # hkf_debye_huckel_params: gradient of A w.r.t. T
    dA_dT = ForwardDiff.derivative(T -> hkf_debye_huckel_params(T, 1.0e5).A, 298.15)
    @test isfinite(dA_dT)

    # activity_model(HKF): gradient of ln(a_Na) w.r.t. mole vector
    cs, h2o, na, cl = _nacl_system()
    model = HKFActivityModel()
    lna = activity_model(cs, model)

    m = 0.3
    n_w = 1.0 / M_W
    n0 = _moles_from_molality(m, n_w)
    p = (ΔₐG⁰overRT = zeros(3), T = 298.15, P = 1.0e5, ϵ = 1.0e-30)

    g = ForwardDiff.gradient(n -> lna(n, p)[2], n0)
    @test all(isfinite, g)

    # Gradient of ln(a_water) w.r.t. moles also finite
    g_w = ForwardDiff.gradient(n -> lna(n, p)[1], n0)
    @test all(isfinite, g_w)
end

# ── Davies model basic verification ──────────────────────────────────────────

@testsection "DaviesActivityModel basic" begin
    cs, h2o, na, cl = _nacl_system()
    model = DaviesActivityModel()
    lna = activity_model(cs, model)

    m = 0.1
    n_w = 1.0 / M_W
    n = _moles_from_molality(m, n_w)
    p = (ΔₐG⁰overRT = zeros(3), T = 298.15, P = 1.0e5, ϵ = 1.0e-30)
    out = lna(n, p)

    # Verify analytical Davies formula for NaCl: z=1, I=m
    A = ChemistryLab._DH_A_25C
    I = m
    sqI = sqrt(I)
    ln10 = log(10.0)
    log10γ = -A * (sqI / (1 + sqI) - 0.3 * I)   # b=0.3, z=1
    @test isapprox(out[2], ln10 * log10γ + log(m); rtol = 1.0e-6)
    @test isapprox(out[3], ln10 * log10γ + log(m); rtol = 1.0e-6)

    # Activity coefficient < 1
    γ = exp(out[2] - log(m))
    @test γ < 1.0

    # AD compatibility
    g = ForwardDiff.gradient(n -> lna(n, p)[2], n)
    @test all(isfinite, g)
end

@testsection "which activity models come from one excess Gibbs energy" begin
    # TWO RELATIONS EVERY MODEL DERIVED FROM A SINGLE G^E MUST SATISFY, measured
    # by automatic differentiation of `ln_activities` alone -- no `G` is built,
    # no solver runs, and nothing is differenced against a second solve.
    #
    #   symmetry     ∂ln aᵢ/∂nⱼ = ∂ln aⱼ/∂nᵢ     (Maxwell: a potential exists)
    #   Gibbs-Duhem  Σᵢ nᵢ ∂ln aᵢ/∂nⱼ = 0        (μᵢ homogeneous of degree 0)
    #
    # Symmetry is the sharper of the two: it fails exactly when no function of
    # `n` has these activities for its gradient, whatever that function is.
    #
    # Both are RELATIVE measures. An absolute residual is meaningless here: the
    # Jacobian holds self-derivatives of trace species of order 1/n, so dividing
    # by its largest entry makes every real inconsistency look like zero. That
    # mistake was made once and reported 1e-8 for a model that is off by 20 %.

    subs = build_species(datapath("slop98-inorganic-thermofun.json"); verbose = false)
    d = Dict(symbol(s) => s for s in subs)
    cs = ChemicalSystem(
        [d[s] for s in ["H2O@", "Na+", "Cl-", "Ca+2", "SO4-2"]],
        ["H2O@", "Na+", "Cl-", "Ca+2", "SO4-2", "Zz"],
    )
    nm = symbol.(cs.species)
    amounts = Dict("H2O@" => 55.5, "Na+" => 0.1, "Cl-" => 0.12, "Ca+2" => 0.02, "SO4-2" => 0.01)
    n0 = [get(amounts, s, 1.0e-10) for s in nm]
    state = ChemicalState(cs, [x * u"mol" for x in n0])
    p = ChemistryLab._build_params(state)
    i_w = findfirst(==("H2O@"), nm)

    # Split ion/ion from solvent/ion: they have different causes, and merging
    # them is what hid the second one.
    function asymmetries(model)
        J = ForwardDiff.jacobian(nn -> activity_model(cs, model)(nn, p), n0)
        ion, solvent = 0.0, 0.0
        for i in eachindex(n0), j in (i + 1):length(n0)
            scale = max(abs(J[i, j]), abs(J[j, i]))
            scale > 1.0e-30 || continue
            r = abs(J[i, j] - J[j, i]) / scale
            (i == i_w || j == i_w) ? (solvent = max(solvent, r)) : (ion = max(ion, r))
        end
        gd = 0.0
        for j in eachindex(n0)
            terms = [n0[i] * J[i, j] for i in eachindex(n0)]
            s = sum(abs, terms)
            s > 1.0e-30 || continue
            gd = max(gd, abs(sum(terms)) / s)
        end
        return (; ion, solvent, gd)
    end

    # THE DEBYE-HÜCKEL LIMITING LAW IS EXACT, which is the classical result and
    # the anchor for everything below: it does come from one excess energy.
    limiting = asymmetries(HKFActivityModel(; å = 0.0, Ḃ = 0.0))
    @test limiting.ion < 1.0e-12
    @test limiting.solvent < 1.0e-10
    @test limiting.gd < 1.0e-10

    # AND SO IS A COMMON ION SIZE WITH NO EXTENDED TERM.
    common = asymmetries(HKFActivityModel(; å = å_NaCl, Ḃ = 0.0))
    @test common.ion < 1.0e-12
    @test common.solvent < 1.0e-10

    # WHAT BREAKS IT, and it takes BOTH to be switched off. Neither alone is
    # enough, and dropping the extended term while keeping ion-specific radii
    # makes the ion/ion asymmetry WORSE, not better.
    #
    #   ∂ln γᵢ/∂nⱼ = ln10 · f′(I; zᵢ, åᵢ) · zⱼ²/(2·kg)
    #
    # so symmetry demands that f′(I; zᵢ, åᵢ)/zᵢ² not depend on i. An ion-specific
    # `å` puts i in the Debye-Hückel denominator, and the extended term `Ḃ·I` is
    # added with the SAME coefficient to every ion, so it contributes `Ḃ` rather
    # than `Ḃ·zᵢ²`. Both are properties of the published extended form, which
    # every geochemical code uses; they are approximations with a stated domain,
    # not defects of this implementation.
    @test asymmetries(HKFActivityModel(; å = å_NaCl)).ion > 1.0e-2       # Ḃ alone
    @test asymmetries(HKFActivityModel(; Ḃ = 0.0)).ion > 1.0e-2        # å alone

    # PITZER IS THE CONTROL, and it holds to machine precision -- γ and the
    # osmotic coefficient are derivatives of one virial expansion, which is
    # exactly what `pitzer.jl` claims for it.
    cs_nacl = ChemicalSystem([d[s] for s in ["H2O@", "Na+", "Cl-"]], ["H2O@", "Na+", "Cl-"])
    nm2 = symbol.(cs_nacl.species)
    m0 = [get(Dict("H2O@" => 55.5, "Na+" => 0.1, "Cl-" => 0.1), s, 1.0e-10) for s in nm2]
    st2 = ChemicalState(cs_nacl, [x * u"mol" for x in m0])
    p2 = ChemistryLab._build_params(st2)
    J = ForwardDiff.jacobian(
        nn -> activity_model(
            cs_nacl,
            PitzerActivityModel(; parameters = build_pitzer_parameters(datapath("pitzer-reardon1990.toml"))),
        )(nn, p2),
        m0,
    )
    pz_sym = maximum(
        abs(J[i, j] - J[j, i]) / max(abs(J[i, j]), abs(J[j, i]))
            for i in eachindex(m0) for j in (i + 1):length(m0)
            if max(abs(J[i, j]), abs(J[j, i])) > 1.0e-30
    )
    @test pz_sym < 1.0e-12
end

@testsection "the written activity kernel and the compiled one are the same kernel" begin
    # The pair `thermo_factories.jl` keeps for thermodynamic models -- a symbolic
    # expression to read and a compiled function to run -- extended to the
    # activity kernels. The compiled path is untouched: `_log10γ_ion` is still
    # what every evaluation calls. What this asserts is that the expression says
    # the same thing, so a formula quoted in the theory page is true by
    # construction rather than by proofreading.
    Ivar = Symbolics.variable(:I)
    cases = (
        (HKFActivityModel(), -2, 4.5),
        (HKFActivityModel(; å = å_NaCl), 1, å_NaCl),
        (HKFActivityModel(; Ḃ = 0.0), 3, 9.0),
        (DaviesActivityModel(), -1, 0.0),
    )
    for (model, z, å) in cases
        expr = log10_gamma_expression(model, z, å)
        @test expr isa Num
        f = Symbolics.build_function(expr, Ivar; expression = Val(false))
        for I in (1.0e-6, 0.01, 0.1, 0.5, 1.0)
            A = model.A
            B = model isa HKFActivityModel ? model.B : zero(A)
            numeric = ChemistryLab._log10γ_ion(model, z, å, I, sqrt(I), A, B)
            @test f(I) ≈ numeric atol = 1.0e-12
        end
    end

    # AND THE DERIVATION OF §5b, SYMBOLICALLY. Maxwell demands that
    # `f'(I; zᵢ, åᵢ) / zᵢ²` not depend on the ion. With a common `å` and no
    # extended term it does not -- the ratio is one expression for every charge,
    # and their difference simplifies to zero, which is stronger than a small
    # residual. With either correction in place it does, which is why the
    # measured asymmetry is what it is.
    D = Symbolics.Differential(Ivar)
    ratio(m, z, å) = Symbolics.expand_derivatives(D(log10_gamma_expression(m, z, å))) / z^2

    # The difference is EVALUATED rather than simplified to zero. `simplify`
    # does not reduce these to `0` -- the parameters are floating point and it
    # declines to cancel them -- and asserting on what a simplifier happens to
    # reach would test the simplifier. Evaluating the derivative expression at
    # several ionic strengths tests the identity.
    function gap(m, z1, å1, z2, å2)
        d = Symbolics.build_function(
            ratio(m, z1, å1) - ratio(m, z2, å2), Ivar; expression = Val(false)
        )
        return maximum(abs(d(I)) for I in (1.0e-4, 0.01, 0.1, 0.5, 1.0))
    end

    exact = HKFActivityModel(; å = å_NaCl, Ḃ = 0.0)
    @test gap(exact, 1, å_NaCl, 2, å_NaCl) < 1.0e-14
    @test gap(exact, 1, å_NaCl, 3, å_NaCl) < 1.0e-14

    # The extended term alone breaks it: it contributes `Ḃ`, not `Ḃ·z²`.
    @test gap(HKFActivityModel(; å = å_NaCl), 1, å_NaCl, 2, å_NaCl) > 1.0e-3
    # And so does an ion-specific size, with no extended term at all.
    @test gap(HKFActivityModel(; Ḃ = 0.0), 1, 3.0, 1, 5.0) > 1.0e-3
end
