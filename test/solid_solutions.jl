using ChemistryLab
using DynamicQuantities
using ForwardDiff
using Symbolics
using LinearAlgebra
using Test

# ── Helpers ────────────────────────────────────────────────────────────────────

# Build a minimal Species for solid solution end-member tests (no thermodynamic data needed)
_ss_em(sym) = Species(sym; aggregate_state = AS_CRYSTAL, class = SC_SSENDMEMBER)

# ── Tests ──────────────────────────────────────────────────────────────────────

@testsection "SC_SSENDMEMBER enum" begin
    @test SC_SSENDMEMBER isa Class
    @test SC_SSENDMEMBER != SC_COMPONENT
    @test SC_SSENDMEMBER != SC_AQSOLUTE
    @test SC_SSENDMEMBER != SC_UNDEF
    sp = _ss_em("AFm1")
    @test class(sp) == SC_SSENDMEMBER
    @test aggregate_state(sp) == AS_CRYSTAL
end

@testsection "IdealSolidSolutionModel constructor" begin
    m = IdealSolidSolutionModel()
    @test m isa IdealSolidSolutionModel
    @test m isa AbstractSolidSolutionModel
end

@testsection "RedlichKisterModel constructor" begin
    # Default all-zero
    m0 = RedlichKisterModel()
    @test m0 isa RedlichKisterModel{Float64}
    @test m0.a0 == 0.0 && m0.a1 == 0.0 && m0.a2 == 0.0

    # Custom parameters
    m1 = RedlichKisterModel(a0 = 4000.0, a1 = 500.0)
    @test m1.a0 ≈ 4000.0
    @test m1.a1 ≈ 500.0
    @test m1.a2 == 0.0

    # Float32 promotion (all three args must be Float32)
    m32 = RedlichKisterModel(a0 = 1000.0f0, a1 = 200.0f0, a2 = 50.0f0)
    @test m32 isa RedlichKisterModel{Float32}

    # Mixed promotes to Float64
    m_mixed = RedlichKisterModel(a0 = 1000.0f0, a1 = 200.0)
    @test m_mixed isa RedlichKisterModel{Float64}
end

@testsection "SolidSolutionPhase constructor" begin
    em1 = _ss_em("Em1")
    em2 = _ss_em("Em2")

    # Ideal (default)
    ss = SolidSolutionPhase("SS", [em1, em2])
    @test ss isa SolidSolutionPhase
    @test name(ss) == "SS"
    @test length(end_members(ss)) == 2
    @test model(ss) isa IdealSolidSolutionModel

    # Explicit model
    rk = RedlichKisterModel(a0 = 3000.0)
    ss_rk = SolidSolutionPhase("SS_RK", [em1, em2]; model = rk)
    @test model(ss_rk) isa RedlichKisterModel

    # Error: wrong aggregate_state
    em_bad = Species("Bad"; aggregate_state = AS_AQUEOUS, class = SC_SSENDMEMBER)
    @test_throws ErrorException SolidSolutionPhase("Bad", [em_bad])

    # Auto-requalification: SC_COMPONENT end-members are silently promoted to SC_SSENDMEMBER
    em_comp = Species("Comp"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)
    ss_auto = @test_nowarn SolidSolutionPhase("Auto", [em_comp])
    @test class(end_members(ss_auto)[1]) == SC_SSENDMEMBER
    @test class(em_comp) == SC_COMPONENT   # original unchanged

    # Error: Redlich-Kister with != 2 end-members
    em3 = _ss_em("Em3")
    @test_throws ErrorException SolidSolutionPhase(
        "SS3", [em1, em2, em3]; model = RedlichKisterModel()
    )
    @test_throws ErrorException SolidSolutionPhase(
        "SS1", [em1]; model = RedlichKisterModel()
    )
end

@testsection "ChemicalSystem without solid solutions" begin
    em1 = _ss_em("Em1")
    cs = ChemicalSystem([em1])
    # Backward compatibility — solid_solutions field is Nothing
    @test cs.solid_solutions === nothing
    @test isempty(cs.ss_groups)
    @test isempty(cs.idx_ssendmembers)
end

@testsection "ChemicalSystem with solid solutions" begin
    em1 = _ss_em("Em1")
    em2 = _ss_em("Em2")
    em3 = _ss_em("Em3")
    ss1 = SolidSolutionPhase("SS1", [em1, em2])
    ss2 = SolidSolutionPhase("SS2", [em3])

    cs = ChemicalSystem([em1, em2, em3]; solid_solutions = [ss1, ss2])

    @test !isnothing(cs.solid_solutions)
    @test length(cs.solid_solutions) == 2
    @test cs.ss_groups == [[1, 2], [3]]
    @test cs.idx_ssendmembers == [1, 2, 3]

    # All end-members are in idx_crystal (AS_CRYSTAL aggregate state)
    @test Set(cs.idx_ssendmembers) ⊆ Set(cs.idx_crystal)

    # Accessor
    @test solid_solutions(cs) === cs.solid_solutions

    # Error: end-member not in species list
    em_extra = _ss_em("Extra")
    ss_bad = SolidSolutionPhase("Bad", [em_extra])
    @test_throws ErrorException ChemicalSystem([em1, em2]; solid_solutions = [ss_bad])
end

@testsection "Ideal SS binary: ln aᵢ = ln xᵢ" begin
    em1 = _ss_em("Em1")
    em2 = _ss_em("Em2")
    ss = SolidSolutionPhase("SS", [em1, em2])
    cs = ChemicalSystem([em1, em2]; solid_solutions = [ss])

    lna = activity_model(cs, DiluteSolutionModel())
    p = (ΔₐG⁰overRT = zeros(2), T = 298.15, P = 1.0e5, ϵ = 1.0e-30)

    # x₁ = 0.7, x₂ = 0.3
    n = [0.7, 0.3]
    out = lna(n, p)
    @test out[1] ≈ log(0.7) rtol = 1.0e-8
    @test out[2] ≈ log(0.3) rtol = 1.0e-8

    # x₁ = 0.5, x₂ = 0.5 (symmetric)
    n2 = [0.5, 0.5]
    out2 = lna(n2, p)
    @test out2[1] ≈ log(0.5) rtol = 1.0e-8
    @test out2[2] ≈ log(0.5) rtol = 1.0e-8
    @test out2[1] ≈ out2[2] rtol = 1.0e-10
end

@testsection "Ideal SS: pure end-member (x₁ → 1)" begin
    em1 = _ss_em("Em1")
    em2 = _ss_em("Em2")
    ss = SolidSolutionPhase("SS", [em1, em2])
    cs = ChemicalSystem([em1, em2]; solid_solutions = [ss])

    lna = activity_model(cs, DiluteSolutionModel())
    p = (ΔₐG⁰overRT = zeros(2), T = 298.15, P = 1.0e5, ϵ = 1.0e-30)

    # n₂ very small → x₁ ≈ 1, ln a₁ ≈ 0
    n_pure = [1.0 - 1.0e-12, 1.0e-12]
    out = lna(n_pure, p)
    @test out[1] ≈ 0.0 atol = 1.0e-8

    # n₁ very small → x₂ ≈ 1, ln a₂ ≈ 0
    n_pure2 = [1.0e-12, 1.0 - 1.0e-12]
    out2 = lna(n_pure2, p)
    @test out2[2] ≈ 0.0 atol = 1.0e-8
end

@testsection "Redlich-Kister binary: analytical values" begin
    em1 = _ss_em("Em1")
    em2 = _ss_em("Em2")
    a0 = 4000.0  # J/mol
    a1 = 800.0   # J/mol
    rk = RedlichKisterModel(a0 = a0, a1 = a1)
    ss = SolidSolutionPhase("SS_RK", [em1, em2]; model = rk)
    cs = ChemicalSystem([em1, em2]; solid_solutions = [ss])

    T = 298.15
    RT = ChemistryLab.R_GAS * T

    lna = activity_model(cs, DiluteSolutionModel())
    p = (ΔₐG⁰overRT = zeros(2), T = T, P = 1.0e5, ϵ = 1.0e-30)

    # x₁ = 0.3, x₂ = 0.7
    x1, x2 = 0.3, 0.7
    n = [x1, x2]
    out = lna(n, p)

    # Analytical Redlich-Kister: ln γ₁
    expected_lng1 = x2^2 * (a0 + a1 * (3 * x1 - x2)) / RT
    expected_lng2 = x1^2 * (a0 - a1 * (3 * x2 - x1)) / RT
    @test out[1] ≈ log(x1) + expected_lng1 rtol = 1.0e-8
    @test out[2] ≈ log(x2) + expected_lng2 rtol = 1.0e-8

    # At x = 0.5 (symmetric), a1=0 case
    rk0 = RedlichKisterModel(a0 = a0)
    ss0 = SolidSolutionPhase("SS0", [em1, em2]; model = rk0)
    cs0 = ChemicalSystem([em1, em2]; solid_solutions = [ss0])
    lna0 = activity_model(cs0, DiluteSolutionModel())
    out0 = lna0([0.5, 0.5], p)
    @test out0[1] ≈ out0[2] rtol = 1.0e-10   # symmetric at x=0.5 with a1=0
end

@testsection "SS in DiluteSolutionModel (mixed aqueous + SS)" begin
    H2O = Species("H2O"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
    NaCl = Species("NaCl"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    em1 = _ss_em("Em1")
    em2 = _ss_em("Em2")
    ss = SolidSolutionPhase("SS", [em1, em2])

    cs = ChemicalSystem([H2O, NaCl, em1, em2]; solid_solutions = [ss])
    lna = activity_model(cs, DiluteSolutionModel())

    n_w = 55.5
    n_Na = 0.1
    n1, n2 = 0.6, 0.4
    n = [n_w, n_Na, n1, n2]
    p = (ΔₐG⁰overRT = zeros(4), T = 298.15, P = 1.0e5, ϵ = 1.0e-30)
    out = lna(n, p)

    # SS end-members: should match ideal mixing formula
    @test out[3] ≈ log(n1 / (n1 + n2)) rtol = 1.0e-8
    @test out[4] ≈ log(n2 / (n1 + n2)) rtol = 1.0e-8

    # Aqueous species: should be non-zero and unaffected by SS
    @test out[1] != 0.0   # solvent
    @test out[2] != 0.0   # solute
end

@testsection "SS in HKFActivityModel (mixed aqueous + SS)" begin
    H2O = Species("H2O"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
    em1 = _ss_em("Em1")
    em2 = _ss_em("Em2")
    ss = SolidSolutionPhase("SS", [em1, em2])

    cs = ChemicalSystem([H2O, em1, em2]; solid_solutions = [ss])
    lna = activity_model(cs, HKFActivityModel())

    n_w = 55.5
    n1, n2 = 0.6, 0.4
    n = [n_w, n1, n2]
    p = (ΔₐG⁰overRT = zeros(3), T = 298.15, P = 1.0e5, ϵ = 1.0e-30)
    out = lna(n, p)

    # SS end-members correct
    @test out[2] ≈ log(n1 / (n1 + n2)) rtol = 1.0e-8
    @test out[3] ≈ log(n2 / (n1 + n2)) rtol = 1.0e-8
end

@testsection "Gibbs-Duhem local SS (ideal)" begin
    # For an ideal SS phase, Σᵢ xᵢ dln(aᵢ)/dξ = 0 where ξ shifts composition
    em1 = _ss_em("Em1")
    em2 = _ss_em("Em2")
    ss = SolidSolutionPhase("SS", [em1, em2])
    cs = ChemicalSystem([em1, em2]; solid_solutions = [ss])

    lna = activity_model(cs, DiluteSolutionModel())
    p = (ΔₐG⁰overRT = zeros(2), T = 298.15, P = 1.0e5, ϵ = 1.0e-30)

    n0 = [0.7, 0.3]
    # Perturbation: shift n₁ ↑ ε, n₂ ↓ ε (internal composition change, total constant)
    ξ = 1.0e-6
    n⁺ = [0.7 + ξ, 0.3 - ξ]
    n⁻ = [0.7 - ξ, 0.3 + ξ]
    dμ = (lna(n⁺, p) - lna(n⁻, p)) / (2 * ξ)   # finite difference

    # Gibbs-Duhem: Σ nᵢ dμᵢ/dξ ≈ 0 (exact for ideal SS)
    residual = abs(sum(n0 .* dμ)) / max(norm(n0 .* abs.(dμ)), 1.0)
    @test residual < 1.0e-6
end

@testsection "ForwardDiff gradient (SS + aqueous)" begin
    H2O = Species("H2O"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
    Na = Species("Na+"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    em1 = _ss_em("Em1")
    em2 = _ss_em("Em2")
    ss = SolidSolutionPhase("SS", [em1, em2])
    cs = ChemicalSystem([H2O, Na, em1, em2]; solid_solutions = [ss])

    lna = activity_model(cs, DiluteSolutionModel())
    p = (ΔₐG⁰overRT = zeros(4), T = 298.15, P = 1.0e5, ϵ = 1.0e-30)
    n0 = [55.5, 0.1, 0.7, 0.3]

    # Gradient of each component
    for i in 1:4
        g = ForwardDiff.gradient(n -> lna(n, p)[i], n0)
        @test all(isfinite, g)
    end

    # Jacobian of the full lna vector
    J = ForwardDiff.jacobian(n -> lna(n, p), n0)
    @test all(isfinite, J)
end

@testsection "ForwardDiff gradient (pure SS, no aqueous)" begin
    em1 = _ss_em("Em1")
    em2 = _ss_em("Em2")
    ss = SolidSolutionPhase("SS", [em1, em2])
    cs = ChemicalSystem([em1, em2]; solid_solutions = [ss])

    lna = activity_model(cs, DiluteSolutionModel())
    p = (ΔₐG⁰overRT = zeros(2), T = 298.15, P = 1.0e5, ϵ = 1.0e-30)
    n0 = [0.7, 0.3]

    J = ForwardDiff.jacobian(n -> lna(n, p), n0)
    @test all(isfinite, J)
end

@testsection "ForwardDiff gradient (Redlich-Kister)" begin
    em1 = _ss_em("Em1")
    em2 = _ss_em("Em2")
    rk = RedlichKisterModel(a0 = 3000.0, a1 = 500.0, a2 = 100.0)
    ss = SolidSolutionPhase("SS", [em1, em2]; model = rk)
    cs = ChemicalSystem([em1, em2]; solid_solutions = [ss])

    lna = activity_model(cs, DiluteSolutionModel())
    p = (ΔₐG⁰overRT = zeros(2), T = 298.15, P = 1.0e5, ϵ = 1.0e-30)
    n0 = [0.6, 0.4]

    J = ForwardDiff.jacobian(n -> lna(n, p), n0)
    @test all(isfinite, J)
end

@testsection "build_potentials with SS" begin
    em1 = _ss_em("Em1")
    em2 = _ss_em("Em2")
    ss = SolidSolutionPhase("SS", [em1, em2])
    cs = ChemicalSystem([em1, em2]; solid_solutions = [ss])

    μ = build_potentials(cs, DiluteSolutionModel())
    n = [0.7, 0.3]
    ΔaGoT = [-100.0, -110.0]
    p = (ΔₐG⁰overRT = ΔaGoT, T = 298.15, P = 1.0e5, ϵ = 1.0e-30)

    out = μ(n, p)
    # μᵢ/RT = ΔₐG⁰ᵢ/RT + ln aᵢ
    lna_out = activity_model(cs, DiluteSolutionModel())(n, p)
    @test out ≈ ΔaGoT .+ lna_out rtol = 1.0e-10
end

@testsection "RegularSolutionModel" begin
    # Construction and its guards.
    m = RegularSolutionModel([0.0 4000.0; 4000.0 0.0])
    @test m isa RegularSolutionModel
    @test m.W[1, 2] == 4000.0
    @test_throws ArgumentError RegularSolutionModel([0.0 1.0 2.0; 1.0 0.0 3.0])
    @test_throws ArgumentError RegularSolutionModel([0.0 4000.0; 1.0 0.0])

    # The binary limit is exactly Redlich-Kister with a0 = W and a1 = a2 = 0.
    # That is what makes this the right generalization rather than a new model:
    # `ln γ₁ = W x₂²/RT`, `ln γ₂ = W x₁²/RT`.
    T = 298.15
    RT = ChemistryLab.R_GAS * T
    W = 4000.0
    reg = RegularSolutionModel([0.0 W; W 0.0])
    rk = RedlichKisterModel(a0 = W)
    for x1 in (0.05, 0.3, 0.5, 0.7, 0.95)
        x = [x1, 1 - x1]
        for k in 1:2
            @test ChemistryLab._excess_ln_gamma(reg, k, x, T) ≈
                ChemistryLab._excess_ln_gamma(rk, k, x, T) rtol = 1.0e-12
        end
        @test ChemistryLab._excess_ln_gamma(reg, 1, x, T) ≈ W * x[2]^2 / RT rtol = 1.0e-12
        @test ChemistryLab._excess_ln_gamma(reg, 2, x, T) ≈ W * x[1]^2 / RT rtol = 1.0e-12
    end

    # W = 0 is ideal mixing, at any arity.
    ideal3 = RegularSolutionModel(zeros(3, 3))
    for k in 1:3
        @test ChemistryLab._excess_ln_gamma(ideal3, k, [0.2, 0.3, 0.5], T) ≈ 0.0 atol = 1.0e-15
    end

    # Ternary: the closed form, checked against the definition
    # `ln γ_k = Σ_{j≠k} W_kj x_j − Σ_{i<j} W_ij x_i x_j`, all over RT.
    W3 = [0.0 3000.0 -1500.0; 3000.0 0.0 2000.0; -1500.0 2000.0 0.0]
    reg3 = RegularSolutionModel(W3)
    x = [0.2, 0.3, 0.5]
    quad = sum(W3[i, j] * x[i] * x[j] for i in 1:3 for j in (i + 1):3) / RT
    for k in 1:3
        lin = sum(W3[k, j] * x[j] for j in 1:3 if j != k) / RT
        @test ChemistryLab._excess_ln_gamma(reg3, k, x, T) ≈ lin - quad rtol = 1.0e-12
    end

    # Gibbs-Duhem at fixed T: Σ x_k d(ln γ_k) = 0, so Σ x_k ln γ_k must equal
    # G^ex/RT. This is the property that makes the partial derivatives mutually
    # consistent, and it is what a hand-written `ln γ` most often gets wrong.
    gex = sum(W3[i, j] * x[i] * x[j] for i in 1:3 for j in (i + 1):3) / RT
    @test sum(x[k] * ChemistryLab._excess_ln_gamma(reg3, k, x, T) for k in 1:3) ≈
        gex rtol = 1.0e-12

    # A phase of any arity accepts it, unlike Redlich-Kister which is binary.
    substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
    dict = Dict(symbol(s) => s for s in substances)
    six = ["CSHQ-TobD", "CSHQ-TobH", "CSHQ-JenH", "CSHQ-JenD", "KSiOH", "NaSiOH"]
    ss = SolidSolutionPhase(
        "CSHQ", [dict[m] for m in six]; model = RegularSolutionModel(zeros(6, 6))
    )
    @test length(end_members(ss)) == 6
    @test model(ss) isa RegularSolutionModel
    @test_throws ErrorException SolidSolutionPhase(
        "CSHQ", [dict[m] for m in six]; model = RedlichKisterModel(a0 = 1.0)
    )

    # And the shipped file loads. Asserted by NAME rather than by count: the
    # file gains phases as the database is exploited further, and a bare count
    # turns every such addition into a spurious failure that says nothing about
    # what broke.
    ss_all = build_solid_solutions(datapath("solid_solutions.toml"), dict)
    names = Set(p.name for p in ss_all)
    for n in (
            "CSHQ", "C3(AF)S0.84H", "Hydrogarnet", "Ettringite_ss",
            "Hydrotalcite",
            # added in 0.15.0
            "Straetlingite_ss", "AFm_SO4_OH", "AFt_SO4_CO3",
            "Hydrotalcite_AlFe", "MSH",
            # the alkali- and aluminum-bearing C-S-H a blended cement needs
            "CNASH_ss",
        )
        @test n in names
    end
    @test length(ss_all) == length(names)   # no phase declared twice
    @test all(length(end_members(p)) >= 2 for p in ss_all)
end

@testset "a mixing energy that is concave is refused, and says where" begin
    # The Gibbs minimum inside a spinodal is two coexisting compositions, and a
    # formulation with one amount per species cannot hold them. Refused at
    # construction rather than discovered as an uncertifiable answer.
    RT = ChemistryLab.R_GAS * 298.15
    em = [
        Species("Ca2SiO4"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT),
        Species("Ca3Si2O7"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT),
    ]

    # The classical symmetric threshold: a regular solution unmixes above 2RT.
    @test spinodal_interval(RegularSolutionModel([0.0 1.9RT; 1.9RT 0.0]), 2) === nothing
    gap = spinodal_interval(RegularSolutionModel([0.0 2.1RT; 2.1RT 0.0]), 2)
    @test gap !== nothing
    @test gap[1] < 0.5 < gap[2]          # symmetric, so it straddles the middle

    # Ideal mixing is convex everywhere, and so is the `AFm` entry this package
    # ships in `data/solid_solutions.toml` — the check does not refuse our own data.
    @test spinodal_interval(IdealSolidSolutionModel(), 2) === nothing
    @test spinodal_interval(RedlichKisterModel(a0 = 3000.0, a1 = 500.0), 2) === nothing

    # The two AFm/AFt parameter sets of the CEM II study, which are concave.
    g1 = spinodal_interval(RedlichKisterModel(a0 = 0.188RT, a1 = 2.49RT), 2)
    g2 = spinodal_interval(RedlichKisterModel(a0 = 1.67RT, a1 = 0.946RT), 2)
    @test g1 !== nothing && isapprox(g1[1], 0.631; atol = 2.0e-3)
    @test g2 !== nothing && isapprox(g2[2], 0.83; atol = 2.0e-3)

    # More than two end-members: a one-dimensional scan is not the right test.
    @test spinodal_interval(RegularSolutionModel([0.0 3RT; 3RT 0.0]), 3) === nothing

    # The refusal names the interval, and can be waived.
    err = try
        SolidSolutionPhase("gap", em; model = RedlichKisterModel(a0 = 0.188RT, a1 = 2.49RT))
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("CONCAVE", err.msg)
    @test occursin("0.631", err.msg)
    @test occursin("check_convexity", err.msg)

    ss = SolidSolutionPhase(
        "gap", em; model = RedlichKisterModel(a0 = 0.188RT, a1 = 2.49RT),
        check_convexity = false,
    )
    @test name(ss) == "gap"

    # Temperature matters: the criterion is a/RT, so the same parameters can be
    # convex hot and concave cold.
    m = RegularSolutionModel([0.0 2.1RT; 2.1RT 0.0])
    @test spinodal_interval(m, 2; T = 298.15) !== nothing
    @test spinodal_interval(m, 2; T = 400.0) === nothing
end

@testset "the common tangent, against an analytic oracle" begin
    # The pair a binary separates into inside a gap, computed from the model
    # alone -- no chemical system, no solver. Checked against an EQUATION rather
    # than a stored number.
    #
    # For a SYMMETRIC model `g(1-x) = g(x)`, so `g'(1-x) = -g'(x)`, and the
    # common-tangent condition `g'(a) = g'(b) = chord` collapses to `g'(x) = 0`:
    #
    #     ln(x/(1-x)) + A(1-2x) = 0,    A = W/RT.
    #
    # That is the oracle. It also fixes the symmetry `b = 1 - a`, which is a
    # second independent check on the same answer.
    RT = ChemistryLab.R_GAS * 298.15

    @testset "symmetric regular solution" begin
        for A in (2.5, 3.0, 4.0)
            m = RegularSolutionModel([0.0 A * RT; A * RT 0.0])
            ct = common_tangent(m, 2)
            @test ct !== nothing
            a, b = ct
            @test 0 < a < b < 1
            @test isapprox(b, 1 - a; atol = 1.0e-8)        # symmetry
            # The oracle, to machine precision.
            @test abs(log(a / (1 - a)) + A * (1 - 2a)) < 1.0e-9
            @test abs(log(b / (1 - b)) + A * (1 - 2b)) < 1.0e-9
            # And the binodal CONTAINS the spinodal, never the other way round.
            sp = spinodal_interval(m, 2)
            @test sp !== nothing
            @test a < sp[1] && sp[2] < b
        end
    end

    @testset "the defining conditions hold for an asymmetric model" begin
        # No closed form here, so the test is the definition itself: equal
        # slopes, and the slope equal to the chord.
        m = RedlichKisterModel(a0 = 0.188RT, a1 = 2.49RT)
        ct = common_tangent(m, 2)
        @test ct !== nothing
        a, b = ct
        A0, A1 = 0.188, 2.49
        g(x) = x * log(x) + (1 - x) * log(1 - x) +
            x * (1 - x) * (A0 + A1 * (2x - 1))
        d(x) = ForwardDiff.derivative(g, x)
        @test isapprox(d(a), d(b); atol = 1.0e-7)
        @test isapprox(d(a), (g(b) - g(a)) / (b - a); atol = 1.0e-7)
        # The tangent line must lie BELOW the curve between the two points --
        # that is what makes the pair the minimum rather than a stationary point.
        line(x) = g(a) + d(a) * (x - a)
        for x in range(a + 1.0e-3, b - 1.0e-3; length = 25)
            @test line(x) <= g(x) + 1.0e-12
        end
        sp = spinodal_interval(m, 2)
        @test a < sp[1] && sp[2] < b
    end

    @testset "nothing where there is no gap" begin
        @test common_tangent(IdealSolidSolutionModel(), 2) === nothing
        @test common_tangent(RegularSolutionModel([0.0 1.9RT; 1.9RT 0.0]), 2) === nothing
        # More than two end-members: a one-dimensional construction is not the
        # right object, and a guess would be worse than a refusal.
        @test common_tangent(RegularSolutionModel([0.0 3RT; 3RT 0.0]), 3) === nothing
    end

    @testset "nothing rather than an unconverged pair" begin
        # Newton is given one iteration, so it cannot reach the root. The
        # function must say so rather than return where it happened to stop: a
        # pair that is not the common tangent is worse than no pair, because
        # everything downstream treats it as exact.
        m = RegularSolutionModel([0.0 3RT; 3RT 0.0])
        @test common_tangent(m, 2; maxit = 1) === nothing
        # And with enough iterations it converges, so the refusal above is the
        # iteration budget and not a broken model.
        @test common_tangent(m, 2; maxit = 100) !== nothing
    end
end

@testset "the lever rule inside the gap" begin
    # Given an overall composition, how the binary separates. Three properties,
    # each checkable without trusting the implementation.
    RT = ChemistryLab.R_GAS * 298.15
    A = 3.0
    m = RegularSolutionModel([0.0 A * RT; A * RT 0.0])
    xa, xb = common_tangent(m, 2)

    @testset "mass balance is exact" begin
        for x̄ in (0.1, 0.2, 0.35, 0.5, 0.65, 0.8, 0.9)
            r = miscibility_split(m, x̄)
            @test r !== nothing
            @test isapprox(r.f_alpha * r.x_alpha + r.f_beta * r.x_beta, x̄; atol = 1.0e-12)
            @test isapprox(r.f_alpha + r.f_beta, 1.0; atol = 1.0e-14)
            @test 0 <= r.f_alpha <= 1
        end
    end

    @testset "the compositions do not depend on the overall one" begin
        # That is the content of the construction: inside a gap only the
        # PROPORTIONS move, never the two compositions.
        for x̄ in (0.2, 0.5, 0.8)
            r = miscibility_split(m, x̄)
            @test isapprox(r.x_alpha, xa; atol = 1.0e-12)
            @test isapprox(r.x_beta, xb; atol = 1.0e-12)
        end
    end

    @testset "the energy released is positive, symmetric, and largest in the middle" begin
        mid = miscibility_split(m, 0.5).Δg
        left = miscibility_split(m, 0.2).Δg
        right = miscibility_split(m, 0.8).Δg
        @test mid > left > 0
        @test isapprox(left, right; atol = 1.0e-9)      # the model is symmetric
        # And zero outside the pair, where the phase is homogeneous.
        out = miscibility_split(m, 0.02)
        @test out.Δg == 0.0
        @test out.f_alpha == 1.0
        @test out.x_alpha == 0.02
    end

    @testset "nothing where there is no gap" begin
        @test miscibility_split(IdealSolidSolutionModel(), 0.5, 2) === nothing
        @test miscibility_split(RegularSolutionModel([0.0 1.9RT; 1.9RT 0.0]), 0.5, 2) === nothing
    end
end

@testset "a miscibility gap can be represented: `instances`" begin
    # Detection was the subject of the test above; this one is about
    # REPRESENTATION. Inside a spinodal the Gibbs minimum is the common-tangent
    # PAIR, and a formulation with one amount per species can only write that
    # down if the substance appears twice.
    RT = ChemistryLab.R_GAS * 298.15
    em = [
        Species("Ca2SiO4"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT),
        Species("Ca3Si2O7"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT),
    ]
    concave = RedlichKisterModel(a0 = 0.188RT, a1 = 2.49RT)

    @testset "refused where it would only add a null direction" begin
        # Two instances of a convex phase are degenerate: every split of the
        # amount between them has the same energy.
        err = try
            SolidSolutionPhase("ideal", em; instances = 2)
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        @test occursin("CONVEX", err.msg)
        @test occursin("degenerate", err.msg)

        @test_throws ErrorException SolidSolutionPhase(
            "gap", em;
            model = concave, instances = 0
        )
    end

    @testset "accepted, and it carries the convexity waiver with it" begin
        # The same declaration that `instances = 1` refuses.
        ss = SolidSolutionPhase("gap", em; model = concave, instances = 2)
        @test ss.instances == 2
        @test ss.declared == "gap"
        @test name(ss) == "gap"
    end

    @testset "ChemicalSystem builds the second composition" begin
        ss1 = SolidSolutionPhase("gap", em; model = concave, check_convexity = false)
        ss2 = SolidSolutionPhase("gap", em; model = concave, instances = 2)

        cs1 = ChemicalSystem(em; solid_solutions = [ss1])
        cs2 = ChemicalSystem(em; solid_solutions = [ss2])

        # One extra copy of each end-member, under a derived symbol.
        @test length(cs2.species) == length(cs1.species) + length(em)
        syms = symbol.(cs2.species)
        @test "Ca2SiO4#2" in syms && "Ca3Si2O7#2" in syms

        # One substance under two labels: byte-identical composition.
        i = findfirst(==("Ca2SiO4"), syms)
        j = findfirst(==("Ca2SiO4#2"), syms)
        @test atoms(cs2.species[i]) == atoms(cs2.species[j])

        # Two phases, and their groups are DISJOINT -- which is what lets the
        # activity assembly, the mole-fraction fill and the certificate stay
        # unchanged.
        @test length(cs2.solid_solutions) == 2
        @test isempty(intersect(cs2.ss_groups[1], cs2.ss_groups[2]))
        @test name.(cs2.solid_solutions) == ["gap", "gap#2"]

        # Conservation is untouched: the new column is a COPY of one already
        # there, so it adds nothing to the row space and the budget `A n` is
        # unchanged as long as the copy starts empty.
        A = Float64.(cs2.CSM.A)
        @test A[:, j] == A[:, i]
    end

    @testset "instances of one declaration are exempt from the overlap refusal" begin
        # Two phases sharing a composition are normally refused -- that is the
        # C-S-H double-count. A miscibility gap is exactly that overlap, on
        # purpose, so the exemption is by provenance and not by composition.
        ss2 = SolidSolutionPhase("gap", em; model = concave, instances = 2)
        @test ChemicalSystem(em; solid_solutions = [ss2]) isa ChemicalSystem

        # ... and a genuine double-count is still refused, instances or not.
        other = SolidSolutionPhase("other", em; model = concave, check_convexity = false)
        one = SolidSolutionPhase("gap", em; model = concave, check_convexity = false)
        @test_throws ErrorException ChemicalSystem(em; solid_solutions = [one, other])
    end
end

# ── C-(N-)A-S-H, and the overlap that must be refused ────────────────────────

@testsection "one gel, three models: the overlap is refused" begin
    substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
    byname = Dict(symbol(s) => s for s in substances)
    mk(n, ms) = SolidSolutionPhase(n, [byname[m] for m in ms])

    CSHQ_MEMBERS = [
        "CSHQ-TobD", "CSHQ-TobH", "CSHQ-JenH", "CSHQ-JenD",
        "KSiOH", "NaSiOH",
    ]
    ECSH_MEMBERS = ["ECSH1-TobCa", "ECSH1-KSH", "ECSH1-NaSH", "ECSH1-SH"]
    CNASH_MEMBERS = [
        "T2C-CNASHss", "T5C-CNASHss", "TobH-CNASHss",
        "5CA", "5CNA", "INFCA", "INFCN", "INFCNA",
    ]

    all_names = vcat(CSHQ_MEMBERS, ECSH_MEMBERS, CNASH_MEMBERS)
    sp = speciation(substances, all_names; aggregate_state = [AS_AQUEOUS])

    @testset "CNASH_ss is shipped and complete" begin
        # All eight end-members are in the public database; the phase was simply
        # not declared before. It is what carries the Al and the alkalis of a
        # blended cement, which `CSHQ` -- having no aluminum end-member at all --
        # cannot.
        ss = build_solid_solutions(datapath("solid_solutions.toml"), byname)
        cnash = findfirst(p -> ChemistryLab.name(p) == "CNASH_ss", ss)
        @test cnash !== nothing
        @test length(end_members(ss[cnash])) == 8
        @test all(haskey(byname, m) for m in CNASH_MEMBERS)
        # It must carry aluminum, which is the whole reason it exists.
        @test any(haskey(atoms(byname[m]), :Al) for m in CNASH_MEMBERS)
        @test !any(haskey(atoms(byname[m]), :Al) for m in CSHQ_MEMBERS)
    end

    @testset "the overlap is exact, not approximate" begin
        # This is what makes it detectable by composition rather than by name.
        @test atoms(byname["KSiOH"]) == atoms(byname["ECSH1-KSH"])
        @test atoms(byname["KSiOH"]) == atoms(byname["ECSH2-KSH"])
    end

    @testset "one at a time builds, two together are refused" begin
        for (nm, members) in (
                ("CSHQ", CSHQ_MEMBERS), ("ECSH1", ECSH_MEMBERS),
                ("CNASH_ss", CNASH_MEMBERS),
            )
            @test ChemicalSystem(
                sp, CEMDATA_PRIMARIES; solid_solutions = [mk(nm, members)]
            ) isa ChemicalSystem
        end

        err = try
            ChemicalSystem(
                sp, CEMDATA_PRIMARIES;
                solid_solutions = [mk("CSHQ", CSHQ_MEMBERS), mk("ECSH1", ECSH_MEMBERS)],
            )
            nothing
        catch e
            sprint(showerror, e)
        end
        @test err !== nothing
        # The message must name both phases and the shared species, or it sends
        # the reader hunting.
        @test occursin("CSHQ", err)
        @test occursin("ECSH1", err)
        @test occursin("KSiOH", err)
    end

    @testset "two models with no composition in common are refused too" begin
        # THE HOLE THIS CLOSES. `CSHQ` and `CNASH_ss` share no composition -- no
        # end-member of one is a substance of the other -- so the composition
        # test above cannot see them, and until 0.25.1 the pair built a system in
        # which the gel was counted twice. The models are now read from
        # data/gel_models.toml and matched by end-member symbol.
        @test isempty(
            [
                (a, b) for a in CSHQ_MEMBERS, b in CNASH_MEMBERS
                    if atoms(byname[a]) == atoms(byname[b])
            ]
        )
        refusal(phases) = try
            ChemicalSystem(sp, CEMDATA_PRIMARIES; solid_solutions = phases)
            nothing
        catch e
            sprint(showerror, e)
        end

        err = refusal([mk("CSHQ", CSHQ_MEMBERS), mk("CNASH", CNASH_MEMBERS)])
        @test err !== nothing
        @test occursin("`CSHQ`", err) && occursin("`CNASH_ss`", err)
        @test occursin("gel_models.toml", err)

        # Matched by the end-members, not by the name the phase is declared
        # under, and whatever subset of a model is declared.
        err = refusal([mk("gel A", CNASH_MEMBERS[1:3]), mk("gel B", CSHQ_MEMBERS[1:4])])
        @test err !== nothing && occursin("`CNASH_ss`", err)

        err = refusal([mk("CNASH_ss", CNASH_MEMBERS), mk("ECSH1", ECSH_MEMBERS)])
        @test err !== nothing && occursin("`ECSH1`", err)

        # One model split over two declared phases is one model, not two.
        @test ChemicalSystem(
            sp, CEMDATA_PRIMARIES;
            solid_solutions = [mk("CSHQ core", CSHQ_MEMBERS[1:4]), mk("CSHQ alkali", CSHQ_MEMBERS[5:6])],
        ) isa ChemicalSystem

        # The shipped registry covers every C-S-H end-member the database has.
        models = ChemistryLab._gel_models()
        for m in vcat(CSHQ_MEMBERS, CNASH_MEMBERS, ECSH_MEMBERS)
            @test haskey(models, m)
        end
        @test models["INFCA"] == (gel = "C-S-H", model = "CNASH_ss")
    end
end

@testsection "an interaction parameter says which convention it is in" begin
    # THE ERROR THIS CATCHES. The literature writes these coefficients two ways.
    # CEMDATA18 gives the AFm sulfate/hydroxide binary as A₀ = 0.188, A₁ = 2.49
    # in RT UNITS; this package takes J/mol and divides by RT internally. The
    # factor between them is RT ≈ 2478 J/mol at 25 °C, and nothing in a bare
    # `Float64` says which was meant -- so both are accepted and one is wrong by
    # three orders of magnitude.
    #
    # A `Quantity` says. RT units have none to give, which is the point: a caller
    # holding 0.188 has to multiply by RT themselves, and that is where the
    # convention becomes visible.

    @test RedlichKisterModel(; a0 = 4000.0u"J/mol").a0 == 4000.0
    @test RedlichKisterModel(; a0 = 4.0u"kJ/mol").a0 == 4000.0      # converted
    @test RedlichKisterModel(; a0 = 4000.0).a0 == 4000.0            # bare: unchanged
    # An energy is not an energy PER MOLE, and the difference is the whole trap.
    @test_throws DimensionError RedlichKisterModel(; a0 = 4000.0u"J")
    @test_throws DimensionError RedlichKisterModel(; a1 = 1.0u"K")
    # Nothing to convert a dimensionless number into: RT units cannot pass as
    # J/mol by accident.
    @test_throws DimensionError RedlichKisterModel(; a0 = 0.188u"1")

    # The diagonal zeros of a `W` literal arrive dimensionless, because Julia
    # promotes the literal before the constructor sees it. A zero carries no
    # convention to get wrong, so it passes; every nonzero value is still checked.
    @test RegularSolutionModel([0.0 4.0u"kJ/mol"; 4.0u"kJ/mol" 0.0]).W[1, 2] == 4000.0
    @test RegularSolutionModel([0.0 4.0u"kJ/mol"; 4.0u"kJ/mol" 0.0]).W[1, 1] == 0.0
    @test_throws DimensionError RegularSolutionModel([0.0 4.0u"J"; 4.0u"J" 0.0])
    @test RedlichKisterModel(; a0 = 0.0u"K").a0 == 0.0

    # Every existing call still means what it meant.
    @test RedlichKisterModel(; a0 = 20_000.0).a0 == 20_000.0
    @test RegularSolutionModel([0.0 4000.0; 4000.0 0.0]).W[1, 2] == 4000.0
end

@testsection "the written formula and the compiled one are the same formula" begin
    # `thermo_factories.jl` keeps a thermodynamic model as a symbolic expression
    # AND a compiled function: the first makes the law readable, the second runs
    # in a solver's inner loop. Solid solutions had only the second, with the
    # formula transcribed into a docstring beside it -- two copies free to drift.
    #
    # `excess_ln_gamma_expression` is the first form, and this asserts that the
    # two agree. That is what makes the written formula true by construction
    # rather than by proofreading. The compiled path is untouched: it is still
    # `_excess_ln_gamma` that every activity evaluation calls.
    models = (
        RedlichKisterModel(; a0 = 4000.0, a1 = 500.0, a2 = -250.0),
        RegularSolutionModel([0.0 4000.0; 4000.0 0.0]),
        IdealSolidSolutionModel(),
    )
    for m in models, k in 1:2
        expr = excess_ln_gamma_expression(m, k)
        @test expr isa Union{Num, Real}
        f = Symbolics.build_function(
            Symbolics.Num(expr),
            Symbolics.variable(:x, 1), Symbolics.variable(:x, 2), Symbolics.variable(:T);
            expression = Val(false),
        )
        for (x1, T) in ((0.3, 298.15), (0.7, 298.15), (0.5, 350.0), (0.05, 273.15))
            numeric = ChemistryLab._excess_ln_gamma(m, k, [x1, 1 - x1], T)
            @test f(x1, 1 - x1, T) ≈ numeric atol = 1.0e-12
        end
    end
end

@testsection "one substance in two phases is said when a mixing is non-ideal" begin
    substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
    byname = Dict(symbol(s) => s for s in substances)
    sp = speciation(
        substances, ["ettringite", "ettringite30", "tricarboalu03", "ettringite03_ss"];
        aggregate_state = [AS_AQUEOUS],
    )
    so4 = SolidSolutionPhase("AFt_SO4", [byname["ettringite"], byname["ettringite30"]])
    aft = [byname["tricarboalu03"], byname["ettringite03_ss"]]
    p = literature_row("Lothenbach2019", "guggenheim_parameters", "AFt SO4/CO3")
    RT = R_GAS * 298.15
    published = RedlichKisterModel(a0 = p.alpha0 * RT, a1 = p.alpha1 * RT)

    # `ettringite03_ss` is ettringite over three, in a formula rounded to seven
    # digits, and the two Gibbs energies agree to that factor within a few J/mol.
    @test ChemistryLab._composition_ratio(byname["ettringite"], byname["ettringite03_ss"]) ≈ 1 / 3 rtol = 1.0e-6
    @test ChemistryLab._composition_ratio(byname["ettringite"], byname["tricarboalu03"]) === nothing
    @test abs(ChemistryLab._g298(byname["ettringite03_ss"]) - ChemistryLab._g298(byname["ettringite"]) / 3) < 10
    # A species built from its formula alone carries no Gibbs energy to compare.
    @test ChemistryLab._g298(Species("CaO")) === nothing

    # Ideal on both sides the double declaration is harmless, and nothing is said.
    @test_logs min_level = Base.CoreLogging.Warn ChemicalSystem(
        sp, CEMDATA_PRIMARIES; solid_solutions = [so4, SolidSolutionPhase("AFt_SO4_CO3", aft)],
    )
    # With the published model on the binary it is said, naming both phases and
    # both species.
    @test_logs (:warn, r"\"AFt_SO4\" and \"AFt_SO4_CO3\" hold one substance twice: \"ettringite03_ss\" is \"ettringite\" times") ChemicalSystem(
        sp, CEMDATA_PRIMARIES; solid_solutions = [
            so4, SolidSolutionPhase("AFt_SO4_CO3", aft; model = published, check_convexity = false),
        ],
    )
end

@testsection "ideal mixing on sublattices" begin
    substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
    byname = Dict(symbol(s) => s for s in substances)
    CNASH = [
        "T2C-CNASHss", "T5C-CNASHss", "TobH-CNASHss",
        "5CA", "5CNA", "INFCA", "INFCN", "INFCNA",
    ]
    CSH3T = ["CSH3T-TobH", "CSH3T-T5C", "CSH3T-T2C"]
    cnash = sublattice_model("Myers2014:cnash", [byname[n] for n in CNASH])
    csh3t = sublattice_model("Kulik2011:csh3t", [byname[n] for n in CSH3T])

    # ln a of every member at the mole fractions x, through the activity path.
    lna(m, x; ϵ = 0.0) = ChemistryLab._ss_log_activities!(zeros(eltype(x), length(x)), eachindex(x), x, m, 298.15, ϵ)
    # A deterministic composition inside the simplex.
    point(K, j) = (v = [1.0 + ((7 * k + 3 * j) % 11) for k in 1:K]; v ./ sum(v))

    @testset "the published models, read for the database records" begin
        @test cnash.multiplicity == [2, 2, 2, 1, 1, 1]
        @test cnash.rank == 8
        # T5C* and 5CA own no species on any site; INFCN owns five of the
        # nine site positions.
        @test cnash.exponents == [4, 0, 1, 0, 2, 1, 5, 4]
        # The Cemdata18 records of CSH3T are half Kulik's formula units.
        @test csh3t.multiplicity == [0.5, 0.5]
        @test csh3t.exponents == [0.5, 0.0, 0.5]
        @test csh3t.rank == 3
    end

    # One transcription against another: the occupancy (Table 1, Eq. 19) and the
    # activities the papers write out term by term (Appendix B, Eq. 19).
    function written_lna(key, prefix, members, x; scale = 1)
        printed = Dict(
            zip(
                literature_table(key, "$(prefix)_end_members").printed,
                literature_table(key, "$(prefix)_end_members").end_member,
            )
        )
        t = literature_table(key, "$(prefix)_site_mixing_terms")
        out = zeros(length(members))
        for (who, c, sum_of) in zip(t.printed, t.coefficient, t.sum_of_mole_fractions)
            k = findfirst(==(printed[who]), members)
            idx = [findfirst(==(printed[strip(p)]), members) for p in split(sum_of, "+")]
            out[k] += scale * ustrip(c) * log(sum(x[idx]))
        end
        return out
    end
    @testset "Myers Appendix B and Table 1 give the same activities" begin
        for j in 1:5
            x = point(8, j)
            @test lna(cnash, x) ≈ written_lna("Myers2014", "cnash", CNASH, x) rtol = 1.0e-13
            # The fictive activity coefficient, as the paper defines it.
            @test [ChemistryLab._excess_ln_gamma(cnash, k, x, 298.15) for k in 1:8] ≈
                written_lna("Myers2014", "cnash", CNASH, x) .- log.(x) rtol = 1.0e-12
        end
    end
    @testset "Kulik Eq. (19) and the occupancy give the same activities" begin
        for j in 1:5
            x = point(3, j)
            @test lna(csh3t, x) ≈ written_lna("Kulik2011", "csh3t", CSH3T, x; scale = 0.5) rtol = 1.0e-13
        end
    end

    @testset "the bridging vacancies of Myers' Table 1 cover the eight end-members" begin
        nu = literature_table("Myers2014", "cnash_bridging_vacancies")
        @test nu.end_member == literature_table("Myers2014", "cnash_end_members").end_member
        @test nu.printed == literature_table("Myers2014", "cnash_end_members").printed
        ν = ustrip.(nu.bridging_vacancies)
        @test all(v -> 0 <= v <= 1, ν)
        # Eq. (11) on a pure end-member: T2C has every bridging site vacant, a gel
        # of dimers; T5C one in two, pentamers.
        chain(v) = 3 / v - 1
        @test chain(ν[findfirst(==("T2C*"), nu.printed)]) == 2
        @test chain(ν[findfirst(==("T5C*"), nu.printed)]) == 5
    end

    @testset "each formula is the sum of its sites" begin
        for (key, prefix, members) in (("Myers2014", "cnash", CNASH), ("Kulik2011", "csh3t", CSH3T))
            sp = Dict(
                zip(
                    literature_table(key, "$(prefix)_species").species,
                    literature_table(key, "$(prefix)_species").formula,
                )
            )
            shared = literature_table(key, "$(prefix)_shared_sites")
            sites = literature_table(key, "$(prefix)_sites")
            occ = literature_table(key, "$(prefix)_occupancy")
            published = literature_table(key, "$(prefix)_end_members")
            for (em, f) in zip(published.end_member, published.formula)
                total = Dict{Symbol, Float64}()
                # A vacancy is written as an empty formula.
                add!(formula, n) = isempty(formula) || for (el, v) in composition(Formula(formula))
                    total[el] = get(total, el, 0.0) + n * Float64(v)
                end
                for (s, n) in zip(shared.species, shared.multiplicity)
                    add!(s, ustrip(n))
                end
                for (who, site, species) in zip(occ.end_member, occ.site, occ.species)
                    who == em || continue
                    add!(sp[species], ustrip(sites.multiplicity[findfirst(==(site), sites.site)]))
                end
                expected = composition(Formula(f))
                @test Set(keys(total)) == Set(k for (k, v) in expected if !iszero(v))
                @test all(isapprox(total[k], Float64(v); atol = 1.0e-12) for (k, v) in expected)
            end
        end
    end

    @testset "one site counted once is ideal mixing, bit for bit" begin
        m = SublatticeModel([1], ["A" "B" "C"])
        for j in 1:3
            x = point(3, j)
            @test lna(m, x; ϵ = 1.0e-30) == lna(IdealSolidSolutionModel(), x; ϵ = 1.0e-30)
        end
    end

    @testset "a site counted twice is the halved pair, per unit" begin
        # (A,B)₂X on its own unit is the pair AX₀.₅/BX₀.₅ with every energy
        # doubled: the activities are the squares.
        x = point(2, 1)
        @test lna(SublatticeModel([2], ["A" "B"]), x) ≈ 2 .* lna(IdealSolidSolutionModel(), x)
        @test lna(SublatticeModel([1, 1], ["A" "B"; "A" "B"]), x) ≈ 2 .* lna(IdealSolidSolutionModel(), x)
    end

    # The Gibbs energy of mixing, coded apart from the model, per amounts n.
    function nG_mix(m, n)
        tot = sum(n)
        g = zero(eltype(n))
        for s in eachindex(m.multiplicity)
            for i in eachindex(m.species[s])
                N = sum(n[k] for k in eachindex(n) if m.occupancy[s, k] == i; init = zero(eltype(n)))
                N > 0 && (g += m.multiplicity[s] * N * log(N / tot))
            end
        end
        return g
    end
    lna_n(m, n) = lna(m, n ./ sum(n))

    @testset "the activities are the gradient of the mixing energy" begin
        for (m, K) in ((cnash, 8), (csh3t, 3)), j in 1:3
            n = 2.5 .* point(K, j)
            @test ForwardDiff.gradient(v -> nG_mix(m, v), n) ≈ lna_n(m, n) rtol = 1.0e-12
            J = ForwardDiff.jacobian(v -> lna_n(m, v), n)
            @test maximum(abs, J - transpose(J)) < 1.0e-12
            # Gibbs-Duhem, and Euler's relation for a homogeneous energy.
            @test maximum(abs, transpose(J) * n) < 1.0e-12
            @test dot(n, lna_n(m, n)) ≈ nG_mix(m, n) rtol = 1.0e-12
        end
    end

    @testset "a pure member has unit activity, and the dilute slope is e_k" begin
        for m in (cnash, csh3t), k in 1:length(m.exponents)
            x = zeros(length(m.exponents)); x[k] = 1.0
            @test lna(m, x)[k] == 0.0
            base = point(length(m.exponents), 2)
            at(t) = (v = copy(base); v[k] = t; v ./ sum(v))
            slope = (lna(m, at(1.0e-9))[k] - lna(m, at(1.0e-8))[k]) / log(at(1.0e-9)[k] / at(1.0e-8)[k])
            @test slope ≈ m.exponents[k] atol = 1.0e-6
        end
    end

    @testset "dual and symbolic numbers go through" begin
        d = ForwardDiff.Dual{Nothing}.(point(8, 1), 1.0)
        @test lna(cnash, d) isa Vector{<:ForwardDiff.Dual}
        # The symbolic form writes Myers' Eq. (B1) for 5CA.
        k = findfirst(==("5CA"), CNASH)
        expr = excess_ln_gamma_expression(cnash, k, 8)
        f = Symbolics.build_function(
            Symbolics.Num(expr), [Symbolics.variable(:x, i) for i in 1:8]...; expression = Val(false),
        )
        x = point(8, 4)
        @test f(x...) ≈ written_lna("Myers2014", "cnash", CNASH, x)[k] - log(x[k]) rtol = 1.0e-12
    end

    @testset "the shipped phase is the published model" begin
        ss = build_solid_solutions(datapath("solid_solutions.toml"), byname)
        ph = ss[findfirst(p -> ChemistryLab.name(p) == "CNASH_ss", ss)]
        @test model(ph) isa SublatticeModel
        @test model(ph).occupancy == cnash.occupancy
        for nm in ("CSH3T", "ECSH1", "ECSH2")
            p = ss[findfirst(q -> ChemistryLab.name(q) == nm, ss)]
            @test model(p) isa IdealSolidSolutionModel
        end
    end

    @testset "what does not describe the phase is refused" begin
        @test_throws ArgumentError SublatticeModel([1, 1], ["A" "A"; "B" "B"])     # two members alike
        @test_throws ArgumentError SublatticeModel([0, 1], ["A" "B"; "A" "B"])     # an empty site
        @test_throws ErrorException SolidSolutionPhase("CSH3T", [byname[n] for n in CSH3T[1:2]]; model = csh3t)
        # A record whose formula is not a multiple of the published one.
        wrong = Species(Formula("(CaO)1.5(SiO2)1(H2O)2.5"); symbol = "CSH3T-TobH", aggregate_state = AS_CRYSTAL)
        err = try
            sublattice_model("Kulik2011:csh3t", [wrong, byname["CSH3T-T5C"], byname["CSH3T-T2C"]])
            nothing
        catch e
            sprint(showerror, e)
        end
        @test err !== nothing && occursin("not a multiple", err)
        @test_throws ErrorException sublattice_model("Kulik2011:csh3t", [byname["CSHQ-TobH"]])
    end
end

# A mixing model with no convexity argument of its own: `mixing_convexity`
# samples it, and a sample proves nothing. Defined at the top level, a struct
# cannot be declared inside a test set.
struct _UnprovedMixing <: ChemistryLab.AbstractSolidSolutionModel end
ChemistryLab._excess_ln_gamma(::_UnprovedMixing, k::Int, x::AbstractVector, T::Real) = zero(eltype(x))

@testsection "convexity beyond two end-members" begin
    RT = ChemistryLab.R_GAS * 298.15
    W3(a, b, c) = [0.0 a b; a 0.0 c; b c 0.0] .* RT

    @testset "the verdicts, and the bound that proves them" begin
        @test mixing_convexity(IdealSolidSolutionModel(), 4).verdict === :convex
        # λmin(QᵀwQ) ≥ −2 on the tangent space proves a regular model convex:
        # a symmetric repulsive ternary up to 2 RT per pair.
        @test mixing_convexity(RegularSolutionModel(W3(1.0, 1.0, 1.0)), 3).verdict === :convex
        @test mixing_convexity(RegularSolutionModel(W3(1.9, 1.9, 1.9)), 3).verdict === :convex
        # One pair past 2 RT: concave at the middle of that edge.
        c = mixing_convexity(RegularSolutionModel(W3(3.0, 0.0, 0.0)), 3)
        @test c.verdict === :nonconvex
        @test c.witness == [0.5, 0.5, 0.0]
        # On two members the bound is the classical threshold W = 2RT.
        for (w, v) in ((1.99, :convex), (2.01, :nonconvex))
            @test mixing_convexity(RegularSolutionModel([0.0 w * RT; w * RT 0.0]), 2).verdict === v
        end
        # A sublattice model is convex by construction.
        @test mixing_convexity(SublatticeModel([2, 1], ["A" "B" "B"; "A" "A" "B"]), 3).verdict === :convex
        # A model with no argument is sampled.
        @test mixing_convexity(_UnprovedMixing(), 3).verdict === :undecided
    end

    sp = Dict(symbol(s) => s for s in build_species(datapath("slop98-inorganic-thermofun.json"); verbose = false))
    members = [sp["Cal"], sp["Mgs"], sp["Str"]]

    @testset "a ternary with a witness is refused, and admitted with two instances" begin
        # REGRESSION: on 0.25.2 only a binary was checked, so this was accepted.
        err = try
            SolidSolutionPhase("carbonate", members; model = RegularSolutionModel(W3(3.0, 0.0, 0.0)))
            nothing
        catch e
            sprint(showerror, e)
        end
        @test err !== nothing && occursin("CONCAVE", err)
        @test SolidSolutionPhase(
            "carbonate", members; model = RegularSolutionModel(W3(3.0, 0.0, 0.0)), instances = 2,
        ) isa SolidSolutionPhase
        # Two instances of a convex ternary are degenerate, and refused.
        @test_throws ErrorException SolidSolutionPhase(
            "carbonate", members; model = RegularSolutionModel(W3(1.0, 1.0, 1.0)), instances = 2,
        )
    end

    @testset "the scope of a certificate follows the verdict" begin
        # REGRESSION: on 0.25.2 a concave ternary was scoped `:global_minimum`.
        names = split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Mg+2 Sr+2 Cal Mgs Str")
        comps = ["H2O@", "H+", "Ca+2", "Mg+2", "Sr+2", "CO3-2", "Zz"]
        exact = HKFActivityModel(å = 4.0, Ḃ = 0.0, Kₙ = 0.0)   # symmetric by construction
        function scope(m)
            ss = [SolidSolutionPhase("carbonate", members; model = m, check_convexity = false)]
            cs = ChemicalSystem([sp[s] for s in names], comps; solid_solutions = ss)
            st = ChemicalState(cs)
            set_quantity!(st, "H2O@", 1.0u"kg")
            for s in ("Cal", "Mgs", "Str")
                set_quantity!(st, s, 0.01u"mol")
            end
            for s in ("H+", "OH-", "Ca+2", "Mg+2", "Sr+2", "CO3-2", "HCO3-", "CO2@")
                set_quantity!(st, s, 1.0e-6u"mol")
            end
            des = DualEquilibriumSolver(cs, exact)
            return ChemistryLab._certificate_scope(
                des, ChemistryLab._build_params(st), ustrip.(us"mol", st.n), FixedTP(),
            )
        end
        @test first(scope(RegularSolutionModel(W3(1.0, 1.0, 1.0)))) === :global_minimum
        s, why = scope(RegularSolutionModel(W3(3.0, 0.0, 0.0)))
        @test s === :kkt_point && occursin("concave", only(why))
        s, why = scope(_UnprovedMixing())
        @test s === :kkt_point && occursin("could not be decided", only(why))
    end
end

@testsection "a member that is a mixture of two others is reported" begin
    substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
    by = Dict(symbol(s) => s for s in substances)
    ldh = [by[n] for n in ("M4A-OH-LDH", "M6A-OH-LDH", "M8A-OH-LDH")]
    # M6A is the average of M4A and M8A, in composition and to 0.3 J/mol in
    # Gibbs energy: the published ternary of Cemdata18 is degenerate.
    @test_logs (:warn, r"M6A-OH-LDH.*is 0.5 ×.*M4A-OH-LDH") SolidSolutionPhase("MgAl_OH_LDH", ldh)
    @test_logs SolidSolutionPhase("MgAl_OH_LDH", ldh; acknowledge_degenerate = true)
    # An ordering energy the model means is not a degeneracy: T5C of CSH3T
    # (−4.35 kJ/mol) and the siliceous hydrogarnet (−10.4 kJ/mol).
    @test_logs SolidSolutionPhase("CSH3T", [by[n] for n in ("CSH3T-TobH", "CSH3T-T5C", "CSH3T-T2C")])
    @test_logs SolidSolutionPhase("Hydrogarnet_Si", [by[n] for n in ("C3AH6", "C3AS0.41H5.18", "C3AS0.84H4.32")])
    # And no shipped phase warns.
    ss = @test_logs min_level = Base.CoreLogging.Warn build_solid_solutions(datapath("solid_solutions.toml"), by)
    @test "MgAl_OH_LDH" in [ChemistryLab.name(p) for p in ss]
end

@testsection "convexity found on the lattice, and flat sublattice directions" begin
    RT = ChemistryLab.R_GAS * 298.15
    # No pair past 2 RT, and not proved convex either: the witness is found on
    # the lattice of compositions.
    w = [0.0 1.95 1.95; 1.95 0.0 -3.0; 1.95 -3.0 0.0]
    Q = ChemistryLab._tangent_basis(3)
    @test minimum(eigvals(Symmetric(transpose(Q) * w * Q))) < -2
    c = mixing_convexity(RegularSolutionModel(w .* RT), 3)
    @test c.verdict === :nonconvex && occursin("sampled", c.how)

    sp = Dict(symbol(s) => s for s in build_species(datapath("slop98-inorganic-thermofun.json"); verbose = false))
    # A reciprocal set of four members on two sites has rank 3: one direction
    # changes no site fraction. With four different carbonates its Gibbs energy
    # moves, and the phase is accepted.
    recip = SublatticeModel([1, 1], ["A" "A" "B" "B"; "A" "B" "A" "B"])
    @test recip.rank == 3
    @test SolidSolutionPhase("reciprocal", [sp["Cal"], sp["Mgs"], sp["Str"], sp["Arg"]]; model = recip) isa SolidSolutionPhase
    # With the same two substances at the ends of that direction, it is flat,
    # and refused: the minimization would have no unique answer.
    err = try
        SolidSolutionPhase("flat", [sp["Cal"], sp["Arg"], sp["Cal"], sp["Arg"]]; model = recip)
        nothing
    catch e
        sprint(showerror, e)
    end
    @test err !== nothing && occursin("no unique answer", err)

    # Records that are not one multiple of the published units.
    db = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-thermofun.json"); verbose = false))
    whole = Species(Formula("(CaO)2(SiO2)3(H2O)5"); symbol = "CSH3T-TobH", aggregate_state = AS_CRYSTAL)
    err = try
        sublattice_model("Kulik2011:csh3t", [whole, db["CSH3T-T5C"], db["CSH3T-T2C"]])
        nothing
    catch e
        sprint(showerror, e)
    end
    @test err !== nothing && occursin("not one multiple", err)
end
