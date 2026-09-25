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
