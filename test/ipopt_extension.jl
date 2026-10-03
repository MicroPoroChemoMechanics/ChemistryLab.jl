# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# The Ipopt back end, an extension. Loading it registers a solver factory, which
# joins the starting points of `equilibrate_certified`: this file runs after every
# other test, whose results must not depend on whether Ipopt is loaded. The
# extension is triggered by `Optimization` and `OptimizationIpopt` together.

using ChemistryLab
using DynamicQuantities
using ForwardDiff
using Optimization
using OptimizationIpopt
using Test

include("reference_species.jl")

const IPOPT_EXT = Base.get_extension(ChemistryLab, :OptimizationIpoptExt)

@testset "the extension registers its back end without taking the default" begin
    @test IPOPT_EXT !== nothing
    @test IPOPT_EXT._default_ipopt_solver in ChemistryLab._SOLVER_FACTORIES
    # OptimaSolver is loaded as well, and its back end stays the default one.
    @test ChemistryLab._DEFAULT_SOLVER_FACTORY[] !== IPOPT_EXT._default_ipopt_solver
    @test IPOPT_EXT._default_ipopt_solver() isa IpoptOptimizer
    # Alone, it would be the default: its `__init__` fills an empty slot.
    saved = ChemistryLab._DEFAULT_SOLVER_FACTORY[]
    try
        ChemistryLab._DEFAULT_SOLVER_FACTORY[] = nothing
        IPOPT_EXT.__init__()
        @test ChemistryLab._DEFAULT_SOLVER_FACTORY[] === IPOPT_EXT._default_ipopt_solver
        @test count(==(IPOPT_EXT._default_ipopt_solver), ChemistryLab._SOLVER_FACTORIES) == 1
    finally
        ChemistryLab._DEFAULT_SOLVER_FACTORY[] = saved
    end
end

@testset "Ipopt reaches the certified equilibrium, in both variable spaces" begin
    # Calcite and carbon dioxide in water, ideal activities.
    db = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-thermofun.json"); verbose = false))
    syms = ["H2O@", "H+", "OH-", "CO2@", "HCO3-", "CO3-2", "Ca+2", "CaOH+", "Ca(CO3)@", "Ca(HCO3)+", "Cal"]
    cs = ChemicalSystem([db[s] for s in syms], ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"])
    st = ChemicalState(cs)
    set_quantity!(st, "H2O@", moles_of_water() * u"mol")
    set_quantity!(st, "Cal", 0.05u"mol")
    set_quantity!(st, "CO2@", 0.01u"mol")
    A = Float64.(cs.SM.A)
    b = A * [ustrip(us"mol", x) for x in st.n]
    amounts(s) = [ustrip(us"mol", x) for x in s.n]

    ipopt = IPOPT_EXT._default_ipopt_solver()
    # Linear variables, element totals taken from the state.
    eq_lin = ChemistryLab.SciMLBase.solve(EquilibriumSolver(cs, DiluteSolutionModel(), ipopt), st)
    # Logarithmic variables, element totals given explicitly, started from the
    # linear answer: from amounts floored at 1e-16, the gradient in log n of every
    # trace is of order 1e-16 μ, and the starting point already passes Ipopt's
    # optimality test.
    eq_log = ChemistryLab.SciMLBase.solve(
        EquilibriumSolver(cs, DiluteSolutionModel(), ipopt; variable_space = Val(:log)), eq_lin; b = b,
    )

    des = DualEquilibriumSolver(cs, DiluteSolutionModel())
    ref = ChemistryLab.SciMLBase.solve(des, eq_lin; b = b)
    @test optimality_certificate(des, ref; b = b).optimal
    for eq in (eq_lin, eq_log)
        n = amounts(eq)
        # Polished by the dual Newton, Ipopt's answer is the certified one, traces
        # included, and conserves matter as it does.
        @test maximum(abs.(A * n .- b)) < 1.0e-10 * maximum(b)
        @test all(isapprox.(n, amounts(ref); rtol = 1.0e-8, atol = 1.0e-14))
    end

    # A composition carrying dual numbers is solved in real arithmetic and
    # differentiated at the answer: here with respect to the calcium of the budget.
    ica = findfirst(p -> symbol(p) == "Ca+2", cs.SM.primaries)
    e = zeros(length(b)); e[ica] = 1.0
    function composition(x, solve)
        seed = ChemicalState(cs, [ustrip(us"mol", v) + zero(x) for v in st.n] .* u"mol")
        return [ustrip(us"mol", v) for v in solve(seed, b .+ x .* e).n]
    end
    via_ipopt(seed, b) = ChemistryLab.SciMLBase.solve(EquilibriumSolver(cs, DiluteSolutionModel(), ipopt), seed; b)
    via_certificate(seed, b) = first(equilibrate_certified(seed; b))
    dn = ForwardDiff.derivative(x -> composition(x, via_ipopt), 0.0)
    dn_ref = ForwardDiff.derivative(x -> composition(x, via_certificate), 0.0)
    # The perturbation of the budget is carried exactly ...
    @test maximum(abs.(A * dn .- e)) < 1.0e-8
    # ... and the derivative, taken at Ipopt's polished answer, is that of the
    # certified one, species by species.
    @test all(isapprox.(dn, dn_ref; rtol = 1.0e-8, atol = 1.0e-12))
end

# The case of the external audit of 2026-10-02: Davies activities, whose neutral
# solutes carry a salting-out term with no partner in the ions' coefficients, so
# that they are not the gradient of a Gibbs energy. `n⋅μ(n)` then has the
# gradient `μ + Jᵀn`, not `μ`, and its minimum is not the equilibrium.
function _audit_davies_case()
    db = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-thermofun.json"); verbose = false))
    syms = ["H2O@", "H+", "OH-", "Na+", "Cl-", "CO2@", "HCO3-", "CO3-2", "Ca+2", "CaOH+", "Ca(CO3)@", "Ca(HCO3)+", "Cal"]
    cs = ChemicalSystem([db[s] for s in syms], ["H2O@", "H+", "Ca+2", "CO3-2", "Na+", "Cl-", "Zz"])
    st = ChemicalState(cs)
    set_quantity!(st, "H2O@", 1.0u"kg")
    set_quantity!(st, "Na+", 0.1u"mol")
    set_quantity!(st, "Cl-", 0.1u"mol")
    set_quantity!(st, "Cal", 0.05u"mol")
    set_quantity!(st, "CO2@", 0.01u"mol")
    return cs, st
end

@testset "every back end returns the equilibrium, not the minimum of n⋅μ(n)" begin
    cs, st = _audit_davies_case()
    model = DaviesActivityModel()
    amounts(s) = [ustrip(us"mol", x) for x in s.n]
    des = DualEquilibriumSolver(cs, model)
    b = des.A * amounts(st)
    ref, cert_ref = equilibrate_certified(st; model)
    @test cert_ref.optimal
    # Not the gradient of a Gibbs energy, and Gibbs–Duhem fails with it.
    J = ForwardDiff.jacobian(n -> des.lna(n, ChemistryLab._build_params(ref)), amounts(ref))
    @test ChemistryLab._gibbs_duhem_defect(J, amounts(ref))[1] > 1.0e-6

    ipopt = IPOPT_EXT._default_ipopt_solver()
    # Ipopt's own answer: matter conserved, the conditions of equilibrium not.
    raw = ChemistryLab._unpolished(() -> ChemistryLab.SciMLBase.solve(EquilibriumSolver(cs, model, ipopt), st))
    c_raw = optimality_certificate(des, raw; b = b)
    @test c_raw.balance < 1.0e-10
    @test c_raw.stationarity_abs > 1.0e-3
    @test !c_raw.optimal

    for (solver, space) in ((ipopt, Val(:linear)), (ipopt, Val(:log)), (OptimaOptimizer(), Val(:linear)), (OptimaOptimizer(), Val(:log)))
        cref = Ref{Any}()
        start = space === Val(:log) ? ref : st      # the log space refines a solved state
        eq = ChemistryLab.SciMLBase.solve(EquilibriumSolver(cs, model, solver; variable_space = space), start; b = b, certificate = cref)
        @test cref[].optimal
        @test optimality_certificate(des, eq; b = b).stationarity_abs < 1.0e-9
        @test all(isapprox.(amounts(eq), amounts(ref); rtol = 1.0e-8, atol = 1.0e-14))
    end
end

@testset "the derivative is that of the composition returned" begin
    cs, st = _audit_davies_case()
    model = DaviesActivityModel()
    ipopt = IPOPT_EXT._default_ipopt_solver()
    des = DualEquilibriumSolver(cs, model)
    A = des.A
    b = A * [ustrip(us"mol", x) for x in st.n]
    # Carbon dioxide added: the budget of its components moves by its column.
    e = A[:, findfirst(==("CO2@"), symbol.(cs.species))]
    function composition(x, solve)
        seed = ChemicalState(cs, [ustrip(us"mol", v) + zero(x) for v in st.n] .* u"mol")
        return [ustrip(us"mol", v) for v in solve(seed, b .+ x .* e).n]
    end
    via_ipopt(seed, bb) = ChemistryLab.SciMLBase.solve(EquilibriumSolver(cs, model, ipopt), seed; b = bb)
    via_certificate(seed, bb) = first(equilibrate_certified(seed; model, b = bb))
    n = composition(0.0, via_ipopt)
    dn = ForwardDiff.derivative(x -> composition(x, via_ipopt), 0.0)
    dn_ref = ForwardDiff.derivative(x -> composition(x, via_certificate), 0.0)
    # The balance carried exactly, `A ṅ = ḃ`.
    @test maximum(abs.(A * dn .- e)) < 1.0e-10
    # The linearized stationarity of the species present, `J ṅ + Aᵀẏ = 0` for
    # some ẏ: the tangent of the conditions the answer satisfies.
    p = ChemistryLab._build_params(ChemicalState(cs, n .* u"mol"))
    F = findall(>(ChemistryLab._CERTIFICATE_FLOOR), n)
    J = ForwardDiff.jacobian(nn -> des.lna(nn, p), n)
    r = J[F, :] * dn
    ẏ = A[:, F]' \ (-r)
    @test maximum(abs.(r .+ A[:, F]' * ẏ)) < 1.0e-8 * max(maximum(abs, r), 1.0)
    # And the derivative of the certified route, species by species.
    @test all(isapprox.(dn, dn_ref; rtol = 1.0e-8, atol = 1.0e-12))
end

@testset "without OptimaSolver, Ipopt refuses a model that breaks Gibbs–Duhem" begin
    cs, st = _audit_davies_case()
    ipopt = IPOPT_EXT._default_ipopt_solver()
    amounts(s) = [ustrip(us"mol", x) for x in s.n]
    ref_dil, _ = equilibrate_certified(st)
    ref_hkf, _ = equilibrate_certified(st; model = HKFActivityModel(; å = 4.0, Ḃ = 0.0, Kₙ = 0.0))
    saved = ChemistryLab._DUAL_AVAILABLE[]
    try
        ChemistryLab._DUAL_AVAILABLE[] = false
        @test_throws ArgumentError ChemistryLab.SciMLBase.solve(EquilibriumSolver(cs, DaviesActivityModel(), ipopt), st)
        err = try
            ChemistryLab.SciMLBase.solve(EquilibriumSolver(cs, HKFActivityModel(), ipopt), st)
        catch ex
            ex
        end
        @test err isa ArgumentError && occursin("Gibbs–Duhem", sprint(showerror, err))
        # The ideal model and a Debye–Hückel form with one ion size and no linear
        # term derive from one Gibbs energy: Ipopt alone minimizes it, and lands
        # on the equilibrium the dual Newton certifies, on every species above a
        # millimole to 2e-6 when this was written. The traces are left to the
        # interior point's absolute tolerance, whatever the model.
        for (model, ref) in ((DiluteSolutionModel(), ref_dil), (HKFActivityModel(; å = 4.0, Ḃ = 0.0, Kₙ = 0.0), ref_hkf))
            eq = ChemistryLab.SciMLBase.solve(EquilibriumSolver(cs, model, ipopt), st)
            major = amounts(ref) .> 1.0e-3
            @test all(isapprox.(amounts(eq)[major], amounts(ref)[major]; rtol = 1.0e-5))
        end
    finally
        ChemistryLab._DUAL_AVAILABLE[] = saved
    end
    # Under Davies, whose neutral solute breaks the relation, the minimum of
    # n⋅μ(n) is elsewhere: the same majors, unpolished, are hundreds of times
    # further off (7e-4 on the dissolved calcium when this was written).
    ref_dav, _ = equilibrate_certified(st; model = DaviesActivityModel())
    raw = ChemistryLab._unpolished(() -> ChemistryLab.SciMLBase.solve(EquilibriumSolver(cs, DaviesActivityModel(), ipopt), st))
    ica = findfirst(==("Ca+2"), symbol.(cs.species))
    @test abs(amounts(raw)[ica] / amounts(ref_dav)[ica] - 1) > 1.0e-4
end

@testset "Ipopt keeps a site budget with its host" begin
    # Hydrous ferric oxide carrying Dzombak and Morel's weak sites, with a budget
    # that follows the amount of oxide, in an acidic chloride solution: the
    # interior point has to be posed on the coupled matrix for the sites to leave
    # with the oxide it dissolves.
    psi = build_species(datapath("psinagra-12-07-thermofun.json"); verbose = false)
    bn = Dict(symbol(s) => s for s in psi)
    M = ustrip(us"kg/mol", bn["Fe(OH)3(am)"][:M])
    ν = literature_value("DzombakMorel1990", "weak_sites_per_mol_Fe")
    aq = speciation(
        psi, ["Fe(OH)3(am)", "Cl-"]; aggregate_state = [AS_AQUEOUS],
        exclude_species = split("H2@ O2@ Fe+2 FeOH+ FeO+ Cl2@ ClO- HClO@ ClO2- ClO3- ClO4- HClO2@"),
    )
    RT = ChemistryLab.R_GAS * 298.15
    prot = reference_oracle("phreeqc_protolysis")
    surface(sym, logK) = begin
        s = Species(sym; aggregate_state = AS_SURFACE, class = SC_SURFCOMPLEX)
        s[:ΔₐG⁰] = SymbolicFunc(-RT * log(10.0) * logK * u"J/mol")
        s
    end
    mem = [surface("XwOH", 0.0), surface("XwOH2+", prot.logK_protonation), surface("XwO-", prot.logK_deprotonation)]
    family = SiteFamily(
        "Xw", mem[1], mem[2:3];
        capacity = MassSiteDensity(ν / M),
        support = SurfaceSupport(
            "hydrous ferric oxide", "Fe(OH)3(am)", FixedSurfaceArea(1.0); coupling = SITES_FOLLOW_HOST,
        ),
    )
    cs = ChemicalSystem(
        AbstractSpecies[vcat(aq, mem)...],
        AbstractSpecies[bn["H2O@"], bn["H+"], bn["Fe+3"], bn["Cl-"], mem[1]];
        site_families = [family],
    )
    st = ChemicalState(cs)
    set_quantity!(st, "H2O@", moles_of_water() * u"mol")
    set_quantity!(st, "Fe(OH)3(am)", 1.0e-3u"mol")
    set_quantity!(st, "H+", 2.0e-3u"mol")
    set_quantity!(st, "Cl-", 2.0e-3u"mol")
    st = host_consistent_state(st)

    Ac = Float64.(conservation_matrix(cs))
    @test Ac != Float64.(cs.SM.A)
    n0 = [ustrip(us"mol", x) for x in st.n]
    eq = ChemistryLab.SciMLBase.solve(EquilibriumSolver(cs, DaviesActivityModel(), IPOPT_EXT._default_ipopt_solver()), st)
    n = [ustrip(us"mol", x) for x in eq.n]
    idx(s) = findfirst(==(s), symbol.(cs.species))
    sites = sum(n[idx(m)] for m in ("XwOH", "XwOH2+", "XwO-"))
    @test n[idx("Fe(OH)3(am)")] < 0.9 * n0[idx("Fe(OH)3(am)")]      # the acid dissolves the host
    @test sites / n[idx("Fe(OH)3(am)")] ≈ ν rtol = 1.0e-8
    @test maximum(abs.(Ac * n .- Ac * n0)) < 1.0e-10 * maximum(abs.(Ac * n0))
end
