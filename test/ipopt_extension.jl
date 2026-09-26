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
        # Matter is conserved to Ipopt's constraint tolerance, relative to the water.
        @test maximum(abs.(A * n .- b)) < 1.0e-10 * maximum(b)
        # Species by species, above a micromole, the interior point stops within
        # a fraction of a percent of the certified answer (2e-3 at most when this
        # was written); traces are not compared.
        major = amounts(ref) .> 1.0e-6
        @test all(abs.(n[major] ./ amounts(ref)[major] .- 1) .< 1.0e-2)
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
    # ... and the derivative, taken at Ipopt's answer, is that of the certified one
    # to within the distance between the two answers.
    icl = findfirst(==("Cal"), symbol.(cs.species))
    @test dn[icl] ≈ dn_ref[icl] rtol = 1.0e-2
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
