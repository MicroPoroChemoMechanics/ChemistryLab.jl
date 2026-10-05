# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# A kinetic run whose partition holds a phase under the compound energy
# formalism, the C-S-H of CASH+ (Kulik et al. 2022). The activities of its
# members depend on their standard Gibbs energies, which the parameters of a
# run did not carry: the first evaluation of the activities stopped on "none
# were passed", before the run had started.

@testsection "a kinetic run with a compound-energy C-S-H" begin
    subs = build_species(datapath("cemdata18-cashplus.json"); verbose = false)
    byname = Dict(symbol(s) => s for s in subs)
    cash = only(filter(p -> ChemistryLab.name(p) == "CASH+", build_solid_solutions(datapath("solid_solutions.toml"), byname)))
    members = ["TSvh", "TSCh", "Tvvh", "TCvh", "TvCh", "TCCh"]
    species = speciation(subs, vcat(["Portlandite", "Amor-Sl"], members); aggregate_state = [AS_AQUEOUS], exclude_species = ["Ca(OH)2@"])
    cs = ChemicalSystem(species, CEMDATA_PRIMARIES; solid_solutions = [cash])
    model = cemdata18_activity_model(:KOH)
    Mw = ustrip(us"g/mol", byname["H2O@"][:M])
    function state(T = 298.15)
        st = ChemicalState(cs; T = T * u"K")
        set_quantity!(st, "H2O@", (1000 / Mw)u"mol")
        set_quantity!(st, "Amor-Sl", 0.05u"mol")
        set_quantity!(st, "Portlandite", 0.06u"mol")
        return st
    end
    # The silica dissolves at a rate of its own amount, into the partition.
    silica = Reaction(OrderedDict(cs["Amor-Sl"] => 1.0), OrderedDict(cs["SiO2@"] => 1.0); symbol = "silica")
    silica[:rate] = KineticFunc((T, P, t, n, lna, n0) -> 1.0e-6 * max(n["Amor-Sl"], zero(eltype(n.data))), NamedTuple(), u"mol/s")
    problem(st) = KineticsProblem(
        cs, [silica], st, (0.0, 1.0e6);
        activity_model = model, equilibrium_solver = EquilibriumSolver(cs, model, OptimaOptimizer()),
    )

    # The activities the rate laws are given are those the equilibrium solver
    # works with, at the composition of an equilibrium holding the gel, at the
    # temperature of the run and at another.
    for T in (298.15, 313.15)
        eq, cert = equilibrate_certified(state(T); model)
        @test cert.optimal
        p = build_kinetics_params(problem(eq))
        n = Float64[ustrip(us"mol", x) for x in eq.n]
        @test p.lna_fn(n, ChemistryLab._lna_params(p, p.T)) ≈ activity_model(cs, model)(n, ChemistryLab._build_params(eq)) rtol = 1.0e-12
    end
    # A cell that warms reads them at its own temperature.
    p = build_kinetics_params(problem(state()))
    q = merge(p, (has_T = true,))
    @test ChemistryLab._lna_params(q, 313.15).ΔₐG⁰overRT ≈ ChemistryLab._build_params(state(313.15)).ΔₐG⁰overRT rtol = 1.0e-12
    # A system with no such phase is handed its parameters as they were.
    plain = ChemicalSystem(speciation(subs, ["Portlandite", "Amor-Sl"]; aggregate_state = [AS_AQUEOUS], exclude_species = ["Ca(OH)2@"]), CEMDATA_PRIMARIES)
    @test !ChemistryLab._reads_standard_g(plain) && ChemistryLab._reads_standard_g(cs)

    # And a run goes through, the gel forming as a silica glass dissolves into
    # the partition. The glass has no standard Gibbs energy, as a glass has
    # none: only the members of the gel are read.
    glass = glass_species(Dict("SiO2" => 1.0); symbol = "SG", M = Species("SiO2")[:M])
    cs_g = ChemicalSystem(vcat(species, [glass]), CEMDATA_PRIMARIES; solid_solutions = [cash])
    dissolution = Reaction(OrderedDict(cs_g["SG"] => 1.0), OrderedDict(cs_g["SiO2@"] => 1.0); symbol = "glass")
    dissolution[:rate] = KineticFunc((T, P, t, n, lna, n0) -> 1.0e-6 * max(n["SG"], zero(eltype(n.data))), NamedTuple(), u"mol/s")
    st = ChemicalState(cs_g)
    set_quantity!(st, "H2O@", (1000 / Mw)u"mol")
    set_quantity!(st, "Portlandite", 0.06u"mol")
    set_quantity!(st, "SG", 0.05u"mol")
    kp = KineticsProblem(
        cs_g, [dissolution], st, (0.0, 1.0e6);
        activity_model = model, equilibrium_solver = EquilibriumSolver(cs_g, model, OptimaOptimizer()),
    )
    sol = integrate(kp, KineticsSolver(; ode_solver = Rodas5P(), reltol = 1.0e-6, abstol = 1.0e-10))
    p = sol.prob.p
    @test SciMLBase.successful_retcode(sol)
    @test !p.rates_read_speciation[]
    @test sol.u[end][p.n_be + 1] ≈ 0.05 * exp(-1.0) rtol = 1.0e-4
    gel = [p.species_index[m] for m in members]
    @test sum(p.n_full[gel]) > 0.01
end
