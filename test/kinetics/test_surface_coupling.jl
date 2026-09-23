# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using ChemistryLab
using DynamicQuantities
using OrderedCollections
using Test

const RT_COUP = ChemistryLab.R_GAS * 298.15
include("../reference_species.jl")

_gc(v) = SymbolicFunc(v * u"J/mol")

# A solid releasing calcium at a constant rate, and an inert sorbent taking some
# of it up. The rate is constant on purpose: the kinetic half is then exactly
# integrable, so anything that disagrees is the **equilibrium map** — which is
# what a surface changes. Same reasoning as `test/coupling_reference.jl`.
function _sorbent_system(; logK_Ca = -3.0, n_sites = 2.0e-4)
    aq(sym, g, cl = SC_AQSOLUTE) = begin
        s = Species(sym; aggregate_state = AS_AQUEOUS, class = cl)
        s[:ΔₐG⁰] = _gc(g)
        s
    end
    sf(sym, g) = begin
        s = Species(sym; aggregate_state = AS_SURFACE, class = SC_SURFCOMPLEX)
        s[:ΔₐG⁰] = _gc(g)
        s
    end

    h2o, hp, oh, ca = reference_species(("H2O@", "H+", "OH-", "Ca+2"))
    # CEMDATA18 names this phase "Portlandite"; the tests below look it up by
    # that symbol, which is the database's, not a spelling chosen here.
    solid = reference_species("Portlandite"; db = :cemdata18)
    G_CA = ustrip(us"J/mol", ca[:ΔₐG⁰](T = 298.15u"K", P = 1.0e5u"Pa"; unit = true))

    free = sf("XsOH", 0.0)
    # XsOH + Ca+2 = XsOCa+ + H+, so ΔrG° = G(XsOCa+) − G(XsOH) − G(Ca+2) = −RT ln K
    bound = sf("XsOCa+", -RT_COUP * log(10.0^logK_Ca) + G_CA)

    support = SurfaceSupport("inert sorbent", nothing, FixedSurfaceArea(1.0))
    family = SiteFamily(
        "Xs", free, [bound]; capacity = TotalSiteAmount(n_sites * u"mol"), support,
    )

    species = [h2o, hp, oh, ca, solid, free, bound]
    cs = ChemicalSystem(species, [h2o, hp, ca, free]; site_families = [family])
    return cs, family, n_sites
end

@testsection "Surface complexation coupled to kinetics" begin

    @testset "the equilibrium partition keeps its site families" begin
        # `_equilibrium_subsystem` rebuilds a `ChemicalSystem` for the partition
        # the equilibrium is solved on, and the `ChemicalSystem` constructor takes
        # its declarations by keyword. Anything not passed is dropped silently —
        # which is exactly what happened to solid solutions until 0.8.2, and the
        # reason this test exists in the same shape for surfaces.
        cs, _, n_sites = _sorbent_system()
        nm = symbol.(cs.species)
        idx_eq = [i for (i, s) in enumerate(nm) if s != "Portlandite"]

        sub = ChemistryLab._equilibrium_subsystem(cs, idx_eq)
        @test sub.site_families !== nothing
        @test length(sub.site_families) == 1
        @test name(first(sub.site_families)) == "Xs"
        @test length(sub.site_groups) == 1
        @test symbol.(sub.species[sub.site_groups[1]]) == ["XsOH", "XsOCa+"]

        # And the site conservation row survives the rebuild: without it the
        # partition has a surface with no budget.
        row = findfirst(p -> symbol(p) == "XsOH", sub.SM.primaries)
        @test row !== nothing
        @test Int.(sub.SM.A[row, sub.site_groups[1]]) == [1, 1]
    end

    @testset "a family split between the partitions is refused, by name" begin
        # Sites equilibrate fast by construction in this release, so a family
        # with one member on each side is a declaration error. Dropping it would
        # leave orphans and a message about orphans; saying which family was
        # split is the useful error.
        cs, _, _ = _sorbent_system()
        nm = symbol.(cs.species)
        split = [i for (i, s) in enumerate(nm) if s != "XsOCa+"]
        err = try
            ChemistryLab._equilibrium_subsystem(cs, split)
            nothing
        catch e
            sprint(showerror, e)
        end
        @test err !== nothing
        @test occursin("Xs", err) && occursin("split", err)
    end

    @testset "no site is classified as kinetic by accident" begin
        # `:auto` puts every non-aqueous participant of a kinetic reaction in the
        # kinetic partition. A surface species is not aqueous, so without an
        # explicit exclusion it would go there — taking its family's conservation
        # row with it and leaving the equilibrium a surface with no sites.
        cs, _, _ = _sorbent_system()
        nm = symbol.(cs.species)
        idx(s) = findfirst(==(s), nm)
        spc = Dict(s => cs.species[idx(s)] for s in nm)

        # A reaction that deliberately involves a surface species alongside the
        # dissolving solid, which is the case the rule has to get right.
        rxn = Reaction(
            OrderedDict(spc["Portlandite"] => 1, spc["XsOH"] => 1),
            OrderedDict(spc["XsOCa+"] => 1, spc["OH-"] => 1, spc["H2O@"] => 1);
            symbol = "release", equal_sign = '→',
        )
        st = zeros(Float64, length(nm))
        st[idx("Portlandite")] = -1.0; st[idx("XsOH")] = -1.0
        st[idx("XsOCa+")] = 1.0; st[idx("OH-")] = 1.0; st[idx("H2O@")] = 1.0
        kr = KineticReaction(rxn, (T, P, t, n, lna, n0) -> 1.0e-9, idx("Portlandite"), st)

        K = ChemistryLab._reactivity_matrix([kr], cs, :auto)
        kinetic_rows = [i for i in 1:length(nm) if any(!iszero, K[i, :])]
        @test idx("Portlandite") in kinetic_rows          # the solid is kinetic
        @test !(idx("XsOH") in kinetic_rows)          # the sites are not
        @test !(idx("XsOCa+") in kinetic_rows)
        @test all(!(i in kinetic_rows) for i in cs.idx_surface)
    end

    @testset "a constant-rate release onto an inert sorbent" begin
        # The kinetic half is exactly integrable — a constant rate — so the
        # amount of calcium released by time t is k·t and nothing else. What is
        # under test is the equilibrium map the coupling calls at every accepted
        # step, now that it has a surface in it.
        logK_Ca, n_sites = -3.0, 2.0e-4
        cs, _, _ = _sorbent_system(; logK_Ca, n_sites)
        nm = symbol.(cs.species)
        idx(s) = findfirst(==(s), nm)
        spc = Dict(s => cs.species[idx(s)] for s in nm)

        k_rate = 2.0e-8                      # mol/s
        rxn = Reaction(
            OrderedDict(spc["Portlandite"] => 1),
            OrderedDict(spc["Ca+2"] => 1, spc["OH-"] => 2);
            symbol = "dissolution", equal_sign = '→',
        )
        st = zeros(Float64, length(nm))
        st[idx("Portlandite")] = -1.0; st[idx("Ca+2")] = 1.0; st[idx("OH-")] = 2.0
        kr = KineticReaction(rxn, (T, P, t, n, lna, n0) -> k_rate, idx("Portlandite"), st)

        n0 = Any[fill(0.0u"mol", length(nm))...]
        n0[idx("H2O@")] = moles_of_water() * u"mol"
        n0[idx("Portlandite")] = 1.0e-3u"mol"
        n0[idx("XsOH")] = n_sites * u"mol"
        state = ChemicalState(cs, n0)

        function run(reltol)
            kp = KineticsProblem(
                cs, [kr], state, (0.0u"s", 3600.0u"s");
                equilibrium_solver = EquilibriumSolver(
                    cs, DiluteSolutionModel(), OptimaOptimizer()
                ),
            )
            ks = KineticsSolver(;
                ode_solver = Rodas5P(), reltol = reltol, abstol = 1.0e-16,
                saveat = [0.0, 900.0, 1800.0, 2700.0, 3600.0],
            )
            return kp, integrate(kp, ks)
        end

        kp, sol = run(1.0e-9)
        @test length(sol.t) == 5

        states = speciated_states(sol, kp)
        bound = Float64[ustrip(us"mol", s.n[idx("XsOCa+")]) for s in states]
        free = Float64[ustrip(us"mol", s.n[idx("XsOH")]) for s in states]

        # The site budget is a conservation row, and it holds at every reported
        # instant — not merely at the end.
        for (f, b) in zip(free, bound)
            @test f + b ≈ n_sites rtol = 1.0e-8
        end

        # The sorbent fills as the calcium arrives, and never empties.
        @test issorted(bound)
        @test bound[end] > bound[1]

        # The released calcium is exactly k·t, and it is somewhere: in solution,
        # on the sorbent, or paired. Nothing is lost between the two partitions.
        released = k_rate * 3600.0
        last = states[end]
        solid_left = ustrip(us"mol", last.n[idx("Portlandite")])
        @test (1.0e-3 - solid_left) ≈ released rtol = 1.0e-6
        ca_accounted = ustrip(us"mol", last.n[idx("Ca+2")]) + bound[end]
        @test ca_accounted ≈ released rtol = 1.0e-5

        # And the law of mass action holds on the integrated state, which is the
        # independent algebraic check: the coupling did not merely conserve
        # matter, it reached the equilibrium it was supposed to.
        n = Float64[ustrip(us"mol", x) for x in last.n]
        des = DualEquilibriumSolver(cs, DiluteSolutionModel())
        p = ChemistryLab._build_params(last)
        lna = des.lna(n, p)
        lhs = lna[idx("XsOCa+")] + lna[idx("H+")]
        rhs = lna[idx("XsOH")] + lna[idx("Ca+2")] + logK_Ca * log(10)
        @test lhs ≈ rhs atol = 1.0e-6

        # Tightening the integrator does not move the answer: the trajectory is
        # converged in the time step, not merely reproducible at one setting.
        _, sol2 = run(1.0e-11)
        states2 = speciated_states(sol2, kp)
        bound2 = Float64[ustrip(us"mol", s.n[idx("XsOCa+")]) for s in states2]
        @test maximum(abs.(bound2 ./ bound .- 1)) < 1.0e-5
    end

end
