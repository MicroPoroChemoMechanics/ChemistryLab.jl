# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# A quaternary cement on the coupled kinetic path: 50 % OPC with blast-furnace
# slag, siliceous fly ash and limestone, from Schöler et al. (2015)
# (`scripts/scholer2015_kinetics.jl`). The materials come from their templates,
# the two glasses become kinetic species through `glass_species` and
# `with_species`, the limestone is at equilibrium. The comparison with the
# thermogravimetry of the article is the page `examples/quaternary_kinetics.md`;
# this file holds a week of one mix to what it must satisfy whatever the data.

isdefined(@__MODULE__, :s15_run) || include(joinpath(pkgdir(ChemistryLab), "scripts", "scholer2015_kinetics.jl"))

@testsection "a quaternary cement on the coupled kinetic path" begin
    day = 86400.0
    setup = s15_setup()
    cs, mats = setup.cs, setup.mats
    mix = s15_mix("30-10-10")

    # The anhydrite brings the SO₃ of the whole to 3 %, counted on the analyses
    # of Table 1 and on the formula of anhydrite.
    so3(m) = get(literature_oxides("Scholer2015", "oxides", m), "SO3", 0.0)
    a = s15_anhydrite(mix)
    f = oxide_content(S15_DB["Anh"], ["SO3"])["SO3"]
    s = 100 * (mix.opc * so3("OPC") + mix.bfs * so3("BFS") + mix.fa * so3("FA") + mix.ls * so3("LS"))
    @test (s + f * a) / (100 + a) ≈ 0.03 rtol = 1.0e-12
    @test a > 0

    # The polymorphs of the Rietveld analysis are one constituent each.
    names = [c.name for c in mats.opc.constituents]
    @test count(==("C2S"), names) == 1 && count(==("C3A"), names) == 1
    @test !any(n -> occursin("C2S", n) && n != "C2S", names)
    # The glasses are species, named by their symbols, the crystals inert.
    @test any(c -> c.name == "BFS" && c isa MineralConstituent, mats.bfs.constituents)
    @test all(c -> c.name == "BFS" || effective_extent(mats.bfs, c, 0.0) == 0, mats.bfs.constituents)

    # The calibration reproduces the degrees the article assumes after a year.
    n = WALLER_PARAMS_FLY_ASH.n
    α(τ, t) = 1 / (1 + (ustrip(us"d", τ) / t)^n)
    year = ustrip(us"d", literature_value("Scholer2015", "assumed_reaction_age"))
    @test α(s15_glass_time(:bfs), year) ≈ literature_value("Scholer2015", "assumed_reaction_slag_glass") / 100 rtol = 1.0e-12
    @test α(s15_glass_time(:fa), year) ≈ literature_value("Scholer2015", "assumed_reaction_fly_ash_glass") / 100 rtol = 1.0e-12

    # A mix without fly ash leaves its glass without a rate law: at equilibrium,
    # where a glass, which has no standard Gibbs energy, cannot be. The system of
    # that mix is built without it; given the full system, the problem is
    # refused by name rather than run with every re-speciation failing.
    nofa = s15_mix("30-0-20")
    err = try
        KineticsProblem(
            s15_recipe(mats, nofa), cs, s15_rates(nofa), (0.0, 86400.0);
            activity_model = cemdata18_activity_model(:KOH),
            equilibrium_solver = EquilibriumSolver(cs, cemdata18_activity_model(:KOH), OptimaOptimizer()),
        )
        nothing
    catch e
        e
    end
    @test err isa ArgumentError && occursin("FA would be at equilibrium", sprint(showerror, err))
    @test !any(s -> symbol(s) == "FA", s15_system(cs, nofa).species)
    @test s15_system(cs, mix) === cs

    r = s15_run(cs, mats, "30-10-10"; days = 7)
    @test SciMLBase.successful_retcode(r.sol)
    # Each glass follows its law: at the reference temperature, the integrated
    # degree is the closed form.
    degree(sym, t) = (
        i = findfirst(s -> symbol(s) == sym, r.kp.system.species);
        k = findfirst(==(i), r.kp.idx_kinetic);
        1 - r.sol(t)[size(r.kp.Ae, 1) + k] / ustrip(us"mol", r.kp.initial_state.n[i])
    )
    for d in (1.0, 7.0)
        @test degree("BFS", d * day) ≈ α(s15_glass_time(:bfs), d) rtol = 1.0e-3
        @test degree("FA", d * day) ≈ α(s15_glass_time(:fa), d) rtol = 1.0e-3
    end
    @test degree("C3S", 7day) > degree("C2S", 7day) > 0

    # The replayed states hold the elements of the mix and are certified.
    A = Float64.(cs.SM.A)
    b0 = A * ustrip.(us"mol", r.kp.initial_state.n)
    states = speciated_states(r.sol, r.kp; times = [1day, 7day])
    for st in states
        @test A * ustrip.(us"mol", st.n) ≈ b0 rtol = 1.0e-9 atol = 1.0e-12
    end
    # The limestone is at equilibrium, and in this system the aluminum the
    # clinker and the glasses release goes mostly to ettringite and to the
    # siliceous hydrogarnet of Cemdata18 rather than to carboaluminates: the
    # calcite is left, and portlandite forms.
    amount(st, s) = ustrip(us"mol", st.n[findfirst(x -> symbol(x) == s, cs.species)])
    @test amount(states[2], "Cal") > 0
    @test amount(states[2], "Portlandite") > amount(states[1], "Portlandite") > 0

    # What the thermobalance weighs: positive, the bound water growing.
    tga = s15_tga(r, [1, 7])
    @test 0 < tga[1].bound_water < tga[2].bound_water
    @test tga[2].portlandite > 0
end
