# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# The aqueous solutes recovered from their potentials through the ionic
# strength (`_aqueous_inverter`), the `invert` of OptimaSolver's aqueous phase.

using ChemistryLab
using DynamicQuantities
using ForwardDiff
using LinearAlgebra
using Test

include("reference_species.jl")

@testsection "the aqueous solutes, recovered through the ionic strength" begin
    CEM = Dict(symbol(x) => x for x in build_species(datapath("cemdata18-thermofun.json")))
    spc = speciation(
        collect(values(CEM)), split("Portlandite Gp H2O@");
        aggregate_state = [AS_AQUEOUS], exclude_species = ["S-2"],
    )
    cs = ChemicalSystem(spc, CEMDATA_PRIMARIES)
    st = ChemicalState(cs)
    set_quantity!(st, "H2O@", 1.0u"kg")
    # A concentrated solution, not an equilibrium: the inversion only has to
    # undo the activity model at whatever composition it is given.
    for (s, v) in ("Ca+2" => 0.4, "SO4-2" => 0.15, "OH-" => 0.5, "CaOH+" => 0.05, "Ca(SO4)@" => 0.01)
        set_quantity!(st, s, v * u"mol")
    end
    n = [ustrip(us"mol", x) for x in st.n]
    p = ChemistryLab._build_params(st; ϵ = 1.0e-16)
    tj = TruesdellJonesActivityModel(; parameters = Dict("Ca+2" => (5.0, 0.165), "SO4-2" => (5.0, -0.04)))
    limiting = HKFActivityModel(å = 0.0, Ḃ = 0.097637, Kₙ = 0.0)
    # The same solutes at a tenth of the amounts, the solvent unchanged.
    dilute = [i in cs.idx_solvent ? n[i] : n[i] / 10 for i in eachindex(n)]

    # Whether the inversion gives back `x` from `c = h(x)`, from a start away
    # from the answer.
    function undone(model, x)
        des = DualEquilibriumSolver(cs, model)
        aq = des.idx_aq
        c = des.lna(x, p)[aq]
        invert = ChemistryLab._aqueous_inverter(des)
        @test invert !== nothing
        ref = x[aq[des.j_solvent]]
        w0 = log.(max.(x[aq], 1.0e-300)) .+ 0.7
        w = invert(c, ref, w0, Float64[], p)
        @test w[des.j_solvent] == w0[des.j_solvent]
        return des, aq, c, w, [t for (t, i) in enumerate(aq) if t != des.j_solvent && x[i] > 1.0e-20]
    end

    # A model whose ionic-strength equation has one root gives the composition
    # back wherever it is.
    models = (HKFActivityModel(), DaviesActivityModel(), tj, DiluteSolutionModel())
    @testset "it undoes the activity model, for $(nameof(typeof(model)))" for model in models
        _, aq, _, w, live = undone(model, n)
        for t in live
            @test w[t] ≈ log(n[aq[t]]) atol = 1.0e-10
        end
    end

    @testset "the limiting law: the dilute branch, which is a solution too" begin
        _, aq, _, w, live = undone(limiting, dilute)
        for t in live
            @test w[t] ≈ log(dilute[aq[t]]) atol = 1.0e-10
        end
        # At 1.4 mol/kg the composition sits on a branch of the limiting law
        # above the first root, past the dip where two roots meet. The inversion
        # returns that first root: another composition, at a lower ionic
        # strength, whose activities meet the same potentials exactly.
        des, aq, c, w, live = undone(limiting, n)
        x = copy(n)
        for (t, i) in enumerate(aq)
            t == des.j_solvent || (x[i] = exp(w[t]))
        end
        @test maximum(t -> abs(des.lna(x, p)[aq[t]] - c[t]), live) < 1.0e-10
        Iof(v) = 0.5 * sum(charge(cs.species[i])^2 * v[i] for i in cs.idx_solutes) / v[only(cs.idx_solvent)] / 0.018015
        @test Iof(x) < Iof(n)
    end

    @testset "no ionic strength solves the limiting law past its range" begin
        aq = DualEquilibriumSolver(cs, DiluteSolutionModel()).idx_aq
        ions = [t for (t, i) in enumerate(aq) if !iszero(charge(cs.species[i]))]
        raise(c) = [t in ions ? c[t] + 10 : c[t] for t in eachindex(c)]
        lim = DualEquilibriumSolver(cs, HKFActivityModel(å = 0.0, Ḃ = 0.0, Kₙ = 0.0))
        ref = n[aq[lim.j_solvent]]
        w0 = log.(max.(n[aq], 1.0e-300))
        c = raise(lim.lna(n, p)[aq])
        @test ChemistryLab._aqueous_inverter(lim)(c, ref, w0, Float64[], p) === nothing
        # The same potentials under an ion size per ion have a solution.
        per = DualEquilibriumSolver(cs, HKFActivityModel())
        @test ChemistryLab._aqueous_inverter(per)(c, ref, w0, Float64[], p) !== nothing
    end

    @testset "its derivative inverts the activity model's" begin
        des = DualEquilibriumSolver(cs, HKFActivityModel())
        aq = des.idx_aq
        sol = [t for (t, i) in enumerate(aq) if t != des.j_solvent && n[i] > 1.0e-20]
        c = des.lna(n, p)[aq]
        ref = n[aq[des.j_solvent]]
        w0 = log.(max.(n[aq], 1.0e-300))
        invert = ChemistryLab._aqueous_inverter(des)
        # ∂w/∂c, by differentiating the inversion; ∂h/∂w, the activity model's.
        dwdc = ForwardDiff.jacobian(cc -> invert(cc, ref, w0, Float64[], p)[sol], c)[:, sol]
        hw(wv) = begin
            x = Vector{eltype(wv)}(n)
            for (k, t) in enumerate(sol)
                x[aq[t]] = exp(wv[k])
            end
            des.lna(x, p)[aq[sol]]
        end
        dhdw = ForwardDiff.jacobian(hw, log.(n[aq[sol]]))
        @test dhdw * dwdc ≈ I(length(sol)) atol = 1.0e-9
    end

    @testset "a narrow dip below zero is found, a dip above it is not a root" begin
        F(s) = (s - 0.3)^2 - 1.0e-4
        dF(s) = 2 * (s - 0.3)
        r = ChemistryLab._ionic_strength_root(F, dF, -5.0)
        @test r ≈ 0.29 atol = 1.0e-12
        G(s) = (s - 0.3)^2 + 1.0e-4
        @test ChemistryLab._ionic_strength_root(G, dF, -5.0; smax = 3.0) === nothing
    end

    @testset "the dilute branch ends at a dip above zero, and at the ceiling" begin
        # A dip that stays above zero, then a rise, and a root far above: the
        # shape of the limiting law with `Ḃ` past its range, where the third root
        # lay at 6500 mol/kg for a paste loaded with sodium chloride. That root is
        # not on the dilute branch, and a solve started from it never recovered.
        H(s) = (s - 0.3)^2 + 0.01 - 0.05 * (s - 0.3)^4
        dH(s) = 2 * (s - 0.3) - 0.2 * (s - 0.3)^3
        @test H(4.8) < 0 < H(0.3) < H(-1.0)       # the far root lies below the ceiling
        @test ChemistryLab._ionic_strength_root(H, dH, -1.0) === nothing
        # A start already past the root is searched downwards, and what it finds
        # has to lie below the ceiling too: potentials that hold more than 1e4
        # mol/kg at zero ionic strength hold no solution.
        L(s) = 3.0 - s
        dL(s) = -1.0
        @test ChemistryLab._ionic_strength_root(L, dL, 5.0) ≈ 3.0 atol = 1.0e-12
        far(s) = 20.0 - s
        @test ChemistryLab._ionic_strength_root(far, dL, 25.0) === nothing
    end

    @testset "a model of more than the ionic strength has an inversion of its own" begin
        # SIT and Pitzer by Newton's method (the testsection below); a model
        # neither inversion covers is left to the sweeps.
        @test ChemistryLab._aqueous_inverter(DualEquilibriumSolver(cs, SITActivityModel())) !== nothing
        @test ChemistryLab._newton_predictor_form(DiluteSolutionModel(), cs, Int[]) === nothing
    end

    @testset "neutral solutes alone: no ionic strength to solve for" begin
        # `ln n = c − ln γ(0) + ln(n_w M_w)`, whatever the start.
        form = (; kind = :ionic, z = [0, 0, 0], M_w = 0.018015, log10γ = (t, z, I, sqrtI, A, B) -> 0.1 * t, AB = p -> (0.51, 0.33))
        c = [0.0, -2.0, -5.0]
        w = ChemistryLab._invert_aqueous(form, c, 55.5, [1.0, 7.0, 7.0], (; ϵ = 1.0e-16), 1)
        @test w[1] == 1.0
        @test w[2:3] ≈ [c[t] - log(10.0) * 0.1 * t + log(55.5 * 0.018015) for t in 2:3]
    end
end

@testsection "SIT and Pitzer: the solutes recovered by Newton's method" begin
    slop = Dict(symbol(s) => s for s in build_species(datapath("slop98-inorganic-thermofun.json"); verbose = false))
    # The three ε the PHREEQC comparison of `test/sit.jl` carries, and the
    # Reardon set of the Pitzer page.
    fixture = reference_oracle("phreeqc_sit")
    sit = SITActivityModel(; parameters = SITParameters([(e.a, e.b) => e.value for e in fixture.epsilon]))
    pitzer = PitzerActivityModel(; parameters = build_pitzer_parameters(datapath("pitzer-reardon1990.toml")))
    # Sodium chloride with the dissociation of water under SIT; without H+ and
    # OH- under Pitzer, whose set has no pair for them, as on the Pitzer page.
    sit_cs = ChemicalSystem([slop[s] for s in split("H2O@ H+ OH- Na+ Cl-")], ["H2O@", "H+", "Na+", "Cl-", "Zz"])
    pz_cs = ChemicalSystem([slop[s] for s in split("H2O@ Na+ Cl- Hl")], ["H2O@", "Na+", "Cl-"])
    function brine(cs, m)
        st = ChemicalState(cs)
        set_quantity!(st, "H2O@", 1.0u"kg")
        set_quantity!(st, "Na+", m * u"mol")
        set_quantity!(st, "Cl-", m * u"mol")
        if haskey(cs.dict_species, "H+")
            set_quantity!(st, "H+", 1.0e-7u"mol")
            set_quantity!(st, "OH-", 1.0e-7u"mol")
        end
        return st
    end
    amounts(st) = [ustrip(us"mol", v) for v in st.n]

    @testset "it undoes the model, for $(nameof(typeof(model))) at $m mol/kg" for (model, cs, m) in
        ((sit, sit_cs, 3.0), (pitzer, pz_cs, 6.0))
        des = DualEquilibriumSolver(cs, model)
        invert = ChemistryLab._aqueous_inverter(des)
        @test invert !== nothing
        st = brine(cs, m)
        x = amounts(st)
        p = ChemistryLab._build_params(st; ϵ = 1.0e-16)
        aq = des.idx_aq
        c = des.lna(x, p)[aq]
        w = invert(c, x[aq[des.j_solvent]], log.(x[aq]) .+ 0.7, Float64[], p)
        for (t, i) in enumerate(aq)
            t == des.j_solvent && continue
            @test w[t] ≈ log(x[i]) atol = 1.0e-10
        end
    end

    @testset "its derivative inverts the model's" begin
        des = DualEquilibriumSolver(sit_cs, sit)
        aq = des.idx_aq
        st = brine(sit_cs, 3.0)
        x = amounts(st)
        p = ChemistryLab._build_params(st; ϵ = 1.0e-16)
        sol = [t for t in eachindex(aq) if t != des.j_solvent]
        c = des.lna(x, p)[aq]
        ref = x[aq[des.j_solvent]]
        w0 = log.(x[aq])
        invert = ChemistryLab._aqueous_inverter(des)
        dwdc = ForwardDiff.jacobian(cc -> invert(cc, ref, w0, Float64[], p)[sol], c)[:, sol]
        hw(wv) = begin
            xx = Vector{eltype(wv)}(x)
            for (k, t) in enumerate(sol)
                xx[aq[t]] = exp(wv[k])
            end
            des.lna(xx, p)[aq[sol]]
        end
        dhdw = ForwardDiff.jacobian(hw, log.(x[aq[sol]]))
        @test dhdw * dwdc ≈ I(length(sol)) atol = 1.0e-9
    end

    @testset "a solute below the reach of a Newton step is placed" begin
        # Hydroxide at 1e-28 mol, near the 1e-30 floor of the activities, where
        # its activity barely moves with it: it is placed from its activity
        # coefficient at the composition found, and its derivative follows.
        des = DualEquilibriumSolver(sit_cs, sit)
        aq = des.idx_aq
        st = brine(sit_cs, 3.0)
        set_quantity!(st, "OH-", 1.0e-28u"mol")
        x = amounts(st)
        p = ChemistryLab._build_params(st; ϵ = 1.0e-16)
        c = des.lna(x, p)[aq]
        ref = x[aq[des.j_solvent]]
        sol = [t for t in eachindex(aq) if t != des.j_solvent]
        invert = ChemistryLab._aqueous_inverter(des)
        w = invert(c, ref, log.(x[aq]) .+ 0.7, Float64[], p)
        for t in sol
            @test w[t] ≈ log(x[aq[t]]) atol = 1.0e-8
        end
        dwdc = ForwardDiff.jacobian(cc -> invert(cc, ref, log.(x[aq]), Float64[], p)[sol], c)[:, sol]
        k = findfirst(t -> symbol(sit_cs.species[aq[t]]) == "OH-", sol)
        @test dwdc[k, k] ≈ 1 atol = 1.0e-6
        @test all(j == k || abs(dwdc[j, k]) < 1.0e-12 for j in eachindex(sol))

        # A tolerance below the rounding of the potentials ends the halving of
        # the step, and the inversion never answers with a wrong composition.
        pred = ChemistryLab._newton_predictor_form(sit, sit_cs, aq)
        out = ChemistryLab._invert_aqueous_newton(
            des.lna, pred, length(sit_cs.species), aq, des.j_solvent, c, ref,
            log.(x[aq]) .+ 0.7, p; tol = 1.0e-300,
        )
        @test out === nothing || all(abs(out[t] - log(x[aq[t]])) < 1.0e-8 for t in sol)

        # A solute started as a trace that is none: placed above the reach, it
        # rejoins the iteration. The start holds hydroxide at e⁻¹⁰⁰ and the
        # Debye–Hückel part is made to find no root, so that start is the one.
        st2 = brine(sit_cs, 3.0)
        x2 = amounts(st2)
        c2 = des.lna(x2, p)[aq]
        w2 = log.(x2[aq]) .+ 0.7
        j = findfirst(t -> symbol(sit_cs.species[aq[t]]) == "OH-", eachindex(aq))
        w2[j] = -100.0
        rootless = merge(pred, (; log10γ = (t, z, I, sqrtI, A, B) -> -1.0e6 * z^2))
        out2 = ChemistryLab._invert_aqueous_newton(
            des.lna, rootless, length(sit_cs.species), aq, des.j_solvent, c2,
            x2[aq[des.j_solvent]], w2, p,
        )
        @test out2 !== nothing
        for t in sol
            @test out2[t] ≈ log(x2[aq[t]]) atol = 1.0e-8
        end
    end

    @testset "a brine under SIT, certified" begin
        eq, cert = equilibrate_certified(brine(sit_cs, 3.0); model = sit)
        @test cert.optimal
        @test cert.balance_relative < 1.0e-10
    end

    @testset "halite in water under Pitzer, against its measured solubility" begin
        st = ChemicalState(pz_cs)
        set_quantity!(st, "H2O@", 1.0u"kg")
        set_quantity!(st, "Hl", 8.0u"mol")
        eq, cert = equilibrate_certified(st; model = pitzer)
        @test cert.optimal
        n = amounts(eq)
        idx = Dict(symbol(s) => i for (i, s) in enumerate(pz_cs.species))
        kg = n[idx["H2O@"]] * ustrip(us"kg/mol", pz_cs.species[idx["H2O@"]][:M])
        m = n[idx["Na+"]] / kg
        measured = ustrip(us"mol/kg", literature_value("HamerWu1972", "nacl_saturated_molality"))
        @info "halite under Pitzer: the saturated molality" computed = m measured
        @test n[idx["Hl"]] > 0
        # 6.1605 mol/kg against the 6.144 Hamer & Wu measured: 0.27 %, from the
        # halite of slop98 and the Na–Cl parameters of the Reardon set, neither
        # fitted here. Pinned at its own value, and against the measurement.
        @test m ≈ 6.1605 rtol = 1.0e-4
        @test m ≈ measured rtol = 5.0e-3
    end
end
