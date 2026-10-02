# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# The aqueous solutes recovered from their potentials through the ionic
# strength (`_aqueous_inverter`), the `invert` of OptimaSolver's aqueous phase.

using ChemistryLab
using DynamicQuantities
using ForwardDiff
using LinearAlgebra
using Test

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

    @testset "a model of more than the ionic strength is left to the sweeps" begin
        @test ChemistryLab._aqueous_inverter(DualEquilibriumSolver(cs, SITActivityModel())) === nothing
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
