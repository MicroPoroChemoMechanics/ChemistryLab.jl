# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using ChemistryLab
using DynamicQuantities
using ForwardDiff
using SciMLBase
using Test

const RT25 = 8.31446261815324 * 298.15

# An amphoteric surface, which is the physics of every oxide and also the one
# shape in which competition has a **closed form under a single constraint**:
# both occupied states are driven by the proton activity, so
#
#     β₁ = K₁ a_H        (XsOH + H⁺ ⇌ XsOH₂⁺)
#     β₂ = K₂ / a_H      (XsOH ⇌ XsO⁻ + H⁺)
#
# and ideal site mixing gives the competitive Langmuir form with its shared
# denominator, without an isotherm being written anywhere:
#
#     n_j / N = β_j / (1 + β₁ + β₂)
#
# The standard energies below are chosen to *be* those constants, so the oracle
# is derived rather than fitted.
_g0(value) = SymbolicFunc(value * u"J/mol")

function _amphoteric_system(; logK1 = 7.0, logK2 = -9.0, n_sites = 1.0e-3)
    h2o = Species("H2O@"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
    h2o[:ΔₐG⁰] = _g0(-237181.0)
    h2o[:M] = 0.018015u"kg/mol"
    h2o[:V⁰] = SymbolicFunc(1.807e-5u"m^3/mol")

    hp = Species("H+"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    hp[:ΔₐG⁰] = _g0(0.0)
    oh = Species("OH-"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
    oh[:ΔₐG⁰] = _g0(-157297.0)

    free = Species("XsOH"; aggregate_state = AS_SURFACE, class = SC_SURFCOMPLEX)
    free[:ΔₐG⁰] = _g0(0.0)
    prot = Species("XsOH2+"; aggregate_state = AS_SURFACE, class = SC_SURFCOMPLEX)
    prot[:ΔₐG⁰] = _g0(-RT25 * log(10.0^logK1))
    depr = Species("XsO-"; aggregate_state = AS_SURFACE, class = SC_SURFCOMPLEX)
    depr[:ΔₐG⁰] = _g0(-RT25 * log(10.0^logK2))

    support = SurfaceSupport("oxide", nothing, FixedSurfaceArea(1.0))
    family = SiteFamily(
        "Xs", free, [prot, depr];
        capacity = TotalSiteAmount(n_sites), support,
    )
    species = [h2o, hp, oh, free, prot, depr]
    cs = ChemicalSystem(species, [h2o, hp, free]; site_families = [family])

    n = Any[fill(1.0e-12u"mol", length(cs.species))...]
    n[1] = 55.5u"mol"
    n[4] = n_sites * u"mol"          # all sites free to start with
    return cs, ChemicalState(cs, n), family
end

# The closed form the plan derives: a shared denominator is the signature of a
# finite capacity, and it is what makes this a test of competition rather than
# of two isotherms side by side.
_langmuir(β) = (1 ./ (1 + sum(β)), β ./ (1 + sum(β)))

@testsection "Surface complexation" begin

    @testset "a site activity is a site fraction, in every activity model" begin
        cs, st, _ = _amphoteric_system()
        n = Float64[ustrip(us"mol", x) for x in st.n]
        # A composition with all three states populated, so the fractions are
        # not degenerate.
        n[4], n[5], n[6] = 5.0e-4, 3.0e-4, 2.0e-4
        p = ChemistryLab._build_params(st)
        N = n[4] + n[5] + n[6]

        for model in (
                DiluteSolutionModel(), HKFActivityModel(), DaviesActivityModel(),
            )
            lna = activity_model(cs, model)(n, p)
            for i in 4:6
                @test lna[i] ≈ log(n[i] / N) rtol = 1.0e-9
            end
        end
    end

    @testset "the mixing differentiates" begin
        cs, st, _ = _amphoteric_system()
        n0 = Float64[ustrip(us"mol", x) for x in st.n]
        n0[4], n0[5], n0[6] = 5.0e-4, 3.0e-4, 2.0e-4
        p = ChemistryLab._build_params(st)
        lna = activity_model(cs, DiluteSolutionModel())

        J = ForwardDiff.jacobian(nn -> lna(nn, p), n0)
        @test all(isfinite, J[4:6, 4:6])
        # d(ln x_i)/dn_i = 1/n_i − 1/N, exactly, and the cross term is −1/N:
        # the free site is coupled to every occupied state, which is what makes
        # saturation a shared effect rather than a per-species one.
        N = n0[4] + n0[5] + n0[6]
        @test J[4, 4] ≈ 1 / n0[4] - 1 / N rtol = 1.0e-6
        @test J[4, 5] ≈ -1 / N rtol = 1.0e-6
        @test J[5, 4] ≈ -1 / N rtol = 1.0e-6
    end

    @testset "competitive Langmuir falls out of the site balance" begin
        logK1, logK2, n_sites = 7.0, -9.0, 1.0e-3
        cs, st, _ = _amphoteric_system(; logK1, logK2, n_sites)
        des = DualEquilibriumSolver(cs, DiluteSolutionModel())
        A = Float64.(cs.SM.A)
        b = A * Float64[ustrip(us"mol", x) for x in st.n]

        for target_pH in (4.0, 6.0, 8.0, 10.0)
            q = Base.RefValue{Any}(nothing)
            eq = SciMLBase.solve(
                des, st; b = b, constraint = FixedpH(target_pH), parameters = q,
            )
            n = Float64[ustrip(us"mol", x) for x in eq.n]

            # The prescribed quantity is the activity, so read it back rather
            # than assume the target was reached.
            p = ChemistryLab._build_params(eq)
            a_H = exp(des.lna(n, p)[2])

            β = [10.0^logK1 * a_H, 10.0^logK2 / a_H]
            x_free, x_occ = _langmuir(β)

            N = n[4] + n[5] + n[6]
            @test N ≈ n_sites rtol = 1.0e-8            # sites are conserved
            @test n[4] / N ≈ x_free rtol = 1.0e-6
            @test n[5] / N ≈ x_occ[1] rtol = 1.0e-6
            @test n[6] / N ≈ x_occ[2] rtol = 1.0e-6
        end
    end

    @testset "saturation, and the surface activity coefficient it makes unnecessary" begin
        # A drive strong enough to fill the sites: at pH 2 with logK₁ = 7 the
        # protonation drive is β₁ = 1e5, so the free sites are nearly gone.
        logK1, logK2, n_sites = 7.0, -9.0, 1.0e-3
        cs, st, _ = _amphoteric_system(; logK1, logK2, n_sites)
        des = DualEquilibriumSolver(cs, DiluteSolutionModel())
        b = Float64.(cs.SM.A) * Float64[ustrip(us"mol", x) for x in st.n]
        q = Base.RefValue{Any}(nothing)
        eq = SciMLBase.solve(des, st; b = b, constraint = FixedpH(2.0), parameters = q)
        n = Float64[ustrip(us"mol", x) for x in eq.n]

        N = n[4] + n[5] + n[6]
        θ = (n[5] + n[6]) / N                      # occupied fraction
        @test θ > 0.999
        @test n[4] / N < 1.0e-3                    # free sites nearly exhausted

        # Kulik's 2006 surface activity coefficient, γ_L = 1/(1 − θ), is what an
        # *eliminated* free site needs in order to reproduce what an explicit one
        # already produces. Here the identity is a check, not a correction:
        # applying γ_L on top of this mixing would count saturation twice.
        @test 1 / (1 - θ) ≈ N / n[4] rtol = 1.0e-8
    end

    @testset "a system without a surface is unchanged" begin
        # The kernel must be a no-op when nothing declares a site, in all four
        # activity models — the guard against a regression reaching every solve.
        h2o = Species("H2O@"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLVENT)
        h2o[:ΔₐG⁰] = _g0(-237181.0)
        h2o[:M] = 0.018015u"kg/mol"
        hp = Species("H+"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
        hp[:ΔₐG⁰] = _g0(0.0)
        oh = Species("OH-"; aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE)
        oh[:ΔₐG⁰] = _g0(-157297.0)
        cs = ChemicalSystem([h2o, hp, oh], [h2o, hp])
        st = ChemicalState(cs, Any[55.5u"mol", 1.0e-7u"mol", 1.0e-7u"mol"])
        n = Float64[ustrip(us"mol", x) for x in st.n]
        p = ChemistryLab._build_params(st)

        @test isempty(cs.site_groups)
        for model in (DiluteSolutionModel(), HKFActivityModel(), DaviesActivityModel())
            lna = activity_model(cs, model)(n, p)
            @test all(isfinite, lna)
            @test length(lna) == 3
        end
    end

end
