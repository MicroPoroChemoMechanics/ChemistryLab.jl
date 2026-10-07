# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# A surface away from 25 °C. A published surface constant carries no enthalpy,
# and PHREEQC then holds its log K constant at every temperature; `site_family`
# does the same, the energy of each complex following the aqueous species it is
# written with. Compared with PHREEQC at 50 °C (test/reference/
# phreeqc_hfo_surface.py --temp 50): the weak sites of hydrous ferric oxide,
# without and with the diffuse layer.

using ChemistryLab
using DynamicQuantities
using SciMLBase
using Test

include("reference_species.jl")

const PHREEQC_PROTOLYSIS_50C = reference_oracle("phreeqc_protolysis_50C")
const PHREEQC_DIFFUSE_LAYER_50C = reference_oracle("phreeqc_diffuse_layer_50C")

# The weak sites of the fixture, from their two reactions as PHREEQC writes them.
function _hfo_weak(f; model = IdealSiteMixing(), T = 298.15u"K", m_na = 0.01, m_cl = 0.01, Tstate = 323.15)
    h2o, hp, oh, na, cl = reference_species(("H2O@", "H+", "OH-", "Na+", "Cl-"))
    reactions = [
        "Hfo_wOH + H+ = Hfo_wOH2+" => f.logK_protonation,
        "Hfo_wOH = Hfo_wO- + H+" => f.logK_deprotonation,
    ]
    area = hasproperty(f, :area) ? f.area : 1.0
    family = site_family(
        "Hfo_w", reactions, [h2o, hp, oh]; master = "Hfo_w", site = "Xw",
        capacity = TotalSiteAmount(f.n_sites * u"mol"),
        support = SurfaceSupport("hfo", nothing, FixedSurfaceArea(area)), model, T,
    )
    members = vcat([family.free_site], family.complexes)
    cs = ChemicalSystem(vcat([h2o, hp, oh, na, cl], members), [h2o, hp, na, cl, family.free_site]; site_families = [family])
    sym = symbol.(cs.species)
    idx(s) = findfirst(==(s), sym)
    n = Any[fill(1.0e-12u"mol", length(cs.species))...]
    n[idx("H2O@")] = moles_of_water() * u"mol"
    n[idx("Na+")] = m_na * u"mol"
    n[idx("Cl-")] = m_cl * u"mol"
    n[idx("XwOH")] = f.n_sites * u"mol"
    return cs, ChemicalState(cs, n; T = Tstate * u"K"), idx
end

function _fractions(cs, st, idx, pH)
    des = DualEquilibriumSolver(cs, DiluteSolutionModel())
    b = Float64.(cs.SM.A) * Float64[ustrip(us"mol", x) for x in st.n]
    eq = SciMLBase.solve(
        des, st; b, constraint = FixedpH(pH), surface_potential = :auto,
        parameters = Base.RefValue{Any}(nothing)
    )
    cert = optimality_certificate(des, eq; b, constraint = FixedpH(pH))
    n(s) = ustrip(us"mol", eq.n[idx(s)])
    N = n("XwOH") + n("XwOH2+") + n("XwO-")
    return (; cert, free = n("XwOH") / N, protonated = n("XwOH2+") / N, deprotonated = n("XwO-") / N)
end

@testsection "a surface away from 25 °C" begin
    f = PHREEQC_PROTOLYSIS_50C

    @testset "the log K of a surface reaction is the same at every temperature" begin
        cs, _, idx = _hfo_weak(f)
        G(s, T) = ustrip(us"J/mol", cs.species[idx(s)][:ΔₐG⁰](T = T * u"K", P = 1.0e5u"Pa"; unit = true))
        logK(T) = -(G("XwOH2+", T) - G("XwOH", T) - G("H+", T)) / (ChemistryLab.R_GAS * T * log(10))
        for T in (278.15, 298.15, 323.15, 353.15)
            @test logK(T) ≈ f.logK_protonation atol = 1.0e-10
        end
        # And the temperature the family is built at no longer matters.
        cs50, _, idx50 = _hfo_weak(f; T = 323.15u"K")
        for T in (298.15, 323.15)
            @test ustrip(us"J/mol", cs50.species[idx50("XwO-")][:ΔₐG⁰](T = T * u"K", P = 1.0e5u"Pa"; unit = true)) ≈
                G("XwO-", T) rtol = 1.0e-12
        end
    end

    @testset "without the diffuse layer, PHREEQC at 50 °C" begin
        # Its fractions at 50 °C are those at 25 °C at the same proton activity,
        # its constants being held: the identity is the oracle's own.
        worst = 0.0
        for pt in f.points
            cs, st, idx = _hfo_weak(f)
            r = _fractions(cs, st, idx, pt.pH)
            @test r.cert.stationarity < 1.0e-8
            for (got, want) in zip((r.free, r.protonated, r.deprotonated), (pt.free, pt.protonated, pt.deprotonated))
                worst = max(worst, abs(got - want))
            end
        end
        @test worst < 1.0e-8
    end

    @testset "with the diffuse layer, PHREEQC at 50 °C" begin
        # The layer screens with the permittivity of water at the temperature of
        # the solve, 69.8 at 50 °C, and agrees with PHREEQC as it does at 25 °C.
        # Held at its 25 °C value, 78.2, as the layer did until 0.35, it misses
        # by a hundred times more: PHREEQC follows the temperature too.
        d = PHREEQC_DIFFUSE_LAYER_50C
        worst_water, worst_held = 0.0, 0.0
        held = DiffuseLayer(; area = d.area, ε_r = water_relative_permittivity(298.15))
        for ser in d.series, pt in ser.points
            mH, mOH = 10.0^(-pt.pH), 10.0^(pt.pH - 14)
            m_cl = 2 * pt.I - ser.nacl - mH - mOH
            for (model, slot) in ((DiffuseLayer(; area = d.area), 1), (held, 2))
                cs, st, idx = _hfo_weak(d; model, m_na = ser.nacl, m_cl)
                r = _fractions(cs, st, idx, pt.pH)
                @test r.cert.stationarity < 1.0e-8
                gap = maximum(abs(g - w) for (g, w) in zip((r.protonated, r.deprotonated), (pt.protonated, pt.deprotonated)))
                slot == 1 ? (worst_water = max(worst_water, gap)) : (worst_held = max(worst_held, gap))
            end
        end
        @test worst_water < 5.0e-4
        @test worst_held > 10 * worst_water
        # The layer that follows the water is the one held, at 25 °C.
        @test ChemistryLab._permittivity(DiffuseLayer(; area = d.area), 298.15) == held.ε_r
        @test ChemistryLab._permittivity(DiffuseLayer(; area = d.area), 323.15) ≈ water_relative_permittivity(323.15)
        @test ChemistryLab._permittivity(held, 323.15) == held.ε_r
    end
end
