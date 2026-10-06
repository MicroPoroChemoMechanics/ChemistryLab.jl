# Standard properties from 0 to 1000 °C and 1 to 5000 bar against the ThermoFun
# library, as Reaktoro 2.13 embeds it, reading the same Cemdata18 file: the
# solvent by the equation of state of water its record declares, aqueous species
# by HKF where water is dense enough for those equations, minerals and a gas by
# their heat-capacity functions. The fixture is written by
# test/reference/thermofun_standard_properties.py.
#
# And the identities the solvent obeys whatever the oracle: one Gibbs energy for
# its five functions, the triple-point convention of SUPCRT92.

include("reference_species.jl")
using JSON

@testsection "Standard properties at high temperature and pressure" begin
    CL = ChemistryLab
    TF = reference_oracle("thermofun_standard_properties")
    cem = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-thermofun.json"); verbose = false))
    rows(name) = [p for p in TF.points if p.species == name]
    raw = JSON.parsefile(datapath("cemdata18-thermofun.json"))
    # The values the solvent's record tabulates at 25 °C and 1 bar, read from
    # the file: ΔfG°, ΔfH° and S°.
    water_record() = let r = only(x for x in raw["substances"] if x["symbol"] == "H2O@")
        Tuple(Float64(r[k]["values"][1]) for k in ("sm_gibbs_energy", "sm_enthalpy", "sm_entropy_abs"))
    end

    @testset "the solvent, the equation of state of water, as SUPCRT92 refers it" begin
        # With the oracle's molar mass the specific properties of the equation of
        # state and the triple-point conversion reproduce it to rounding. The
        # library writes 273.15 K where SUPCRT92 writes the triple point,
        # 273.16 K, in the constant term of the Gibbs energy (Ttr S_tr), which
        # puts its Gibbs energy of water 0.633 J/mol below the published
        # conversion's at every state; its enthalpy and entropy are unaffected.
        tp = CL._WATER_TRIPLE_POINT
        M = TF.water_molar_mass
        library_offset = (tp.T - 273.15) * tp.S
        @test library_offset ≈ 0.633 atol = 5.0e-4
        for p in rows("H2O@")
            h = CL._hgk_specific_enthalpy(p.T, p.P)
            s = CL._hgk_specific_entropy(p.T, p.P)
            @test M * h + tp.H ≈ p.H rtol = 1.0e-10
            @test M * h - p.T * (tp.S + M * s) + tp.T * tp.S + tp.G - library_offset ≈ p.G rtol = 1.0e-10
            @test M / CL._hgk_density(p.T, p.P) ≈ p.V rtol = 1.0e-9
            @test M * ForwardDiff.derivative(t -> CL._hgk_specific_enthalpy(t, p.P), p.T) ≈ p.Cp rtol = 1.0e-8
        end
    end

    @testset "the solvent: the record's values at 25 °C, the equation's increments" begin
        # The package's solvent is SUPCRT92's water shifted by three constants,
        # so that it takes the record's tabulated values at 25 °C and 1 bar, and
        # its molar mass is the one its atomic masses give, 18.015 g/mol against
        # the oracle's 18.015268. Both differences are recomposed exactly.
        tp = CL._WATER_TRIPLE_POINT
        M = TF.water_molar_mass
        library_offset = (tp.T - 273.15) * tp.S
        w = cem["H2O@"]
        Mw = ustrip(us"kg/mol", w[:M])
        Tr, Pr = 298.15, 1.0e5
        G_tab, H_tab, S_tab = water_record()
        hr, sr = CL._hgk_specific_enthalpy(Tr, Pr), CL._hgk_specific_entropy(Tr, Pr)
        δS = S_tab - (tp.S + Mw * sr)
        δH = H_tab - (tp.H + Mw * hr)
        δG = G_tab - (Mw * hr - Tr * (tp.S + Mw * sr) + tp.T * tp.S + tp.G)
        # The constants SUPCRT92 takes from Helgeson and Kirkham (1974) put the
        # Gibbs energy 1.35 J/mol and the enthalpy 48.6 J/mol above the record's.
        @test δG ≈ -1.352 atol = 1.0e-3
        @test δH ≈ -48.607 atol = 1.0e-3
        @test abs(δS) < 1.0e-3
        for p in rows("H2O@")
            h = CL._hgk_specific_enthalpy(p.T, p.P)
            s = CL._hgk_specific_entropy(p.T, p.P)
            @test w[:ΔₐH⁰](T = p.T, P = p.P) ≈ p.H + (Mw - M) * h + δH rtol = 1.0e-10
            @test w[:ΔₐG⁰](T = p.T, P = p.P) ≈
                p.G + library_offset + (Mw - M) * (h - p.T * s) + δG - δS * (p.T - Tr) rtol = 1.0e-10
            @test w[:V⁰](T = p.T, P = p.P) ≈ p.V * Mw / M rtol = 1.0e-9
            @test w[:Cp⁰](T = p.T, P = p.P) ≈ p.Cp * Mw / M rtol = 1.0e-8
        end
        # At the reference, the record's values, to the last bit.
        @test w[:ΔₐG⁰](T = Tr, P = Pr) == G_tab
        @test w[:ΔₐH⁰](T = Tr, P = Pr) == H_tab
        @test w[:S⁰](T = Tr, P = Pr) == S_tab
        # The element entropies of the file require G − H + T S = Tr Σ S_el of
        # every species, so that a reaction's enthalpy is the one its Gibbs
        # energy implies: the record's values meet it within 2.4 J/mol, the
        # triple-point constants miss it by 45.
        S_el = Dict(e["symbol"] => Float64(e["entropy"]["values"][1]) for e in raw["elements"])
        Sel = 2 * S_el["H"] + S_el["O"]
        @test abs(G_tab - H_tab + Tr * S_tab - Tr * Sel) < 2.5
        @test tp.T * tp.S + tp.G - tp.H - Tr * Sel ≈ -45.0 atol = 0.5
    end

    @testset "the solvent: one Gibbs energy for its five functions" begin
        w = cem["H2O@"]
        G_tab, H_tab, S_tab = water_record()
        for (T, P) in ((283.15, 1.0e5), (298.15, 1.0e5), (353.15, 2.0e6), (573.15, 3.0e7), (873.15, 2.0e8))
            G(T, P) = w[:ΔₐG⁰](T = T, P = P)
            @test -ForwardDiff.derivative(t -> G(t, P), T) ≈ w[:S⁰](T = T, P = P) rtol = 1.0e-10
            @test ForwardDiff.derivative(p -> G(T, p), P) ≈ w[:V⁰](T = T, P = P) rtol = 1.0e-10
            @test ForwardDiff.derivative(t -> w[:ΔₐH⁰](T = t, P = P), T) ≈ w[:Cp⁰](T = T, P = P) rtol = 1.0e-10
            # G − H + T S is the constant of the record's values at 25 °C.
            @test G(T, P) - w[:ΔₐH⁰](T = T, P = P) + T * w[:S⁰](T = T, P = P) ≈
                G_tab - H_tab + 298.15 * S_tab rtol = 1.0e-12
        end
        # The heat capacity is the equation's own second derivative, and the
        # derivative of the enthalpy.
        for (T, P) in ((278.15, 1.0e5), (473.15, 5.0e7))
            @test CL._hgk_specific_state(T, P).cp ≈ ForwardDiff.derivative(t -> CL._hgk_specific_enthalpy(t, P), T) rtol = 1.0e-10
        end
        # Untagged dual numbers, as a caller differentiating by hand passes them,
        # and a derivative of a derivative: no derivative is nested inside.
        d = w[:S⁰](T = ForwardDiff.Dual(298.15, 1.0), P = ForwardDiff.Dual(1.0e5, 0.0))
        @test ForwardDiff.partials(d)[1] ≈ w[:Cp⁰](T = 298.15, P = 1.0e5) / 298.15 rtol = 1.0e-10
        # A second derivative goes through the density lifted into the dual
        # numbers by two Newton steps, exact to the first order and good to
        # 1e-7 at the second.
        @test ForwardDiff.derivative(t -> ForwardDiff.derivative(u -> w[:ΔₐH⁰](T = u, P = 1.0e5), t), 320.0) ≈
            ForwardDiff.derivative(t -> w[:Cp⁰](T = t, P = 1.0e5), 320.0) rtol = 1.0e-6
        # Away from the reference, the same function as at a differentiated
        # reference.
        @test w[:ΔₐG⁰](T = 298.15, P = 1.0e5) ≈ ForwardDiff.value(w[:ΔₐG⁰](T = ForwardDiff.Dual(298.15, 1.0), P = 1.0e5)) rtol = 1.0e-14
    end

    @testset "aqueous species by HKF, where water is dense" begin
        for name in ("CO2@", "Ca+2", "OH-", "HCO3-")
            s = cem[name]
            for p in rows(name)
                @test s[:ΔₐG⁰](T = p.T, P = p.P) ≈ p.G atol = 0.5
                @test s[:ΔₐH⁰](T = p.T, P = p.P) ≈ p.H atol = 2.0
                @test s[:V⁰](T = p.T, P = p.P) ≈ p.V atol = 1.0e-9 rtol = 1.0e-4
            end
        end
    end

    @testset "minerals and a gas, by their heat-capacity functions" begin
        for name in ("Portlandite", "Cal", "CO2")
            s = cem[name]
            for p in rows(name)
                @test s[:ΔₐG⁰](T = p.T, P = p.P) ≈ p.G rtol = 1.0e-12
                @test s[:ΔₐH⁰](T = p.T, P = p.P) ≈ p.H rtol = 1.0e-12
                @test s[:Cp⁰](T = p.T, P = p.P) ≈ p.Cp rtol = 1.0e-12
            end
        end
    end

    @testset "the density of water: every root, every start" begin
        # The liquid where it exists, the vapor or the supercritical fluid
        # otherwise; at 0 °C and 5 kbar the earlier iteration returned 150.8.
        for (T, P, ρ_range) in (
                (273.15, 5.0e8, (1150, 1155)), (298.15, 1.0e5, (997.0, 997.1)),
                (423.15, 1.0e5, (916, 918)), (623.15, 1.0e5, (0.34, 0.36)),
                (673.15, 1.0e5, (0.32, 0.33)), (1273.15, 5.0e8, (609, 610)),
            )
            ρ = [CL.water_density_hgk(T, P; D0 = D0) for D0 in (1000.0, 1200.0, 1.0)]
            @test ρ_range[1] < ρ[1] < ρ_range[2]
            @test ρ[2] ≈ ρ[1] rtol = 1.0e-11
            @test ρ[3] ≈ ρ[1] rtol = 1.0e-11
        end
        # The lift: dρ/dP = 1/(∂P/∂ρ), the implicit-function theorem.
        ρ = CL.water_density_hgk(298.15, 1.0e7)
        @test ForwardDiff.derivative(p -> CL.water_density_hgk(298.15, p), 1.0e7) ≈
            1 / CL._hgk_pressure(298.15, ρ)[2] rtol = 1.0e-10
    end

    @testset "HKF refuses water too dilute for its equations" begin
        # At 350 °C and 1 bar water is a vapor: no aqueous solution to describe.
        @test_throws DomainError cem["Ca+2"][:ΔₐG⁰](T = 623.15, P = 1.0e5)
        @test_throws DomainError cem["OH-"][:V⁰](T = 1073.15, P = 1.0e6)
        @test cem["H2O@"][:V⁰](T = 623.15, P = 1.0e5) > 1.0e-2     # the vapor's volume
    end
end
