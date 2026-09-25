# The C-S-H surface of Guo et al. (2018), against PHREEQC on the same model.
#
# Silanol sites on the C-S-H, with the four reactions of Guo's Table 1 that
# these solutions exercise — deprotonation, calcium, chloride and sodium — and a
# diffuse layer. The fixture, test/reference/phreeqc_csh_surface.json, is written
# by test/reference/phreeqc_csh_surface.py from the same table of
# data/literature/Guo2018.json that `site_family` reads here: one source for the
# reactions and their constants, on both sides.
#
# What is compared is the surface. The aqueous side is held identical — the free
# ions and the ion product of water, Davies activity coefficients on both sides —
# and the systems are closed and of the same total composition, so the pH is an
# output of both codes rather than an input to either.

include("reference_species.jl")

const PHREEQC_CSH = reference_oracle("phreeqc_csh_surface")

@testsection "C-S-H surface against PHREEQC" begin
    f = PHREEQC_CSH
    aq = Dict(
        symbol(s) => s for s in reference_species(("H2O@", "H+", "OH-", "Na+", "Ca+2", "Cl-"))
    )
    G(s) = ustrip(us"J/mol", s[:ΔₐG⁰](T = 298.15u"K", P = 1.0e5u"Pa"; unit = true))

    @testset "the two codes share their inputs" begin
        # The ion product of water PHREEQC was given is the one these species
        # carry: the generator asks the package for the energies.
        logkw = -(G(aq["OH-"]) + G(aq["H+"]) - G(aq["H2O@"])) / (R_GAS * 298.15 * log(10))
        @test logkw ≈ f.log_kw atol = 1.0e-10
        # And its reactions are the table this test reads.
        t = literature_table("Guo2018", "surface_reactions_phreeqc")
        mine = [eq for eq in t.reaction if !occursin("K+", eq)]
        @test [r.equation for r in f.reactions] == mine
        # The site density the fixture used is the one the file gives.
        @test f.n_sites ≈ ustrip(u"mol/g", literature_value("Guo2018", "csh_site_density")) *
            f.csh_grams_per_kg_water
    end

    t = literature_table("Guo2018", "surface_reactions_phreeqc")
    reactions = [eq => lk for (eq, lk) in zip(t.reaction, t.log_K) if !occursin("K+", eq)]
    area = ustrip(u"m^2/g", literature_value("Guo2018", "csh_specific_surface_area")) *
        f.csh_grams_per_kg_water
    dl = DiffuseLayer(; area)
    family = site_family(
        "Csh_w", reactions, collect(values(aq)); master = "Csh_w", site = "Xw",
        capacity = TotalSiteAmount(f.n_sites * u"mol"),
        support = SurfaceSupport("C-S-H", nothing, FixedSurfaceArea(area)), model = dl,
    )
    aqueous = [aq[s] for s in ("H2O@", "H+", "OH-", "Na+", "Ca+2", "Cl-")]
    members = vcat([family.free_site], family.complexes)
    cs = ChemicalSystem(
        vcat(aqueous, members),
        [aq["H2O@"], aq["H+"], aq["Na+"], aq["Ca+2"], aq["Cl-"], family.free_site];
        site_families = [family],
    )
    idx = Dict(symbol(sp) => k for (k, sp) in enumerate(cs.species))
    model = DaviesActivityModel()
    z = Float64[charge(sp) for sp in members]
    i_H = idx["H+"]

    worst = (fraction = 0.0, la_H = 0.0, I = 0.0, psi = 0.0)
    certified = 0
    for pt in f.points
        n0 = Any[fill(1.0e-14u"mol", length(cs.species))...]
        n0[idx["H2O@"]] = moles_of_water() * u"mol"
        n0[idx["Na+"]] = (pt.naoh + pt.nacl) * u"mol"
        n0[idx["OH-"]] = pt.naoh * u"mol"
        n0[idx["Ca+2"]] = pt.cacl2 * u"mol"
        n0[idx["Cl-"]] = (2 * pt.cacl2 + pt.nacl) * u"mol"
        n0[idx["XwOH"]] = f.n_sites * u"mol"
        st = ChemicalState(cs, n0)
        eq, cert = equilibrate_certified(st; model)
        cert.optimal && (certified += 1)
        @test cert.optimal

        n = Float64[ustrip(us"mol", x) for x in eq.n]
        N = sum(n[idx[symbol(sp)]] for sp in members)
        for (name, frac) in pairs(pt.fractions)
            sym = replace(String(name), "Csh_w" => "Xw"; count = 1)
            worst = merge(worst, (fraction = max(worst.fraction, abs(n[idx[sym]] / N - frac)),))
        end
        p = ChemistryLab._build_params(eq)
        la_H = activity_model(cs, model)(n, p)[i_H] / log(10)
        I = ionic_strength(eq)
        psi = diffuse_layer_potential(dl, z, [n[idx[symbol(sp)]] for sp in members], I, 298.15)
        psi_phreeqc = pt.psi_V * FARADAY / (R_GAS * 298.15)
        worst = merge(
            worst, (
                la_H = max(worst.la_H, abs(la_H - pt.la_H)),
                I = max(worst.I, abs(I / pt.I - 1)),
                psi = max(worst.psi, abs(psi - psi_phreeqc)),
            ),
        )
    end
    @info "C-S-H surface against PHREEQC" certified worst
    @test certified == length(f.points)
    # Measured at 4.4e-4 on a site fraction, 2.2e-4 on log a(H+), 3e-5 relative
    # on the ionic strength and 1.5e-3 on FΨ/RT (0.04 mV); bounded here with a
    # factor of about two to five, the margin a change of solver version needs.
    @test worst.fraction < 1.0e-3
    @test worst.la_H < 1.0e-3
    @test worst.I < 1.0e-4
    @test worst.psi < 5.0e-3
end

const PHREEQC_CSH_PASTE = reference_oracle("phreeqc_csh_paste")

@testsection "Guo's paste, surface and hydrates, against PHREEQC" begin
    # The same surface in Guo's hydrated paste: portlandite, the two AFm
    # hydrates, ettringite, Friedel's and Kuzel's salts, over a NaCl sweep, all
    # from CEMDATA18 on both sides. The chloride the paste binds is split between
    # the salts and the surface, and that split is what is compared.
    f = PHREEQC_CSH_PASTE
    subs = Dict(
        symbol(s) => s for s in build_species(datapath("cemdata18-thermofun.json"); verbose = false)
    )
    ions = [subs[s] for s in ("H2O@", "H+", "OH-", "Na+", "Cl-", "Ca+2", "AlO2-", "SO4-2")]
    solids = [subs[String(ph)] for ph in keys(f.initial_mol)]

    @testset "the two codes share their inputs" begin
        # The log K PHREEQC was given are the ones these species carry, and the
        # moles it started from are Guo's grams over the molar masses the
        # package computed for these solids.
        G(s) = ustrip(us"J/mol", s[:ΔₐG⁰](T = 298.15u"K", P = 1.0e5u"Pa"; unit = true))
        for (ph, prod) in pairs(f.phase_reactions)
            drg = sum(ν * G(subs[String(sp)]) for (sp, ν) in pairs(prod)) - G(subs[String(ph)])
            @test -drg / (R_GAS * 298.15 * log(10)) ≈ f.phase_log_k[ph] atol = 1.0e-10
            # Each reaction is balanced in atoms and in charge.
            lhs = atoms(subs[String(ph)])
            rhs = Dict{Symbol, Float64}()
            for (sp, ν) in pairs(prod), (el, k) in atoms(subs[String(sp)])
                rhs[el] = get(rhs, el, 0.0) + ν * k
            end
            @test all(get(rhs, el, 0.0) ≈ k for (el, k) in lhs)
            @test sum(ν * charge(subs[String(sp)]) for (sp, ν) in pairs(prod)) == 0
        end
        # Guo's inventory is per liter of concrete; a liter holds 10 × porosity
        # grams of pore water, and the fixture is per kilogram of it.
        scale = 1000 / (10 * ustrip(literature_value("Guo2018", "porosity_percent")))
        guo = Dict(
            :Portlandite => "ch_per_liter", :monosulphate12 => "afm_per_liter",
            :ettringite => "aft_per_liter",
        )
        for (ph, g) in pairs(f.grams_per_kg_water)
            expected = haskey(guo, ph) ? ustrip(us"g", literature_value("Guo2018", guo[ph])) * scale : 0.0
            @test g ≈ expected rtol = 1.0e-12
            M = ustrip(us"g/mol", subs[String(ph)][:M])
            @test f.initial_mol[ph] ≈ g / M rtol = 1.0e-12
        end
        @test f.n_sites ≈ ustrip(u"mol/g", literature_value("Guo2018", "csh_site_density")) *
            ustrip(us"g", literature_value("Guo2018", "csh_per_liter")) * scale rtol = 1.0e-12
    end

    t = literature_table("Guo2018", "surface_reactions_phreeqc")
    reactions = [eq => lk for (eq, lk) in zip(t.reaction, t.log_K) if !occursin("K+", eq)]
    dl = DiffuseLayer(; area = f.area_m2)
    family = site_family(
        "Csh_w", reactions, ions; master = "Csh_w", site = "Xw",
        capacity = TotalSiteAmount(f.n_sites * u"mol"),
        support = SurfaceSupport("C-S-H", nothing, FixedSurfaceArea(f.area_m2)), model = dl,
    )
    members = vcat([family.free_site], family.complexes)
    cs = ChemicalSystem(
        vcat(ions, solids, members),
        [
            subs["H2O@"], subs["H+"], subs["Na+"], subs["Ca+2"], subs["Cl-"], subs["AlO2-"],
            subs["SO4-2"], family.free_site,
        ];
        site_families = [family],
    )
    idx = Dict(symbol(sp) => k for (k, sp) in enumerate(cs.species))
    model = DaviesActivityModel()

    function charged(nacl)
        n0 = Any[fill(1.0e-14u"mol", length(cs.species))...]
        n0[idx["H2O@"]] = moles_of_water() * u"mol"
        for (ph, g) in pairs(f.grams_per_kg_water)
            g > 0 && (n0[idx[String(ph)]] = g * u"g" / subs[String(ph)][:M])
        end
        n0[idx["XwOH"]] = f.n_sites * u"mol"
        s = max(nacl, 1.0e-10)
        n0[idx["Na+"]] = s * u"mol"
        n0[idx["Cl-"]] = s * u"mol"
        return ChemicalState(cs, n0)
    end

    # Walked down from the saltiest point, where the assemblage is simplest.
    pts = sort(collect(f.points); by = p -> -p.nacl)
    A = Float64.(cs.SM.A)
    budgets = [A * ustrip.(us"mol", charged(p.nacl).n) for p in pts]
    states, certs = equilibrate_path(charged(first(pts).nacl), budgets; model)
    @test all(c.optimal for c in certs)

    M_w = ustrip(us"kg/mol", subs["H2O@"][:M])
    worst = (phase = 0.0, surface = 0.0, la_H = 0.0, psi = 0.0, cl_bound = 0.0, water = 0.0)
    for (pt, eq) in zip(pts, states)
        n = Float64[ustrip(us"mol", x) for x in eq.n]
        # The water the hydrates and the surface traded with the solution: an
        # amount compared as a molality would be off by this ratio, 9 % here.
        worst = merge(worst, (water = max(worst.water, abs(n[idx["H2O@"]] * M_w / pt.water_kg - 1)),))
        for (ph, m) in pairs(pt.phases)
            worst = merge(worst, (phase = max(worst.phase, abs(n[idx[String(ph)]] - m)),))
        end
        for (sp, m) in pairs(pt.surface)
            sym = replace(String(sp), "Csh_w" => "Xw"; count = 1)
            worst = merge(worst, (surface = max(worst.surface, abs(n[idx[sym]] - m)),))
        end
        p = ChemistryLab._build_params(eq)
        la_H = activity_model(cs, model)(n, p)[idx["H+"]] / log(10)
        z = Float64[charge(sp) for sp in members]
        psi = diffuse_layer_potential(
            dl, z, [n[idx[symbol(sp)]] for sp in members], ionic_strength(eq), 298.15,
        )
        # The chloride the paste holds out of solution: in the two salts and on
        # the surface, against PHREEQC's.
        bound = 2n[idx["C4AClH10"]] + n[idx["C4AsClH12"]] + n[idx["XwOHCl-"]]
        bound_ref = 2pt.phases.C4AClH10 + pt.phases.C4AsClH12 + pt.surface[Symbol("Csh_wOHCl-")]
        worst = merge(
            worst, (
                la_H = max(worst.la_H, abs(la_H + pt.pH)),
                psi = max(worst.psi, abs(psi - pt.psi_V * FARADAY / (R_GAS * 298.15))),
                cl_bound = max(worst.cl_bound, abs(bound - bound_ref)),
            ),
        )
    end
    @info "Guo's paste against PHREEQC" worst
    # Measured at 2.9e-4 mol on a phase, 4.5e-4 mol on a surface species,
    # 1.0e-4 on log a(H+), 3.9e-4 on FΨ/RT, 5e-5 mol on the bound chloride and
    # 8.5e-6 relative on the water; bounded with the same margin as above.
    @test worst.phase < 1.0e-3
    @test worst.surface < 1.0e-3
    @test worst.la_H < 1.0e-3
    @test worst.psi < 5.0e-3
    @test worst.cl_bound < 1.0e-3
    @test worst.water < 1.0e-4
end
