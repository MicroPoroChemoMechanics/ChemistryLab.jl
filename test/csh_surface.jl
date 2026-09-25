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
        # carry. The generator reads the tabulated Gibbs energies; the package
        # rebuilds them from ΔfH° and S°, which moves water by 0.15 J/mol, 3e-5
        # in log K — an order of magnitude below the agreement asserted below.
        logkw = -(G(aq["OH-"]) + G(aq["H+"]) - G(aq["H2O@"])) / (R_GAS * 298.15 * log(10))
        @test logkw ≈ f.log_kw atol = 1.0e-4
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
