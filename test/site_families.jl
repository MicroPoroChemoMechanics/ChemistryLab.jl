# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using ChemistryLab
using DynamicQuantities
using Test

# Helpers — a surface species is an ordinary species in AS_SURFACE, and the
# family requalifies it anyway; building it plainly is the point.
_surf(sym) = Species(sym; aggregate_state = AS_SURFACE, class = SC_SURFCOMPLEX)
_aq(sym, cl = SC_AQSOLUTE) = Species(sym; aggregate_state = AS_AQUEOUS, class = cl)

function _hfo_system()
    free, prot, depr = _surf("XsOH"), _surf("XsOH2+"), _surf("XsO-")
    species = [
        _aq("H2O@", SC_AQSOLVENT), _aq("H+"), _aq("OH-"), _aq("Ca+2"),
        free, prot, depr,
    ]
    support = SurfaceSupport("hydrous ferric oxide", nothing, FixedSurfaceArea(600.0))
    family = SiteFamily(
        "Hfo_s", free, [prot, depr];
        capacity = TotalSiteAmount(5.0e-6), support,
    )
    # The primary is the **free site species**, never a bare `Species("Xs")`.
    primaries = [species[1], species[2], species[4], free]
    return ChemicalSystem(species, primaries; site_families = [family]), family
end

@testsection "Surface site families" begin

    @testset "site symbols" begin
        @test length(SITE_SYMBOLS) == 24
        @test is_site_symbol(:Xs) && is_site_symbol(:Xw) && is_site_symbol(:Xv)
        # Xenon is a real element, and `Xx` is this package's own placeholder
        # for "not an element" — a typo must not become a valid site family.
        @test !is_site_symbol(:Xe)
        @test !is_site_symbol(:Xx)
        @test :Xe ∉ SITE_SYMBOLS && :Xx ∉ SITE_SYMBOLS
        @test !is_site_symbol(:Ca) && !is_site_symbol(:Zz) && !is_site_symbol(:X)

        # Every site symbol has a place in the canonical ordering, and they all
        # come before the charge placeholder, which `speciation` relies on
        # staying last.
        for s in SITE_SYMBOLS
            @test s in ATOMIC_ORDER
            @test findfirst(==(s), ATOMIC_ORDER) < findfirst(==(:Zz), ATOMIC_ORDER)
        end
    end

    @testset "a site in a formula" begin
        f = Formula("XsOH")
        @test ChemistryLab.composition(f) == Dict(:Xs => 1, :O => 1, :H => 1)
        @test check_mendeleev(f)

        # Denticity is the coefficient, read from the formula.
        @test ChemistryLab.composition(Formula("Xs2OCa"))[:Xs] == 2

        # A symbol that is neither an element nor a site is still refused, so a
        # typo does not quietly become a site.
        @test !check_mendeleev(Formula("QqOH"))

        # The mass is the adsorbed part alone: the support is weighed once, by
        # its own mineral species, and counting it on every occupied site would
        # be double counting.
        sp = _surf("XsOH")
        m_site = ustrip(us"kg/mol", sp[:M])
        m_water = ustrip(us"kg/mol", Species("H2O")[:M])
        m_h = ustrip(us"kg/mol", Species("H2")[:M]) / 2
        @test m_site ≈ m_water - m_h rtol = 1.0e-10

        # A formula carrying a symbol in no ordering names it, rather than
        # raising `isless(::Int64, ::Nothing)` from inside a sort.
        @test_throws ArgumentError Formula(Dict(:Zq => 1, :O => 1))
    end

    @testset "with_aggregate_state" begin
        sp = Species("XsOH"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)
        moved = with_aggregate_state(sp, AS_SURFACE)
        @test aggregate_state(moved) == AS_SURFACE
        @test class(moved) == SC_COMPONENT      # everything else preserved
        @test symbol(moved) == symbol(sp)
        @test atoms(moved) == atoms(sp)
    end

    @testset "capacities" begin
        @test site_capacity(_hfo_system()[2]) isa TotalSiteAmount
        @test TotalSiteAmount(5.0e-6u"mol").N ≈ 5.0e-6 rtol = 1.0e-12
        @test MassSiteDensity(2.0u"mol/kg").q ≈ 2.0 rtol = 1.0e-12
        @test AreaSiteDensity(3.84e-6u"mol/m^2").Γ_C ≈ 3.84e-6 rtol = 1.0e-12

        # Units are converted or refused, never stripped in silence.
        @test_throws ArgumentError TotalSiteAmount(5.0u"kg")
        @test_throws ArgumentError AreaSiteDensity(1.0u"mol/kg")
        # A negative number of sites is not a small one.
        @test_throws ArgumentError TotalSiteAmount(-1.0)

        support = SurfaceSupport("s", nothing, BETSurfaceArea(600.0))
        # A prescribed total ignores the host; the other two do not, which is
        # what will let a support that precipitates carry its sites with it.
        @test site_moles(TotalSiteAmount(5.0e-6), support, 0.0, 0.0, 0.0) == 5.0e-6
        @test site_moles(MassSiteDensity(2.0), support, 0.01, 0.01, 0.1) ≈ 2.0 * 0.01 * 0.1
        @test site_moles(AreaSiteDensity(1.0e-5), support, 0.01, 0.01, 0.1) ≈
            1.0e-5 * 600.0 * 0.01 * 0.1
    end

    @testset "a family, and what it refuses" begin
        _, family = _hfo_system()
        @test name(family) == "Hfo_s"
        @test family.site === :Xs
        @test symbol.(site_members(family)) == ["XsOH", "XsOH2+", "XsO-"]
        @test symbol(first(site_members(family))) == "XsOH"     # free site first
        @test all(s -> aggregate_state(s) == AS_SURFACE, site_members(family))
        @test all(s -> class(s) == SC_SURFCOMPLEX, site_members(family))
        @test denticity(family, family.free_site) == 1
        @test denticity(family, _aq("H+")) == 0
        @test site_moles(family, 0.0, 0.0, 0.0) == 5.0e-6

        support = SurfaceSupport("s", nothing, FixedSurfaceArea(1.0))
        cap = TotalSiteAmount(1.0e-6)

        # No site symbol at all: nothing says which budget these share.
        @test_throws ArgumentError SiteFamily(
            "bad", _surf("SiOH"), [_surf("SiO-")]; capacity = cap, support
        )
        # Two families in one: a bridge across budgets is a different model.
        @test_throws ArgumentError SiteFamily(
            "bad", _surf("XsOH"), [_surf("XsXwO")]; capacity = cap, support
        )
        # Members disagreeing on the symbol.
        @test_throws ArgumentError SiteFamily(
            "bad", _surf("XsOH"), [_surf("XwOH2+")]; capacity = cap, support
        )
        # A free site occupying two sites is not a free site.
        @test_throws ArgumentError SiteFamily(
            "bad", _surf("Xs2OH"), [_surf("XsO-")]; capacity = cap, support
        )
        # Multidentate: the balance would hold, the ideal mixing would not.
        @test_throws ArgumentError SiteFamily(
            "bad", _surf("XsOH"), [_surf("Xs2OCa")]; capacity = cap, support
        )
        # The same species twice would be counted twice in the balance.
        @test_throws ArgumentError SiteFamily(
            "bad", _surf("XsOH"), [_surf("XsO-"), _surf("XsO-")]; capacity = cap, support
        )
    end

    @testset "the site balance is a row of the conservation matrix" begin
        cs, _ = _hfo_system()
        @test symbol.(cs.SM.primaries) == ["H2O@", "H+", "Ca+2", "XsOH"]
        @test cs.idx_surface == [5, 6, 7]
        @test cs.site_groups == [[5, 6, 7]]
        @test symbol.(surface(cs)) == ["XsOH", "XsOH2+", "XsO-"]
        @test site_families(cs) !== nothing
        @test length(site_families(cs)) == 1

        row = cs.SM.A[findfirst(p -> symbol(p) == "XsOH", cs.SM.primaries), :]
        # One per monodentate surface species, nothing on the aqueous ones:
        # that row *is* n_free + n_prot + n_depr = N_t.
        @test Int.(row) == [0, 0, 0, 0, 1, 1, 1]

        # And the aqueous decomposition is the protolysis: XsOH2+ = XsOH + H+,
        # XsO- = XsOH - H+.
        h_row = cs.SM.A[findfirst(p -> symbol(p) == "H+", cs.SM.primaries), :]
        @test Int(h_row[6]) == 1
        @test Int(h_row[7]) == -1

        # No spurious charge component. The free site is the primary, and a bare
        # neutral `Species("Xs")` in its place would flip the exact rank test
        # and add one -- silently.
        @test :Zz ∉ Symbol.(symbol.(cs.SM.primaries))
        @test size(cs.SM.A, 1) == 4
    end

    @testset "what a system refuses" begin
        support = SurfaceSupport("s", nothing, FixedSurfaceArea(1.0))
        cap = TotalSiteAmount(1.0e-6)
        free, prot = _surf("XsOH"), _surf("XsOH2+")
        base = [_aq("H2O@", SC_AQSOLVENT), _aq("H+")]
        fam = SiteFamily("f", free, [prot]; capacity = cap, support)

        # An AS_SURFACE species belonging to no family would carry a site
        # pseudo-element into the matrix, and so a row, with nothing mixing on it.
        @test_throws ArgumentError ChemicalSystem(
            vcat(base, [free, prot]), vcat(base, [free])
        )

        # A member missing from the species list.
        @test_throws ArgumentError ChemicalSystem(
            vcat(base, [free]), vcat(base, [free]); site_families = [fam]
        )

        # Two families on one pseudo-element share a budget without saying so.
        other = SiteFamily("g", _surf("XsO-"), AbstractSpecies[]; capacity = cap, support)
        @test_throws ArgumentError ChemicalSystem(
            vcat(base, [free, prot, _surf("XsO-")]), vcat(base, [free]);
            site_families = [fam, other],
        )
    end

    @testset "a system without a surface is what it always was" begin
        species = [_aq("H2O@", SC_AQSOLVENT), _aq("H+"), _aq("OH-"), _aq("Ca+2")]
        cs = ChemicalSystem(species, [species[1], species[2], species[4]])
        @test site_families(cs) === nothing
        @test isempty(cs.idx_surface)
        @test isempty(cs.site_groups)
        @test isempty(surface(cs))
        # And the matrices are untouched by the existence of the machinery.
        @test size(cs.SM.A, 1) == 3
    end

    @testset "surface matter is counted with the solid" begin
        cs, _ = _hfo_system()
        n = fill(1.0e-8u"mol", length(cs.species))
        n[1] = 55.5u"mol"
        state = ChemicalState(cs, n)

        # Three surface species at 1e-8 mol each land in the solid compartment,
        # not the liquid one, and not nowhere.
        @test ustrip(us"mol", moles(state).solid) ≈ 3.0e-8 rtol = 1.0e-10
        m_expected = sum(1.0e-8 * ustrip(us"kg/mol", cs.species[i][:M]) for i in cs.idx_surface)
        @test ustrip(us"kg", mass(state).solid) ≈ m_expected rtol = 1.0e-10
        # They carry no standard molar volume, so they add none: the volume they
        # occupy is the host's.
        @test ustrip(us"m^3", volume(state).solid) == 0.0
    end

end

@testsection "a member is matched by identity, not by its label" begin

    # `_resolve_site_families` matched members to the system by `symbol` alone.
    # A symbol is a label, and the species a label lands on need not be the one
    # the family validated. The code carried a comment asserting that a shared
    # member was unreachable; these are the inputs that reach it.

    support = SurfaceSupport("oxide", nothing, FixedSurfaceArea(600.0))
    cap = TotalSiteAmount(5.0e-6)
    aq = [_aq("H2O@", SC_AQSOLVENT), _aq("H+")]

    @testset "two families, different site symbols, one set of labels" begin
        # Same label `S1`/`S2` on both sides, different pseudo-elements behind
        # them. Resolving by name gives BOTH families the same indices, so the
        # second family's conservation row is simply absent from the matrix —
        # and nothing says so.
        s_free = Species("XsOH"; symbol = "S1", aggregate_state = AS_SURFACE, class = SC_SURFCOMPLEX)
        s_occ = Species("XsONa"; symbol = "S2", aggregate_state = AS_SURFACE, class = SC_SURFCOMPLEX)
        w_free = Species("XwOH"; symbol = "S1", aggregate_state = AS_SURFACE, class = SC_SURFCOMPLEX)
        w_occ = Species("XwONa"; symbol = "S2", aggregate_state = AS_SURFACE, class = SC_SURFCOMPLEX)

        fam_s = SiteFamily("Xs", s_free, [s_occ]; capacity = cap, support)
        fam_w = SiteFamily("Xw", w_free, [w_occ]; capacity = cap, support)

        # The two families are legitimate on their own: different pseudo-elements,
        # so the "one symbol, one family" guard does not fire. It is the label
        # collision that has to be caught, and it was not.
        @test fam_s.site === :Xs
        @test fam_w.site === :Xw

        species = vcat(aq, [s_free, s_occ])
        @test_throws ArgumentError ChemicalSystem(
            species, [aq[1], aq[2], s_free]; site_families = [fam_s, fam_w],
        )
    end

    @testset "a family's qualified copy is not the system's species" begin
        # `SiteFamily` requalifies copies of its members as AS_SURFACE. The
        # system keeps what the caller passed. Pass unqualified species to both
        # and the family believes itself on a surface while `idx_surface` is
        # empty — the site mixing runs and the phase accounting disagrees.
        raw_free = Species("XsOH")                       # AS_UNDEF, deliberately
        raw_occ = Species("XsONa")
        fam = SiteFamily("Xs", raw_free, [raw_occ]; capacity = cap, support)

        # The family did qualify its own copies…
        @test all(sp -> aggregate_state(sp) == AS_SURFACE, ChemistryLab.site_members(fam))
        # …so the mismatch with the caller's originals is real, and refused.
        species = vcat(aq, [raw_free, raw_occ])
        @test_throws ArgumentError ChemicalSystem(
            species, [aq[1], aq[2], raw_free]; site_families = [fam],
        )
    end

    @testset "same label, different formula" begin
        # The narrowest form: one label, two chemistries.
        declared = Species("XsOH"; symbol = "T", aggregate_state = AS_SURFACE, class = SC_SURFCOMPLEX)
        present = Species("XsOH2+"; symbol = "T", aggregate_state = AS_SURFACE, class = SC_SURFCOMPLEX)
        other = Species("XsO-"; aggregate_state = AS_SURFACE, class = SC_SURFCOMPLEX)
        fam = SiteFamily("Xs", declared, [other]; capacity = cap, support)
        species = vcat(aq, [present, other])
        @test_throws ArgumentError ChemicalSystem(
            species, [aq[1], aq[2], present]; site_families = [fam],
        )
    end

    @testset "the ordinary declaration still builds" begin
        # The guard must not cost anything to a system declared properly.
        cs, _ = _hfo_system()
        @test cs isa ChemicalSystem
        @test length(cs.site_groups) == 1
        @test !isempty(cs.idx_surface)
    end
end

@testsection "the site-density scale a constant refers to" begin

    # An intrinsic adsorption constant is not a property of a surface alone: it
    # is fitted at some total site density, and its value depends on that
    # choice. Kulik (2002) eq. 21 is the conversion, and he makes the point on
    # Dzombak & Morel's own two densities — which is what makes this checkable
    # against a published number rather than against itself.

    @testset "Kulik's own worked example" begin
        # Published: +log(2.254/12.05) = -0.73 for the weak sites and
        # log(0.056/12.05) = -2.33 for the strong ones.
        @test round(convert_logk_site_density(0.0, 2.254); digits = 2) == -0.73
        @test round(convert_logk_site_density(0.0, 0.056); digits = 2) == -2.33

        # The two shifts differ by 1.6 log units, which is the substance of the
        # remark: correlating one of their constants against the other without
        # converting compares two different scales.
        weak = convert_logk_site_density(0.0, 2.254)
        strong = convert_logk_site_density(0.0, 0.056)
        @test abs(weak - strong) ≈ log10(2.254 / 0.056) rtol = 1.0e-12
        @test abs(weak - strong) > 1.6
    end

    @testset "the conversion is a change of scale, with the properties of one" begin
        # Identity at the reference density: nothing to convert.
        @test convert_logk_site_density(3.7, REFERENCE_SITE_DENSITY_NM2) == 3.7
        # Additive in logK, since it only shifts.
        @test convert_logk_site_density(3.7, 2.254) - convert_logk_site_density(0.0, 2.254) ≈ 3.7
        # Reversible: converting to Γ° and back gives the original.
        there = convert_logk_site_density(2.5, 2.254)
        back = convert_logk_site_density(there, REFERENCE_SITE_DENSITY_NM2^2 / 2.254)
        @test back ≈ 2.5 rtol = 1.0e-12
        # The side the neutral group is written on flips the sign, and nothing
        # else — that is the whole content of the second half of eq. 21.
        @test convert_logk_site_density(0.0, 2.254; free_site_side = :product) ≈
            -convert_logk_site_density(0.0, 2.254)
        # Only the RATIO enters, so any consistent unit works.
        @test convert_logk_site_density(0.0, 2.254) ≈
            convert_logk_site_density(0.0, 2.254e18 / AVOGADRO; Γ0 = REFERENCE_SITE_DENSITY)
    end

    @testset "the reference density is derived, not transcribed twice" begin
        # 12.05 nm⁻² is what Kulik writes; the mol/m² form comes from it through
        # the library's Avogadro constant, so the two cannot drift apart.
        @test REFERENCE_SITE_DENSITY ≈ REFERENCE_SITE_DENSITY_NM2 * 1.0e18 / AVOGADRO
        @test REFERENCE_SITE_DENSITY ≈ 2.0e-5 rtol = 1.0e-3
        @test AVOGADRO == ustrip(us"1/mol", AVOGADRO_Q)
    end

    @testset "refusals" begin
        @test_throws ArgumentError convert_logk_site_density(0.0, 0.0)
        @test_throws ArgumentError convert_logk_site_density(0.0, -1.0)
        @test_throws ArgumentError convert_logk_site_density(0.0, 1.0; Γ0 = 0.0)
        @test_throws ArgumentError convert_logk_site_density(0.0, 1.0; free_site_side = :both)
    end
end

@testsection "a site budget that follows its host" begin

    fs = _surf("XsOH")
    mk(cap, sup) = SiteFamily("Xs", fs, AbstractSpecies[]; capacity = cap, support = sup)
    coupled(area) = SurfaceSupport("s", "Host", area; coupling = SITES_FOLLOW_HOST)
    M = 0.1                                   # kg/mol, the host's molar mass

    @testset "the coupling is asked for, never inferred" begin
        # Naming a host does not couple anything: the kinetics has named one
        # since long before, to find the amount a rate law scales with.
        @test SurfaceSupport("s", "Host", FixedSurfaceArea(1.0)).coupling === SITES_FIXED
        @test SurfaceSupport("s", FixedSurfaceArea(1.0)).coupling === SITES_FIXED
        @test coupled(BETSurfaceArea(90.0)).coupling === SITES_FOLLOW_HOST
        # And following a host with no host named is refused at construction.
        @test_throws ArgumentError SurfaceSupport(
            "s", nothing, FixedSurfaceArea(1.0); coupling = SITES_FOLLOW_HOST,
        )
        @test_throws ArgumentError SurfaceSupport(
            "s", FixedSurfaceArea(1.0); coupling = SITES_FOLLOW_HOST,
        )
    end

    @testset "ν is the coefficient, and it is the obvious product" begin
        # q·M for a mass density, Γ·a·M for an area density over a specific
        # area — both checkable by hand, which is why they are checked that way.
        @test ChemistryLab.sites_per_host(
            mk(MassSiteDensity(2.0), coupled(BETSurfaceArea(90.0))), M,
        ) ≈ 2.0 * M
        @test ChemistryLab.sites_per_host(
            mk(AreaSiteDensity(1.0e-5), coupled(BETSurfaceArea(90.0))), M,
        ) ≈ 1.0e-5 * 90.0 * M
    end

    @testset "what is refused is refused on evidence, not on a type list" begin
        # A total amount and an area over a FIXED area are constants: they do
        # not follow anything, and the probe finds that by scaling the host.
        @test_throws ArgumentError ChemistryLab.sites_per_host(
            mk(TotalSiteAmount(5.0e-6), coupled(BETSurfaceArea(90.0))), M,
        )
        @test_throws ArgumentError ChemistryLab.sites_per_host(
            mk(AreaSiteDensity(1.0e-5), coupled(FixedSurfaceArea(600.0))), M,
        )
        # A shrinking core is linear exactly at p = 1 and nowhere else, which is
        # the case the probe exists for: the nonlinearity is invisible if `n₀`
        # is scaled along with `n`, because the ratio is then always one.
        @test ChemistryLab.sites_per_host(
            mk(AreaSiteDensity(1.0e-5), coupled(ShrinkingCoreArea(BETSurfaceArea(90.0); exponent = 1))), M,
        ) ≈ 1.0e-5 * 90.0 * M rtol = 1.0e-6
        @test_throws ArgumentError ChemistryLab.sites_per_host(
            mk(AreaSiteDensity(1.0e-5), coupled(ShrinkingCoreArea(BETSurfaceArea(90.0); exponent = 2 // 3))), M,
        )
        # A capacity of zero passes proportionality emptily, and is refused for
        # that reason rather than admitted as a family with no sites.
        @test_throws ArgumentError ChemistryLab.sites_per_host(
            mk(MassSiteDensity(0.0), coupled(BETSurfaceArea(90.0))), M,
        )
    end
end

@testsection "the coupling row, built from the declaration" begin

    # `site_coupling_rows` states `Σ dₖ nₖ − ν n_host = 0` for a coupled family.
    # It is built and checked here; it is deliberately NOT imposed by the solver
    # yet, for the reason its docstring measures — `SM.A` already carries a site
    # row pinning the same total, and two of them forbid the host to move.

    h2o = _aq("H2O@", SC_AQSOLVENT)
    hp = _aq("H+")
    ca = _aq("Ca+2")
    host = Species(
        "Ca(OH)2"; symbol = "Portlandite",
        aggregate_state = AS_CRYSTAL, class = SC_COMPONENT,
    )
    free, occ = _surf("XsOH"), _surf("XsOH2+")

    function build(coupling)
        sup = SurfaceSupport("sorbent", "Portlandite", BETSurfaceArea(90.0); coupling)
        cap = coupling === SITES_FOLLOW_HOST ?
            AreaSiteDensity(1.0e-5) : TotalSiteAmount(1.0e-3)
        fam = SiteFamily("Xs", free, [occ]; capacity = cap, support = sup)
        cs = ChemicalSystem(
            [h2o, hp, ca, host, free, occ], [h2o, hp, ca, free]; site_families = [fam],
        )
        return cs, fam
    end

    @testset "uncoupled systems get nothing at all" begin
        cs, _ = build(SITES_FIXED)
        rows, labels = site_coupling_rows(cs)
        @test size(rows) == (0, length(cs.species))
        @test isempty(labels)
    end

    @testset "a coupled family gets one row, and it is the obvious one" begin
        cs, _ = build(SITES_FOLLOW_HOST)
        rows, labels = site_coupling_rows(cs)
        @test labels == ["Xs"]
        @test size(rows) == (1, length(cs.species))

        idx(s) = findfirst(==(s), symbol.(cs.species))
        M = ustrip(us"kg/mol", host[:M])
        ν = 1.0e-5 * 90.0 * M                 # Γ · a · M, checkable by hand
        @test rows[1, idx("XsOH")] == 1.0     # denticity, read from the formula
        @test rows[1, idx("XsOH2+")] == 1.0
        @test rows[1, idx("Portlandite")] ≈ -ν
        # Everything else is untouched: the row says nothing about the aqueous
        # species, which is what makes it a site balance and not a mass balance.
        for s in ("H2O@", "H+", "Ca+2")
            @test rows[1, idx(s)] == 0.0
        end

        # And the row is exactly `Σ d n − ν n_host` evaluated on any composition.
        n = [55.5, 1.0e-6, 1.0e-3, 0.1, 5.0e-6, 2.0e-6]
        @test rows[1, :]' * n ≈ 5.0e-6 + 2.0e-6 - ν * 0.1
    end

    @testset "a host that is not in the system is refused by name" begin
        sup = SurfaceSupport("sorbent", "Ghost", BETSurfaceArea(90.0); coupling = SITES_FOLLOW_HOST)
        fam = SiteFamily("Xs", free, [occ]; capacity = AreaSiteDensity(1.0e-5), support = sup)
        cs = ChemicalSystem(
            [h2o, hp, ca, host, free, occ], [h2o, hp, ca, free]; site_families = [fam],
        )
        @test_throws ArgumentError site_coupling_rows(cs)
    end
end
