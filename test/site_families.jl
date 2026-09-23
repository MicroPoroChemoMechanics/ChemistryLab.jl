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
