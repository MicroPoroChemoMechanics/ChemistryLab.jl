# Seawater (docs/src/examples/seawater.md): the paste of De Weerdt et al.
# (2014) titrated by their seawater. The numbers the page states are pinned.

isdefined(@__MODULE__, :sw_titration) || include(joinpath(pkgdir(ChemistryLab), "scripts", "seawater.jl"))

@testsection "durability: seawater" begin
    g = sw_charge_gap()
    @test g.anions - g.cations ≈ 0.04 atol = 5.0e-4
    @test g.Na_closed ≈ 0.518 atol = 5.0e-4

    volumes = 10.0 .^ range(0, 4; length = 41)
    rows = sw_titration(volumes)
    @test all(r -> r.certified, rows)
    tr = sw_transitions(rows)
    near(x, y) = isapprox(x, y; rtol = 1.0e-6)
    @test near(tr[:friedel].first, 10^1.5)
    @test near(tr[:monocarbonate].last, 100.0)
    @test near(tr[:friedel].last, 10^2.7)
    @test near(tr[:thaumasite].first, 10^2.8)
    @test near(tr[:brucite].first, 10^3.2)
    @test near(tr[:portlandite].last, 10^3.6)
    @test near(tr[:calcite].last, 10^3.4)
    @test isnan(tr[:kuzel].first)
    @test tr[:hydrotalcite].first == 1.0
    # The order of the sequence, which the core shows from the surface inward.
    @test tr[:friedel].first < tr[:thaumasite].first < tr[:brucite].first < tr[:portlandite].last
end
