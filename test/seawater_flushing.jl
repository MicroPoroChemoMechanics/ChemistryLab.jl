# Seawater flushed through a ground paste (docs/src/examples/seawater_flushing.md):
# the paste of De Weerdt and Justnes (2015) renewed two thousand times by their
# seawater. The numbers the page states are pinned.

isdefined(@__MODULE__, :dj_path) || include(joinpath(pkgdir(ChemistryLab), "scripts", "seawater_flushing.jl"))

@testsection "durability: seawater flushed through a paste" begin
    s = dj_sample()
    @test s.cement ≈ 25.62 atol = 0.01
    @test s.cement_moist ≈ 20.5 atol = 0.01
    @test s.filling ≈ 195.2 atol = 0.1
    @test s.total ≈ 3.903e5 rtol = 1.0e-3
    @test s.retained ≈ 74.16 atol = 0.01
    g = dj_charge_gap()
    @test g.anions - g.cations ≈ 0.0707 atol = 5.0e-4
    @test g.Na_closed ≈ 0.484 atol = 5.0e-4

    # The budget of the seawater is neutral once closed on sodium.
    cs = sw_system()
    @test abs(ChemistryLab._budget_charge(cs, dj_seawater(cs))) < 1.0e-12

    path = dj_path(2000)
    rows = path.rows
    @test all(r -> r.certified, rows)
    m = dj_measured()
    @test m.retained ≈ 0.196 atol = 5.0e-4
    # Before the exposure the ratios are those of the cement.
    for e in (:Mg, :S, :Na, :K)
        @test path.initial[e] / path.initial[:Ca] ≈ m.original[e] rtol = 0.05
    end
    @test path.initial[:Al] / path.initial[:Ca] ≈ m.original[:Al] rtol = 0.06

    # The calcium left by the end of the exposure, and the point where it is
    # the measured one.
    @test rows[end].retained < 0.01
    x = dj_match(rows, m.retained)
    @test x.V ≈ 15511 rtol = 2.0e-3
    @test x.ratios[:Mg] ≈ 2.393 rtol = 2.0e-3
    @test x.ratios[:S] ≈ 0.494 rtol = 2.0e-3
    @test x.ratios[:Al] ≈ 0.301 rtol = 2.0e-3
    @test x.ratios[:Cl] < 1.0e-6
    @test x.ratios_pw[:Cl] ≈ 0.18 rtol = 5.0e-3
    @test x.ratios_pw[:Na] ≈ 0.103 rtol = 5.0e-3
    @test !haskey(x.after, "Portlandite") && !haskey(x.after, "ettringite")
    @test all(k -> haskey(x.after, k), ("M15SH", "Brc", "Gp"))

    # The sequence of the page.
    V = [r.V for r in rows]
    first_with(sp) = V[findfirst(r -> get(r.solids, sp, 0.0) > 1.0e-6, rows)]
    last_with(sp) = V[findlast(r -> get(r.solids, sp, 0.0) > 1.0e-6, rows)]
    @test last_with("C4AClH10") ≈ 585 rtol = 2.0e-3
    @test first_with("Brc") ≈ 1366 rtol = 2.0e-3
    @test last_with("Portlandite") ≈ 4489 rtol = 2.0e-3
    @test first_with("Gp") ≈ 6440 rtol = 2.0e-3
    @test first_with("thaumasite") ≈ 781 rtol = 2.0e-3
    @test first_with("M15SH") ≈ 9368 rtol = 2.0e-3

    # The solids take up the chloride first, the sulfur next, the magnesium
    # last: the volume at which each holds the most.
    peak(e) = V[argmax([r.elements[e] for r in rows])]
    @test peak(:Cl) ≈ 195 rtol = 2.0e-3
    @test peak(:S) ≈ 5074 rtol = 2.0e-3
    @test peak(:Mg) ≈ 16394 rtol = 2.0e-3
end
