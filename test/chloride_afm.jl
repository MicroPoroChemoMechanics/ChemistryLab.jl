# Chloride and the AFm phases (docs/src/examples/chloride_afm.md): the mixture
# of Balonis et al. (2010) titrated by CaCl2. The numbers the page states are
# pinned.

isdefined(@__MODULE__, :ba10_titration) || include(joinpath(pkgdir(ChemistryLab), "scripts", "balonis2010_chloride.jl"))

@testsection "durability: chloride and the AFm phases" begin
    ratios = 0.0:0.02:1.1
    free = ba10_titration(ratios)
    carb = ba10_titration(ratios; calcite = true)
    @test all(r -> r.certified, free) && all(r -> r.certified, carb)
    b = ba10_boundaries(free)
    @test b.monosulfate_gone ≈ 0.44
    @test b.kuzel_gone ≈ 0.78
    @test isnan(ba10_boundaries(carb).kuzel_first)
    row(rows, r) = rows[argmin([abs(x.ratio - r) for x in rows])]
    @test ba10_phases(row(free, 0.0)) == "E, Ms, P"
    @test ba10_phases(row(free, 0.3)) == "E, Ms, Ks, P"
    @test ba10_phases(row(free, 0.5)) == "E, Ks, P"
    @test ba10_phases(row(free, 1.0)) == "E, Fs, P"
    @test ba10_phases(row(carb, 0.0)) == "E, Mc, P, Cc"
    @test ba10_phases(row(carb, 1.0)) == "E, Fs, P, Cc"
    x = row(free, 0.3)
    @test x.Ca ≈ 24.6 atol = 0.05
    @test x.Cl ≈ 10.3 atol = 0.05
    @test x.pH ≈ 12.44 atol = 0.005
    @test row(free, 0.5).Cl / 18.3 ≈ 1.9 atol = 0.05
    @test row(free, 1.0).Cl / 143.0 ≈ 0.82 atol = 0.01
    @test row(free, 1.0).solids / row(free, 0.0).solids ≈ 1.26 atol = 0.005
    @test row(carb, 1.0).solids / row(carb, 0.0).solids ≈ 1.059 atol = 0.005
end
