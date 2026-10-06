# Hemicarbonate and monocarbonate (docs/src/examples/hemicarbonate.md): the
# carbonate series of Georget et al. (2022) against their own calculation and
# their samples. The numbers the page states are pinned.

isdefined(@__MODULE__, :ge22_series) || include(joinpath(pkgdir(ChemistryLab), "scripts", "georget2022_hemicarbonate.jl"))

@testsection "durability: hemicarbonate and monocarbonate" begin
    rows = ge22_series()
    t = ge22_steps()
    @test length(rows) == 101
    @test all(r -> r.certified, rows)
    g(f, k) = ustrip(us"g", getproperty(t, f)[k])
    gap = maximum(maximum(abs(getproperty(rows[k], f) - g(f, k)) for f in (:katoite, :hemicarbonate, :monocarbonate, :calcite, :portlandite)) for k in eachindex(rows))
    @test gap < 0.03
    b = ge22_breakpoints(rows)
    @test b.katoite_gone ≈ 0.4515 atol = 1.0e-4
    @test b.hemicarbonate_gone ≈ 0.8565 atol = 1.0e-4
    @test ge22_phases(rows, 0.25) == ["Hc", "katoite", "CH"]
    @test ge22_phases(rows, 0.5) == ["Hc", "Mc", "CH"]
    @test ge22_phases(rows, 0.75) == ["Hc", "Mc", "CH"]
    @test ge22_phases(rows, 1.0) == ["Mc", "CH", "CC"]
    @test rows[1].pH ≈ 12.66 atol = 0.005
end
