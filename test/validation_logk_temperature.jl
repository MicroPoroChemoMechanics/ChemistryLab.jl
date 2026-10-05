# The equilibrium constants of Cemdata18 from 5 to 90 °C against the analytical
# expressions of phreeqc.dat (`scripts/logk_temperature_check.jl`), which the
# page `tutorials/validation_logk_temperature.md` reads. What is asserted: the
# records are read as PHREEQC writes them, and the numbers the page states are
# pinned.

isdefined(@__MODULE__, :lkt_compare) || include(joinpath(pkgdir(ChemistryLab), "scripts", "logk_temperature_check.jl"))

@testsection "validation: equilibrium constants from 5 to 90 °C against phreeqc.dat" begin
    # A phase, its reaction on the line after its name; a species, headed by
    # its reaction; gypsum, two expressions, the second kept.
    @test lkt_analytic("Calcite") == [17.118, -0.046528, -3496.0, 0.0, 0.0, 0.0]
    @test lkt_analytic("H2O = OH- + H+") == [293.29227, 0.1360833, -10576.913, -123.73158, 0.0, -6.996455e-5]
    @test lkt_analytic("Gypsum") == [93.7, 5.99e-3, -4.0e3, -35.019, 0.0, 0.0]
    @test lkt_phreeqc(lkt_analytic("Calcite"), 298.15) ≈ -8.48 atol = 0.001
    @test_throws ErrorException lkt_analytic("Gibbsite")

    rows = Dict(r.name => r for r in lkt_compare())
    @test length(rows) == 12
    @test count(r -> abs(r.phreeqc25 - r.cemdata25) < 0.01, values(rows)) == 10
    @test rows["quartz"].cemdata25 - rows["quartz"].phreeqc25 ≈ 0.234 atol = 0.001
    @test maximum(r -> abs(r.at50), values(rows)) ≈ 0.063 atol = 0.001
    @test rows["calcite"].worst_T == 90
    @test rows["calcite"].worst ≈ 0.235 atol = 0.001
    @test maximum(r -> abs(r.worst), values(rows)) == abs(rows["calcite"].worst)
end
