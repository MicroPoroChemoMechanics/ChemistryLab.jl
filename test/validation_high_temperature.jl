# Water and quartz from 0 to 1000 °C (docs/src/tutorials/validation_high_temperature.md):
# the ionization constant of water against the formulation of Bandura and Lvov
# (2006), the solubility of quartz against the equation of Manning (1994), with
# the standard states of slop98. The numbers the page states are pinned.

isdefined(@__MODULE__, :hw_pkw) || include(joinpath(pkgdir(ChemistryLab), "scripts", "water_high_temperature.jl"))

@testsection "validation: water and quartz from 0 to 1000 °C" begin
    r = hw_bandura_lvov()
    Δ = r.ours .- r.ref
    @test length(Δ) == 175
    largest(sel) = maximum(abs, Δ[sel])
    mean_at(T_C) = (i = findall(==(T_C), r.T_C); sum(Δ[i]) / length(i))
    @test largest(25 .<= r.T_C .<= 400) ≈ 0.203 atol = 5.0e-4
    @test largest(r.T_C .== 0) ≈ 0.298 atol = 5.0e-4
    @test mean_at(600.0) ≈ 0.377 atol = 5.0e-4
    @test mean_at(800.0) ≈ 0.826 atol = 5.0e-4
    @test mean_at(1000.0) ≈ 1.364 atol = 5.0e-4
    @test hw_density(1273.15, 2.5e8) ≈ 0.397 atol = 5.0e-4
    @test hw_density(1273.15, 5.0e8) ≈ 0.609 atol = 5.0e-4

    gap(T_C, P_bar) = hw_quartz_log_m(T_C + 273.15, P_bar * 1.0e5) - hw_manning(T_C + 273.15, P_bar * 1.0e5)
    low = [abs(gap(T, P)) for (T, P) in ((25, 1), (100, 1), (200, 500), (300, 1000), (400, 2000))]
    @test minimum(low) ≈ 0.004 atol = 5.0e-4
    @test maximum(low) ≈ 0.022 atol = 5.0e-4
    @test abs(gap(500, 2000)) ≈ 0.058 atol = 5.0e-4
    @test gap(800, 5000) ≈ -0.12 atol = 5.0e-4
    # The α-β transition of quartz, at the top of its first interval.
    @test temperature_range(HW_DB["Qtz"])[1] == 273.15
end
