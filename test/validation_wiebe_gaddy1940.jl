# Carbon dioxide in water under pressure, Wiebe and Gaddy (1940), which the page
# `tutorials/validation_co2_solubility.md` computes with the gas ideal and with
# the equation of state of Peng and Robinson. What is asserted: the
# transcription of Table I gives, converted, the mole fractions Spycher et al.
# (2003) tabulate from it; every point certifies; the four terms of ln m add up
# to the equilibrium; and the numbers the page states are pinned.

isdefined(@__MODULE__, :wg40_compute) || include(joinpath(pkgdir(ChemistryLab), "scripts", "wiebe_gaddy1940_co2.jl"))

@testsection "validation: CO2 solubility under pressure (Wiebe and Gaddy 1940)" begin
    d = wg40_measured()
    @test length(d.T) == 42

    # Converted as Spycher et al. did, point by point: the same temperatures,
    # the pressures from atmospheres to bar, and the mole fractions to their
    # last printed digit, 5e-4 %, within the molar mass of water: 18.015 g/mol
    # from the atomic masses of the library, 18.0153 for them (55.508 mol/kg),
    # 2e-5 of the fraction.
    sp = literature_table("Spycher2003", "appendix_a_wiebe_gaddy1940")
    @test length(sp.temperature_C) == 42
    per_kg = 1 / ustrip(us"kg/mol", WG40_DB["H2O@"][:M])
    for i in eachindex(d.T)
        @test d.T_C[i] == sp.temperature_C[i]
        @test d.P[i] ≈ ustrip(us"Pa", sp.pressure[i]) rtol = 2.0e-3   # bar to one decimal
        @test 100 * d.m[i] / (d.m[i] + per_kg) ≈ sp.x_CO2_percent[i] atol = 6.0e-4
    end

    r = wg40_compute()
    dev_real = 100 .* (r.real ./ r.m .- 1)
    dev_ideal = 100 .* (r.ideal ./ r.m .- 1)

    # The four terms add up to the equilibrium molality, short of the
    # bicarbonate, which the equilibrium counts and they do not.
    for i in eachindex(r.T)
        t = r.terms[i]
        δ = t.henry + t.pressure + t.fugacity + t.volume - log(r.real[i])
        @test -1.0e-3 < δ < 0
    end

    # The figures of the page.
    @test round(sum(dev_real) / 42; digits = 1) == -2.1
    @test round(sum(abs, dev_real) / 42; digits = 1) == 3.0
    @test round(maximum(abs, dev_real); digits = 1) == 7.3
    i_worst = argmax(abs.(dev_real))
    @test (r.T_C[i_worst], round(Int, r.P[i_worst] / 101325)) == (40.0, 200)
    @test round(Int, minimum(dev_ideal)) == 12
    @test round(Int, maximum(dev_ideal)) == 464
    i12 = findall(==(12.0), r.T_C)
    @test all(1.5 < x < 4 for x in dev_real[i12])
    above = [i for i in eachindex(r.T) if r.T_C[i] > 31 && r.P[i] > 40 * 101325]
    @test all(-7.5 < x < -2.8 for x in dev_real[above])
    i500 = findfirst(i -> r.T_C[i] == 40.0 && round(Int, r.P[i] / 101325) == 500, eachindex(r.T))
    @test round(r.terms[i500].fugacity; digits = 2) == -1.41
    @test round(r.terms[i500].volume; digits = 2) == -0.65
end
