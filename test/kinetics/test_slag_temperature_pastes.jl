# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# The slag-limestone cement of Snellings et al. (2022) integrated in time
# (`scripts/snellings2022_pastes.jl`): its three materials built from Table 1,
# whose phase names some tables print as minerals, and two days of the paste at
# 20 °C, what the page `examples/slag_temperature_pastes.md` reports at one and
# two days pinned.

isdefined(@__MODULE__, :sn22p_run) || include(joinpath(pkgdir(ChemistryLab), "scripts", "snellings2022_pastes.jl"))

@testsection "a slag-limestone cement integrated at its curing temperature" begin
    # The Rietveld names of Table 1: alite and belite are the database's C3S
    # and C2S, the two C3A one constituent, aphthitalite entered by its oxides.
    pc = material_template("PC (Snellings 2022)", SN22_DB)
    species(name) = symbol(only(c for c in pc.constituents if c.name == name).species)
    @test (species("Alite"), species("Belite"), species("C3A")) == ("C3S", "C2S", "C3A")
    c3a = only(c for c in pc.constituents if c.name == "C3A")
    @test c3a.mass_fraction ≈ only(c for c in pc.constituents if c.name == "Alite").mass_fraction * (4.2 + 6.5) / 62.2
    aph = only(c for c in pc.constituents if c.name == "Aphthitalite")
    @test aph isa OxideConstituent && sum(values(aph.oxides)) ≈ 1
    @test sum(c.mass_fraction for c in pc.constituents) ≈ 1
    slag = material_template("slag (Snellings 2022)", SN22_DB)
    @test only(c.mass_fraction for c in slag.constituents if c.name == "glass") ≈ 0.955

    setup = sn22p_setup()
    @test haskey(setup.cs.dict_species, "SLAG")
    θ = (τ = 5.339, n = 0.5906, Ea = 67.4, α_max = 0.653)
    @test sort!(collect(keys(sn22p_rates(setup, θ)))) == ["Alite", "Belite", "C3A", "C4AF", "SLAG"]

    run = sn22p_run(setup, 20.0, θ; days = 2)
    @test SciMLBase.successful_retcode(run.sol)
    obs = sn22p_observables(run, [1, 2])
    # All the sulfate of the cement is in ettringite from the first day.
    @test [o.ettringite for o in obs] ≈ [9.85, 9.86] atol = 0.01
    @test [o.bound_water for o in obs] ≈ [14.03, 17.04] atol = 0.02
    @test [o.portlandite for o in obs] ≈ [4.45, 5.34] atol = 0.02
    # The measurements the page compares with, read back.
    @test sn22p_measured("ettringite", "XRD", 5.0)[3] == 14.66
    @test length(sn22p_measured("bound water", "TGA", 40.0)) == 6
end
