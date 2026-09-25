# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using ChemistryLab
using DynamicQuantities
using Test

include("reference_species.jl")

@testsection "what a solid assemblage loses on heating" begin

    c18 = Dict(
        symbol(s) => s for s in
            build_species(datapath("cemdata18-thermofun.json"); verbose = false)
    )
    aq = Dict(
        symbol(s) => s for s in
            build_species(datapath("slop98-inorganic-thermofun.json"); verbose = false)
    )

    # Portlandite, gypsum and calcite: one loses a water per formula unit, one
    # loses two, one loses a CO₂ and no water. Three arithmetics, one system.
    solids = [c18[k] for k in ("Portlandite", "Gp", "Cal")]
    aqueous = [aq[k] for k in ("H2O@", "H+", "Ca+2", "SO4-2", "CO3-2")]
    cs = ChemicalSystem(
        vcat(aqueous, solids),
        [aq["H2O@"], aq["H+"], aq["Ca+2"], aq["SO4-2"], aq["CO3-2"]],
    )
    idx = Dict(symbol(sp) => k for (k, sp) in enumerate(cs.species))

    n = Any[fill(1.0e-12u"mol", length(cs.species))...]
    n[idx["H2O@"]] = moles_of_water() * u"mol"
    n[idx["Portlandite"]] = 1.0u"mol"
    n[idx["Gp"]] = 2.0u"mol"
    n[idx["Cal"]] = 3.0u"mol"
    st = ChemicalState(cs, n)

    M_H2O = ustrip(us"kg/mol", aq["H2O@"][:M])

    @testset "hydrogen leaves as water, whatever it is written as" begin
        loss = ignition_loss(st)
        # Ca(OH)₂ has no H₂O in its formula and loses one: 1 mol × H/2 = 1 mol.
        # CaSO₄·2H₂O loses two, so 2 mol × 2 = 4 mol. Total 5 mol of water.
        @test ustrip(us"kg", loss.water) ≈ 5.0 * M_H2O rtol = 1.0e-10
        @test ustrip(us"kg", bound_water(st)) ≈ ustrip(us"kg", loss.water)

        # A rule that looked for `H2O` in the formula would report 4, not 5 —
        # which is why this counts hydrogen.
        @test get(atoms(c18["Portlandite"]), :H, 0) == 2
    end

    @testset "carbon leaves as carbon dioxide" begin
        loss = ignition_loss(st)
        # The system declares no CO₂, so the gas is weighed from its formula,
        # with the library's atomic masses.
        M_CO2 = ustrip(us"kg/mol", calculate_molar_mass(Dict(:C => 1, :O => 2)))
        @test ustrip(us"kg", loss.carbon_dioxide) ≈ 3.0 * M_CO2 rtol = 1.0e-10
        @test ustrip(us"kg", loss.total) ≈
            ustrip(us"kg", loss.water) + ustrip(us"kg", loss.carbon_dioxide)
    end

    @testset "the pore solution is not bound water" begin
        # Fifty-five moles of solvent are in the state and none of it counts:
        # evaporable water is not what a hydrate bound, and a rule that summed
        # all the hydrogen would be reporting the beaker.
        loss = ignition_loss(st)
        @test ustrip(us"kg", loss.water) < 0.1          # 5 mol, not 55
        dry = ChemicalState(cs, Any[fill(1.0e-12u"mol", length(cs.species))...])
        @test ustrip(us"kg", ignition_loss(dry).total) < 1.0e-12
    end

    @testset "per phase, which is where a thermogram would attach" begin
        per = bound_water_per_phase(st)
        names = first.(per)
        @test "Gp" in names && "Portlandite" in names
        @test !("Cal" in names)                          # calcite has no hydrogen
        # Largest first: gypsum carries 4 mol of water, portlandite 1.
        @test first(names) == "Gp"
        @test ustrip(us"kg", per[1].second) ≈ 4.0 * M_H2O rtol = 1.0e-10
        @test sum(ustrip(us"kg", p.second) for p in per) ≈
            ustrip(us"kg", bound_water(st)) rtol = 1.0e-10
    end
end
