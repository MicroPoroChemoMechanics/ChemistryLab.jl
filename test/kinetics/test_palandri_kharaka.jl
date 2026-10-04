# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# The dissolution mechanisms of Palandri and Kharaka (2004), built from the
# tables transcribed in data/literature/PalandriKharaka2004.json: each mechanism
# the report gives, its Arrhenius constant and its catalyst, and the refusals
# that keep a missing parameter from becoming a number.

@testsection "Palandri and Kharaka (2004): the mechanisms of a mineral" begin
    minerals = palandri_kharaka_minerals()
    @test length(minerals) == 34 && allunique(minerals)
    @test issubset(["calcite", "quartz, BET surface area", "gibbsite", "anhydrite", "gypsum", "wollastonite", "hematite"], minerals)

    k25(m) = m.k(; T = 298.15)
    row(t, m) = literature_row("PalandriKharaka2004", t, m)

    # Gibbsite has the three mechanisms of Table 32, with their orders in H+,
    # negative for the base mechanism, and their rate constants at 25 °C.
    g = palandri_kharaka("gibbsite")
    r = row("hydroxide_rates", "gibbsite")
    @test length(g) == 3
    @test [only(m.catalysts).n for m in g[[1, 3]]] == [0.992, -0.784]
    @test all(only(m.catalysts).species == "H+" for m in g[[1, 3]]) && isempty(g[2].catalysts)
    @test k25.(g) ≈ 10.0 .^ [r.acid_log_k, r.neutral_log_k, r.base_log_k] rtol = 1.0e-12
    # The activation energy is the table's: the Arrhenius ratio between two
    # temperatures is exactly its own.
    E = ustrip(us"J/mol", r.base_E)
    @test g[3].k(; T = 310.0) / k25(g[3]) ≈ exp(-E / R_GAS * (1 / 310.0 - 1 / 298.15)) rtol = 1.0e-12
    @test all(m -> m.p == 1 && m.q == 1, g)

    # A dash in the table is no mechanism: anhydrite dissolves by the neutral one
    # alone, and quartz in the row of Knauss and Wolery by the base one alone.
    @test length(palandri_kharaka("anhydrite")) == 1 && isempty(only(palandri_kharaka("anhydrite")).catalysts)
    qb = only(palandri_kharaka("quartz, base mechanism (Knauss and Wolery 1988)"))
    @test only(qb.catalysts).n == -0.5
    @test_throws ArgumentError palandri_kharaka("anhydrite"; mechanisms = (:acid,))

    # Where the table gives the pre-exponential factor, log k is its value at
    # 25 °C to the rounding of the table. The base row of quartz is the
    # exception the notes record: its log k is that of the factor 491 of Knauss
    # and Wolery, not of the 10 the table prints.
    lk(A, E) = log10(A) - E / (log(10) * R_GAS * 298.15)
    s = literature_table("PalandriKharaka2004", "silica_rates")
    for k in eachindex(s.mineral)
        ismissing(s.neutral_A[k]) && continue
        @test abs(lk(ustrip(s.neutral_A[k]), ustrip(us"J/mol", s.neutral_E[k])) - s.neutral_log_k[k]) <= 0.011
    end
    q = row("silica_rates", "quartz, base mechanism (Knauss and Wolery 1988)")
    @test lk(491, ustrip(us"J/mol", q.base_E)) ≈ q.base_log_k atol = 0.005
    @test !(lk(ustrip(q.base_A), ustrip(us"J/mol", q.base_E)) ≈ q.base_log_k)

    # The neutral mechanism of gypsum has no activation energy in the report:
    # refused, unless the caller states one.
    @test ismissing(row("sulfate_rates", "gypsum").neutral_E)
    @test_throws ArgumentError palandri_kharaka("gypsum")
    gy = only(palandri_kharaka("gypsum"; assume_Ea = 0.0))
    @test k25(gy) ≈ 10.0^-2.79 rtol = 1.0e-12
    @test gy.k(; T = 310.0) ≈ k25(gy) rtol = 1.0e-12

    # The carbonate mechanism reads P(CO2): asked for, it needs the species.
    @test length(palandri_kharaka("calcite")) == 2
    @test_throws ArgumentError palandri_kharaka("calcite"; mechanisms = (:acid, :neutral, :carbonate))
    c3 = palandri_kharaka("calcite"; mechanisms = (:acid, :neutral, :carbonate), pco2 = "CO2")
    @test only(c3[3].catalysts).species == "CO2" && only(c3[3].catalysts).n == 1.0
    @test_throws ArgumentError palandri_kharaka("not a mineral")
    @test_throws ArgumentError palandri_kharaka("calcite"; mechanisms = (:alkaline,))

    # The mechanisms on a reaction: calcite in water. Its carbonate mechanism
    # reads a gas the system does not hold, which `transition_state` refuses by
    # name rather than reading its activity as one.
    subs = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-thermofun.json"); verbose = false))
    sp = speciation(collect(values(subs)), ["Cal", "H2O@"]; aggregate_state = [AS_AQUEOUS], exclude_species = split("H2@ O2@ CH4@"))
    cs = ChemicalSystem(sp, CEMDATA_PRIMARIES)
    @test !haskey(cs.dict_species, "CO2")
    rxn = Reaction(OrderedDict(cs["Cal"] => 1.0), OrderedDict(cs["Ca+2"] => 1.0, cs["CO3-2"] => 1.0); symbol = "calcite dissolution")
    err = try
        transition_state(c3, cs, rxn, BETSurfaceArea(90.0))
        nothing
    catch e
        e
    end
    @test err isa ArgumentError && occursin("CO2", sprint(showerror, err))
    # Without it, the law far from equilibrium is the sum of the mechanisms, the
    # acid one through the activity of H+.
    tst = transition_state(palandri_kharaka("calcite"), cs, rxn, BETSurfaceArea(90.0))
    @test tst isa KineticFunc
end
