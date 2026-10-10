# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# The dissolution mechanisms of Palandri and Kharaka (2004), built from the
# tables transcribed in data/literature/PalandriKharaka2004.json: each mechanism
# the report gives, its Arrhenius constant and its catalyst, and the refusals
# that keep a missing parameter from becoming a number.

@testsection "Palandri and Kharaka (2004): the mechanisms of a mineral" begin
    minerals = palandri_kharaka_minerals()
    @test length(minerals) == 82 && allunique(minerals)
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
    # Its saturation terms are read by symbol, as those of `saturation_ratio`:
    # by formula, calcite and aragonite would read one activity.
    @test first.(ChemistryLab._stoich_named(cs, rxn)) == ["Cal", "Ca+2", "CO3-2"]
    # A participant the system does not hold, or one without a standard Gibbs
    # energy, is passed over by `transition_state` and refused by
    # `saturation_ratio`.
    absent = Species("CO2"; symbol = "CO2(g)", aggregate_state = AS_GAS, class = SC_GASFLUID)
    bare = Species("CaO"; symbol = "Lim_bare", aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)
    cs2 = ChemicalSystem(vcat(sp, [bare]), CEMDATA_PRIMARIES)
    r2 = Reaction(
        OrderedDict(cs2["Cal"] => 1.0, cs2["Lim_bare"] => 1.0, absent => 1.0),
        OrderedDict(cs2["Ca+2"] => 1.0, cs2["CO3-2"] => 1.0); symbol = "with two passed over",
    )
    @test first.(ChemistryLab._stoich_named(cs2, r2)) == ["Cal", "Ca+2", "CO3-2"]
    @test_throws ArgumentError saturation_ratio(cs2, r2)

    # The sulfides of Table 35 carry orders in Fe3+ and in dissolved O2, the
    # activities of the report's Eq. (3a); pyrite's acid mechanism is its Eq. (3b).
    py = palandri_kharaka("pyrite")
    @test length(py) == 2
    @test [(c.species, c.n) for c in py[1].catalysts] == [("H+", -0.5), ("Fe+3", 0.5)]
    @test [(c.species, c.n) for c in py[2].catalysts] == [("O2@", 0.5)]
    @test k25(py[2]) ≈ 10.0^-4.55 rtol = 1.0e-12
    py2 = palandri_kharaka("pyrite"; fe3 = "Fe[3+]", o2 = "O2(aq)")
    @test py2[1].catalysts[2].species == "Fe[3+]" && only(py2[2].catalysts).species == "O2(aq)"

    # K-feldspar has the three mechanisms of Table 15, the base one up to the
    # alkaline range of a pore solution.
    kf = palandri_kharaka("K-feldspar")
    @test length(kf) == 3 && only(kf[3].catalysts).n == -0.823
    @test k25.(kf) ≈ 10.0 .^ [-10.06, -12.41, -21.2] rtol = 1.0e-12

    # Printed as -53.9 kJ/mol, the acid activation energy of kyanite is stored as
    # printed and used only as the caller decides; a blank order in H+ (the acid
    # mechanism of fayalite) is not an order of zero.
    @test ustrip(us"J/mol", row("orthosilicate_rates", "kyanite").acid_E) == -53900.0
    @test_throws ArgumentError palandri_kharaka("kyanite")
    ky = palandri_kharaka("kyanite"; assume_Ea = 53.9e3)
    @test ky[1].k(; T = 310.0) / k25(ky[1]) ≈ exp(-53.9e3 / R_GAS * (1 / 310.0 - 1 / 298.15)) rtol = 1.0e-12
    @test ismissing(row("orthosilicate_rates", "fayalite").acid_n_H)
    @test_throws ArgumentError palandri_kharaka("fayalite")
    @test only(palandri_kharaka("fayalite"; mechanisms = (:neutral,))).catalysts == RateModelCatalyst{Float64}[]

    # Quartz dissolving into water, through `transition_state` and a run: with
    # the neutral mechanism of the BET row and one silica species, the amount
    # dissolved follows n(t) = n_eq (1 − e^{−A k t / n_eq}), n_eq the solubility
    # in the kilogram of water.
    qsp = [subs["H2O@"], subs["SiO2@"], subs["Qtz"]]
    csq = ChemicalSystem(qsp, [subs["H2O@"], subs["SiO2@"]])
    G(s) = ustrip(us"J/mol", subs[s][:ΔₐG⁰](T = 298.15u"K", P = 1.0e5u"Pa"; unit = true))
    m_eq = exp(-(G("SiO2@") - G("Qtz")) / (R_GAS * 298.15))      # mol/kg, ideal
    nw = ustrip(us"mol", 1.0u"kg" / subs["H2O@"][:M])
    n_eq = m_eq * 1.0                                              # one kilogram of water
    area = 1.0e4                                                    # m², a test value
    rq = Reaction(OrderedDict(csq["Qtz"] => 1.0), OrderedDict(csq["SiO2@"] => 1.0); symbol = "quartz dissolution")
    law = transition_state(palandri_kharaka("quartz, BET surface area"), csq, rq, FixedSurfaceArea(area))
    iq(sym) = findfirst(s -> symbol(s) == sym, csq.species)
    amounts(w, si, qz) = (n = fill(0.0u"mol", 3); n[iq("H2O@")] = w * u"mol"; n[iq("SiO2@")] = si * u"mol"; n[iq("Qtz")] = qz * u"mol"; n)
    nq0 = amounts(nw, 0.0, 1.0)
    kq = KineticsProblem(csq, [KineticReaction(csq, rq, law)], ChemicalState(csq, nq0), (0.0, 3.0e6); equilibrium_solver = nothing)
    tq = [0.0, 3.0e5, 1.0e6, 3.0e6]
    solq = integrate(kq, KineticsSolver(; ode_solver = Rodas5P(), reltol = 1.0e-10, abstol = 1.0e-20, saveat = tq))
    kqz = 10.0^row("silica_rates", "quartz, BET surface area").neutral_log_k
    dissolved = [1.0 - u[1] for u in solq.u]
    @test dissolved ≈ [n_eq * (1 - exp(-area * kqz * t / n_eq)) for t in tq] rtol = 1.0e-6 atol = 1.0e-14
    @test area * kqz * tq[end] / n_eq > 1                           # the run reaches the saturation

    # The saturation ratio of the reaction, as a rate law reads it: one at the
    # solubility, and the value of the method over the whole system.
    sat = saturation_ratio(csq, rq)
    st_eq = ChemicalState(csq, amounts(nw, n_eq, 1.0))
    @test sat(298.15, 1.0e5, log_activities(st_eq, DiluteSolutionModel())) ≈ 1 rtol = 1.0e-10
    st_half = ChemicalState(csq, amounts(nw, 0.5n_eq, 1.0))
    la = log_activities(st_half, DiluteSolutionModel())
    ν = [symbol(s) == "SiO2@" ? 1.0 : symbol(s) == "Qtz" ? -1.0 : 0.0 for s in csq.species]
    g = [ustrip(us"J/mol", s[:ΔₐG⁰](T = 298.15u"K", P = 1.0e5u"Pa"; unit = true)) / (R_GAS * 298.15) for s in csq.species]
    @test sat(298.15, 1.0e5, la) ≈ saturation_ratio(ν, [la[symbol(s)] for s in csq.species], g) rtol = 1.0e-12
    @test sat(298.15, 1.0e5, la) ≈ 0.5 rtol = 1.0e-10
    @test_throws ArgumentError saturation_ratio(csq, Reaction(OrderedDict(csq["Qtz"] => 1.0), OrderedDict(subs["Portlandite"] => 1.0); symbol = "not in the system"))
end
