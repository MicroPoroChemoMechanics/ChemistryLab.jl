# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

# Inhibitors of a dissolution mechanism, and the dissolution of the glasses of
# supplementary cementitious materials.

@testsection "an inhibitor divides a mechanism by (1 + K a)^m" begin
    inh = RateModelInhibitor("Ca+2", 2000.0)
    @test inh.species == "Ca+2" && inh.K === 2000.0 && inh.m === 1.0
    @test RateModelInhibitor("Ca+2", 1; m = 2).m === 2.0
    @test_throws ArgumentError RateModelInhibitor("Ca+2", -1.0)
    @test_throws ArgumentError RateModelInhibitor("Ca+2", 1.0; m = -1)

    # The four-argument constructors, typed or not, still build a mechanism
    # without inhibitors, and the five-argument one promotes its numbers.
    k = arrhenius_rate_constant(1.0e-8, 40.0e3)
    @test isempty(RateMechanism(k, 1.0, 1.0, RateModelCatalyst{Float64}[]).inhibitors)
    @test isempty(RateMechanism{typeof(k), Float64}(k, 1, 1, RateModelCatalyst{Float64}[]).inhibitors)
    m = RateMechanism(k, 1, 1, [RateModelCatalyst("H+", 0.5)], [RateModelInhibitor("Ca+2", 2)])
    @test m.p === 1.0 && only(m.catalysts).n === 0.5 && only(m.inhibitors).K === 2.0

    # On calcite in water, far from equilibrium.
    subs = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-thermofun.json"); verbose = false))
    sp = speciation(collect(values(subs)), ["Cal", "H2O@"]; aggregate_state = [AS_AQUEOUS], exclude_species = split("H2@ O2@ CH4@"))
    cs = ChemicalSystem(sp, CEMDATA_PRIMARIES)
    rxn = Reaction(OrderedDict(cs["Cal"] => 1.0), OrderedDict(cs["Ca+2"] => 1.0, cs["CO3-2"] => 1.0); symbol = "calcite dissolution")
    # Activities and amounts are read by formula and by symbol, as in a run.
    index = Dict{String, Int}()
    for (i, s) in enumerate(cs.species)
        index[phreeqc(formula(s))] = i
        index[symbol(s)] = i
    end

    # A rate law built on a fresh system, before any of its species was asked
    # for its Gibbs energy, has the saturation ratio of the whole reaction: it
    # vanishes where the ion activity product is the solubility product.
    tst = transition_state([RateMechanism(k, 1.0, 1.0)], cs, rxn, BETSurfaceArea(90.0))
    g(s) = cs[s][:ΔₐG⁰](; T = 298.15, unit = false)
    lnK = -(g("Ca+2") + g("CO3-2") - g("Cal")) / (ChemistryLab.R_GAS * 298.15)
    leq = zeros(length(cs.species))
    leq[index["Ca+2"]] = leq[index["CO3-2"]] = lnK / 2
    lna = fill(-7.0, length(cs.species))
    lna[index["Cal"]] = 0.0
    lna[index["Ca+2"]] = log(3.0e-4)
    lna[index["CO3-2"]] = log(1.0e-8)                       # undersaturated
    n = fill(1.0e-3, length(cs.species))
    nv, lv = StateView(n, index), StateView(lna, index)
    far = tst(298.15, 1.0e5, 0.0, nv, lv, nv)
    @test far > 0
    @test abs(tst(298.15, 1.0e5, 0.0, nv, StateView(leq, index), nv)) < 1.0e-12 * far
    rate(mech) = transition_state([mech], cs, rxn, BETSurfaceArea(90.0))(298.15, 1.0e5, 0.0, nv, lv, nv)

    plain = rate(RateMechanism(k, 1.0, 1.0))
    @test plain > 0
    # Without the species, or with K = 0, the rate is the one without the
    # inhibitor, to the last bit.
    @test rate(RateMechanism(k, 1.0, 1.0, RateModelCatalyst{Float64}[], [RateModelInhibitor("Ca+2", 0.0)])) === plain
    # Otherwise it is divided by (1 + K a)^m, a the activity of Ca²⁺.
    a = 3.0e-4
    for (K, mo) in ((2000.0, 1.0), (500.0, 2.5))
        r = rate(RateMechanism(k, 1.0, 1.0, RateModelCatalyst{Float64}[], [RateModelInhibitor("Ca+2", K; m = mo)]))
        @test r ≈ plain / (1 + K * a)^mo rtol = 1.0e-14
    end
    # Its derivative with respect to K, through the rate law, is the one of the
    # factor: −m a (1 + K a)^(−m−1) times the plain rate.
    K0 = 2000.0
    dr = ForwardDiff.derivative(K -> rate(RateMechanism(k, 1.0, 1.0, RateModelCatalyst{Float64}[], [RateModelInhibitor("Ca+2", K)])), K0)
    @test dr ≈ -a * (1 + K0 * a)^(-2) * plain rtol = 1.0e-12

    # An inhibitor the system does not hold is refused by name.
    err = try
        transition_state([RateMechanism(k, 1.0, 1.0, RateModelCatalyst{Float64}[], [RateModelInhibitor("Al(OH)4-", 1.0)])], cs, rxn, BETSurfaceArea(90.0))
        nothing
    catch e
        e
    end
    @test err isa ArgumentError && occursin("Al(OH)4-", sprint(showerror, err))
end

@testsection "the dissolution of a glass at pH 13 after Snellings (2013)" begin
    key = "Snellings2013"
    glasses = literature_table(key, "glasses")
    markers = literature_table(key, "fig8_markers")
    rates = literature_table(key, "initial_rates")
    a = ustrip(literature_value(key, "fig8_slope"))
    b = ustrip(literature_value(key, "fig8_intercept"))
    M(ox) = ustrip(us"g/mol", Species(ox)[:M])
    analysis(i) = Dict(
        "CaO" => ustrip(glasses.CaO_percent[i]) / 100, "Al2O3" => ustrip(glasses.Al2O3_percent[i]) / 100,
        "SiO2" => ustrip(glasses.SiO2_percent[i]) / 100,
    )

    # One mole of the glass is one mole of its cations.
    g1 = analysis(1)
    Mc = cation_molar_mass(g1)
    @test ustrip(us"g/mol", Mc) ≈ 1 / (g1["CaO"] / M("CaO") + 2 * g1["Al2O3"] / M("Al2O3") + g1["SiO2"] / M("SiO2")) rtol = 1.0e-14
    sp = glass_species(g1; symbol = "G1", M = Mc)
    @test sum(v for (el, v) in atoms(sp) if el != :O) ≈ 1 rtol = 1.0e-14
    @test ustrip(us"g/mol", cation_molar_mass(merge(g1, Dict("H2O" => 0.05)))) ≈ ustrip(us"g/mol", Mc) rtol = 1.0e-14
    @test_throws ArgumentError cation_molar_mass(Dict("H2O" => 1.0))

    for i in eachindex(glasses.glass)
        g = glasses.glass[i]
        mech = snellings2013_glass(analysis(i); Ea = 60.0e3)
        logk = log10(mech.k(; T = 293.15))
        # The abscissa is the one read on the figure, to its reading.
        x = ChemistryLab._snellings2013_abscissa(analysis(i), ("CaO",))
        @test x ≈ ustrip(markers.abscissa[findfirst(==(g), markers.glass)]) atol = 1.0e-3
        # The constant is the printed regression, exactly, and within 0.1 of the
        # rate the paper measured in NaOH alone (its error is 0.15).
        @test logk ≈ a * x + b rtol = 1.0e-14
        j = findfirst(k -> rates.glass[k] == g && all(iszero ∘ ustrip, (rates.Al_initial[k], rates.Ca_initial[k], rates.Si_initial[k])), eachindex(rates.glass))
        @test abs(logk - ustrip(rates.log_rate[j])) < 0.1
        @test isempty(mech.catalysts) && isempty(mech.inhibitors)
    end

    # The derivative with respect to the lime content is the slope times that of
    # the abscissa, n_Ca over the denominator.
    den = 2 * 2 * g1["Al2O3"] / M("Al2O3") + g1["SiO2"] / M("SiO2")
    d = ForwardDiff.derivative(f -> log10(snellings2013_glass(merge(g1, Dict("CaO" => f)); Ea = 60.0e3).k(; T = 293.15)), g1["CaO"])
    @test d ≈ a / M("CaO") / den rtol = 1.0e-12

    # Arrhenius about 20 °C, and no activation energy taken for granted.
    mech = snellings2013_glass(g1; Ea = 60.0e3)
    @test mech.k(; T = 313.15) ≈ mech.k(; T = 293.15) * exp(-60.0e3 / ChemistryLab.R_GAS * (1 / 313.15 - 1 / 293.15)) rtol = 1.0e-13
    @test_throws UndefKeywordError snellings2013_glass(g1)

    # A slag analysis lies beyond the most calcic glass: refused, unless asked
    # for; a modifier named counts in the abscissa.
    slag = Dict("CaO" => 0.4164, "Al2O3" => 0.1127, "SiO2" => 0.3521, "MgO" => 0.0596)
    @test_throws DomainError snellings2013_glass(slag; Ea = 60.0e3)
    x_ca = ChemistryLab._snellings2013_abscissa(slag, ("CaO",))
    x_camg = ChemistryLab._snellings2013_abscissa(slag, ("CaO", "MgO"))
    @test x_camg - x_ca ≈ slag["MgO"] / M("MgO") / (2 * 2 * slag["Al2O3"] / M("Al2O3") + slag["SiO2"] / M("SiO2")) rtol = 1.0e-12
    @test log10(snellings2013_glass(slag; Ea = 60.0e3, extrapolate = true).k(; T = 293.15)) ≈ a * x_ca + b rtol = 1.0e-14
    @test_throws ArgumentError ChemistryLab._snellings2013_abscissa(Dict("CaO" => 1.0), ("CaO",))
end

isdefined(@__MODULE__, :sn13_rows) || include(joinpath(pkgdir(ChemistryLab), "scripts", "snellings2013_glass.jl"))

@testsection "calcium and aluminum slow the glasses down (fitted on Snellings 2013, Table II)" begin
    rows = sn13_rows()
    # Every solution of Table II at pH 13 certifies, with the NaOH that holds it.
    @test length(rows) == 51 && all(r -> r.certified, rows)
    @test all(r -> 0.08 < r.NaOH < 0.12, rows)
    # The activities follow what was added: none where nothing was.
    @test all(r -> r.ca > 0 || r.a_Ca < 1.0e-20, rows) && all(r -> r.al > 0 || r.a_Al < 1.0e-20, rows)
    fit = sn13_fit(rows)
    # Calcium: one factor on all six glasses, within the error of the paper.
    @test fit.ca.n == 12 && fit.ca.rms < ustrip(literature_value("Snellings2013", "log_rate_uncertainty"))
    @test fit.ca.K ≈ 8617.6 rtol = 1.0e-4
    # Aluminum on the tectosilicate glasses: one factor cannot follow G3 and G6
    # both, and the residual says so.
    @test fit.al.n == 12 && 0.2 < fit.al.rms < 0.3
    @test fit.al.K ≈ 1310.5 rtol = 1.0e-4
    # The fit is a minimum: its derivative in ln K vanishes.
    d(r) = r.log_rate - r.log_base
    sse(u) = sum((d(r) + log10(1 + exp(u) * r.a_Ca))^2 for r in rows if r.ca > 0)
    @test abs(ForwardDiff.derivative(sse, log(fit.ca.K))) < 1.0e-10
    # The inhibitors of a glass: calcium always, aluminum on a tectosilicate one.
    @test [i.species for i in sn13_inhibitors(fit, "G1")] == ["Ca+2"]
    @test [i.species for i in sn13_inhibitors(fit, :tectosilicate)] == ["Ca+2", "AlO2-"]
end
