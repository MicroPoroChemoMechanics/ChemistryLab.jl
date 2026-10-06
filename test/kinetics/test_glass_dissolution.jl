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
