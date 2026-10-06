# The extension of CEMDATA18 to the products of the alkali-silica reaction.
#
# `cemdata18-asr.json` is BUILT by ChemistryLab on first use
# (src/databases/derived.jl) from the Cemdata18 file and Table 1 of Jin et al.
# (2023). Its promises are asserted on the database `datapath` returns: the
# published rows, nothing of the base overwritten, the log Ksp round trip
# through the Cemdata18 aqueous species, an enthalpy that agrees with the Gibbs
# energy and the entropy, and the solubility products at 80 °C the authors
# refined.

using JSON

@testsection "ASR extension of CEMDATA18" begin
    base = JSON.parsefile(datapath("cemdata18-thermofun.json"); dicttype = Dict{String, Any})
    ext = JSON.parsefile(datapath("cemdata18-asr.json"); dicttype = Dict{String, Any})
    by_symbol(db) = Dict(String(s["symbol"]) => s for s in db["substances"])
    B, E = by_symbol(base), by_symbol(ext)
    t = literature_table("Jin2023", "products")
    corr = literature_table("Jin2024", "corrected_enthalpies")

    # The added entries are the published rows, with their provenance.
    @test setdiff(keys(E), keys(B)) == Set(t.symbol)
    for (i, sym) in enumerate(t.symbol)
        e = E[sym]
        @test only(e["sm_gibbs_energy"]["values"]) ≈ ustrip(u"J/mol", t.dfG[i])
        @test only(e["sm_entropy_abs"]["values"]) ≈ ustrip(u"J/(mol*K)", t.S[i])
        @test only(e["sm_heat_capacity_p"]["values"]) ≈ ustrip(u"J/(mol*K)", t.Cp[i])
        @test only(e["sm_volume"]["values"]) ≈ 0.1 * ustrip(us"cm^3/mol", t.V[i])
        p = e["asr_provenance"]
        @test p["doi"] == literature("Jin2023").source["doi"] && p["corrigendum_doi"] == literature("Jin2024").source["doi"]
        @test p["log_Ksp_298K"] == ustrip(t.log_K[i])
    end
    # Nothing of the base is overwritten.
    gibbs(s) = get(get(s, "sm_gibbs_energy", Dict()), "values", nothing)
    @test all(gibbs(E[k]) == gibbs(B[k]) for k in keys(B))

    # The log Ksp round trip: the published Gibbs energies give back the
    # published solubility products through the Cemdata18 aqueous species.
    check = ChemistryLab.asr_logK_check(datapath("cemdata18-thermofun.json"))
    @test all(abs(c.difference) < 0.05 for c in check)

    # The enthalpy of each entry is the one its Gibbs energy and entropy give
    # with the element entropies of the base; the article's agrees with it to
    # 0.5 kJ/mol, the corrigendum's does not (recorded in the data file).
    S_el = Dict(Symbol(e["symbol"]) => Float64(e["entropy"]["values"][1]) for e in base["elements"])
    for (i, sym) in enumerate(t.symbol)
        s = sum(Float64(n) * S_el[el] for (el, n) in parse_formula(t.formula[i]))
        H = ustrip(u"kJ/mol", t.dfG[i]) + 298.15 * (ustrip(u"J/(mol*K)", t.S[i]) - s) / 1000
        @test only(E[sym]["sm_enthalpy"]["values"]) / 1000 ≈ H rtol = 1.0e-12
        @test abs(ustrip(u"kJ/mol", t.dfH[i]) - H) < 0.5
        @test abs(ustrip(u"kJ/mol", corr.dfH[findfirst(==(sym), corr.symbol)]) - H) > 500
    end

    # The entropy and the heat capacity follow from the volume through the
    # correlations of Eqs. (3)-(4), to the printed digit.
    NA = 6.02214076e23
    for i in eachindex(t.symbol)
        Vm = ustrip(us"cm^3/mol", t.V[i]) * 1.0e21 / NA
        @test 1579 * Vm + 6 ≈ ustrip(u"J/(mol*K)", t.S[i]) atol = 0.15
        @test 1388.1 * Vm + 4.58 ≈ ustrip(u"J/(mol*K)", t.Cp[i]) atol = 0.15
    end

    # At 80 °C, the solubility products the database gives through Cemdata18
    # fall within the uncertainties of the ones the authors refined (Fig. 1).
    sp = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-asr.json"); verbose = false))
    g(x, T) = ustrip(us"J/mol", sp[x][:ΔₐG⁰](T = T * u"K", P = 1.0e5u"Pa"; unit = true))
    for z in ChemistryLab.asr_records()
        logK(T) = -(sum(ν * g(p, T) for (p, ν) in z.products) - g(z.symbol, T)) / (ChemistryLab.R_GAS * T * log(10))
        @test logK(298.15) ≈ z.logKsp atol = 0.05
        name = z.symbol == "K-shlykovite" ? "log_K_80C_K_shlykovite_refined" : "log_K_80C_Na_shlykovite_refined"
        q = literature("Jin2023")[name]
        @test abs(logK(353.15) - ustrip(literature_value("Jin2023", name))) <= ustrip(ChemistryLab.uncertainty(q))
    end
end
