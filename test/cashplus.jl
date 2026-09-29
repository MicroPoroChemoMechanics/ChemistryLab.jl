# The CASH+ core model of C-S-H (Kulik, Miron & Lothenbach 2022): the compound
# energy formalism it is written in, its data, the database that carries it, and
# the equilibria the paper reports for the Ca-Si-H2O system at 25 °C.

using JSON, LinearAlgebra

const CASHPLUS = ["TSvh", "TSCh", "Tvvh", "TCvh", "TvCh", "TCCh"]

# n G / RT of a compound-energy model, coded from the Gibbs energy of its
# docstring and nothing else, so that ForwardDiff gives the chemical potentials
# the model must return. With `split = true`, plus the divergence D(x) the solver
# minimizes to choose the split of the members.
_xlogx(v) = v > 0 ? v * log(v) : zero(v)
function _cef_nG(mdl, g, T; split = true)
    RT = ChemistryLab.R_GAS * T
    o = mdl.lattice.occupancy
    return function (n)
        N = sum(n)
        x = n ./ N
        y = site_fractions(mdl, x)
        G = sum(g[j] * prod(y[s][o[s, j]] for s in axes(o, 1)) for j in axes(o, 2))
        G += sum(mdl.lattice.multiplicity[s] * sum(_xlogx, y[s]) for s in eachindex(y))
        for (s, i, l, W) in mdl.interactions
            G += W / RT * y[s][i] * y[s][l]
        end
        split && (G += sum(_xlogx, x) - sum(sum(_xlogx, ys) for ys in y))
        return N * G
    end
end

# The amounts with the given site fractions that the model selects: the product
# of the fractions.
_product_split(mdl, y) = [prod(y[s][mdl.lattice.occupancy[s, k]] for s in eachindex(y)) for k in axes(mdl.lattice.occupancy, 2)]

_cef_lna(mdl, x, g; T = 298.15, ϵ = 0.0) =
    ChemistryLab._ss_log_activities!(zeros(promote_type(eltype(x), eltype(g)), length(x)), eachindex(x), x, mdl, T, ϵ, g)

@testsection "CASH+: the compound energy formalism" begin
    lat = SublatticeModel([1.0, 2.0], ["A" "A" "A" "B" "B" "B"; "X" "Y" "Z" "X" "Y" "Z"]; sites = ["s1", "s2"])
    mdl = CompoundEnergyModel(lat; interactions = [("s1", "A", "B", -4000.0), ("s2", "X", "Z", 2500.0), ("s2", "Y", "Z", -1500.0)])
    g = [-3.1, -1.7, 0.4, -2.2, 1.9, -0.6] .* 10
    n = [0.3, 0.1, 0.05, 0.2, 0.15, 0.2]
    T = 310.0

    # The activities are the gradient of n G: the chemical potentials minus the
    # standard ones. Euler's relation follows, and the Jacobian is symmetric.
    nG = _cef_nG(mdl, g, T)
    μ = ForwardDiff.gradient(nG, n)
    @test _cef_lna(mdl, n ./ sum(n), g; T) .+ g ≈ μ rtol = 1.0e-12
    @test dot(n, μ) ≈ nG(n) rtol = 1.0e-12
    J = ForwardDiff.jacobian(m -> _cef_lna(mdl, m ./ sum(m), g; T), n)
    @test J ≈ J' atol = 1.0e-10

    # The divergence that selects the split is not negative, and it vanishes, with
    # its gradient, at the product of the site fractions: there the activities
    # are those of the compound energy formalism alone.
    y = site_fractions(mdl, n ./ sum(n))
    xp = _product_split(mdl, y)
    @test sum(xp) ≈ 1
    @test site_fractions(mdl, xp) ≈ y
    @test nG(n) > _cef_nG(mdl, g, T; split = false)(n)
    @test nG(xp) ≈ _cef_nG(mdl, g, T; split = false)(xp) rtol = 1.0e-13
    @test _cef_lna(mdl, xp, g; T) .+ g ≈ ForwardDiff.gradient(_cef_nG(mdl, g, T; split = false), xp) rtol = 1.0e-12

    # A pure end-member has unit activity.
    for k in eachindex(n)
        x = zeros(6)
        x[k] = 1
        @test _cef_lna(mdl, x, g; T, ϵ = 1.0e-300)[k] ≈ 0 atol = 1.0e-12
    end

    # Shifting each standard energy by those of its elements changes nothing: the
    # formalism does not depend on the convention of the database.
    shift = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0] .* 7.3 .+ repeat([0.0, -2.1, 4.4], 2) .+ 11.0
    @test _cef_lna(mdl, n ./ sum(n), g .+ shift; T) ≈ _cef_lna(mdl, n ./ sum(n), g; T) rtol = 1.0e-12

    # With no interaction and no reciprocal energy (standard energies additive
    # over the sites), it is ideal mixing on the sites, at the split selected.
    ideal = CompoundEnergyModel(lat)
    x = xp
    @test _cef_lna(ideal, x, shift; T) ≈ ChemistryLab._ss_log_activities!(zeros(6), 1:6, x, lat, T, 0.0) rtol = 1.0e-12
    # Any reciprocal energy makes it differ.
    @test !(_cef_lna(ideal, x, g; T) ≈ ChemistryLab._ss_log_activities!(zeros(6), 1:6, x, lat, T, 0.0))

    # Integer multiplicities and real interaction energies make one real model,
    # with the same occupancy.
    whole = CompoundEnergyModel(
        SublatticeModel([1, 1], ["A" "A" "A" "B" "B" "B"; "X" "Y" "Z" "X" "Y" "Z"]; sites = ["s1", "s2"]);
        interactions = [("s1", "A", "B", -4000.0)]
    )
    @test whole isa CompoundEnergyModel{Float64}
    @test whole.lattice.occupancy == lat.occupancy
    @test whole.lattice.species == lat.species

    # The end-members must be every compound of the sites, each once.
    partial = SublatticeModel([1.0, 1.0], ["A" "A" "B"; "X" "Y" "X"]; sites = ["s1", "s2"])
    @test_throws ErrorException CompoundEnergyModel(partial)
    @test_throws ErrorException CompoundEnergyModel(lat; interactions = [("s1", "A", "C", 1.0)])
    @test_throws ErrorException CompoundEnergyModel(lat; interactions = [("s1", "A", "B", 1.0), ("s1", "B", "A", 2.0)])
    @test_throws ErrorException CompoundEnergyModel(lat; interactions = [("s3", "A", "B", 1.0)])
    # Its activities cannot be had from the mole fractions alone.
    @test_throws ErrorException ChemistryLab._ss_log_activities!(zeros(6), 1:6, x, mdl, T, 0.0)
    @test mixing_convexity(mdl, 6).verdict === :undecided
    # Given the energies, convexity is decided in the site fractions. A reciprocal
    # energy of 30 RT on two sites of two species each makes the energy concave at
    # the center, which the lattice finds; the same sites with none are convex by
    # the bound.
    sq = SublatticeModel([1.0, 1.0], ["A" "A" "B" "B"; "X" "Y" "X" "Y"]; sites = ["s1", "s2"])
    concave = mixing_convexity(CompoundEnergyModel(sq), 4; g = [0.0, 0.0, 0.0, 30.0])
    @test concave.verdict === :nonconvex && sum(concave.witness) ≈ 1
    @test mixing_convexity(CompoundEnergyModel(sq), 4; g = [0.0, 1.0, 2.0, 3.0]).verdict === :convex
    @test isempty(ChemistryLab._bounded_members(mdl))
end

@testsection "CASH+: the data of Kulik et al. (2022)" begin
    props = literature_table("Kulik2022", "cashplus_standard_properties")
    formulas = literature_table("Kulik2022", "cashplus_end_members")
    species = literature_table("Kulik2022", "cashplus_species")
    shared = literature_table("Kulik2022", "cashplus_shared_sites")
    occ = literature_table("Kulik2022", "cashplus_occupancy")
    @test props.end_member == formulas.end_member == CASHPLUS

    # Each bulk formula (Table 10) is the sum of the moieties its name puts on
    # the sites (Table 2), the two transcribed independently; and each formula
    # unit is neutral.
    moiety(site, sp) = Formula(species.formula[findfirst(r -> species.site[r] == site && species.species[r] == sp, eachindex(species.site))])
    for (k, em) in enumerate(CASHPLUS)
        total = Dict{Symbol, Float64}()
        charge = 0
        parts = [Formula(f) for f in shared.formula]
        for r in eachindex(occ.end_member)
            occ.end_member[r] == em || continue
            push!(parts, moiety(occ.site[r], occ.species[r]))
            charge += species.charge[findfirst(i -> species.site[i] == occ.site[r] && species.species[i] == occ.species[r], eachindex(species.site))]
        end
        for f in parts, (el, c) in composition(f)
            total[el] = get(total, el, 0.0) + c
        end
        @test total == Dict(el => Float64(c) for (el, c) in composition(Formula(formulas.formula[k])))
        @test charge == 0
    end

    # G = H - T (S - the entropies of the elements) holds to 0.01 kJ/mol for five
    # end-members and for CaSiO3@; for TSvh it holds to 0.28 kJ/mol, a
    # transposition of digits in its printed enthalpy (-5032.03 for -5032.30),
    # which the notes of the file record. The element entropies are those of the
    # database.
    db = JSON.parsefile(datapath("cemdata18-cashplus.json"); dicttype = Dict{String, Any})
    S_el = Dict(Symbol(e["symbol"]) => Float64(e["entropy"]["values"][1]) for e in db["elements"])
    mismatch(f, G, H, S) = H - 298.15 * (S - sum(c * S_el[el] for (el, c) in composition(Formula(f)))) / 1000 - G
    kJ(q) = ustrip(us"kJ/mol", q)
    JK(q) = ustrip(us"J/(mol*K)", q)
    for (k, em) in enumerate(CASHPLUS)
        d = mismatch(formulas.formula[k], kJ(props.G[k]), kJ(props.H[k]), JK(props.S[k]))
        em == "TSvh" ? (@test d ≈ 0.28 atol = 0.01) : (@test abs(d) < 0.011)
    end
    q(name) = literature_value("Kulik2022", name)
    @test abs(mismatch("CaSiO3", kJ(q("CaSiO3_aq_G")), kJ(q("CaSiO3_aq_H")), JK(q("CaSiO3_aq_S")))) < 0.01
end

@testsection "CASH+: the database cemdata18-cashplus.json" begin
    base = JSON.parsefile(datapath("cemdata18-thermofun.json"); dicttype = Dict{String, Any})
    db = JSON.parsefile(datapath("cemdata18-cashplus.json"); dicttype = Dict{String, Any})
    bysym(d) = Dict(s["symbol"] => s for s in d["substances"])
    b, c = bysym(base), bysym(db)
    # The base is copied through, CaSiO3@ excepted, and the twelve end-members of
    # CASH+NK added.
    nk = literature_table("Miron2022a", "cashplus_nk_end_members").end_member
    @test sort(collect(setdiff(keys(c), keys(b)))) == sort(nk)
    @test CASHPLUS ⊆ nk
    @test all(c[k] == b[k] for k in keys(b) if k != "CaSiO3@")
    @test c["CaSiO3@"]["sm_gibbs_energy"]["values"][1] ≈ -1514140 rtol = 1.0e-12
    @test b["CaSiO3@"]["sm_gibbs_energy"]["values"][1] ≈ -1517556.9 rtol = 1.0e-12
    subs = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-cashplus.json"); verbose = false))
    at25(sp, key, unit) = ustrip(unit, subs[sp][key](T = 298.15u"K", P = 1.0e5u"Pa"; unit = true))
    # The CaSiO3@ complex holds the reaction properties of Table 9 against the
    # Ca+2 and SiO3-2 of the database, which is what sets its temperature trend,
    # each within the precision the table prints: its H to 0.1 kJ/mol, its S and
    # the reaction heat capacity to the unit.
    Δr(key, unit) = at25("CaSiO3@", key, unit) - at25("Ca+2", key, unit) - at25("SiO3-2", key, unit)
    @test Δr(:ΔₐH⁰, us"kJ/mol") ≈ ustrip(us"kJ/mol", literature_value("Kulik2022", "CaSiO3_aq_dH_reaction")) atol = 0.05
    @test Δr(:S⁰, us"J/(mol*K)") ≈ ustrip(us"J/(mol*K)", literature_value("Kulik2022", "CaSiO3_aq_dS_reaction")) atol = 0.5
    @test Δr(:Cp⁰, us"J/(mol*K)") ≈ ustrip(us"J/(mol*K)", literature_value("Kulik2022", "CaSiO3_aq_dCp_reaction")) atol = 0.5
    props = literature_table("Kulik2022", "cashplus_standard_properties")
    for (k, em) in enumerate(CASHPLUS)
        @test at25(em, :ΔₐG⁰, us"kJ/mol") ≈ ustrip(us"kJ/mol", props.G[k]) rtol = 1.0e-12
        @test at25(em, :V⁰, us"cm^3/mol") ≈ ustrip(us"cm^3/mol", props.V[k]) rtol = 1.0e-12
    end
    # The phase the TOML declares reads the model from the same file.
    shipped = build_solid_solutions(datapath("solid_solutions.toml"), subs)
    cash = only(filter(p -> ChemistryLab.name(p) == "CASH+", shipped))
    m = ChemistryLab.model(cash)
    @test m isa CompoundEnergyModel
    @test m.lattice.species == [["S", "v", "C"], ["v", "C"]]
    @test length(m.interactions) == 4
    # The core model is convex, by the bound on its reference surface.
    g25 = [at25(n, :ΔₐG⁰, us"J/mol") / (ChemistryLab.R_GAS * 298.15) for n in CASHPLUS]
    c = mixing_convexity(m, 6; T = 298.15, g = g25)
    @test c.verdict === :convex && occursin("two sites", c.how)
end

@testsection "CASH+: the Ca-Si-H2O system of Kulik et al. (2022)" begin
    subs = build_species(datapath("cemdata18-cashplus.json"); verbose = false)
    byname = Dict(symbol(s) => s for s in subs)
    cash = only(filter(p -> ChemistryLab.name(p) == "CASH+", build_solid_solutions(datapath("solid_solutions.toml"), byname)))
    species = speciation(subs, vcat(["Portlandite", "Amor-Sl"], CASHPLUS); aggregate_state = [AS_AQUEOUS])
    cs = ChemicalSystem(species, CEMDATA_PRIMARIES; solid_solutions = [cash])
    # The activity model the paper fits with (its Eq. 13): Cemdata18's, for KOH.
    model = cemdata18_activity_model(:KOH)
    mdl = ChemistryLab.model(cash)
    Mw = ustrip(us"g/mol", byname["H2O@"][:M])
    at(eq, sp) = ustrip(us"mol", eq.n[findfirst(s -> symbol(s) == sp, cs.species)])
    function paste(ca_si)
        st = ChemicalState(cs)
        set_quantity!(st, "H2O@", (1000 / Mw)u"mol")
        set_quantity!(st, "Amor-Sl", 0.05u"mol")
        set_quantity!(st, "Portlandite", (0.05 * ca_si)u"mol")
        eq, cert = equilibrate_certified(st; model)
        x = [at(eq, m) for m in CASHPLUS]
        el(e) = sum(x[k] * ChemistryLab.atoms(byname[CASHPLUS[k]])[e] for k in eachindex(x))
        return (; eq, cert, y = site_fractions(mdl, x ./ sum(x)), ca_si = el(:Ca) / el(:Si))
    end
    q(name) = literature_value("Kulik2022", name)
    printed = literature_table("Kulik2022", "portlandite_saturated_site_fractions")

    # Beside portlandite the C-S-H is fixed (three phases, three components): its
    # Ca/Si and site fractions are the paper's, to the digits it prints.
    for ca_si in (2.4, 4.0)
        r = paste(ca_si)
        @test r.cert.optimal
        @test at(r.eq, "Portlandite") > 0
        @test r.ca_si ≈ q("portlandite_saturation_Ca_Si") atol = 0.005
        for (site, sp, f) in zip(printed.site, printed.species, printed.fraction)
            s = findfirst(==(site), mdl.lattice.sites)
            @test r.y[s][findfirst(==(sp), mdl.lattice.species[s])] ≈ f atol = 0.005
        end
    end
    # Beside amorphous silica, the C-S-H with the lowest Ca/Si.
    r = paste(0.6)
    @test r.cert.optimal
    @test at(r.eq, "Amor-Sl") > 0
    @test r.ca_si ≈ q("silica_saturation_Ca_Si") atol = 0.005
    # Between the two, the C-S-H alone.
    r = paste(1.2)
    @test r.cert.optimal
    @test at(r.eq, "Amor-Sl") == 0 && at(r.eq, "Portlandite") == 0
    @test r.ca_si ≈ 1.09 atol = 0.01
end

@testsection "CASH+NK: the alkali extension of Miron et al. (2022a, b)" begin
    names = literature_table("Miron2022a", "cashplus_nk_end_members")
    species = literature_table("Miron2022a", "cashplus_nk_species")
    shared = literature_table("Miron2022a", "cashplus_nk_shared_sites")
    occ = literature_table("Miron2022a", "cashplus_nk_occupancy")
    # Each formula is the sum of its moieties, and neutral.
    row(site, sp) = findfirst(r -> species.site[r] == site && species.species[r] == sp, eachindex(species.site))
    for (k, em) in enumerate(names.end_member)
        total = Dict{Symbol, Float64}()
        charge = 0
        parts = [Formula(f) for f in shared.formula]
        for r in eachindex(occ.end_member)
            occ.end_member[r] == em || continue
            push!(parts, Formula(species.formula[row(occ.site[r], occ.species[r])]))
            charge += species.charge[row(occ.site[r], occ.species[r])]
        end
        for f in parts, (el, c) in composition(f)
            total[el] = get(total, el, 0.0) + c
        end
        @test total == Dict(el => Float64(c) for (el, c) in composition(Formula(names.formula[k])))
        @test charge == 0
    end

    # The database carries the twelve end-members: the core of Kulik2022 with the
    # enthalpy of TSvh Miron2022a reprints, the alkali end-members of Miron2022a,
    # and TCNh and TCKh as Miron2022b fine-tunes them.
    subs = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-cashplus.json"); verbose = false))
    at25(sp, key, unit) = ustrip(unit, subs[sp][key](T = 298.15u"K", P = 1.0e5u"Pa"; unit = true))
    @test all(haskey(subs, n) for n in names.end_member)
    @test at25("TCNh", :ΔₐG⁰, us"kJ/mol") ≈ -4873.07 rtol = 1.0e-12
    @test at25("TCKh", :ΔₐG⁰, us"kJ/mol") ≈ -4891.65 rtol = 1.0e-12
    @test at25("TSNh", :ΔₐG⁰, us"kJ/mol") ≈ -5087.35 rtol = 1.0e-12
    @test at25("TSvh", :ΔₐH⁰, us"kJ/mol") ≈ -5032.3 rtol = 1.0e-12
    shipped = build_solid_solutions(datapath("solid_solutions.toml"), subs)
    nk = only(filter(p -> ChemistryLab.name(p) == "CASH+NK", shipped))
    @test ChemistryLab.model(nk).lattice.species == [["S", "v", "C"], ["v", "C", "N", "K"]]
    # CASH+NK is not decided: the bound fails and no concave point is found.
    gnk = [at25(n, :ΔₐG⁰, us"J/mol") / (ChemistryLab.R_GAS * 298.15) for n in names.end_member]
    @test mixing_convexity(ChemistryLab.model(nk), 12; T = 298.15, g = gnk).verdict === :undecided

    # The discretized model of Miron2022a: 110 pseudocompounds, each a C-S-H of
    # fixed composition whose Gibbs energy the authors computed with their
    # implementation of the model. The formula of each fixes the site fractions,
    # so the Gibbs energy of our model there must be the one its dissolution
    # reaction gives, with the aqueous species of the database. The formulas
    # are printed to four decimals, which moves the energy by up to 0.07 kJ/mol;
    # the tolerance is 0.1 kJ/mol, three orders of magnitude below the energies
    # of mixing and of the reciprocal reactions.
    mdl = compound_energy_model("Miron2022a:cashplus_nk", names.end_member)
    alkali = literature_table("Miron2022a", "alkali_standard_properties")
    core = literature_table("Kulik2022", "cashplus_standard_properties")
    G°(n) = ustrip(
        us"J/mol", n in core.end_member ? core.G[findfirst(==(n), core.end_member)] :
            alkali.G[findfirst(==(n), alkali.end_member)]
    )
    RT = ChemistryLab.R_GAS * 298.15
    g = [G°(n) / RT for n in names.end_member]
    nG = _cef_nG(mdl, g, 298.15; split = false)
    dsp = literature_table("Miron2022a", "dsp_pseudocompounds")
    @test length(dsp.name) == 110
    lat = mdl.lattice
    worst = 0.0
    for r in eachindex(dsp.name)
        lhs, rhs = split(dsp.reaction[r], " = ")
        # The coefficients by their digits: `Formula` would read one within 1e-3
        # of a simple fraction as that fraction (the 2.599 Si of CSH081+K0239 as
        # 13/5), as a database formula means it, and these are computed ones.
        comp = Dict(Symbol(m[1]) => parse(Float64, m[2]) for m in eachmatch(r"([A-Z][a-z]?)([0-9.]+)", first(split(lhs, " + "))))
        Gp = sum(split(rhs, " + ")) do term
            m = match(r"^\s*([0-9.]*)(\S+)$", term)
            (isempty(m[1]) ? 1.0 : parse(Float64, m[1])) * at25(String(m[2]), :ΔₐG⁰, us"J/mol")
        end
        Gx = Gp - ustrip(us"J/mol", dsp.dG[r])
        # The site fractions the formula fixes, from its Si, alkali, H and Ca.
        el(e) = Float64(get(comp, e, 0))
        yA = el(:Na) + el(:K)
        yCIC = (el(:H) - 6 - yA) / 2
        yS = el(:Si) - 2
        yCBT = el(:Ca) - 2 - yCIC
        yBT = Dict("S" => yS, "C" => yCBT, "v" => 1 - yS - yCBT)
        yIC = Dict("v" => 1 - yCIC - yA, "C" => yCIC, "N" => el(:Na), "K" => el(:K))
        # Any amounts with these site fractions: the product of the fractions,
        # which sum to one. The rounding of the formula makes a fraction of 1e-4
        # come out as -1e-4 in two pseudocompounds; it is kept, and its
        # configurational term taken as zero.
        x = [yBT[lat.species[1][lat.occupancy[1, k]]] * yIC[lat.species[2][lat.occupancy[2, k]]] for k in axes(lat.occupancy, 2)]
        d = abs(RT * nG(x) - Gx)
        worst = max(worst, d)
    end
    @test worst < 100

    # The gel in a solution of both hydroxides, three times as much potassium as
    # sodium, as in a cement pore solution: it certifies at the Ca/Si of a gel
    # beside portlandite and takes up both alkalis.
    sp = speciation(collect(values(subs)), vcat(["Portlandite", "Amor-Sl"], names.end_member); aggregate_state = [AS_AQUEOUS])
    complexes(s) = charge(s) == 0 && any(el -> haskey(atoms(s), el), (:Na, :K)) && aggregate_state(s) == AS_AQUEOUS
    cs = ChemicalSystem(filter(!complexes, sp), CEMDATA_PRIMARIES; solid_solutions = [nk])
    st = ChemicalState(cs)
    Mw = ustrip(us"g/mol", subs["H2O@"][:M])
    set_quantity!(st, "H2O@", (1000 / Mw)u"mol")
    set_quantity!(st, "Amor-Sl", 0.05u"mol")
    set_quantity!(st, "Portlandite", 0.08u"mol")
    set_quantity!(st, "Na+", 0.05u"mol")
    set_quantity!(st, "K+", 0.15u"mol")
    set_quantity!(st, "OH-", 0.2u"mol")
    eq, cert = equilibrate_certified(st; model = cemdata18_activity_model(:KOH))
    @test cert.optimal
    gel = solid_solution_totals(eq, "CASH+NK").elements
    @test gel[:Na] > 0 && gel[:K] > 0
    # The members are the product of the site fractions, the split selected.
    x = [ustrip(us"mol", eq.n[findfirst(s -> symbol(s) == m, cs.species)]) for m in names.end_member]
    x ./= sum(x)
    @test x ≈ _product_split(ChemistryLab.model(nk), site_fractions(ChemistryLab.model(nk), x)) rtol = 1.0e-6
end
