# Chloride in a paste whose C-S-H an equilibrium has decided, in two stages.
#
# The surface model of Guo et al. (2018), silanol sites on a C-S-H that binds
# calcium, sodium and chloride, cannot be put on the CSHQ solid solution: its
# end members hold the same calcium and alkalis, and would count them twice.
# The model therefore needs a C-S-H of fixed composition, and leaves open which
# composition. A first equilibrium with CSHQ gives it; `freeze_solid_solution`
# sets the gel aside, returns its alkalis to the solution, and the second
# equilibrium puts the sites on the frozen gel.
#
# Two parts. The helpers on the gel of Hong & Glasser (1999) at Ca/Si 1.8 in
# sodium and potassium hydroxide, a published system in which the alkali end
# members of CSHQ hold something and portlandite forms: the totals, the exact
# conservation across the stages, dual numbers, and what is refused. Then Guo's
# paste through both stages against PHREEQC, test/reference/phreeqc_csh_frozen.json,
# which declares CSHQ as an ideal SOLID_SOLUTIONS of its own at the first stage
# and carries no solid solution at the second: the first stage is checked as
# well as the second.

include("reference_species.jl")

const PHREEQC_CSH_FROZEN = reference_oracle("phreeqc_csh_frozen")

@testsection "A solid solution, read and then frozen" begin
    subs = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
    byname = Dict(symbol(s) => s for s in subs)
    cshq = ["CSHQ-TobD", "CSHQ-TobH", "CSHQ-JenH", "CSHQ-JenD", "NaSiOH", "KSiOH"]
    model = HKFActivityModel(å = 0.0, Ḃ = gems_bdot(), Kₙ = 0.0)
    M(s) = ustrip(us"g/mol", Species(s)[:M])
    excluded = split("H2@ O2@ CH4@")

    sp1 = speciation(
        subs, vcat(["Portlandite", "Amor-Sl", "Lim"], cshq);
        aggregate_state = [AS_AQUEOUS], exclude_species = excluded
    )
    cs1 = ChemicalSystem(
        sp1, CEMDATA_PRIMARIES;
        solid_solutions = [SolidSolutionPhase("CSHQ", [byname[m] for m in cshq])],
    )
    # Hong & Glasser's experiment at Ca/Si 1.8, as test/hong_glasser1999_reference.jl
    # sets it up, in 100 mM NaOH with 10 mM KOH so that both alkali end members
    # hold something.
    hg(q) = literature_value("HongGlasser1999", q)
    scale = 1000 / ustrip(u"mL", hg("solution_volume"))
    solid_g = ustrip(us"g", hg("solid_mass"))
    w = ustrip(hg("water_content_high_CaSi"))
    n_Si = solid_g * (1 - w) / (1.8 * M("CaO") + M("SiO2")) * scale
    st = ChemicalState(cs1; T = (273.15 + ustrip(hg("temperature_C"))) * u"K")
    set_quantity!(st, "Lim", 1.8 * n_Si * u"mol")
    set_quantity!(st, "Amor-Sl", n_Si * u"mol")
    set_quantity!(st, "H2O@", (1000 + solid_g * w * scale) / M("H2O") * u"mol")
    set_quantity!(st, "Na+", 0.1u"mol")
    set_quantity!(st, "K+", 0.01u"mol")
    set_quantity!(st, "OH-", 0.11u"mol")
    eq1, cert1 = equilibrate_certified(st; model)
    @test cert1.optimal
    n1 = Float64[ustrip(us"mol", x) for x in eq1.n]
    idx1 = Dict(symbol(s) => i for (i, s) in enumerate(cs1.species))
    @test n1[idx1["Portlandite"]] > 0

    # The second system: the same elements and chloride, no C-S-H.
    sp2 = speciation(
        subs, ["Portlandite", "HSiO3-", "Na+", "K+", "Cl-"];
        aggregate_state = [AS_AQUEOUS], exclude_species = excluded
    )
    cs2 = ChemicalSystem(sp2, CEMDATA_PRIMARIES)
    elements(cs, n) = let d = Dict{Symbol, Float64}()
        for (s, x) in zip(cs.species, n), (e, v) in atoms(s)
            d[e] = get(d, e, 0.0) + v * x
        end
        d
    end

    @testset "the totals, by hand" begin
        t = solid_solution_totals(eq1, "CSHQ")
        @test collect(keys(t.members)) == cshq
        @test t.amount ≈ sum(n1[idx1[m]] for m in cshq) rtol = 1.0e-14
        for e in (:Ca, :Si, :Na, :K, :O, :H)
            @test t.elements[e] ≈ sum(Float64(get(atoms(byname[m]), e, 0)) * n1[idx1[m]] for m in cshq) rtol = 1.0e-14
        end
        @test t.mass ≈ sum(ustrip(us"kg/mol", byname[m][:M]) * n1[idx1[m]] for m in cshq) rtol = 1.0e-14
        # The gel at portlandite saturation, 1.618, below the 1.8 it was made
        # at: the rest of the calcium is in the portlandite.
        @test t.elements[:Ca] / t.elements[:Si] ≈ 1.618 atol = 5.0e-4
        @test t.elements[:Na] > 0 && t.elements[:K] > 0
    end

    @testset "freezing conserves every element" begin
        res = freeze_solid_solution(eq1, "CSHQ", cs2; buffer = "Portlandite")
        n2 = Float64[ustrip(us"mol", x) for x in res.state.n]
        b1, b2 = elements(cs1, n1), elements(cs2, n2)
        for e in keys(b1)
            @test b2[e] + get(res.frozen.elements, e, 0.0) ≈ b1[e] atol = 1.0e-12 * max(1.0, b1[e])
        end
        # The alkalis are all released, each with one hydroxide: NaSiOH gives
        # back 0.5 NaOH per formula unit and keeps 0.2 SiO2 and 0.2 H2O.
        t = solid_solution_totals(eq1, "CSHQ")
        @test res.frozen.released[:Na] == t.elements[:Na]
        @test res.frozen.released[:K] == t.elements[:K]
        @test res.frozen.elements[:Na] == 0 && res.frozen.elements[:K] == 0
        idx2 = Dict(symbol(s) => j for (j, s) in enumerate(cs2.species))
        @test n2[idx2["Na+"]] ≈ n1[idx1["Na+"]] + t.elements[:Na] rtol = 1.0e-14
        @test n2[idx2["OH-"]] ≈ n1[idx1["OH-"]] + t.elements[:Na] + t.elements[:K] rtol = 1.0e-14
        @test res.frozen.mass ≈ t.mass - t.elements[:Na] * (M("Na") + M("OH")) / 1000 -
            t.elements[:K] * (M("K") + M("OH")) / 1000 rtol = 1.0e-12
        # And what stays frozen is neutral C-S-H: calcium, silicon, oxygen and
        # hydrogen, oxygen and hydrogen in the proportions of oxides and water.
        fe = res.frozen.elements
        @test fe[:O] ≈ fe[:Ca] + 2fe[:Si] + fe[:H] / 2 rtol = 1.0e-12
        # With nothing released, the whole gel is frozen, alkalis included.
        all_in = freeze_solid_solution(eq1, "CSHQ", cs2; release = ())
        @test isempty(all_in.frozen.released)
        @test all_in.frozen.elements[:Na] == t.elements[:Na]
        @test all_in.frozen.mass == t.mass
    end

    @testset "dual numbers pass through" begin
        # d/ds of everything at n(s) = s n₁, s = 1: each total is its own
        # derivative.
        nd = [ForwardDiff.Dual(x, x) for x in n1]
        sd = ChemicalState(cs1, nd .* u"mol"; T = temperature(eq1))
        t = solid_solution_totals(sd, "CSHQ")
        @test ForwardDiff.partials(t.elements[:Ca])[1] ≈ ForwardDiff.value(t.elements[:Ca]) rtol = 1.0e-14
        @test ForwardDiff.partials(t.mass)[1] ≈ ForwardDiff.value(t.mass) rtol = 1.0e-14
        res = freeze_solid_solution(sd, "CSHQ", cs2)
        x = ustrip(us"mol", moles(res.state, byname["Na+"]))
        @test ForwardDiff.partials(x)[1] ≈ ForwardDiff.value(x) rtol = 1.0e-14
        @test ForwardDiff.partials(res.frozen.elements[:Si])[1] ≈
            ForwardDiff.value(res.frozen.elements[:Si]) rtol = 1.0e-14
    end

    @testset "what is refused" begin
        # A second system that cannot hold the first stage's alkalis.
        bare = ChemicalSystem(
            speciation(subs, ["Portlandite", "HSiO3-"]; aggregate_state = [AS_AQUEOUS], exclude_species = excluded),
            CEMDATA_PRIMARIES,
        )
        @test_throws ArgumentError freeze_solid_solution(eq1, "CSHQ", bare)
        # One that could make the gel again.
        again = ChemicalSystem(vcat(sp2, [byname["CSHQ-JenD"]]), CEMDATA_PRIMARIES)
        @test_throws ArgumentError freeze_solid_solution(eq1, "CSHQ", again)
        # A symbol is a label: "Na+" that is potassium.
        fake = [symbol(s) == "Na+" ? Species("K+"; symbol = "Na+", aggregate_state = AS_AQUEOUS, class = SC_AQSOLUTE) : s for s in sp2]
        @test_throws ArgumentError freeze_solid_solution(eq1, "CSHQ", ChemicalSystem(fake))
        # A released element with no monovalent cation to go to.
        @test_throws ArgumentError freeze_solid_solution(eq1, "CSHQ", cs2; release = (:Na, :Ca))
        # No buffer, an unknown buffer, an unknown phase, a system with none.
        gone = ChemicalState(cs1, [symbol(s) == "Portlandite" ? 0.0u"mol" : x for (s, x) in zip(cs1.species, eq1.n)])
        @test_throws ArgumentError freeze_solid_solution(gone, "CSHQ", cs2; buffer = "Portlandite")
        @test_throws ArgumentError freeze_solid_solution(eq1, "CSHQ", cs2; buffer = "Cal")
        @test_throws ArgumentError solid_solution_totals(eq1, "CNASH_ss")
        @test_throws ArgumentError solid_solution_totals(ChemicalState(cs2), "CSHQ")
    end

    @testset "sites hosted on the solid solution are refused" begin
        ions = [byname[s] for s in ("H2O@", "H+", "OH-", "Na+", "Cl-", "Ca+2")]
        t = literature_table("Guo2018", "surface_reactions_phreeqc")
        support = SurfaceSupport("C-S-H", "CSHQ-TobD", FixedSurfaceArea(1.0e3))
        on_gel(reactions) = site_family(
            "Csh_w", reactions, ions; master = "Csh_w", site = "Xw",
            capacity = TotalSiteAmount(1.0u"mol"), support, model = DiffuseLayer(; area = 1.0e3),
        )
        members = [byname[m] for m in cshq]
        prim = [byname["H2O@"], byname["H+"], byname["Na+"], byname["Ca+2"], byname["Cl-"], byname["HSiO3-"], byname["K+"]]
        build(f) = ChemicalSystem(
            vcat(ions, [byname["HSiO3-"], byname["K+"]], members, [f.free_site], f.complexes),
            vcat(prim, [f.free_site]);
            solid_solutions = [SolidSolutionPhase("CSHQ", members)], site_families = [f],
        )
        # Guo's sites bind calcium, sodium and chloride; CSHQ holds the first two.
        guo = on_gel([eq => lk for (eq, lk) in zip(t.reaction, t.log_K) if !occursin("K+", eq)])
        err = try
            build(guo)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("Ca, Na", err.msg)
        # Deprotonation alone binds nothing the gel holds.
        @test build(on_gel([first(t.reaction) => first(t.log_K)])) isa ChemicalSystem
    end
end

@testsection "Guo's paste on CSHQ, then frozen, against PHREEQC" begin
    f = PHREEQC_CSH_FROZEN
    subs = Dict(
        symbol(s) => s for s in build_species(datapath("cemdata18-thermofun.json"); verbose = false)
    )
    G(s) = ustrip(us"J/mol", s[:ΔₐG⁰](T = 298.15u"K", P = 1.0e5u"Pa"; unit = true))
    cshq = String.(collect(keys(f.stage1.cshq)))
    guo(q) = literature_value("Guo2018", q)
    molar(s) = ustrip(us"g/mol", subs[s][:M])

    @testset "the two codes share their inputs" begin
        for (ph, prod) in pairs(f.phase_reactions)
            drg = sum(ν * G(subs[String(sp)]) for (sp, ν) in pairs(prod)) - G(subs[String(ph)])
            @test -drg / (R_GAS * 298.15 * log(10)) ≈ f.phase_log_k[ph] atol = 1.0e-10
            rhs = Dict{Symbol, Float64}()
            for (sp, ν) in pairs(prod), (el, k) in atoms(subs[String(sp)])
                rhs[el] = get(rhs, el, 0.0) + ν * k
            end
            @test all(get(rhs, el, 0.0) ≈ k for (el, k) in atoms(subs[String(ph)]))
            @test abs(sum(ν * charge(subs[String(sp)]) for (sp, ν) in pairs(prod))) < 1.0e-12
        end
        # Guo's C-S-H per silicon and its molar mass, from the library.
        t = literature_table("Guo2018", "csh_formula")
        per = Dict(zip(t.oxide, t.coefficient ./ t.coefficient[findfirst(==("SiO2"), t.oxide)]))
        @test f.csh_molar_mass_per_si ≈ sum(k * ustrip(us"g/mol", Species(ox)[:M]) for (ox, k) in per) rtol = 1.0e-12
    end

    ions1 = [subs[s] for s in ("H2O@", "H+", "OH-", "Ca+2", "AlO2-", "SO4-2", "HSiO3-")]
    hyd1 = [subs[s] for s in ("Portlandite", "monosulphate12", "monosulphate14", "ettringite")]
    gel = [subs[m] for m in cshq]
    cs1 = ChemicalSystem(
        vcat(ions1, hyd1, gel), [subs[s] for s in ("H2O@", "H+", "Ca+2", "AlO2-", "SO4-2", "HSiO3-")];
        solid_solutions = [SolidSolutionPhase("CSHQ", gel)],
    )
    # Guo's inventory brought to a kilogram of pore water. The C-S-H enters as
    # its elements: its lime with the portlandite, its silica as HSiO3- + H+,
    # its water with the solution's, as PHREEQC's REACTION adds them.
    scale = 1000 / (10 * ustrip(guo("porosity_percent")))
    grams(q) = ustrip(us"g", guo(q)) * scale
    per = f.csh_per_si
    n_csh = grams("csh_per_liter") / f.csh_molar_mass_per_si
    @test n_csh ≈ f.n_csh_mol rtol = 1.0e-12
    st = ChemicalState(cs1)
    set_quantity!(st, "Portlandite", (grams("ch_per_liter") / molar("Portlandite") + per.CaO * n_csh) * u"mol")
    set_quantity!(st, "monosulphate12", grams("afm_per_liter") / molar("monosulphate12") * u"mol")
    set_quantity!(st, "ettringite", grams("aft_per_liter") / molar("ettringite") * u"mol")
    set_quantity!(st, "HSiO3-", n_csh * u"mol")
    set_quantity!(st, "H+", n_csh * u"mol")
    set_quantity!(st, "H2O@", (1000 / molar("H2O@") + (per.H2O - per.CaO - 1) * n_csh) * u"mol")
    for (ph, n) in pairs(f.initial_mol)
        @test grams(Dict(:Portlandite => "ch_per_liter", :monosulphate12 => "afm_per_liter", :ettringite => "aft_per_liter")[ph]) / molar(String(ph)) ≈ n rtol = 1.0e-12
    end
    model = DaviesActivityModel()
    b1 = Float64.(cs1.SM.A) * ustrip.(us"mol", st.n)
    eq1, cert1 = equilibrate_certified(st; model, b = b1)
    @test cert1.optimal
    idx1 = Dict(symbol(s) => i for (i, s) in enumerate(cs1.species))
    n1 = Float64[ustrip(us"mol", x) for x in eq1.n]
    gel1 = solid_solution_totals(eq1, "CSHQ")

    @testset "the first stage: CSHQ against PHREEQC's solid solution" begin
        worst_phase = maximum(abs(n1[idx1[String(ph)]] - m) for (ph, m) in pairs(f.stage1.phases))
        worst_member = maximum(abs(n1[idx1[String(em)]] - m) for (em, m) in pairs(f.stage1.cshq))
        @info "stage 1 against PHREEQC" worst_phase worst_member
        # Measured at 1.1e-5 mol on a phase and 5.5e-5 mol on an end member,
        # out of 8.6 and 4.3 mol, 2.3e-6 relative on the calcium of the gel and
        # 3.0e-4 on the pH; bounded with the margins of test/csh_surface.jl.
        @test worst_phase < 1.0e-4
        @test worst_member < 1.0e-3
        @test gel1.elements[:Si] ≈ f.frozen.Si_mol rtol = 1.0e-5
        @test gel1.elements[:Ca] ≈ f.frozen.Ca_mol rtol = 1.0e-5
        @test pH(eq1, model) ≈ f.stage1.pH atol = 1.0e-3
        # Portlandite is there to buffer the gel: 8.62 mol, Ca/Si 1.628.
        @test n1[idx1["Portlandite"]] ≈ 8.62 atol = 5.0e-4
        @test gel1.elements[:Ca] / gel1.elements[:Si] ≈ 1.6276 atol = 5.0e-5
    end

    # The second stage: the sites on the frozen gel, referred to its silicon
    # through Guo's molar mass per silicon.
    ions2 = [subs[s] for s in ("H2O@", "H+", "OH-", "Na+", "Cl-", "Ca+2", "AlO2-", "SO4-2", "HSiO3-")]
    hyd2 = [subs[String(ph)] for ph in keys(f.points[1].phases)]
    n_sites = ustrip(u"mol/g", guo("csh_site_density")) * f.csh_molar_mass_per_si * gel1.elements[:Si]
    area = ustrip(u"m^2/g", guo("csh_specific_surface_area")) * f.csh_molar_mass_per_si * gel1.elements[:Si]
    @test n_sites ≈ f.frozen.n_sites rtol = 1.0e-6
    @test area ≈ f.frozen.area_m2 rtol = 1.0e-6
    t = literature_table("Guo2018", "surface_reactions_phreeqc")
    reactions = [eq => lk for (eq, lk) in zip(t.reaction, t.log_K) if !occursin("K+", eq)]
    dl = DiffuseLayer(; area)
    family = site_family(
        "Csh_w", reactions, ions2; master = "Csh_w", site = "Xw",
        capacity = TotalSiteAmount(n_sites * u"mol"),
        support = SurfaceSupport("C-S-H", nothing, FixedSurfaceArea(area)), model = dl,
    )
    members = vcat([family.free_site], family.complexes)
    cs2 = ChemicalSystem(
        vcat(ions2, hyd2, members),
        vcat([subs[s] for s in ("H2O@", "H+", "Na+", "Ca+2", "Cl-", "AlO2-", "SO4-2", "HSiO3-")], [family.free_site]);
        site_families = [family],
    )
    res = freeze_solid_solution(eq1, "CSHQ", cs2; buffer = "Portlandite")
    idx = Dict(symbol(sp) => k for (k, sp) in enumerate(cs2.species))
    function charged(nacl)
        n = Float64[ustrip(us"mol", x) for x in res.state.n]
        n[idx["XwOH"]] = n_sites
        n[idx["Na+"]] += max(nacl, 1.0e-10)
        n[idx["Cl-"]] += max(nacl, 1.0e-10)
        return ChemicalState(cs2, n .* u"mol")
    end
    pts = sort(collect(f.points); by = p -> -p.nacl)
    A = Float64.(cs2.SM.A)
    budgets = [A * ustrip.(us"mol", charged(p.nacl).n) for p in pts]
    states, certs = equilibrate_path(charged(first(pts).nacl), budgets; model)
    @test all(c.optimal for c in certs)

    @testset "the second stage against PHREEQC" begin
        worst = (phase = 0.0, surface = 0.0, pH = 0.0, psi = 0.0, cl_bound = 0.0)
        for (pt, eq) in zip(pts, states)
            n = Float64[ustrip(us"mol", x) for x in eq.n]
            for (ph, m) in pairs(pt.phases)
                worst = merge(worst, (phase = max(worst.phase, abs(n[idx[String(ph)]] - m)),))
            end
            for (sp, m) in pairs(pt.surface)
                sym = replace(String(sp), "Csh_w" => "Xw"; count = 1)
                worst = merge(worst, (surface = max(worst.surface, abs(n[idx[sym]] - m)),))
            end
            z = Float64[charge(sp) for sp in members]
            psi = diffuse_layer_potential(dl, z, [n[idx[symbol(sp)]] for sp in members], ionic_strength(eq), 298.15)
            bound = 2n[idx["C4AClH10"]] + n[idx["C4AsClH12"]] + n[idx["XwOHCl-"]]
            bound_ref = 2pt.phases.C4AClH10 + pt.phases.C4AsClH12 + pt.surface[Symbol("Csh_wOHCl-")]
            worst = merge(
                worst, (
                    pH = max(worst.pH, abs(pH(eq, model) - pt.pH)),
                    psi = max(worst.psi, abs(psi - pt.psi_V * FARADAY / (R_GAS * 298.15))),
                    cl_bound = max(worst.cl_bound, abs(bound - bound_ref)),
                ),
            )
        end
        @info "stage 2 against PHREEQC" worst
        # Measured at 2.8e-4 mol on a phase, 4.8e-4 mol on a surface species,
        # 1.3e-4 on the pH, 3.8e-4 on FΨ/RT and 5.5e-5 mol on the bound
        # chloride: the one-stage paste's figures, bounded with its margins.
        @test worst.phase < 1.0e-3
        @test worst.surface < 1.0e-3
        @test worst.pH < 1.0e-3
        @test worst.psi < 5.0e-3
        @test worst.cl_bound < 1.0e-3
    end

    @testset "the portlandite still buffers the gel" begin
        # The frozen composition is a model only while portlandite coexists
        # with the gel. The surface takes 2.5 mol of calcium from it, which is
        # the calcium a surface on CSHQ would have counted twice, and 70 % of it
        # is left at every point.
        ch = [ustrip(us"mol", eq.n[idx["Portlandite"]]) for eq in states]
        @test minimum(ch) / n1[idx1["Portlandite"]] ≈ 0.703 atol = 5.0e-4
        @test ustrip(us"mol", states[end].n[idx["XwOCa+"]]) ≈ 2.551 atol = 5.0e-4
    end
end
