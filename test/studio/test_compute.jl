# A studio started as a user starts it, on the environment CLSTUDIO_WORKER_PROJECT
# names: its boot and warm-up, the catalog, a calculation run through the page's
# route and compared with the same one written by hand against ChemistryLab in
# the same process, a failing script, a stop, the refusal for lack of memory;
# then the watchdog and the diagnosis of `check = true`, each on a process of
# its own.

project = ENV["CLSTUDIO_WORKER_PROJECT"]
julia = first(Base.julia_cmd().exec)

h = redirect_stdout(devnull) do
    clstudio(; port = 0, no_browser = true, wait = false, project)
end
cp = h.session.compute
token = last(split(CLS.studio_url(h), "token="))
base = "http://127.0.0.1:$(h.session.port)"
api(path) = JSON.parse(String(HTTP.get(base * path, ["X-Studio-Token" => token]).body); dicttype = Dict{String, Any})
function post(path, body)
    r = HTTP.post(base * path, ["X-Studio-Token" => token, "Content-Type" => "application/json"], JSON.json(body))
    return JSON.parse(String(r.body); dicttype = Dict{String, Any})
end

# the boot prepares nothing (the environment is complete), loads the packages,
# reads the catalog and compiles the default calculation
for _ in 1:1200
    local s = CLS.status(cp)
    (s.state == "ready" && isempty(s.message) || s.state == "failed") && break
    sleep(0.25)
end
@test CLS.status(cp).state == "ready"

@testset "catalog" begin
    cat = api("/api/catalog")
    @test "Portland paste (Lothenbach and Winnefeld 2006)" in [l["name"] for l in cat["phase_lists"]]
    pc = only(e for e in cat["templates"] if e["name"] == "PC (Durdzinski 2017), Bogue")
    @test [x["name"] for x in pc["constituents"]][1:4] == ["C3S", "C2S", "C3A", "C4AF"]
    # the fractions add up to what the published analysis totals, a little under 1
    @test 0.95 < sum(x["mass_fraction"] for x in pc["constituents"]) <= 1 + 1.0e-9
    @test pc["kind"] == "cement" && pc["source"] == "Durdzinski2017"
    cem = only(e for e in cat["templates"] if e["name"] == "CEM I 52.5 N (Gruyaert 2010), Bogue")
    @test cem["source"] == "Gruyaert2010"
    @test all(e -> haskey(e, "constituents") || occursin("species", e["unavailable"]), cat["templates"])
    @test any(d -> d["name"] == "cemdata18-thermofun.json", cat["databases"])
    @test cat["versions"]["julia"] == string(VERSION)
    # the published works, read from ChemistryLab's own files
    @test cat["references"]["Gruyaert2010"]["doi"] == "10.1007/s10973-010-0841-6"
    @test startswith(cat["references"]["Gruyaert2010"]["citation"], "Gruyaert, E.")
end

@testset "diagnosis and species" begin
    d = post("/api/diagnose", Dict())
    @test d["ok"] && [x["code"] for x in d["diagnostics"]] == ["reactant_temperature_range"]
    species = JSON.parse(String(HTTP.get(base * "/api/species?database=cemdata18-thermofun.json", ["X-Studio-Token" => token]).body); dicttype = Dict{String, Any})
    portlandite = only(x for x in species if x["symbol"] == "Portlandite")
    @test portlandite["state"] == "crystal" && portlandite["cement"] == "CH" && sort(portlandite["elements"]) == ["Ca", "H", "O"]
    @test portlandite["molar_mass"] ≈ 74.09 rtol = 1.0e-3
    # a second read is served from the server's copy
    @test length(api("/api/species?database=cemdata18-thermofun.json")) == length(species)
end

@testset "a path against the same one written by hand" begin
    model = api("/api/state")["model"]
    model["calculation"] = Dict("kind" => "carbonation", "amounts" => [0.0, 0.3], "age" => 28.0)
    post("/api/model", Dict("model" => model))
    run = post("/api/run", Dict())
    @test run["ok"] && run["results"]["kind"] == "path"
    res = run["results"]
    # the packages first: `u"K"` below is expanded before the block runs
    @test CLS.call!(cp, :(using ChemistryLab, OptimaSolver, DynamicQuantities)).ok
    reference = CLS.call!(
        cp, quote
            let substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
                cs = phase_list_system("Portland paste (Lothenbach and Winnefeld 2006)", substances)
                cement = with_extents(
                    material_template("CEM I 52.5 N (Gruyaert 2010), Bogue", cs),
                    Dict("C3S" => 0.8, "C2S" => 0.45, "C3A" => 0.9, "C4AF" => 0.6),
                )
                model = cemdata18_activity_model(:KOH)
                rs, _ = equilibrate_certified(Recipe(cement => 1.0; w_b = 0.5, T = T_ZERO_CELSIUS_Q + 20u"K"), cs; model)
                p = carbonate(rs, [0.0, 0.3])
                (; pH = [pore_solution(x).pH for x in p], calcite = [get(phase_masses(x), "Cal", 0.0) for x in p])
            end
        end
    )
    @test reference.ok
    @test res["parameter"] == [0.0, 0.3]
    @test res["pH"] ≈ reference.value.pH rtol = 1.0e-9
    calcite = res["phase_masses"]["series"][findfirst(==("Cal"), res["phase_masses"]["labels"])]
    @test calcite ≈ reference.value.calcite rtol = 1.0e-9
    @test all(res["certified"])
    @test length(res["porosity"]) == 2 && length(res["bound_water"]) == 2
    @test haskey(res["volume_fractions"], "series")
    @test "Ca" in res["molalities"]["labels"]
    post("/api/library", Dict("id" => "cem1-52.5n-gruyaert2010"))
end

@testset "the library" begin
    entries = Dict(e["id"] => e for e in CLS.library_entries())
    # the proportions a formulation takes from its paper, against ChemistryLab's data
    published = CLS.call!(
        cp, quote
            (;
                slag = ustrip(literature_value("Durdzinski2017", "slag_replacement")) / 100,
                fly_ash = ustrip(literature_value("Durdzinski2017", "fly_ash_replacement")) / 100,
                w_b = ustrip(literature_value("Durdzinski2017", "water_binder_ratio")),
                gruyaert_w_b = ustrip(literature_value("Gruyaert2010", "water_binder_ratio")),
                gruyaert_slag = sort(unique(Float64.(literature_table("Gruyaert2010", "bound_water").slag_to_binder))),
                # the mortar C100 of Lavergne et al. in their semi-adiabatic cell:
                # the paste once the water the sand absorbs is set aside, and the
                # heat capacity of the vessel (its 380, printed in kJ/K, is in J/K:
                # scripts/ionic_hydration.jl), of the sand and of that water
                lavergne = let m = literature_row("Lavergne2018", "semi_adiabatic_mixes_wb050", "C100")
                    absorbed = ustrip(literature_value("Lavergne2018", "sand_water_absorption")) * m.dry_sand
                    db = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-thermofun.json"); verbose = false))
                    heat(s, mass) = ustrip(us"J/(mol*K)", db[s][:Cp⁰](T = 293.15, P = 1.0e5, unit = true)) /
                        ustrip(us"kg/mol", db[s][:M]) * ustrip(us"kg", mass)
                    (;
                        binder = ustrip(us"g", m.binder), w_b = ustrip((m.water - absorbed) / m.binder),
                        Cp = ustrip(us"kJ/K", literature_value("Lavergne2018", "calorimeter_heat_capacity")) +
                            heat("Qtz", m.dry_sand) + heat("H2O@", absorbed),
                        a = ustrip(us"W/K", literature_value("Lavergne2018", "heat_loss_a")),
                        b = ustrip(us"W/K^2", literature_value("Lavergne2018", "heat_loss_b")),
                    )
                end,
            )
        end
    )
    @test published.ok
    pub = published.value
    fractions(id) = [m["fraction"] for m in entries[id]["model"]["materials"]]
    @test fractions("pc-s1-slag-rilem") ≈ [1 - pub.slag, pub.slag]
    @test fractions("pc-siliceous-fly-ash-rilem") ≈ [1 - pub.fly_ash, pub.fly_ash]
    @test entries["pc-s1-slag-rilem"]["model"]["mix"]["w_b"] == entries["pc-siliceous-fly-ash-rilem"]["model"]["mix"]["w_b"] == pub.w_b
    @test entries["slag-replacement-gruyaert2010"]["model"]["calculation"]["fractions"] == pub.gruyaert_slag
    @test all(entries[id]["model"]["mix"]["w_b"] == pub.gruyaert_w_b for id in keys(entries) if (occursin("cem1", id) || occursin("gruyaert", id)) && !occursin("lavergne", id))
    # the two water-to-binder ratios of the mortars of Lavergne et al., as their Table 11 names them
    @test entries["cem1-lavergne2018-kinetics-wb050"]["model"]["mix"]["w_b"] == 0.5
    @test entries["cem1-lavergne2018-kinetics-wb032"]["model"]["mix"]["w_b"] == 0.32
    lav = entries["cem1-lavergne2018-semiadiabatic"]["model"]
    cell = lav["calculation"]["calorimeter"]
    @test lav["mix"]["binder_mass"] == pub.lavergne.binder
    @test lav["mix"]["w_b"] ≈ pub.lavergne.w_b rtol = 1.0e-4
    @test cell["heat_capacity"] ≈ pub.lavergne.Cp rtol = 1.0e-4
    @test cell["loss"] ≈ pub.lavergne.a rtol = 1.0e-5
    @test cell["loss_quadratic"] ≈ pub.lavergne.b rtol = 1.0e-4
    # every entry runs, but the kinetic one (its own test below integrates a shorter run)
    for (id, e) in entries
        e["process"] == "kinetics" && continue
        post("/api/library", Dict("id" => id))
        r = post("/api/run", Dict())
        r["ok"] || @error "library entry $id failed" r["error"] Base.get(r, "output", "")
        @test r["ok"] && r["results"]["kind"] == (e["process"] == "equilibrium" ? "equilibrium" : "path")
    end
    post("/api/library", Dict("id" => "cem1-52.5n-gruyaert2010"))
end

@testset "the thermogram and the surface of an equilibrium, against the calls themselves" begin
    model = api("/api/state")["model"]
    model["surface"] = Dict("kind" => "csh_surface")
    # the water of portlandite read between 400 and 520 °C, every other release
    # as the literature reads it
    model["tga"] = Dict(
        "kind" => "tga_windows",
        "windows" => [Dict("key" => "Portlandite (water)", "value" => Dict("kind" => "window_interval", "from" => 400.0, "to" => 520.0))],
    )
    post("/api/model", Dict("model" => model))
    run = post("/api/run", Dict())
    run["ok"] || @error "the run with a surface failed" run["error"] Base.get(run, "output", "")
    @test run["ok"]
    r = run["results"]
    reference = CLS.call!(
        cp, quote
            let substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
                cs = phase_list_system("Portland paste (Lothenbach and Winnefeld 2006)", substances)
                cement = with_extents(
                    material_template("CEM I 52.5 N (Gruyaert 2010), Bogue", cs),
                    Dict("C3S" => 0.8, "C2S" => 0.45, "C3A" => 0.9, "C4AF" => 0.6),
                )
                model = cemdata18_activity_model(:KOH)
                recipe = Recipe(cement => 1.0; w_b = 0.5, T = T_ZERO_CELSIUS_Q + 20u"K")
                rs, _ = equilibrate_certified(recipe, cs; model)
                own = [DecompositionWindow("Portlandite"; between = (673.15, 793.15))]
                tg = thermogram(rs.state, vcat(own, filter(w -> w.phase != "Portlandite", literature_windows(rs.state))); temperatures = collect(303.15:5.0:1253.15))
                portlandite = window_interval(only(w for w in literature_windows(rs.state) if w.phase == "Portlandite"))
                s = with_csh_surface(rs)
                (; mass = tg.mass, portlandite, potential = s.potential, film = s.film, Ca = s.bound[:Ca])
            end
        end
    )
    @test reference.ok
    ref = reference.value
    g = r["tga"]
    @test length(g["temperature"]) == length(ref.mass) && g["temperature"][1] ≈ 30.0 && g["temperature"][end] ≈ 980.0
    @test g["states"][1]["mass"] ≈ 100 .* 1000 .* ref.mass ./ r["binder_mass_g"] rtol = 1.0e-9
    @test haskey(g["states"][1], "by_phase") && sort(g["sources"]) == sort(["LHopital2016b", "Moschner2009", "Shi2016"])
    # the window set, and the literature's beside it, in °C
    ch = only(w for w in g["windows"] if w["key"] == "Portlandite (water)")
    @test ch["own"] && ch["from"] ≈ 400 && ch["to"] ≈ 520 && ch["shape"] == "interval"
    @test [ch["literature"]["from"], ch["literature"]["to"]] ≈ collect(ref.portlandite) .- 273.15 && ch["literature"]["source"] == "DeWeerdt2011"
    @test all(!w["own"] && w["source"] == w["literature"]["source"] for w in g["windows"] if w["key"] != "Portlandite (water)")
    s = r["surface"]
    @test s["potential_mV"] ≈ 1000 * ref.potential rtol = 1.0e-9
    @test s["film_nm"] ≈ 1.0e9 * ref.film rtol = 1.0e-9
    @test only(v for (k, v) in zip(s["bound"]["labels"], s["bound"]["values"]) if k == "Ca") ≈ 100 * ref.Ca / r["binder_mass_g"] rtol = 1.0e-9
    @test s["certified"] && !haskey(s, "donnan")
    # a Donnan layer on this paste is refused, with the reason: its pore water
    # is a film thinner than the Debye length
    model["surface"] = Dict("kind" => "csh_surface", "donnan" => true, "thickness" => 0.1)
    post("/api/model", Dict("model" => model))
    run = post("/api/run", Dict())
    @test !run["ok"] && occursin("Debye length", run["error"])
    # in a suspension, one Debye length thick: its water is its area times its
    # thickness, and its excess is reported
    model["mix"]["w_b"] = 20.0
    model["surface"] = Dict("kind" => "csh_surface", "donnan" => true, "thickness" => 0.0)
    post("/api/model", Dict("model" => model))
    run = post("/api/run", Dict())
    run["ok"] || @error "the suspension with a Donnan layer failed" run["error"] Base.get(run, "output", "")
    s = run["results"]["surface"]
    @test run["ok"] && s["film_nm"] > s["debye_nm"]
    @test s["donnan"]["water_g"] ≈ 100 / run["results"]["binder_mass_g"] * 1.0e-3 * s["area_m2"] * s["debye_nm"] rtol = 1.0e-6
    @test !isempty(s["donnan"]["excess"]["labels"])
    post("/api/library", Dict("id" => "cem1-52.5n-gruyaert2010"))
end

@testset "a glass reacting by its law, with its heat" begin
    post("/api/library", Dict("id" => "pc-siliceous-fly-ash-rilem-kinetics"))
    model = api("/api/state")["model"]
    model["calculation"]["horizon"] = 1.0
    model["calculation"]["ages"] = [1.0]
    post("/api/model", Dict("model" => model))
    run = post("/api/run", Dict())
    run["ok"] || @error "the run with a glass failed" run["error"] Base.get(run, "output", "")
    @test run["ok"] && run["results"]["retcode"] == "Success"
    res = run["results"]
    k = findfirst(==("glass"), res["degrees"]["labels"])
    @test k !== nothing && 0 < res["degrees"]["series"][k][end] < 1
    # the glass has an enthalpy, estimated from its oxides: the heat is known
    @test haskey(res["heat"], "heat") && res["heat"]["heat"][end] > 0
    @test length(res["tga"]["states"]) == 1
    post("/api/library", Dict("id" => "cem1-52.5n-gruyaert2010"))
end

@testset "a kinetic run against the same one written by hand" begin
    model = api("/api/state")["model"]
    rates = [Dict("key" => k, "value" => Dict("kind" => "parrott_killoh_avrami")) for k in ("C3S", "C2S")]
    model["calculation"] = Dict("kind" => "kinetics", "horizon" => 2.0, "ages" => [1.0, 2.0], "rates" => rates)
    post("/api/model", Dict("model" => model))
    run = post("/api/run", Dict())
    run["ok"] || @error "the kinetic run failed" run["error"] Base.get(run, "output", "")
    @test run["ok"] && run["results"]["kind"] == "kinetics"
    res = run["results"]
    @test CLS.call!(cp, :(using OrdinaryDiffEq)).ok
    reference = CLS.call!(
        cp, quote
            let substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
                cs = phase_list_system("Portland paste (Lothenbach and Winnefeld 2006)", substances)
                cement = with_extents(
                    material_template("CEM I 52.5 N (Gruyaert 2010), Bogue", cs),
                    Dict("C3S" => 0.8, "C2S" => 0.45, "C3A" => 0.9, "C4AF" => 0.6),
                )
                model = cemdata18_activity_model(:KOH)
                recipe = Recipe(cement => 1.0; w_b = 0.5, T = T_ZERO_CELSIUS_Q + 20u"K")
                law(p, q) = parrott_killoh_avrami(q, p; blaine = 385.0, w_c = 0.5, α_max = powers_alpha_max(0.5))
                rates = Dict{String, Any}("C3S" => law("C3S", PK84_PARAMS_C3S), "C2S" => law("C2S", PK84_PARAMS_C2S))
                kp = KineticsProblem(
                    recipe, cs, rates, (0.0, 2 * 86400.0); activity_model = model,
                    equilibrium_solver = EquilibriumSolver(cs, model, OptimaOptimizer()),
                )
                sol = integrate(kp, KineticsSolver(; ode_solver = Rodas5P(), reltol = 1.0e-6, abstol = 1.0e-10))
                times = [1.0, 2.0] .* 86400
                states = speciated_states(sol, kp; times)
                _, Q, _ = heat_release(sol, kp; times, states, reference = kp.initial_state)
                # a Dict: this process does not load OrderedCollections to read one back
                (; α = Dict(degrees_of_hydration(sol, kp; times)), pH = [pore_solution(s, model).pH for s in states], Q)
            end
        end
    )
    @test reference.ok
    ref = reference.value
    @test res["parameter"] == [1.0, 2.0] && res["retcode"] == "Success"
    @test res["time"][end] ≈ 2.0
    for (k, name) in enumerate(res["degrees"]["labels"])
        @test res["degrees"]["series"][k][end] ≈ ref.α[name][end] rtol = 1.0e-9
    end
    @test sort(res["degrees"]["labels"]) == ["C2S", "C3S"]
    # a long run is drawn from at most 300 of its instants, the last kept
    t = CLS.call!(cp, :(Main.StudioWorker.sampled(collect(0.0:1.0:1000.0)))).value
    @test length(t) <= 301 && last(t) == 1000.0 && first(t) > 0 && issorted(t)
    @test res["pH"] ≈ ref.pH rtol = 1.0e-9
    @test res["heat"]["heat"] ≈ ref.Q ./ 100 rtol = 1.0e-9
    @test all(x -> 0 < x < 1, res["porosity"]) && length(res["bound_water"]) == 2
    @test "Portlandite" in res["phase_masses"]["labels"] && "pore solution" in res["volume_fractions"]["labels"]
    post("/api/library", Dict("id" => "cem1-52.5n-gruyaert2010"))
end

# Half a day in a semi-adiabatic cell, short of the change of assemblage where
# ChemistryLab's integration stops: the heat and the temperature of the cell.
@testset "a semi-adiabatic run against the same one written by hand" begin
    model = api("/api/state")["model"]
    model["calculation"] = Dict(
        "kind" => "kinetics", "horizon" => 0.5, "ages" => [0.25, 0.5],
        "rates" => [Dict("key" => "C3S", "value" => Dict("kind" => "parrott_killoh_avrami"))],
        "calorimeter" => Dict("kind" => "semi_adiabatic", "heat_capacity" => 900.0, "loss" => 0.05, "ambient" => 20.0),
    )
    post("/api/model", Dict("model" => model))
    run = post("/api/run", Dict())
    run["ok"] || @error "the semi-adiabatic run failed" run["error"] Base.get(run, "output", "")
    @test run["ok"] && run["results"]["retcode"] == "Success"
    h = run["results"]["heat"]
    reference = CLS.call!(
        cp, quote
            let substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
                cs = phase_list_system("Portland paste (Lothenbach and Winnefeld 2006)", substances)
                cement = with_extents(
                    material_template("CEM I 52.5 N (Gruyaert 2010), Bogue", cs),
                    Dict("C3S" => 0.8, "C2S" => 0.45, "C3A" => 0.9, "C4AF" => 0.6),
                )
                model = cemdata18_activity_model(:KOH)
                T = T_ZERO_CELSIUS_Q + 20u"K"
                recipe = Recipe(cement => 1.0; w_b = 0.5, T)
                cal = SemiAdiabaticCalorimeter(; Cp = 900.0u"J/K", L = 0.05u"W/K", T_env = T, T0 = T)
                rates = Dict{String, Any}("C3S" => parrott_killoh_avrami(PK84_PARAMS_C3S, "C3S"; blaine = 385.0, w_c = 0.5, α_max = powers_alpha_max(0.5)))
                kp = KineticsProblem(
                    recipe, cs, rates, (0.0, 0.5 * 86400.0); activity_model = model,
                    equilibrium_solver = EquilibriumSolver(cs, model, OptimaOptimizer()), calorimeter = cal,
                )
                sol = integrate(kp, KineticsSolver(; ode_solver = Rodas5P(), reltol = 1.0e-6, abstol = 1.0e-10))
                _, Q = cumulative_heat(sol, cal)
                _, θ = temperature_profile(sol, cal; times = [last(sol.t)])
                (; Q = last(Q), T = only(θ))
            end
        end
    )
    @test reference.ok
    @test last(h["heat"]) ≈ reference.value.Q / 100 rtol = 1.0e-9
    @test last(h["temperature"]) ≈ reference.value.T - 273.15 rtol = 1.0e-9
    # the paste warms the cell above the room, and its heat only grows
    @test maximum(h["temperature"]) > 20 && issorted(h["heat"])
    @test length(h["time"]) == length(h["heat"]) == length(h["temperature"])
    post("/api/library", Dict("id" => "cem1-52.5n-gruyaert2010"))
end

@testset "a calculation against the same one written by hand" begin
    run = post("/api/run", Dict())
    @test run["ok"]
    res = run["results"]
    reference = CLS.call!(
        cp, quote
            let substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
                cs = phase_list_system("Portland paste (Lothenbach and Winnefeld 2006)", substances)
                cement = with_extents(
                    material_template("CEM I 52.5 N (Gruyaert 2010), Bogue", cs),
                    Dict("C3S" => 0.8, "C2S" => 0.45, "C3A" => 0.9, "C4AF" => 0.6),
                )
                rs, _ = equilibrate_certified(
                    Recipe(cement => 1.0; w_b = 0.5, T = T_ZERO_CELSIUS_Q + 20u"K"), cs;
                    model = cemdata18_activity_model(:KOH),
                )
                # the aqueous species counted as one, the pore solution
                liquid = Set(symbol(s) for s in cs.species if aggregate_state(s) == AS_AQUEOUS)
                fractions = Dict{String, Float64}()
                for (k, f) in volume_fractions(rs)
                    key = k in liquid ? "pore solution" : k
                    fractions[key] = get(fractions, key, 0.0) + f
                end
                (;
                    pH = pore_solution(rs).pH, bound_water = bound_water(rs), masses = Dict(phase_masses(rs)),
                    fractions, porosity = porosity(rs).total,
                )
            end
        end
    )
    @test reference.ok
    ref = reference.value
    @test res["pore_solution"]["pH"] ≈ ref.pH rtol = 1.0e-9
    @test res["bound_water"] ≈ ref.bound_water rtol = 1.0e-9
    masses = Dict(zip(res["phase_masses"]["labels"], res["phase_masses"]["values"]))
    @test keys(masses) == keys(ref.masses)
    @test all(isapprox(masses[k], ref.masses[k]; rtol = 1.0e-9) for k in keys(masses))
    @test res["certificate"]["optimal"] === true
    fractions = Dict(zip(res["volume_fractions"]["labels"], res["volume_fractions"]["values"]))
    @test keys(fractions) == keys(ref.fractions)
    @test all(isapprox(fractions[k], ref.fractions[k]; rtol = 1.0e-9) for k in keys(fractions))
    @test sum(values(fractions)) ≈ 1 rtol = 1.0e-9          # the void closes the sum
    @test res["porosity"]["total"] ≈ ref.porosity rtol = 1.0e-9
    # a cement whose minor oxides have no sourced density: the volumes are unknown, and said so
    post("/api/library", Dict("id" => "portland-rilem-pc"))
    other = post("/api/run", Dict())["results"]
    @test occursin("No sourced density", other["volume_fractions"]["reason"])
    @test other["porosity"]["total"] === nothing
    # and along a path
    model = api("/api/state")["model"]
    model["calculation"] = Dict("kind" => "carbonation", "amounts" => [0.0, 0.3], "age" => 28.0)
    post("/api/model", Dict("model" => model))
    path = post("/api/run", Dict())["results"]
    @test occursin("No sourced density", path["volume_fractions"]["reason"])
    post("/api/library", Dict("id" => "cem1-52.5n-gruyaert2010"))
    @test res["saturation_indices"]["values"] == sort(res["saturation_indices"]["values"]; rev = true)
end

@testset "what the process prints, a script that fails" begin
    @test CLS.call!(cp, :(println("from the compute process"); flush(stdout))).ok
    for _ in 1:40
        "from the compute process" in CLS.status(cp).log && break
        sleep(0.05)
    end
    @test "from the compute process" in CLS.status(cp).log
    r = CLS.call!(cp, :(Main.StudioWorker.run_compute("x = 1\nerror(\"boom\")\n")))
    @test r.ok && !r.value["ok"]
    @test occursin("boom", r.value["error"]) && occursin("studio_script.jl:2", r.value["error"])
    r = CLS.call!(cp, :(Main.StudioWorker.run_compute("println(\"hello\")\n")))
    @test r.value["ok"] && r.value["output"] == "hello\n" && r.value["results"] === nothing
    r = CLS.call!(cp, :(error("internal")))
    @test !r.ok && occursin("internal", r.error)
end

@testset "busy, stop, refusal" begin
    # A stop may kill the process (on Windows, where the interrupt does not
    # reach it), and its coverage with it: it exits normally first, writing what
    # the tests above measured.
    CLS.shutdown!(cp)
    CLS.start!(cp)
    @test !CLS.stop!(cp)
    t = @async CLS.call!(cp, :(Main.StudioWorker.run_compute("sleep(60)\n")))
    for _ in 1:100
        CLS.status(cp).state == "busy" && break
        sleep(0.05)
    end
    @test CLS.call!(cp, :(1 + 1)).error == "A calculation is already running."
    sleep(1)
    @test CLS.stop!(cp)
    r = fetch(t)
    @test !r.ok
    if r.stopped === :user
        # the interrupt did not take (on Windows, where the console event does
        # not reach the process): the watchdog killed it after 5 s
        @test occursin("restarting", r.error) && CLS.status(cp).state == "stopped"
        CLS.start!(cp)
    else
        # the interrupt took: no restart, the compilation is kept
        @test r.error == "Stopped." && CLS.status(cp).state == "ready"
    end
    Sys.iswindows() || @test r.stopped === :none
    r = CLS.call!(cp, :(1 + 1); floor = typemax(Int))
    @test !r.ok && occursin("GB of memory are free", r.error)
    @test CLS.call!(cp, :(1 + 1)).value == 2
end

close(h)          # a normal exit of the compute process
@test CLS.status(cp).state == "stopped"

# The compute process ends when the process that started it is killed: its
# connection closes, and Malt's loop ends it.
Sys.iswindows() || @testset "the compute process follows its studio" begin
    code = """
    using ChemistryLab
    const CLS = ChemistryLab.Studio
    cp = CLS.ComputeProcess(; project = $(repr(project)), julia = $(repr(julia)), boot = joinpath(CLS.WORKER_DIR, "boot.jl"))
    CLS.start!(cp)
    println("worker ", cp.worker.proc_pid)
    flush(stdout)
    sleep(600)
    """
    log = tempname()
    driver = run(pipeline(`$(Base.julia_cmd()) $(CLS.coverage_flags()) --project=$(Base.active_project()) -e $code`; stdout = log, stderr = log); wait = false)
    found = nothing
    for _ in 1:2400
        found = match(r"worker (\d+)", read(log, String))
        (found !== nothing || !process_running(driver)) && break
        sleep(0.1)
    end
    @test found !== nothing
    pid = parse(Cint, found[1])
    alive(p) = ccall(:kill, Cint, (Cint, Cint), p, 0) == 0
    @test alive(pid)
    kill(driver, Base.SIGKILL)
    for _ in 1:300
        alive(pid) || break
        sleep(0.1)
    end
    @test !alive(pid)
    rm(log; force = true)
end

@testset "watchdog" begin
    cp2 = CLS.ComputeProcess(; project, julia, boot = joinpath(CLS.WORKER_DIR, "boot.jl"))
    CLS.start!(cp2)
    rss = something(CLS.process_rss(cp2.worker.proc_pid), 0)
    r = CLS.call!(cp2, :(x = ones(UInt8, 600_000_000); sleep(30); sum(x)); cap = rss + 200 * 2^20)
    @test !r.ok && r.stopped === :memory && occursin("Make the calculation smaller", r.error)
    @test CLS.status(cp2).state == "stopped"
    @test occursin("not ready", CLS.call!(cp2, :(1 + 1)).error)
end

@testset "check = true" begin
    out = Pipe()
    ok = redirect_stdout(out) do
        withenv(() -> clstudio(; check = true, project), "CLSTUDIO_PROJECT" => nothing)
    end
    close(out.in)
    text = read(out, String)
    @test ok
    @test occursin("ChemistryLab", text) && occursin("The compute environment is ready.", text)
end

@testset "a check that fails" begin
    a = CLS.APP
    app = CLS.StudioApp(;
        name = a.name, version = a.version, default_port = a.default_port,
        project_variable = a.project_variable, default_project = a.default_project,
        environment = a.environment, boot = a.boot, web_dir = a.web_dir, new_state = a.new_state,
        routes! = a.routes!, warm_up! = a.warm_up!, check = :(error("no check today")),
    )
    cp3 = CLS.ComputeProcess(; project, julia, boot = a.boot)
    out = Pipe()
    ok = redirect_stdout(out) do
        CLS.run_check(app, cp3; managed = false, update = false)
    end
    close(out.in)
    @test !ok && occursin("no check today", read(out, String))
end
