# A studio started as a user starts it, on the environment CLSTUDIO_WORKER_PROJECT
# names: its boot and warm-up, the catalog, a calculation run through the page's
# route and compared with the same one written by hand against ChemistryLab in
# the same process, a failing script, a stop, the refusal for lack of memory;
# then the watchdog and the diagnosis of `check = true`, each on a process of
# its own.

using HTTP, JSON

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
    @test all(entries[id]["model"]["mix"]["w_b"] == pub.gruyaert_w_b for id in keys(entries) if occursin("cem1", id) || occursin("gruyaert", id))
    # every entry runs, but the kinetic one (its own test below integrates a shorter run)
    for (id, e) in entries
        e["process"] == "kinetics" && continue
        post("/api/library", Dict("id" => id))
        r = post("/api/run", Dict())
        r["ok"] || @error "library entry $id failed" r["error"] get(r, "output", "")
        @test r["ok"] && r["results"]["kind"] == (e["process"] == "equilibrium" ? "equilibrium" : "path")
    end
    post("/api/library", Dict("id" => "cem1-52.5n-gruyaert2010"))
end

@testset "a kinetic run against the same one written by hand" begin
    model = api("/api/state")["model"]
    rates = [Dict("key" => k, "value" => Dict("kind" => "parrott_killoh_avrami")) for k in ("C3S", "C2S")]
    model["calculation"] = Dict("kind" => "kinetics", "horizon" => 2.0, "ages" => [1.0, 2.0], "rates" => rates)
    post("/api/model", Dict("model" => model))
    run = post("/api/run", Dict())
    run["ok"] || @error "the kinetic run failed" run["error"] get(run, "output", "")
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
    using ChemistryLabStudio
    const CLS = ChemistryLabStudio
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
