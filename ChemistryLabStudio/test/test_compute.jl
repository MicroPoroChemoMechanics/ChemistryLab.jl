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
end

@testset "a calculation against the same one written by hand" begin
    run = post("/api/run", Dict())
    @test run["ok"]
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
                rs, _ = equilibrate_certified(
                    Recipe(cement => 1.0; w_b = 0.5, T = T_ZERO_CELSIUS_Q + 20u"K"), cs;
                    model = cemdata18_activity_model(:KOH),
                )
                (;
                    pH = pore_solution(rs).pH, bound_water = bound_water(rs), masses = Dict(phase_masses(rs)),
                    fractions = Dict(volume_fractions(rs)), porosity = porosity(rs).total,
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
    @test !r.ok && r.error == "Stopped." && r.stopped === :none
    @test CLS.status(cp).state == "ready"          # the interrupt took: no restart
    r = CLS.call!(cp, :(1 + 1); floor = typemax(Int))
    @test !r.ok && occursin("GB of memory are free", r.error)
    @test CLS.call!(cp, :(1 + 1)).value == 2
end

close(h)          # a normal exit of the compute process
@test CLS.status(cp).state == "stopped"

# The watchdog reads the resident memory in /proc; elsewhere the floor checked
# before each run is the protection.
Sys.islinux() && @testset "watchdog" begin
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
