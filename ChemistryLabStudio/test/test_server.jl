# The HTTP side, with a session whose compute process is never started: the
# guard of every `/api` route, the static files, and the routes that need no
# calculation. Requests go through a real HTTP client.

using HTTP, JSON

cp0 = CLS.ComputeProcess(; project = joinpath(tempdir(), "no-such-environment"), julia = first(Base.julia_cmd().exec), boot = "none.jl")
session = CLS.Session(CLS.APP, "s"^32, "127.0.0.1", 0, cp0, CLS.StudioState(), false)
router = HTTP.Router()
CLS.core_routes!(router, session)
CLS.chemistrylab_routes!(router, session)
server = CLS.serve(session, router)
base = "http://127.0.0.1:$(session.port)"
token = ["X-Studio-Token" => session.token]
json_headers = vcat(token, ["Content-Type" => "application/json"])
get(path; headers = token) = HTTP.get(base * path, headers; status_exception = false)
post(path, body; headers = json_headers) = HTTP.post(base * path, headers, JSON.json(body); status_exception = false)
payload(r) = JSON.parse(String(r.body); dicttype = Dict{String, Any})

@testset "static files" begin
    r = get("/"; headers = Pair{String, String}[])
    @test r.status == 200
    @test occursin("ChemistryLab Studio", String(r.body))
    @test occursin("script-src 'self'", HTTP.header(r, "Content-Security-Policy"))
    @test get("/app.js"; headers = Pair{String, String}[]).status == 200
    @test get("/i18n/fr.json"; headers = Pair{String, String}[]).status == 200
    @test get("/nothing.js"; headers = Pair{String, String}[]).status == 404
    @test get("/%2e%2e/%2e%2e/etc/passwd"; headers = Pair{String, String}[]).status == 404
    @test CLS.static_file(CLS.WEB_DIR, "/../src/api.jl") === nothing
    @test CLS.static_file(CLS.WEB_DIR, "/a%00b") === nothing
    @test HTTP.post(base * "/", [], ""; status_exception = false).status == 405
end

@testset "guard" begin
    @test get("/api/state"; headers = Pair{String, String}[]).status == 403
    @test get("/api/state"; headers = ["X-Studio-Token" => "t"^32]).status == 403
    @test get("/api/state"; headers = vcat(token, ["Origin" => "http://evil.example"])).status == 403
    @test get("/api/state"; headers = vcat(token, ["Origin" => "http://127.0.0.1:$(session.port)"])).status == 200
    @test get("/api/state"; headers = vcat(token, ["Host" => "evil.example:$(session.port)"])).status == 403
    @test HTTP.post(base * "/api/run", vcat(token, ["Content-Type" => "text/plain"]), "{}"; status_exception = false).status == 403
    @test CLS.secure_equal("abc", "abc")
    @test !CLS.secure_equal("abc", "abd")
    @test !CLS.secure_equal("abc", "ab")
end

@testset "model and script" begin
    r = get("/api/state")
    @test r.status == 200
    st = payload(r)
    @test st["path"] === nothing
    @test isempty(st["problems"])
    @test startswith(st["script"], "# ===")
    model = st["model"]
    model["mix"]["w_b"] = 0.4
    st = payload(post("/api/model", Dict("model" => model)))
    @test occursin("w_b = 0.4 ", st["script"])
    r = post("/api/model", Dict("model" => Dict("mix" => Dict("kind" => "mix", "w_b" => "x"))))
    @test r.status == 400 && occursin("w_b", payload(r)["error"])
    @test post("/api/model", Dict("nothing" => 1)).status == 400
    model["mix"]["temperature"] = 150
    st = payload(post("/api/model", Dict("model" => model)))
    @test [p["code"] for p in st["problems"]] == ["out_of_range"]
    # a model that cannot run is neither run nor diagnosed
    r = post("/api/run", Dict())
    @test r.status == 422 && payload(r)["problems"][1]["code"] == "out_of_range"
    @test payload(post("/api/diagnose", Dict()))["pending"]
    model["mix"]["temperature"] = 20
    post("/api/model", Dict("model" => model))
    @test occursin("not ready", payload(post("/api/diagnose", Dict()))["error"])
    @test get("/api/species?database=cemdata18-thermofun.json").status == 503
    r = payload(post("/api/run", Dict()))
    @test r["ok"] == false && occursin("not ready", r["error"])
end

@testset "library and catalog" begin
    entries = payload(get("/api/library"))
    @test [e["id"] for e in entries][1:2] == ["cem1-52.5n-gruyaert2010", "portland-rilem-pc"]
    @test all(e -> haskey(e, "sources") && !isempty(e["sources"]), entries)
    st = payload(post("/api/library", Dict("id" => "portland-rilem-pc")))
    @test st["model"]["mix"]["w_b"] == 0.45
    @test post("/api/library", Dict("id" => "nothing")).status == 400
    @test payload(get("/api/catalog")) == Dict{String, Any}("pending" => true)
end

@testset "files" begin
    mktempdir() do dir
        @test post("/api/save", Dict()).status == 400
        @test post("/api/save", Dict("path" => joinpath(dir, "missing", "x"))).status == 400
        st = payload(post("/api/save", Dict("path" => joinpath(dir, "paste"))))
        file = joinpath(dir, "paste.jl")
        @test st["path"] == file && isfile(file)
        @test !isfile(file * ".part")
        @test post("/api/save", Dict()).status == 200          # the bound file
        other = joinpath(dir, "other.jl")
        write(other, "x = 1\n")
        r = post("/api/save", Dict("path" => other))
        @test r.status == 409 && payload(r)["exists"]
        @test read(other, String) == "x = 1\n"
        @test post("/api/save", Dict("path" => other, "overwrite" => true)).status == 200
        @test CLS.read_script(read(other, String)).model !== nothing

        st = payload(post("/api/open", Dict("path" => file)))
        @test st["path"] == file && st["message"] == ""
        write(file, replace(read(file, String), "w_b = 0.45 " => "w_b = 0.6 "))
        st = payload(post("/api/open", Dict("path" => file)))
        @test st["path"] === nothing && occursin("edited by hand", st["message"])
        plain = joinpath(dir, "plain.jl")
        write(plain, "x = 1\n")
        r = post("/api/open", Dict("path" => plain))
        @test r.status == 422 && occursin("not written", payload(r)["error"])
        damaged = joinpath(dir, "damaged.jl")
        write(damaged, "x = 1\n" * CLS.MODEL_TAG * " sha256=0\n{\"mix\": {\"kind\": \"mix\", \"w_b\": \"x\"}}\n=#\n")
        @test post("/api/open", Dict("path" => damaged)).status == 422
        @test post("/api/open", Dict("path" => joinpath(dir, "absent.jl"))).status == 400
        @test post("/api/open", Dict()).status == 400

        mkdir(joinpath(dir, "sub"))
        write(joinpath(dir, ".hidden.jl"), "")
        write(joinpath(dir, "notes.txt"), "")
        listing = payload(get("/api/browse?path=" * HTTP.escapeuri(dir)))
        @test listing["path"] == dir
        @test listing["folders"] == ["sub"]
        @test sort(listing["files"]) == ["damaged.jl", "other.jl", "paste.jl", "plain.jl"]
        @test payload(get("/api/browse"))["path"] == homedir()
        # a folder the user may not read lists as empty
        if Sys.isunix() && !(Base.Libc.getuid() == 0)
            locked = joinpath(dir, "locked")
            mkdir(locked)
            chmod(locked, 0o000)
            listing = payload(get("/api/browse?path=" * HTTP.escapeuri(locked)))
            @test isempty(listing["folders"]) && isempty(listing["files"])
            chmod(locked, 0o700)
        end
        @test payload(get("/api/browse?path=" * HTTP.escapeuri(joinpath(dir, "absent"))))["path"] == homedir()
    end
end

@testset "pasted rows and the schema" begin
    r = payload(post("/api/rows/parse", Dict("text" => "CaO;63,7\nSiO2;19,3\n", "unit" => "%")))
    @test [x["formula"] for x in r["rows"]] == ["CaO", "SiO2"] && r["decimal_comma"]
    @test post("/api/rows/parse", Dict("text" => 3)).status == 400
    @test post("/api/rows/parse", Dict("text" => "x 1", "unit" => "furlong")).status == 400
    @test payload(get("/api/schema")) isa Dict
end

@testset "status, stop, restart" begin
    s = payload(get("/api/status"))
    @test s["state"] == "stopped"
    @test s["app"] == "ChemistryLab Studio"
    @test s["memory_gb"] > 0
    @test payload(post("/api/compute/stop", Dict()))["stopping"] == false
    # restarting boots in the background; this environment lacks the packages
    @test payload(post("/api/compute/restart", Dict()))["restarting"]
    for _ in 1:200
        payload(get("/api/status"))["state"] == "failed" && break
        sleep(0.05)
    end
    s = payload(get("/api/status"))
    @test s["state"] == "failed" && occursin("lacks ChemistryLab", s["message"])
end

@testset "an app is built by keyword" begin
    app = CLS.StudioApp(;
        name = "Test Studio", version = v"0.0.1", default_port = 9999, project_variable = "X",
        default_project = "@x", environment = CLS.EnvironmentSpec(Pair{String, String}[]),
        boot = "none.jl", web_dir = CLS.WEB_DIR, new_state = () -> nothing,
        routes! = (r, s) -> nothing, warm_up! = s -> nothing, check = :(String[]),
    )
    @test app.name == "Test Studio" && app.default_port == 9999
end

@testset "a failing route answers 500" begin
    r = HTTP.Router()
    HTTP.register!(r, "GET", "/api/boom", _ -> error("boom"))
    resp = @test_logs (:error,) match_mode = :any CLS.answer(r, HTTP.Request("GET", "/api/boom"))
    @test resp.status == 500 && occursin("boom", String(resp.body))
    @test_throws ArgumentError CLS.read_json(HTTP.Request("POST", "/api/x", [], "[1, 2]"))
    @test CLS.read_json(HTTP.Request("POST", "/api/x", [], "")) == Dict{String, Any}()
end

HTTP.forceclose(server)
