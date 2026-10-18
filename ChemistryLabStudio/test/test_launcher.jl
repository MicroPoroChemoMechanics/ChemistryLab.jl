# Starting a studio: the address and the remote-session hint, the browser
# command, the port fallback, a handle that stops cleanly, and a blocking run
# stopped by an interrupt.

using HTTP

@testset "remote session and browser" begin
    withenv("SSH_CONNECTION" => "1 2 3 4", "VSCODE_IPC_HOOK_CLI" => nothing, "BROWSER" => nothing) do
        @test CLS.is_remote_session()
        lines = CLS.banner("Studio", "http://127.0.0.1:9000/", 9000)
        @test length(lines) == 2 && occursin("ssh -L 9000:127.0.0.1:9000", lines[2])
        @test CLS.browser_command("http://x/") === nothing
        @test !CLS.open_browser("http://x/")
    end
    withenv("SSH_CONNECTION" => nothing, "VSCODE_IPC_HOOK_CLI" => nothing, "BROWSER" => "mybrowser --new-tab") do
        @test !CLS.is_remote_session()
        @test length(CLS.banner("Studio", "http://127.0.0.1:9000/", 9000)) == 1
        @test CLS.browser_command("http://x/").exec == ["mybrowser", "--new-tab", "http://x/"]
    end
    withenv("SSH_CONNECTION" => nothing, "VSCODE_IPC_HOOK_CLI" => nothing, "BROWSER" => nothing) do
        cmd = CLS.browser_command("http://x/")
        @test last(cmd.exec) == "http://x/"
        @test first(cmd.exec) == (Sys.isapple() ? "open" : Sys.iswindows() ? "rundll32" : "xdg-open")
    end
    if Sys.isunix()
        withenv("BROWSER" => "true") do
            @test CLS.open_browser("http://x/")
        end
        withenv("BROWSER" => "/no/such/browser") do
            @test !CLS.open_browser("http://x/")
        end
    end
    @test CLS.url_host("0.0.0.0") == "localhost"
    @test CLS.url_host("::1") == "[::1]"
    @test CLS.url_host("127.0.0.1") == "127.0.0.1"
end

lacking = mktempdir()     # an environment without the packages: the boot fails at once
newsession(port) = CLS.Session(
    CLS.APP, "p"^32, "127.0.0.1", port,
    CLS.ComputeProcess(; project = lacking, julia = first(Base.julia_cmd().exec), boot = "none.jl"),
    CLS.StudioState(), false,
)

@testset "port taken" begin
    blocker = HTTP.serve!(_ -> HTTP.Response(200), "127.0.0.1", 0)
    taken = parse(Int, last(split(HTTP.server_addr(blocker), ':')))
    s = newsession(taken)
    e = try
        CLS.serve(s, HTTP.Router())
    catch err
        err
    end
    @test e isa ArgumentError && occursin("`port = $(taken + 1)`", e.msg)
    s = newsession(taken)
    server = redirect_stdout(devnull) do
        CLS.serve(s, HTTP.Router(); tries = 5)
    end
    @test s.port > taken
    HTTP.forceclose(server)
    HTTP.forceclose(blocker)
end

@testset "a studio that cannot boot" begin
    h = redirect_stdout(devnull) do
        withenv("CLSTUDIO_PROJECT" => nothing) do
            clstudio(; port = 0, no_browser = true, wait = false, project = lacking)
        end
    end
    @test occursin("?token=", CLS.studio_url(h))
    shown = sprint(show, h)
    @test occursin("ChemistryLab Studio", shown) && !occursin("token", shown)
    for _ in 1:200
        CLS.status(h.session.compute).state == "failed" && break
        sleep(0.05)
    end
    @test CLS.status(h.session.compute).state == "failed"
    @test occursin("lacks", CLS.status(h.session.compute).message)
    # the page loads, with the token of the address
    port = h.session.port
    token = last(split(CLS.studio_url(h), "token="))
    @test HTTP.get("http://127.0.0.1:$port/api/state", ["X-Studio-Token" => token]).status == 200
    close(h)
    close(h)                        # a second close does nothing
    @test h.closed
    @test occursin("closed", sprint(show, h))
end

@testset "a user's environment is checked, not modified" begin
    e = try
        redirect_stdout(devnull) do
            clstudio(; check = true, project = lacking)
        end
    catch err
        err
    end
    @test e isa ArgumentError && occursin("lacks ChemistryLab", e.msg)
    # the environment variable names the environment just as the keyword does
    e = try
        redirect_stdout(devnull) do
            withenv(() -> clstudio(; check = true), "CLSTUDIO_PROJECT" => lacking)
        end
    catch err
        err
    end
    @test e isa ArgumentError && occursin(lacking, e.msg)
end

@testset "a blocking run stops on an interrupt" begin
    result = redirect_stdout(devnull) do
        task = @async clstudio(; port = 0, no_browser = true, project = lacking)
        sleep(2)
        schedule(task, InterruptException(); error = true)
        fetch(task)
    end
    @test result === nothing
end

@testset "a boot that prepares its environment" begin
    # the whole boot on a managed environment holding a standard library only,
    # and a boot file that needs nothing: no download, no ChemistryLab
    mktempdir() do dir
        env = joinpath(dir, "env")
        boot = joinpath(dir, "boot.jl")
        write(boot, "const BOOTED = true\n")
        app = CLS.StudioApp(;
            name = "Test Studio", version = v"0.0.1", default_port = 9999, project_variable = "X",
            default_project = "@x", environment = CLS.EnvironmentSpec(["Printf" => "1"]), boot,
            web_dir = CLS.WEB_DIR, new_state = () -> nothing, routes! = (r, s) -> nothing,
            warm_up! = s -> CLS.set_state!(s.compute, :ready, "warmed"), check = :(String[]),
        )
        cp = CLS.ComputeProcess(; project = env, julia = first(Base.julia_cmd().exec), boot)
        session = CLS.Session(app, "b"^32, "127.0.0.1", 0, cp, nothing, true)
        CLS.boot!(session)
        s = CLS.status(cp)
        @test s.state == "ready" && s.message == "warmed"
        @test any(occursin("Preparing the compute environment", l) for l in s.log)
        @test CLS.call!(cp, :(Main.BOOTED)).value === true
        CLS.shutdown!(cp)
        @test CLS.status(cp).state == "stopped"
    end
end
