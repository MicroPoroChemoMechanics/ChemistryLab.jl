# The compute environment: where it lives, what it lacks, the script that
# builds it, and what is done (or refused) with a user's own environment.

spec = CLS.EnvironmentSpec(["Alpha" => "0.3", "Beta" => "1"])

@test CLS.environment_dir("@chemistrylab-studio") == joinpath(first(DEPOT_PATH), "environments", "chemistrylab-studio")
@test CLS.environment_dir("some/dir") == abspath("some/dir")

script = CLS.environment_script("@studio-test", spec; update = false)
@test !any(a -> a isa Expr && a.head in (:error, :incomplete), Meta.parseall(script).args)
@test occursin("Pkg.activate(\"studio-test\"; shared = true)", script)
@test occursin("Pkg.PackageSpec(name = \"Alpha\", version = \"0.3\")", script)
@test occursin("Pkg.compat(\"Beta\", \"1\")", script)
@test !occursin("Pkg.update()", script)
@test occursin("Pkg.update()", CLS.environment_script("/tmp/x", spec; update = true))
@test occursin("shared = false", CLS.environment_script("/tmp/x", spec; update = false))

mktempdir() do dir
    @test CLS.missing_packages(dir, spec) == ["Alpha", "Beta"]
    write(joinpath(dir, "Project.toml"), "[deps]\nAlpha = \"00000000-0000-0000-0000-000000000001\"\n")
    @test CLS.missing_packages(dir, spec) == ["Beta"]
    lines = String[]
    log = line -> push!(lines, line)
    julia = first(Base.julia_cmd().exec)
    # a user's environment is checked, never modified
    e = try
        CLS.prepare_environment(dir, spec; managed = false, update = false, julia, log)
    catch err
        err
    end
    @test e isa ArgumentError && occursin("lacks Beta", e.msg)
    write(joinpath(dir, "Project.toml"), "[deps]\nAlpha = \"00000000-0000-0000-0000-000000000001\"\nBeta = \"00000000-0000-0000-0000-000000000002\"\n")
    @test CLS.prepare_environment(dir, spec; managed = false, update = false, julia, log) == dir
    # a managed environment that is complete is left alone
    write(joinpath(dir, "Manifest.toml"), "")
    @test CLS.prepare_environment(dir, spec; managed = true, update = false, julia, log) == dir
    @test isempty(lines)
end

# a managed environment is created in a process of its own: a standard
# library stands for the packages, so that no download is needed
mktempdir() do dir
    lines = String[]
    julia = first(Base.julia_cmd().exec)
    spec1 = CLS.EnvironmentSpec(["Printf" => "1"])
    @test CLS.prepare_environment(dir, spec1; managed = true, update = false, julia, log = l -> push!(lines, l)) == dir
    @test occursin("Preparing the compute environment", first(lines))
    @test isempty(CLS.missing_packages(dir, spec1))
    @test occursin("Printf = \"1\"", read(joinpath(dir, "Project.toml"), String))
    empty!(lines)
    CLS.prepare_environment(dir, spec1; managed = true, update = true, julia, log = l -> push!(lines, l))
    @test occursin("Updating the compute environment", first(lines))
end

# the output of a command reaches the log line by line; a failure names its last lines
lines = String[]
julia = first(Base.julia_cmd().exec)
CLS.run_logged(`$julia --startup-file=no -e 'println("one"); println("two")'`, l -> push!(lines, l))
@test lines == ["one", "two"]
e = try
    CLS.run_logged(`$julia --startup-file=no -e 'println("last words"); exit(3)'`, l -> nothing)
catch err
    err
end
@test e isa ErrorException && occursin("status 3", e.msg) && occursin("last words", e.msg)
