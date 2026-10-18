# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2026 Jean-François Barthélémy (Cerema, UMR MCD)

# The Julia environment the compute process runs in. By default the studio
# manages a shared environment of its own (`@chemistrylab-studio`): created at
# the first launch, with the packages a calculation needs and the `[compat]`
# bounds the studio is tested against. The same environment runs the scripts
# the studio writes, outside it. A project the user names instead is never
# modified: a missing package is reported with the command that adds it.

"""
    EnvironmentSpec(packages)

The packages the compute process needs, as `name => compat` pairs, the bound
being the range the studio is tested against.
"""
struct EnvironmentSpec
    packages::Vector{Pair{String, String}}
end

"""
    environment_dir(project) -> String

The directory of an environment named as for `julia --project`: a path, or
`@name` for a shared environment of the first depot.
"""
function environment_dir(project::AbstractString)
    startswith(project, "@") || return abspath(project)
    return joinpath(first(DEPOT_PATH), "environments", project[2:end])
end

"The packages of `spec` that the environment at `dir` does not list."
function missing_packages(dir::AbstractString, spec::EnvironmentSpec)
    names = first.(spec.packages)
    file = joinpath(dir, "Project.toml")
    isfile(file) || return names
    deps = get(TOML.parsefile(file), "deps", Dict{String, Any}())
    return filter(p -> !haskey(deps, p), names)
end

"""
    environment_script(project, spec; update) -> String

The code a separate Julia process runs to create, complete or update the
environment `project`, then precompile it.
"""
function environment_script(project::AbstractString, spec::EnvironmentSpec; update::Bool)
    shared = startswith(project, "@")
    target = shared ? project[2:end] : environment_dir(project)
    adds = join(("Pkg.PackageSpec(name = $(repr(p)), version = $(repr(c)))" for (p, c) in spec.packages), ", ")
    bounds = join(("Pkg.compat($(repr(p)), $(repr(c)))" for (p, c) in spec.packages), "\n")
    return """
    import Pkg
    Pkg.activate($(repr(target)); shared = $shared)
    Pkg.add([$adds])
    $bounds
    $(update ? "Pkg.update()" : "")
    Pkg.precompile()
    """
end

"""
    run_logged(cmd, log) -> Nothing

Run `cmd`, passing each line it prints to `log`; an error naming the last lines
when it fails.
"""
function run_logged(cmd::Cmd, log::F) where {F}
    out = Pipe()
    proc = run(pipeline(cmd; stdout = out, stderr = out); wait = false)
    close(out.in)
    tail = String[]
    for line in eachline(out)
        log(line)
        push!(tail, line)
        length(tail) > 20 && popfirst!(tail)
    end
    wait(proc)
    success(proc) || error("the command exited with status $(proc.exitcode):\n" * join(tail, "\n"))
    return nothing
end

"""
    prepare_environment(project, spec; managed, update, julia, log) -> String

Make sure the environment `project` holds the packages of `spec`, and return
its directory. A `managed` environment (the studio's own) is created or
completed in a separate process, its output passed to `log`; it is left as it
is when complete, unless `update`. An environment the user named is only
checked.
"""
function prepare_environment(
        project::AbstractString, spec::EnvironmentSpec;
        managed::Bool, update::Bool, julia::AbstractString, log::F,
    ) where {F}
    dir = environment_dir(project)
    lacking = missing_packages(dir, spec)
    if !managed
        isempty(lacking) && return dir
        throw(
            ArgumentError(
                "the environment $project lacks $(join(lacking, ", ")). Add them " *
                    "(`] activate $project`, then `] add $(join(lacking, ' '))`), or start " *
                    "without `project` to let the studio manage an environment of its own."
            )
        )
    end
    complete = isempty(lacking) && isfile(joinpath(dir, "Manifest.toml"))
    (complete && !update) && return dir
    log(
        update ? "Updating the compute environment $project." :
            "Preparing the compute environment $project: this happens once and takes several minutes."
    )
    script = environment_script(project, spec; update)
    cmd = addenv(
        `$julia --startup-file=no -e $script`,
        "JULIA_NUM_PRECOMPILE_TASKS" => "1", "OPENBLAS_NUM_THREADS" => "1",
    )
    run_logged(cmd, log)
    return dir
end
