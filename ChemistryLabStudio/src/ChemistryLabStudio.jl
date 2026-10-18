# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2026 Jean-François Barthélémy (Cerema, UMR MCD)

"""
    ChemistryLabStudio

A local web interface for ChemistryLab, for cementitious materials: build a
calculation with forms, run it, read its plots and tables, and get the Julia
script that does the same on its own. Start it with [`clstudio`](@ref).

The server is this package, light; the calculations run in a separate Julia
process, in an environment of their own (`@chemistrylab-studio`) that the first
launch prepares. The files under `src/core/` hold what any studio needs and
never refer to ChemistryLab.
"""
module ChemistryLabStudio

using Dates: Dates
using HTTP: HTTP
using JSON: JSON
using Malt: Malt
using Random: Random
using SHA: SHA
using TOML: TOML

export clstudio

const PACKAGE_ROOT = dirname(@__DIR__)
const WEB_DIR = joinpath(PACKAGE_ROOT, "web")
const WORKER_DIR = joinpath(PACKAGE_ROOT, "worker")
const LIBRARY_DIR = joinpath(PACKAGE_ROOT, "library")

include_dependency(joinpath(PACKAGE_ROOT, "Project.toml"))
const STUDIO_VERSION = VersionNumber(TOML.parsefile(joinpath(PACKAGE_ROOT, "Project.toml"))["version"])

include("core/memory.jl")
include("core/browser.jl")
include("core/environment.jl")
include("core/compute.jl")
include("core/server.jl")
include("core/launcher.jl")

include("rows.jl")
include("schema.jl")
include("features/extents.jl")
include("features/materials.jl")
include("features/database.jl")
include("features/system.jl")
include("features/activity.jl")
include("features/mix.jl")
include("features/constraint.jl")
include("features/calculation.jl")
include("features/kinetics.jl")
include("model.jl")
include("readback.jl")
include("codegen.jl")
include("api.jl")

"""
The packages the compute process needs, with the bounds this studio is tested
against. OrdinaryDiffEq integrates the kinetic runs (loaded by the first one);
Plots draws the figures of the scripts the studio writes, which run in the same
environment.
"""
const ENVIRONMENT = EnvironmentSpec(
    [
        "ChemistryLab" => "0.39.3",
        "OptimaSolver" => "0.8",
        "DynamicQuantities" => "1",
        "OrdinaryDiffEq" => "7",
        "Plots" => "1",
    ]
)

const APP = StudioApp(;
    name = "ChemistryLab Studio",
    version = STUDIO_VERSION,
    default_port = 8766,
    project_variable = "CLSTUDIO_PROJECT",
    default_project = "@chemistrylab-studio",
    environment = ENVIRONMENT,
    boot = joinpath(WORKER_DIR, "boot.jl"),
    web_dir = WEB_DIR,
    new_state = StudioState,
    routes! = chemistrylab_routes!,
    warm_up! = warm_up!,
    check = :(Main.StudioWorker.check_report()),
)

"""
    clstudio(; host = "127.0.0.1", port = nothing, no_browser = false,
               project = nothing, julia = nothing, check = false, wait = true,
               update = false)

Start ChemistryLab Studio, the graphical interface for cementitious materials,
and block until it stops (press Ctrl-C in the REPL). With `wait = false`, return
a handle at once and keep the REPL free: `close(handle)` stops the studio.

The first launch prepares the environment the calculations run in,
`@chemistrylab-studio`, which takes several minutes; the page opens at once and
shows the progress. The first calculation of a session then compiles, which
takes a while depending on the machine; later ones take a fraction of a second.

# Arguments
- `host`/`port` — where to listen; the address printed is
  `http://host:port/?token=…`. Without `port`, the studio takes 8766, or the
  next free port when 8766 is taken (MFH Studio uses 8765, so both run
  together); a `port` given explicitly is used or refused.
- `no_browser` — start without opening the browser (the address is printed).
- `project` — the environment the calculations run in: a folder, or `@name`
  for a shared one. By default the studio's own `@chemistrylab-studio`, which
  it creates and keeps; the `CLSTUDIO_PROJECT` environment variable is honored.
  An environment named this way is only checked, never modified.
- `julia` — the Julia executable of the compute process (default: the running
  one; the `JULIA` environment variable is honored).
- `check` — prepare and check the compute environment, print what it holds, and
  stop instead of serving.
- `wait` — block until the studio stops (default), or return its handle.
- `update` — update the packages of the studio's environment before starting.

On a remote machine (SSH, a VS Code remote window) the studio prints its
address: open it on your own machine, after forwarding the port
(`ssh -L 8766:127.0.0.1:8766 …`; VS Code forwards it by itself).

# Examples
```julia
using ChemistryLabStudio
clstudio()                     # opens the browser on http://127.0.0.1:8766
clstudio(port = 9000)          # another port
clstudio(no_browser = true)    # stay in the terminal
h = clstudio(wait = false)     # keep the REPL; later `close(h)`
clstudio(check = true)         # diagnose the compute environment
```
"""
function clstudio(;
        host::AbstractString = "127.0.0.1",
        port::Union{Nothing, Integer} = nothing,
        no_browser::Bool = false,
        project::Union{Nothing, AbstractString} = nothing,
        julia::Union{Nothing, AbstractString} = nothing,
        check::Bool = false,
        wait::Bool = true,
        update::Bool = false,
    )
    return run_studio(APP; host, port, no_browser, project, julia, check, wait, update)
end

end
