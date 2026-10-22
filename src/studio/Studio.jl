# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2026 Jean-François Barthélémy (Cerema, UMR MCD)

"""
    ChemistryLab.Studio

ChemistryLab Studio, a local web interface for cementitious materials: build a
calculation with forms, run it, read its plots and tables, and get the Julia
script that does the same on its own. Start it with [`clstudio`](@ref).

The server runs in the user's session; the calculations run in a separate Julia
process, in an environment of their own (`@chemistrylab-studio`) holding this
very version of ChemistryLab, which the first launch prepares. The files under
`core/` hold what any studio needs and never refer to ChemistryLab's
chemistry; `worker/` is the code the compute process includes; `web/` is the
page.
"""
module Studio

using ..ChemistryLab: ChemistryLab
using Dates: Dates
using HTTP: HTTP
using JSON: JSON
using Malt: Malt
using Random: Random
using SHA: SHA
using TOML: TOML

export clstudio

const STUDIO_ROOT = @__DIR__
const WEB_DIR = joinpath(STUDIO_ROOT, "web")
const WORKER_DIR = joinpath(STUDIO_ROOT, "worker")
const LIBRARY_DIR = joinpath(STUDIO_ROOT, "library")

"The version of the studio: that of ChemistryLab, which ships it."
const STUDIO_VERSION = pkgversion(ChemistryLab)

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
include("features/surface.jl")
include("shared/releases.jl")
include("features/tga.jl")
include("model.jl")
include("readback.jl")
include("codegen.jl")
include("api.jl")

# The bounds ChemistryLab declares, read from its own Project.toml. Copied by
# hand, the studio's `"OptimaSolver" => "0.8"` outlived the 0.40 it was written
# for: ChemistryLab 0.41.0, which requires 0.9.3, could not be added beside it,
# and the compute process never started.
const OWN_COMPAT = TOML.parsefile(joinpath(@__DIR__, "..", "..", "Project.toml"))["compat"]

"""
The packages the compute process needs: ChemistryLab itself, at the version
running (the studio and the chemistry it drives are one), and the others with
ChemistryLab's own bounds where it declares one. OrdinaryDiffEq integrates the
kinetic runs (loaded by the first one); Plots draws the figures of the scripts
the studio writes, which run in the same environment, and is no dependency of
ChemistryLab: its bound is the studio's.
"""
const ENVIRONMENT = EnvironmentSpec(
    PackageSource(ChemistryLab),
    [
        "OptimaSolver" => OWN_COMPAT["OptimaSolver"],
        "DynamicQuantities" => OWN_COMPAT["DynamicQuantities"],
        "OrdinaryDiffEq" => OWN_COMPAT["OrdinaryDiffEq"],
        "Plots" => "1",
    ],
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

# The docstring of `clstudio` is in `src/ChemistryLab.jl`, where it is exported:
# the documentation leaves this module out of its docstring check
# (`checkdocs_ignored_modules` in docs/make.jl), which also drops its docstrings
# from every `@docs` block, so that one written here showed as missing on the
# studio's page.
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
