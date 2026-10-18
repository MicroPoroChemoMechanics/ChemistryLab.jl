# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2026 Jean-François Barthélémy (Cerema, UMR MCD)

# The tests of ChemistryLab Studio (src/studio/), in a module of their own so
# that their helpers stay apart from those of the package's other tests. A unit
# of the `core` part of the suite; the Studio workflow runs them alone too, with
# a compute process.

module StudioTests

using ChemistryLab, Test
using ChemistryLab.Studio: HTTP, JSON, Malt, Dates

const CLS = ChemistryLab.Studio

@testset "Studio" begin
    @testset "rows" include("test_rows.jl")
    @testset "schema" include("test_schema.jl")
    @testset "words of the page" include("test_i18n.jl")
    @testset "model" include("test_model.jl")
    @testset "codegen" include("test_codegen.jl")
    @testset "readback" include("test_readback.jl")
    @testset "memory" include("test_memory.jl")
    @testset "environment" include("test_environment.jl")
    @testset "server" include("test_server.jl")
    @testset "launcher" include("test_launcher.jl")
    # The compute process needs an environment holding ChemistryLab and its
    # solver; CI builds one, and a local run names it.
    if haskey(ENV, "CLSTUDIO_WORKER_PROJECT")
        @testset "compute" include("test_compute.jl")
    else
        @info "Compute-process tests skipped: set CLSTUDIO_WORKER_PROJECT to an environment holding ChemistryLab."
    end
end

end
