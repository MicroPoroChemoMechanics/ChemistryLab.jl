# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2026 Jean-François Barthélémy (Cerema, UMR MCD)

using ChemistryLabStudio
using Test

const CLS = ChemistryLabStudio

@testset "ChemistryLabStudio" begin
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
    @testset "Aqua" include("test_aqua.jl")
end
