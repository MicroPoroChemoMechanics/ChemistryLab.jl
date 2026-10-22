# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2026 Jean-François Barthélémy (Cerema, UMR MCD)

# Included by the compute process at start, in its `Main`. The server never
# loads ChemistryLab: everything that reads ChemistryLab's objects lives here.

module StudioWorker

using ChemistryLab, OptimaSolver, DynamicQuantities
using Logging: ConsoleLogger, with_logger
using TOML: TOML

include(joinpath(@__DIR__, "..", "shared", "releases.jl"))
include(joinpath(@__DIR__, "catalog.jl"))
include(joinpath(@__DIR__, "run.jl"))
include(joinpath(@__DIR__, "extract.jl"))

end
