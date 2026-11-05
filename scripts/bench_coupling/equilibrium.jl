# SPDX-License-Identifier: LGPL-2.1-or-later
#
# What one equilibrium costs at a node of a coupled computation, by entry point.
#
#   julia --project=scripts scripts/bench_coupling/equilibrium.jl [outdir]
#
# A Portland paste (CEM I 52.5 N, Bogue, w/b = 0.5) with and without the C-S-H
# surface of Guo et al., solved for the budget of a neighboring node (every
# element 1 % more), as a transport step hands it over:
#
#   cold: certified             `equilibrate_certified(state; b)`
#   cold: one back end          `solve(EquilibriumSolver)`, the solver built once
#   cold: solver rebuilt        the same, the solver built at every call
#   warm                        `solve_certified(des, (previous,); b)`
#   warm + multipliers          the same, resumed from the previous multipliers
#   warm + multipliers - report the same, `report = false`
#
# Reported: the minimum time of `REPEATS` calls after one warm-up, and what one
# call allocates.

using ChemistryLab, DynamicQuantities, OptimaSolver, Printf
const SciMLBase = ChemistryLab.SciMLBase

const REPEATS = 10
const OUTDIR = get(ARGS, 1, joinpath(@__DIR__, "..", "..", "docs", "src", "assets", "bench_coupling"))

amounts(st) = Float64[ustrip(us"mol", x) for x in st.n]

substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
cs = phase_list_system("Portland paste (Lothenbach and Winnefeld 2006)", substances)
cement = with_extents(
    material_template("CEM I 52.5 N (Gruyaert 2010), Bogue", cs),
    Dict("C3S" => 0.8, "C2S" => 0.45, "C3A" => 0.9, "C4AF" => 0.6),
)
model = cemdata18_activity_model(:KOH)
paste, _ = equilibrate_certified(Recipe(cement => 1.0; w_b = 0.5), cs; model)
surface = with_csh_surface(paste)

function entry_points(st)
    sys = st.system
    des = DualEquilibriumSolver(sys, model)
    es = EquilibriumSolver(sys, model, OptimaOptimizer())
    b = des.A * amounts(st) .* 1.01
    y = Ref{Any}(nothing)
    first(solve_certified(des, (st,); b = des.A * amounts(st), dual_state = y))
    yprev = y[].y
    return [
        "cold: certified" => () -> equilibrate_certified(st; model, b),
        "cold: one back end" => () -> SciMLBase.solve(es, st; b),
        "cold: solver rebuilt" => () -> SciMLBase.solve(EquilibriumSolver(sys, model, OptimaOptimizer()), st; b),
        "warm" => () -> solve_certified(des, (st,); b),
        "warm + multipliers" => () -> solve_certified(des, (st,); b, y0 = yprev),
        "warm + multipliers - report" => () -> solve_certified(des, (st,); b, y0 = yprev, report = false),
    ]
end

mkpath(OUTDIR)
open(joinpath(OUTDIR, "equilibrium.csv"), "w") do io
    println(io, "system,entry_point,time_ms,alloc_MB")
    for (name, st) in (("paste", paste.state), ("paste with C-S-H surface", surface.state))
        for (label, f) in entry_points(st)
            f()
            t = minimum(@elapsed(f()) for _ in 1:REPEATS)
            a = @allocated f()
            @printf(io, "%s,%s,%.3f,%.2f\n", name, label, 1.0e3t, a / 2^20)
            @printf("%-26s %-22s %8.3f ms %8.2f MB\n", name, label, 1.0e3t, a / 2^20)
        end
    end
end
