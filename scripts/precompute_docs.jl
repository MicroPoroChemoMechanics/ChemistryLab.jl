# SPDX-License-Identifier: LGPL-2.1-or-later
# =============================================================================
#  precompute_docs.jl — the expensive trajectories, computed once
#
#  A coupled hydration run costs a Gibbs minimization per accepted ODE step,
#  and replaying it certifies one equilibrium per reported instant. Measured on
#  a 2-core machine: 198 s for one `run_ionic_hydration` to 28 days, and 364 s
#  for one coupled forward solve replayed on twenty instants. The documentation
#  called for a dozen of those, which put the build past three hours.
#
#  So they are computed HERE, once, and the pages plot the result. Nothing is
#  approximated to make that possible: the runs below are at a FINER sampling
#  than the pages used when they computed it themselves (80 instants against
#  40), because the cost is paid once rather than at every build.
#
#  Every number the documentation shows therefore remains a computed number,
#  traceable to this script, to the package version recorded in the header of
#  each file, and to the commit it ran at. Nothing here is fitted by hand,
#  adjusted, or transcribed from a paper.
#
#  Usage:
#      julia --project=docs scripts/precompute_docs.jl            # everything
#      julia --project=docs scripts/precompute_docs.jl --quick    # coarse, for
#                                                                 # checking the
#                                                                 # plumbing only
# =============================================================================

using ChemistryLab
using DynamicQuantities
using OptimaSolver
using OrdinaryDiffEq
using Optimization, OptimizationOptimJL
using OrderedCollections
using Printf
using Dates

const QUICK = "--quick" in ARGS
const OUT = joinpath(@__DIR__, "..", "docs", "src", "assets", "precomputed")

isdefined(Main, :run_ionic_hydration) ||
    include(joinpath(@__DIR__, "ionic_hydration.jl"))

"""
Sampling of the reported instants.

80 log-spaced points over 28 days, against the 40 the pages used when they
computed this themselves. Log-spaced because hydration is a decade-wise process:
a uniform grid spends its resolution on the plateau and none on the acceleration
peak.
"""
const N_INSTANTS = QUICK ? 8 : 80

const TEND = 28 * 86400.0
const TIMES = 10 .^ range(log10(0.05 * 86400), log10(TEND); length = N_INSTANTS)

# The CEM I 52.5 N of Lavergne et al. (2018), Table 9 — the formulation the
# documentation already runs, unchanged here.
const CLINKER = (C3S = 0.65, C2S = 0.11, C3A = 0.11, C4AF = 0.08)

"""
    provenance(io, what)

The header every file carries: what produced it, when, with which package, and
at which commit. A precomputed number is only as good as its traceability.
"""
function provenance(io, what)
    commit = try
        readchomp(`git -C $(pkgdir(ChemistryLab)) rev-parse --short HEAD`)
    catch
        "unknown"
    end
    ver = try
        string(pkgversion(ChemistryLab))
    catch
        "unknown"
    end
    println(io, "# ", what)
    println(io, "# produced by scripts/precompute_docs.jl -- do not edit by hand")
    println(io, "# ChemistryLab version: ", ver)
    println(io, "# commit: ", commit)
    println(io, "# generated: ", Dates.format(now(), "yyyy-mm-dd HH:MM"))
    println(io, "# instants: ", N_INSTANTS, " log-spaced over ", TEND / 86400, " days")
    println(io, "# w/b 0.50, Blaine default, 20 C, 1 kg of binder")
    println(io, "# clinker (mass fractions): ", CLINKER)
    QUICK && println(io, "# !! QUICK MODE -- coarse sampling, not for publication !!")
    return nothing
end

"""
    write_case(tag, filler, label)

Run one paste, replay it on `TIMES`, and write two files: the phase history
(volume fractions of each family, pH, porosity) and the calorimetry.
"""
function write_case(tag, filler, label)
    mkpath(OUT)          # cheap, and survives the directory being moved under us
    @info "running" tag filler
    t0 = time()
    run = run_ionic_hydration(;
        wb = 0.5, clinker = CLINKER, gypsum = 0.046, filler = filler, tend = TEND,
    )
    @printf("  integrated in %.0f s, %d accepted steps, retcode %s\n",
        time() - t0, length(run.sol.t), run.sol.retcode)
    flush(stdout)

    t1 = time()
    states = speciated_states(run.sol, run.kp; times = TIMES)
    times, fracs, pore_pH, poro = ionic_phase_history(run, TIMES; states = states)
    @printf("  %d instants certified in %.0f s\n", length(TIMES), time() - t1)
    flush(stdout)

    groups = sort(collect(union(keys.(fracs)...)))

    open(joinpath(OUT, "$(tag)_phases.csv"), "w") do io
        provenance(io, "$label -- volume fractions of the phase families, pH, porosity")
        println(io, "time_s,", join(groups, ","), ",pore_pH,poro_liquid,poro_void,poro_total")
        for (i, t) in enumerate(times)
            vals = [get(fracs[i], g, 0.0) for g in groups]
            @printf(io, "%.6e,%s,%.6f,%.6f,%.6f,%.6f\n", t,
                join((@sprintf("%.6e", v) for v in vals), ","),
                pore_pH[i], poro[i].liquid, poro[i].void, poro[i].total)
        end
    end

    t2 = time()
    tcal, Q, qd = heat_release(run.sol, run.kp; times = TIMES, states = states)
    @printf("  calorimetry in %.0f s\n", time() - t2)

    open(joinpath(OUT, "$(tag)_heat.csv"), "w") do io
        provenance(io, "$label -- isothermal calorimetry at 20 C, per gram of binder")
        println(io, "time_s,Q_J_per_g,heat_flow_W_per_g")
        for i in eachindex(tcal)
            @printf(io, "%.6e,%.6e,%.6e\n", tcal[i], Q[i] / 1000, qd[i] / 1000)
        end
    end

    @printf("  wrote %s_phases.csv and %s_heat.csv  (total %.0f s)\n",
        tag, tag, time() - t0)
    flush(stdout)
    return nothing
end

function main()
    mkpath(OUT)
    @info "precomputing the documentation's coupled runs" N_INSTANTS QUICK
    t0 = time()
    write_case("ionic_opc", 0.035, "CEM I with 3.5 % limestone filler")
    write_case("ionic_nolimestone", 0.0, "the same paste with the limestone removed")
    @printf("\nall cases written to %s in %.0f s\n", OUT, time() - t0)
    return nothing
end

main()
