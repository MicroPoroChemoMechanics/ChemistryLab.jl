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
using Logging

const QUICK = "--quick" in ARGS
const OUT = joinpath(@__DIR__, "..", "docs", "src", "assets", "precomputed")

# Both at top level, and NOT inside the functions that use them: new methods
# defined by an `include` during a call are not visible to that call (Julia's
# world age), so a nested include compiles fine and fails at run time with a
# `no method matching` for something the file plainly defines.
isdefined(Main, :run_ionic_hydration) ||
    include(joinpath(@__DIR__, "ionic_hydration.jl"))
isdefined(Main, :forward_Q) ||
    include(joinpath(@__DIR__, "hydration_calibration.jl"))

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
function provenance(io, what; instants = N_INSTANTS, window = TEND)
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
    println(io, "# instants: ", instants, " log-spaced over ",
        round(window / 86400; digits = 2), " days")
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

    # How many instants are PROVED optimal, not merely converged. The number
    # belongs in the file: "80 instants, 79 of them proved" is a different claim
    # from "80 instants", and a stored result is worth what its provenance is.
    #
    # Taken from `speciated_states` itself, which already decides it per instant
    # and says so in a warning. Re-deriving it here would be a second opinion
    # with no authority over the first, and an earlier attempt at exactly that
    # reported 0 of 80 where the replay had certified 79.
    # Captured by listening to `speciated_states`, which already decides it per
    # instant and warns. The warning goes to a buffer so it does not also land in
    # the console twice; the text is then echoed, because a silenced warning is
    # worse than a noisy one.
    t1 = time()
    logbuf = IOBuffer()
    states = Logging.with_logger(Logging.SimpleLogger(logbuf, Logging.Info)) do
        speciated_states(run.sol, run.kp; times = TIMES)
    end
    logged = String(take!(logbuf))
    m = match(r"(\d+) of (\d+) replayed instants could not be certified", logged)
    n_uncertified = m === nothing ? 0 : parse(Int, m.captures[1])
    n_certified = length(TIMES) - n_uncertified
    @printf("  %d of %d instants proved optimal\n", n_certified, length(TIMES))
    n_uncertified > 0 && println("  (", strip(first(split(logged, "\n"))), ")")
    flush(stdout)
    times, fracs, pore_pH, poro = ionic_phase_history(run, TIMES; states = states)
    @printf("  %d instants certified in %.0f s\n", length(TIMES), time() - t1)
    flush(stdout)

    groups = sort(collect(union(keys.(fracs)...)))


    open(joinpath(OUT, "$(tag)_phases.csv"), "w") do io
        provenance(io, "$label -- volume fractions of the phase families, pH, porosity")
        println(io, "# certified: ", n_certified, " of ", length(TIMES),
            " replayed instants proved optimal against the KKT conditions")
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

"""
    write_calibration()

The coupled curves of `examples/hydration_calibration.md`: the published
Parrott-Killoh parameters and the calibrated ones, each on the calibration
target and on the holdout record.

Five coupled forward solves, some six minutes each. The fitted parameter vector
itself is NOT refitted here — it is `CALIBRATED_THETA`, already a stored
constant of `hydration_calibration.jl`, obtained by the optimization that script
performs and that the page has never run at build time either.
"""
function write_calibration()
    mkpath(OUT)
    θ0 = prior_vector()
    θ̂ = CALIBRATED_THETA

    for (tag, record, label) in (
            ("calibration_target", CEM_I_TARGET,
                "CEM I 52.5 R Cizkovice, w/b 0.50 -- the calibration target"),
            ("calibration_holdout", CEM_I_HOLDOUT,
                "CEM I 52.5 R Ladce, w/b 0.45 -- the holdout, never fitted"),
        )
        @info "coupled forward solves" tag
        t0 = time()
        data = resample_log(record, N_RESIDUALS_COUPLED)
        Q_prior = forward_Q(θ0, data; mode = :coupled)
        Q_fit = forward_Q(θ̂, data; mode = :coupled)
        @printf("  two coupled solves in %.0f s\n", time() - t0)
        flush(stdout)

        open(joinpath(OUT, "$(tag).csv"), "w") do io
            provenance(
                io, "$label -- measured against published and calibrated";
                instants = length(data.t), window = data.t[end],
            )
            println(io, "# measured data: Smilauer & Reiterman (2025), Zenodo")
            println(io, "#   10.5281/zenodo.15212785, CC-BY-4.0")
            println(io, "# Q_prior: published Parrott-Killoh parameters, untouched")
            println(io, "# Q_fit:   CALIBRATED_THETA of scripts/hydration_calibration.jl")
            println(io, "time_s,Q_measured,Q_prior,Q_fit,Q_depositors_fit")
            for i in eachindex(data.t)
                @printf(io, "%.6e,%.6e,%.6e,%.6e,%.6e\n",
                    data.t[i], data.Q[i], Q_prior[i], Q_fit[i], data.Qref[i])
            end
        end
        @printf("  wrote %s.csv\n", tag)
        flush(stdout)
    end

    # The clinker-sensitivity of section 5: the same fit with the alite content
    # moved by plus and minus twenty percent, far more than a Bogue or a QXRD
    # analysis is uncertain by. Two more full coupled runs.
    @info "clinker sensitivity"
    t0 = time()
    data = resample_log(CEM_I_TARGET, N_RESIDUALS_COUPLED)
    open(joinpath(OUT, "calibration_sensitivity.csv"), "w") do io
        provenance(
            io, "sensitivity of the fit to the alite content";
            instants = length(data.t), window = data.t[end],
        )
        println(io, "# the fitted rate constants are held at CALIBRATED_THETA;")
        println(io, "# only the clinker composition moves, the rest rescaled to close")
        println(io, "delta_C3S,C3S,Q_end_J_per_g,RMSE_J_per_g")
        for δ in (0.0, -0.20, 0.20)
            c = CALIB_CLINKER.C3S * (1 + δ)
            scale = (1 - c) / (1 - CALIB_CLINKER.C3S)
            clinker = (
                C3S = c, C2S = CALIB_CLINKER.C2S * scale,
                C3A = CALIB_CLINKER.C3A * scale, C4AF = CALIB_CLINKER.C4AF * scale,
            )
            run = run_ionic_hydration(;
                wb = data.meta.wb, clinker, gypsum = CALIB_GYPSUM,
                filler = CALIB_FILLER, blaine = data.meta.blaine * u"m^2/kg",
                tend = data.t[end], pk_params = apply_parameters(θ̂),
            )
            _, Q, _ = heat_release(run.sol, run.kp; times = data.t)
            rmse = sqrt(sum(abs2, Q ./ 1000 .- data.Q) / length(data.Q))
            @printf(io, "%.2f,%.4f,%.4f,%.4f\n", δ, c, Q[end] / 1000, rmse)
            @printf("  delta %+.0f %% done\n", 100δ)
            flush(stdout)
        end
    end
    @printf("  wrote calibration_sensitivity.csv in %.0f s\n", time() - t0)
    return nothing
end

"""
    main()

Produce every file the documentation reads, or the subset named on the command
line: `--only=ionic`, `--only=calibration`. Selecting is not an optimization for
its own sake — each group costs minutes, and re-running the ones already in hand
wastes them for nothing.
"""
function main()
    mkpath(OUT)
    only = ""
    for a in ARGS
        startswith(a, "--only=") && (only = split(a, "=")[2])
    end
    want(group) = isempty(only) || only == group

    @info "precomputing the documentation's coupled runs" N_INSTANTS QUICK only
    t0 = time()
    if want("ionic")
        write_case("ionic_opc", 0.035, "CEM I with 3.5 % limestone filler")
        write_case("ionic_nolimestone", 0.0, "the same paste with the limestone removed")
    end
    want("calibration") && write_calibration()
    @printf("\nwritten to %s in %.0f s\n", OUT, time() - t0)
    return nothing
end

main()
