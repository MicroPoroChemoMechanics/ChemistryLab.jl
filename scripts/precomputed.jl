# SPDX-License-Identifier: LGPL-2.1-or-later
#
# The coupled trajectories the documentation plots, COMPUTED BY THE BUILD.
#
# They were stored in CSV files for a while, because a coupled hydration run cost
# tens of minutes and the site called for several. Two solver fixes later a warm
# equilibrium costs 17 ms instead of 583 ms, and computing them is affordable
# again -- which is worth the build time, because a stored result is a claim
# about code that may since have changed, and keeping the two in step needed a
# guard, a procedure and a documented list of traps.
#
# `read_precomputed` keeps its name and its return shape, so the pages that read
# it did not change. What changed is where the numbers come from.
#
# MEMOIZED, and that is not an optimization detail: Documenter runs every
# `@example` block of the whole site in ONE process, and a page asks for the
# phase history and the calorimetry of the same run as two separate tables. Two
# calls, one integration.

using Printf
using DynamicQuantities
using Logging

isdefined(Main, :run_ionic_hydration) ||
    include(joinpath(@__DIR__, "ionic_hydration.jl"))
isdefined(Main, :forward_Q) ||
    include(joinpath(@__DIR__, "hydration_calibration.jl"))

"""
Sampling of the reported instants.

Log-spaced because hydration is a decade-wise process: a uniform grid spends its
resolution on the plateau and none on the acceleration peak.
"""
const N_INSTANTS = 40

const TEND_DOC = 28 * 86400.0
const TIMES_DOC = 10 .^ range(log10(0.05 * 86400), log10(TEND_DOC); length = N_INSTANTS)

# The CEM I 52.5 N of Lavergne et al. (2018), Table 9.
const DOC_CLINKER = (C3S = 0.65, C2S = 0.11, C3A = 0.11, C4AF = 0.08)

# One entry per coupled case the documentation plots.
const IONIC_CASES = Dict(
    "ionic_opc" => (filler = 0.035,
        label = "CEM I with 3.5 % limestone filler"),
    "ionic_nolimestone" => (filler = 0.0,
        label = "the same paste with the limestone removed"),
)

const _CACHE = Dict{String, Any}()

"""
    _provenance(what; instants, window, extra) -> Vector{String}

What a figure can say about where its numbers came from: the package version,
the commit, the sampling and the composition. The same block the stored files
used to carry, now describing a run this build performed.
"""
function _provenance(what; instants, window, extra = String[])
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
    osver = try
        string(pkgversion(OptimaSolver))
    catch
        "unknown"
    end
    return vcat(
        [
            what,
            "computed by this documentation build -- no stored result is read",
            "ChemistryLab version: $ver, OptimaSolver version: $osver",
            "commit: $commit",
            @sprintf("%d log-spaced instants over %.2f days", instants, window / 86400),
            "w/b 0.50, Blaine default, 20 C, 1 kg of binder",
            "clinker (mass fractions): $(DOC_CLINKER)",
        ],
        extra,
    )
end

"""
    _ionic_case(tag) -> NamedTuple

Integrate one paste and replay it on `TIMES_DOC`, returning the phase history
and the calorimetry together — they come from one trajectory and one certified
replay, and splitting them into two calls would integrate twice.
"""
function _ionic_case(tag::AbstractString)
    return get!(_CACHE, "case:" * tag) do
        case = IONIC_CASES[tag]
        run = run_ionic_hydration(;
            wb = 0.5, clinker = DOC_CLINKER, gypsum = 0.046,
            filler = case.filler, tend = TEND_DOC,
        )

        # How many instants are PROVED optimal, not merely converged. Taken from
        # `speciated_states` itself, which decides it per instant and warns;
        # re-deriving it here would be a second opinion with no authority over
        # the first, and an earlier attempt at exactly that reported 0 of 80
        # where the replay had certified 79. The warning goes to a buffer so it
        # does not also reach the page twice, and its text is kept.
        logbuf = IOBuffer()
        states = Logging.with_logger(Logging.SimpleLogger(logbuf, Logging.Info)) do
            speciated_states(run.sol, run.kp; times = TIMES_DOC)
        end
        logged = String(take!(logbuf))
        m = match(r"(\d+) of (\d+) replayed instants could not be certified", logged)
        n_uncertified = m === nothing ? 0 : parse(Int, m.captures[1])
        n_certified = length(TIMES_DOC) - n_uncertified

        times, fracs, pore_pH, poro = ionic_phase_history(run, TIMES_DOC; states = states)
        tcal, Q, qd = heat_release(run.sol, run.kp; times = TIMES_DOC, states = states)
        T_semi = langavant_temperature(tcal, qd ./ 1000, states)

        return (;
            run, states, times, fracs, pore_pH, poro, tcal, Q, qd, T_semi,
            n_certified, label = case.label,
        )
    end
end

_table(columns, header, provenance) = (; columns, provenance, header)

function _ionic_phases_table(tag)
    c = _ionic_case(tag)
    groups = sort(collect(union(keys.(c.fracs)...)))
    columns = Dict{String, Vector{Float64}}("time_s" => collect(c.times))
    for g in groups
        columns[g] = [get(c.fracs[i], g, 0.0) for i in eachindex(c.times)]
    end
    columns["pore_pH"] = collect(c.pore_pH)
    columns["poro_liquid"] = [p.liquid for p in c.poro]
    columns["poro_void"] = [p.void for p in c.poro]
    columns["poro_total"] = [p.total for p in c.poro]
    header = vcat(
        ["time_s"], groups, ["pore_pH", "poro_liquid", "poro_void", "poro_total"]
    )
    prov = _provenance(
        "$(c.label) -- volume fractions of the phase families, pH, porosity";
        instants = length(c.times), window = TEND_DOC,
        extra = [
            "certified: $(c.n_certified) of $(length(TIMES_DOC)) replayed " *
                "instants proved optimal against the KKT conditions",
        ],
    )
    return _table(columns, header, prov)
end

function _ionic_heat_table(tag)
    c = _ionic_case(tag)
    columns = Dict{String, Vector{Float64}}(
        "time_s" => collect(c.tcal),
        "Q_J_per_g" => c.Q ./ 1000,
        "heat_flow_W_per_g" => c.qd ./ 1000,
        "T_semiadiabatic_K" => collect(c.T_semi),
    )
    prov = _provenance(
        "$(c.label) -- isothermal calorimetry at 20 C, per gram of binder";
        instants = length(c.tcal), window = TEND_DOC,
        extra = ["T_semiadiabatic_K: the Langavant cell of NF EN 196-9, same run"],
    )
    return _table(
        columns,
        ["time_s", "Q_J_per_g", "heat_flow_W_per_g", "T_semiadiabatic_K"],
        prov,
    )
end

const CALIBRATION_CASES = Dict(
    "calibration_target" => (record = :target,
        label = "CEM I 52.5 R Cizkovice, w/b 0.50 -- the calibration target"),
    "calibration_holdout" => (record = :holdout,
        label = "CEM I 52.5 R Ladce, w/b 0.45 -- the holdout, never fitted"),
)

function _calibration_table(tag)
    return get!(_CACHE, tag) do
        case = CALIBRATION_CASES[tag]
        record = case.record === :target ? CEM_I_TARGET : CEM_I_HOLDOUT
        data = resample_log(record, N_RESIDUALS_COUPLED)
        Q_prior = forward_Q(prior_vector(), data; mode = :coupled)
        Q_fit = forward_Q(CALIBRATED_THETA, data; mode = :coupled)
        columns = Dict{String, Vector{Float64}}(
            "time_s" => collect(data.t),
            "Q_measured" => collect(data.Q),
            "Q_prior" => collect(Q_prior),
            "Q_fit" => collect(Q_fit),
            "Q_depositors_fit" => collect(data.Qref),
        )
        prov = _provenance(
            "$(case.label) -- measured against published and calibrated";
            instants = length(data.t), window = data.t[end],
            extra = [
                "measured data: Smilauer & Reiterman (2025), Zenodo",
                "  10.5281/zenodo.15212785, CC-BY-4.0",
                "Q_prior: published Parrott-Killoh parameters, untouched",
                "Q_fit:   CALIBRATED_THETA of scripts/hydration_calibration.jl",
            ],
        )
        return _table(
            columns,
            ["time_s", "Q_measured", "Q_prior", "Q_fit", "Q_depositors_fit"],
            prov,
        )
    end
end

"""
    _sensitivity_table() -> table

The clinker sensitivity of the calibration page: the same fitted rate constants
with the alite content moved by ±20 %, far more than a Bogue or a QXRD analysis
is uncertain by. Three coupled runs.
"""
function _sensitivity_table()
    return get!(_CACHE, "calibration_sensitivity") do
        data = resample_log(CEM_I_TARGET, N_RESIDUALS_COUPLED)
        δs = Float64[]
        c3s = Float64[]
        qend = Float64[]
        rmse = Float64[]
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
                tend = data.t[end], pk_params = apply_parameters(CALIBRATED_THETA),
            )
            _, Q, _ = heat_release(run.sol, run.kp; times = data.t)
            push!(δs, δ)
            push!(c3s, c)
            push!(qend, Q[end] / 1000)
            push!(rmse, sqrt(sum(abs2, Q ./ 1000 .- data.Q) / length(data.Q)))
        end
        prov = _provenance(
            "sensitivity of the fit to the alite content";
            instants = length(data.t), window = data.t[end],
            extra = [
                "the fitted rate constants are held at CALIBRATED_THETA;",
                "only the clinker composition moves, the rest rescaled to close",
            ],
        )
        return _table(
            Dict{String, Vector{Float64}}(
                "delta_C3S" => δs, "C3S" => c3s,
                "Q_end_J_per_g" => qend, "RMSE_J_per_g" => rmse,
            ),
            ["delta_C3S", "C3S", "Q_end_J_per_g", "RMSE_J_per_g"],
            prov,
        )
    end
end

"""
    read_precomputed(name) -> (; columns, provenance, header)

One coupled trajectory of the documentation, computed on demand.

`columns` maps a column name to its vector, `provenance` is a block of lines
saying what produced the numbers — package and solver versions, commit, sampling,
composition — and `header` is the column order. The name and the shape are what
they were when these tables were read from stored CSV files, so a page written
against that interface is unchanged.

The result is **cached for the process**, which matters because Documenter runs
the whole site in one: a page asking for a run's phase history and its
calorimetry gets one integration, not two.
"""
function read_precomputed(name::AbstractString)
    tag = endswith(name, ".csv") ? name[1:(end - 4)] : String(name)
    endswith(tag, "_phases") && return _ionic_phases_table(tag[1:(end - 7)])
    endswith(tag, "_heat") && return _ionic_heat_table(tag[1:(end - 5)])
    haskey(CALIBRATION_CASES, tag) && return _calibration_table(tag)
    tag == "calibration_sensitivity" && return _sensitivity_table()
    return error(
        "unknown coupled trajectory \"$name\". Known: " *
            join(
            sort(
                vcat(
                    [k * "_phases" for k in keys(IONIC_CASES)],
                    [k * "_heat" for k in keys(IONIC_CASES)],
                    collect(keys(CALIBRATION_CASES)),
                    ["calibration_sensitivity"],
                ),
            ), ", ",
        ),
    )
end

"""
    phase_families(tbl) -> Vector{String}

The phase-family columns of a phase-history table, in a fixed order chosen so a
stacked plot reads from the anhydrous grains at the bottom to the pore water at
the top — the order the classical hydration diagram uses.

Families absent from the table are skipped, so the same order serves a paste
with limestone and one without.
"""
function phase_families(tbl)
    order = [
        "anhydrous", "gypsum", "calcite",
        "C-S-H", "CH", "AFt", "AFm", "hydrogarnet", "FH3",
        "water", "void",
    ]
    return [f for f in order if haskey(tbl.columns, f)]
end

"""
    family_colors() -> Dict{String,Symbol}

One color per phase family, fixed here so that every figure in the
documentation uses the same one for the same phase — a reader should not have to
re-read the legend between two plots.
"""
function family_colors()
    return Dict(
        "anhydrous" => :grey40, "gypsum" => :khaki, "calcite" => :tan,
        "C-S-H" => :steelblue, "CH" => :seagreen, "AFt" => :orchid,
        "AFm" => :mediumpurple, "hydrogarnet" => :peru, "FH3" => :indianred,
        "water" => :lightskyblue, "void" => :white,
    )
end
