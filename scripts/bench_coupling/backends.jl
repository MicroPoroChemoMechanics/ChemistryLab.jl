# SPDX-License-Identifier: LGPL-2.1-or-later
#
# What each back end returns on its own, before the dual Newton polishes it, and
# what the certified route costs.
#
#   julia --project=scripts scripts/bench_coupling/backends.jl [outdir]
#
# Cases: calcite in water at 25 °C, with 1 and 50 mmol of carbon dioxide, pure
# water, and a Portland paste (CEM I 52.5 N, Bogue) at w/b = 0.30, 0.45, 0.60.
# Back ends: the active set of OptimaSolver (its default), its interior point
# (`algorithm = :interior_point`) and Ipopt, each run unpolished; then the
# certified route, `equilibrate_certified`.
#
# Reported per run, in `backends.csv`: the worst element balance relative to
# what its row holds, the stationarity residual (scaled as the certificate
# scales it), whether the certificate proves the answer optimal, and the time,
# the minimum of `REPEATS` runs after a warm-up. A back end leaves an absent
# species on its lower bound, which the state records at the amount floor; it is
# set to zero before the certificate reads the answer, as the polish does, so
# that the certificate judges it absent rather than present and unstable.

using ChemistryLab, DynamicQuantities, OptimaSolver, Optimization, OptimizationIpopt, Printf, Logging
const SciMLBase = ChemistryLab.SciMLBase

const REPEATS = 3
const OUTDIR = get(ARGS, 1, joinpath(@__DIR__, "..", "..", "docs", "src", "assets", "bench_coupling"))

function calcite_state(co2)
    sp = Dict(symbol(x) => x for x in build_species(datapath("slop98-inorganic-thermofun.json"); verbose = false))
    cs = ChemicalSystem([sp[x] for x in split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Cal")], ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"])
    st = ChemicalState(cs)
    set_quantity!(st, "H2O@", 1.0u"kg")
    set_quantity!(st, "H+", 1.0e-7u"mol")
    set_quantity!(st, "OH-", 1.0e-7u"mol")
    co2 === nothing && return st, DiluteSolutionModel()
    set_quantity!(st, "Cal", 0.05u"mol")
    iszero(co2) || set_quantity!(st, "CO2@", co2 * u"mol")
    return st, DiluteSolutionModel()
end

const SUBSTANCES = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
const PASTE_SYSTEM = phase_list_system("Portland paste (Lothenbach and Winnefeld 2006)", SUBSTANCES)
const CEMENT = with_extents(
    material_template("CEM I 52.5 N (Gruyaert 2010), Bogue", PASTE_SYSTEM),
    Dict("C3S" => 0.8, "C2S" => 0.45, "C3A" => 0.9, "C4AF" => 0.6),
)
const CEMENT_MODEL = cemdata18_activity_model(:KOH)

# the reactants of the recipe, the start `equilibrate_certified(recipe, …)` takes
paste_state(w_b) = (budget(Recipe(CEMENT => 1.0; w_b), PASTE_SYSTEM).state, CEMENT_MODEL)

const CASES = [
    ("calcite in water", () -> calcite_state(0.0)),
    ("calcite + 1 mmol CO2", () -> calcite_state(1.0e-3)),
    ("calcite + 50 mmol CO2", () -> calcite_state(5.0e-2)),
    ("pure water", () -> calcite_state(nothing)),
    ("CEM I paste w/b 0.30", () -> paste_state(0.3)),
    ("CEM I paste w/b 0.45", () -> paste_state(0.45)),
    ("CEM I paste w/b 0.60", () -> paste_state(0.6)),
]

const BACKENDS = [
    ("OptimaSolver active set", () -> OptimaOptimizer()),
    ("OptimaSolver interior point", () -> OptimaOptimizer(algorithm = :interior_point)),
    ("Ipopt", () -> IpoptOptimizer()),
]

amounts(st) = Float64[ustrip(us"mol", x) for x in st.n]

function on_bounds_absent(st)
    out = copy(st)
    for (i, x) in enumerate(amounts(st))
        x <= ChemistryLab._AMOUNT_FLOOR && (out.n[i] = 0.0u"mol")
    end
    return out
end

function timed(f)
    f()
    return minimum(@elapsed(f()) for _ in 1:REPEATS)
end

mkpath(OUTDIR)
open(joinpath(OUTDIR, "backends.csv"), "w") do io
    println(io, "case,route,balance_relative,stationarity,optimal,time_s")
    for (cname, build) in CASES
        st0, model = build()
        des = DualEquilibriumSolver(st0.system, model)
        b = des.A * amounts(st0)
        record(route, st, t) = begin
            c = optimality_certificate(des, st; b)
            line = @sprintf("%s,%s,%.2e,%.2e,%s,%.4g", cname, route, c.balance_relative, c.stationarity, c.optimal, t)
            println(io, line)
            flush(io)
            println(line)
        end
        for (bname, mk) in BACKENDS
            solve_raw() = with_logger(NullLogger()) do
                ChemistryLab._unpolished(() -> equilibrate(st0, mk(); model))
            end
            record(bname, on_bounds_absent(solve_raw()), timed(solve_raw))
        end
        certified() = with_logger(NullLogger()) do
            equilibrate_certified(st0; model)
        end
        record("certified route", first(certified()), timed(certified))
    end
end
