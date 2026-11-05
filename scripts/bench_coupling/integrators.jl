# SPDX-License-Identifier: LGPL-2.1-or-later
#
# What a kinetic run costs, by integrator and tolerance, and what a host code
# that imposes its own steps (a point of a finite-element mesh) pays per step.
#
#   julia --project=scripts scripts/bench_coupling/integrators.jl [outdir]
#
# `PARTS=integrators` or `PARTS=march` in the environment writes one of the two
# files only.
#
# Three cases, the first two those of the clinker page, the third that of the
# coupling tutorial:
#
#   clinker      four Parrott–Killoh laws, products imposed by the reactions,
#                semi-adiabatic cell, no equilibrium partition, 7 days;
#   clinker_eq   the same clinker with the products left to thermodynamics,
#                isothermal: the laws read no activity, so the partition is
#                solved once per accepted step, 2 days;
#   calcite      calcite dissolving under r = k(1 − Ω): the law reads the
#                partition, which is solved at every evaluation of the
#                right-hand side, 10⁵ s.
#
# `integrators.csv`: each case integrated over its interval by every algorithm at
# three tolerances; the time is the minimum of `REPEATS` runs after a warm-up (one
# run when it exceeds `SLOW` seconds), and the deviation is the largest one from
# a Rodas5P run at reltol = 1e-10 over its instants, relative to the largest
# magnitude of each component.
#
# `march.csv`: the two cases with an equilibrium partition advanced in `NSTEPS`
# steps imposed by the host, either by an ODE restarted at every step from the
# speciated state the last one ended at, or by one implicit step
# (`kinetic_step`) per host step, the rate laws reading the first state as their
# `n_initial`. Reported: the time of the whole march, per step, and the largest
# deviation of the states at the host's steps from the reference.

using ChemistryLab, DynamicQuantities, OptimaSolver, OrderedCollections, OrdinaryDiffEq, Printf, Logging
const SciMLBase = ChemistryLab.SciMLBase

const REPEATS = 3
const SLOW = 20.0
const NSTEPS = 50
const PARTS = split(get(ENV, "PARTS", "integrators,march"), ',')
const OUTDIR = get(ARGS, 1, joinpath(@__DIR__, "..", "..", "docs", "src", "assets", "bench_coupling"))

# What `using OrdinaryDiffEq` exports: three stiff methods (Rosenbrock of orders
# 4 and 2, variable-order BDF), one explicit Runge–Kutta, and the explicit one
# that switches to a stiff one when it detects stiffness.
const ALGORITHMS = [
    ("Rodas5P", Rodas5P),
    ("Rosenbrock23", Rosenbrock23),
    ("FBDF", FBDF),
    ("Tsit5", Tsit5),
    ("AutoTsit5(Rosenbrock23())", () -> AutoTsit5(Rosenbrock23())),
]
const TOLERANCES = [1.0e-4, 1.0e-6, 1.0e-8]

# ── the cases ────────────────────────────────────────────────────────────────

const SUBSTANCES = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
const CLINKER = (C3S = 0.619, C2S = 0.165, C3A = 0.08, C4AF = 0.087)   # kg per kg of cement

function clinker_laws()
    α_max = powers_alpha_max(0.4)
    blaine = 380.0u"m^2/kg"
    return Dict(
        "C3S" => parrott_killoh_avrami(PK84_PARAMS_C3S, "C3S"; α_max, blaine),
        "C2S" => parrott_killoh_avrami(PK84_PARAMS_C2S, "C2S"; α_max, blaine),
        "C3A" => parrott_killoh_avrami(PK84_PARAMS_C3A, "C3A"; α_max, blaine),
        "C4AF" => parrott_killoh_avrami(PK84_PARAMS_C4AF, "C4AF"; α_max, blaine),
    )
end

clinker_species() = speciation(
    SUBSTANCES, split("C3S C2S C3A C4AF Portlandite Jennite ettringite monosulphate12 C3AH6 C3FH6 H2O@");
    aggregate_state = [AS_AQUEOUS],
)

function clinker_state(cs)
    st = ChemicalState(cs)
    for (name, frac) in pairs(CLINKER)
        set_quantity!(st, string(name), frac * u"kg")
    end
    set_quantity!(st, "H2O@", 0.4u"kg")
    return st
end

function clinker_case()
    cs = ChemicalSystem(clinker_species(), CEMDATA_PRIMARIES)
    pk = clinker_laws()
    function rxn(reactants, products, name)
        r = Reaction(
            OrderedDict(cs[k] => v for (k, v) in reactants), OrderedDict(cs[k] => v for (k, v) in products);
            symbol = name,
        )
        r[:rate] = pk[name]
        return r
    end
    reactions = [
        rxn(("C3S" => 1.0, "H2O@" => 103 / 30), ("Jennite" => 1.0, "Portlandite" => 4 / 3), "C3S"),
        rxn(("C2S" => 1.0, "H2O@" => 73 / 30), ("Jennite" => 1.0, "Portlandite" => 1 / 3), "C2S"),
        rxn(("C3A" => 1.0, "H2O@" => 6.0), ("C3AH6" => 1.0,), "C3A"),
        rxn(("C4AF" => 1.0, "Portlandite" => 2.0, "H2O@" => 10.0), ("C3AH6" => 1.0, "C3FH6" => 1.0), "C4AF"),
    ]
    cal = SemiAdiabaticCalorimeter(;
        Cp = 900.0u"J/K", T_env = 293.15u"K", heat_loss = ΔT -> 0.3 * ΔT + 0.003 * ΔT^2, T0 = 293.15u"K",
    )
    build(st, tspan; n_initial = nothing) =
        KineticsProblem(cs, reactions, st, tspan; calorimeter = cal, equilibrium_solver = nothing, n_initial)
    return (; build, state = clinker_state(cs), tspan = (0.0, 7.0 * 86400.0), stepper = nothing)
end

function clinker_eq_case()
    cs = ChemicalSystem(clinker_species(), CEMDATA_PRIMARIES; kinetic_species = clinker_laws())
    model = DiluteSolutionModel()
    es = EquilibriumSolver(cs, model, OptimaOptimizer())
    build(st, tspan; n_initial = nothing) = KineticsProblem(cs, st, tspan; equilibrium_solver = es, n_initial)
    st = clinker_state(cs)
    stepper = KineticStepSolver(cs, model, build(st, (0.0, 1.0)).kinetic_reactions)
    return (; build, state = st, tspan = (0.0, 2.0 * 86400.0), stepper)
end

function calcite_case()
    sp = Dict(symbol(x) => x for x in build_species(datapath("slop98-inorganic-thermofun.json"); verbose = false))
    cs = ChemicalSystem([sp[x] for x in split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Cal")], ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"])
    model = DiluteSolutionModel()
    st = ChemicalState(cs)
    set_quantity!(st, "Cal", 0.05u"mol")
    set_quantity!(st, "H2O@", 1.0u"kg")
    set_quantity!(st, "H+", 1.0e-7u"mol")
    set_quantity!(st, "OH-", 1.0e-7u"mol")
    dissolution = Reaction(OrderedDict(cs["Cal"] => 1.0), OrderedDict(cs["Ca+2"] => 1.0, cs["CO3-2"] => 1.0); symbol = "calcite")
    g = [ustrip(us"J/mol", x[:ΔₐG⁰](T = 298.15u"K", P = 1.0e5u"Pa"; unit = true)) / (R_GAS * 298.15) for x in cs.species]
    ν = Float64.(KineticReaction(cs, dissolution, KineticFunc((T, P, t, n, lna, n0) -> 0.0, NamedTuple(), u"mol/s")).stoich)
    names = [symbol(x) for x in cs.species]
    law = KineticFunc(
        (T, P, t, n, lna, n0) -> 1.0e-4 * (1 - saturation_ratio(ν, [lna[s] for s in names], g)), NamedTuple(), u"mol/s",
    )
    dissolution[:rate] = law
    es = EquilibriumSolver(cs, model, OptimaOptimizer())
    build(st, tspan; n_initial = nothing) =
        KineticsProblem(cs, [dissolution], st, tspan; activity_model = model, equilibrium_solver = es, n_initial)
    stepper = KineticStepSolver(cs, model, [KineticReaction(cs, dissolution, law)])
    return (; build, state = st, tspan = (0.0, 1.0e5), stepper)
end

const CASES = [("clinker", clinker_case), ("clinker_eq", clinker_eq_case), ("calcite", calcite_case)]

# ── the runs ─────────────────────────────────────────────────────────────────

run_once(kp, alg, tol) = with_logger(NullLogger()) do
    integrate(kp, KineticsSolver(; ode_solver = alg, reltol = tol, abstol = tol * 1.0e-3, maxiters = 100_000))
end

# The largest deviation of `sol` from `ref` at the reference's instants, relative
# to the largest magnitude of each component along the reference.
function deviation(sol, ref)
    scale = [max(maximum(abs, getindex.(ref.u, j)), eps()) for j in eachindex(ref.u[1])]
    return maximum(maximum(abs.(sol(t) .- u) ./ scale) for (t, u) in zip(ref.t, ref.u))
end

function timed(f)
    first = @elapsed f()
    first > SLOW && return @elapsed f()
    return minimum(@elapsed(f()) for _ in 1:REPEATS)
end

amounts(st) = Float64[ustrip(us"mol", x) for x in st.n]
state_deviation(st, ref) = maximum(abs.(amounts(st) .- amounts(ref))) / maximum(abs, amounts(ref))

# The host's march: `NSTEPS` steps over the case's interval, each from the
# state the last one ended at, the rate laws reading the amounts of the first
# state as their `n_initial`. Returns the state at the end of every step.
host_times(case) = range(case.tspan[1], case.tspan[2]; length = NSTEPS + 1)[2:end]

function march_ode(case, alg, tol)
    st = case.state
    t_prev = case.tspan[1]
    states = typeof(st)[]
    for t in host_times(case)
        kp = case.build(st, (t_prev, t); n_initial = case.state)
        st = only(speciated_states(run_once(kp, alg, tol), kp; times = [t], start = st))
        push!(states, st)
        t_prev = t
    end
    return states
end

function march_implicit(case)
    st = case.state
    t_prev = case.tspan[1]
    states = typeof(st)[]
    with_logger(NullLogger()) do
        for t in host_times(case)
            st = kinetic_step(case.stepper, st, (t - t_prev) * u"s"; t = t_prev, n_initial = case.state)
            push!(states, st)
            t_prev = t
        end
    end
    return states
end

# The largest deviation over the march, relative to the largest amount.
march_deviation(states, refs) = maximum(state_deviation(s, r) for (s, r) in zip(states, refs))

mkpath(OUTDIR)
cases = [(name, build()) for (name, build) in CASES]

"integrators" in PARTS && open(joinpath(OUTDIR, "integrators.csv"), "w") do io
    println(io, "case,algorithm,reltol,time_s,naccept,nreject,nf,njacs,deviation,retcode")
    for (cname, case) in cases
        kp = case.build(case.state, case.tspan)
        ref = run_once(kp, Rodas5P(), 1.0e-10)
        for (aname, mk) in ALGORITHMS, tol in TOLERANCES
            sol = run_once(kp, mk(), tol)
            t = timed(() -> run_once(kp, mk(), tol))
            s = sol.stats
            dev = SciMLBase.successful_retcode(sol) ? deviation(sol, ref) : NaN
            line = @sprintf(
                "%s,%s,%.0e,%.4g,%d,%d,%d,%d,%.3e,%s", cname, aname, tol, t, s.naccept, s.nreject, s.nf, s.njacs, dev, sol.retcode,
            )
            println(io, line)
            flush(io)
            println(line)
        end
    end
end

"march" in PARTS && open(joinpath(OUTDIR, "march.csv"), "w") do io
    println(io, "case,method,reltol,steps,time_s,time_per_step_ms,deviation")
    for (cname, case) in cases
        case.stepper === nothing && continue
        kp = case.build(case.state, case.tspan)
        ref_sol = run_once(kp, Rodas5P(), 1.0e-10)
        refs = speciated_states(ref_sol, kp; times = collect(host_times(case)))
        # Tsit5 left out: on the stiff case it takes a hundred seconds per interval.
        runs = Any[("ODE restarted: $aname", tol, () -> march_ode(case, mk(), tol)) for (aname, mk) in ALGORITHMS[[1, 2, 5]] for tol in (1.0e-4, 1.0e-6)]
        push!(runs, ("implicit step", NaN, () -> march_implicit(case)))
        for (label, tol, f) in runs
            states = f()
            t = timed(f)
            line = @sprintf(
                "%s,%s,%.0e,%d,%.4g,%.4g,%.3e", cname, label, tol, NSTEPS, t, 1.0e3 * t / NSTEPS, march_deviation(states, refs),
            )
            println(io, line)
            flush(io)
            println(line)
        end
    end
end
