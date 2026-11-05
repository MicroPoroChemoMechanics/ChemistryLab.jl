# [Inside a host code: equilibria and kinetics at every node](@id sec-host-coupling)

!!! info "Before this page"
    [Solving an equilibrium](@ref sec-solving) and
    [Writing a kinetic model](@ref sec-kinetics-syntax). The figures below come
    from three scripts of `scripts/bench_coupling/`, named in
    [Reproducing the figures](@ref sec-host-coupling-reproduce).

A transport or a structural code, such as a finite-volume or finite-element
solver of a porous medium, calls the chemistry at each of its nodes (a Gauss
point, a cell) at every one of its steps. Between two calls it changes the
element budget `b` of the node, or it advances the time by its own step `Δt`.
The chemistry answers with the composition at that budget, or after that step,
and the host reads from it what its own equations need: a porosity, a pH, the
amount of a solid, a flux.

The cost of the coupled computation is then

```
cost = (number of nodes) × (number of host steps) × (cost of one call)
```

and only the last factor is the chemistry's. This page measures it for each way
of calling, and says which call to make when.

```@setup host
using ChemistryLab, Markdown, Printf

# The tables of this page are the CSV files the scripts of
# `scripts/bench_coupling/` write; nothing is typed by hand.
const BENCH = joinpath(pkgdir(ChemistryLab), "docs", "src", "assets", "bench_coupling")

function bench_rows(name)
    lines = readlines(joinpath(BENCH, name))
    header = split(lines[1], ',')
    return [Dict(zip(header, split(l, ','))) for l in lines[2:end] if !isempty(l)]
end

function bench_table(rows, columns; filter = r -> true)
    head = "| " * join(first.(columns), " | ") * " |"
    rule = "|" * join((i == 1 ? ":--" : "--:" for i in eachindex(columns)), "|") * "|"
    body = ["| " * join((f(r) for (_, f) in columns), " | ") * " |" for r in rows if filter(r)]
    return Markdown.parse(join(vcat(head, rule, body), "\n"))
end

ms(x) = @sprintf("%.1f", 1.0e3 * parse(Float64, x))
fixed(x; d = 1) = @sprintf("%.*f", d, parse(Float64, x))
sci(x) = (v = parse(Float64, x); isnan(v) ? "—" : @sprintf("%.0e", v))
nothing
```

## [Equilibrium at a node](@id sec-host-coupling-equilibrium)

### The entry points, and what each one costs

A node whose budget changed by a little between two calls is the common case,
and the one this section measures: a Portland paste (CEM I 52.5 N, Bogue, w/b =
0.5), with and without the C-S-H surface of [Guo2018](@citet), solved for a budget 1 %
richer in every element than the one of its previous answer.

```@example host
bench_table(bench_rows("equilibrium.csv"), [
    "system" => r -> r["system"],
    "call" => r -> r["entry_point"],
    "time (ms)" => r -> fixed(r["time_ms"]),
    "allocated (MB)" => r -> fixed(r["alloc_MB"]),
])
```

The calls, from the most to the least general:

- **cold: certified** is [`equilibrate_certified`](@ref) on the state given: it
  solves from every back end that is loaded and returns the answer its
  certificate proves.
- **cold: one back end** is `solve` of an [`EquilibriumSolver`](@ref) built once
  for the system: OptimaSolver's answer, polished by the dual Newton method.
  **cold: solver rebuilt** builds that solver at every call.
- **warm** is [`solve_certified`](@ref) of a [`DualEquilibriumSolver`](@ref)
  built once, started from the node's previous answer: the dual Newton alone.
- **warm + multipliers** resumes it also from the multipliers of that answer,
  `y0`, which a `dual_state` reference returns with every solve.
- **warm + multipliers - report** skips what the certificate covers (`scope`,
  the curvature of the answer, the range of the activity model) with
  `report = false`, for a caller that reads only `optimal` and the residuals.

What decides the cost is where the solve starts. Building the solver is not: the
times with the solver rebuilt at every call are those with it built once. From
the node's previous answer and its multipliers, a call takes under a third of the
time of a cold one on the paste and half of it with the surface, and allocates
five and two and a half times less; `report = false` takes a further 5 to 15 %
off.

### What to do at a node

A node keeps its last answer and the multipliers of that answer, and each thread
uses a solver of its own, since a solver keeps the start of its next solve. Below, four nodes of a column of calcite receive carbon dioxide at
every step of their host, each more than the one before it; each call is the warm
one, with the certified cold route as the fallback of an answer the certificate
does not prove.

```@example host
using DynamicQuantities, OptimaSolver

sp = Dict(symbol(x) => x for x in build_species(datapath("slop98-inorganic-thermofun.json"); verbose = false))
cs = ChemicalSystem([sp[x] for x in split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Cal")],
                    ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"])
model = DiluteSolutionModel()
st0 = ChemicalState(cs)
set_quantity!(st0, "Cal", 0.05u"mol")
set_quantity!(st0, "H2O@", 1.0u"kg")
set_quantity!(st0, "H+", 1.0e-7u"mol")
set_quantity!(st0, "OH-", 1.0e-7u"mol")

des = DualEquilibriumSolver(cs, model)               # once
first_answer, _ = equilibrate_certified(st0; model)  # the nodes start equal: solved once
states = fill(first_answer, 4)                       # the last answer of each node
ys = Vector{Union{Nothing, Vector{Float64}}}(nothing, 4)    # and its multipliers
co2 = des.A[:, findfirst(x -> symbol(x) == "CO2@", cs.species)]   # the budget of 1 mol of CO₂

fallbacks = 0
for step in 1:5, k in eachindex(states)
    b = des.A * [ustrip(us"mol", x) for x in states[k].n] .+ 1.0e-4 * k .* co2
    y = Ref{Any}(nothing)
    st, cert = solve_certified(des, (states[k],); b, y0 = ys[k], dual_state = y, report = false)
    if cert.optimal
        ys[k] = y[].y
    else                                             # the cold route, where the warm one failed
        global fallbacks += 1
        st, cert = equilibrate_certified(states[k]; model, b)
        ys[k] = nothing
    end
    states[k] = st
end
@printf("fallbacks: %d of 20 calls; calcite left: %s mol\n", fallbacks,
        join((@sprintf("%.5f", ustrip(us"mol", moles(st, "Cal"))) for st in states), ", "))
```

Two things a call should not do:

- **start from the host's initial state** at every step: the cold route then
  solves again what the previous step already solved;
- **solve two identical nodes twice**: nodes that share a budget, at the start
  of a computation or where the transport has not yet arrived, share their
  answer, as the four nodes above share their first one.

A derivative with respect to `b`, which an implicit host needs for its Jacobian,
comes from the conditions of equilibrium at the answer
([Differentiating an equilibrium](@ref)): one factorization serves every partial
derivative, and no second equilibrium is solved.

## [Kinetics at a node](@id sec-host-coupling-kinetics)

A kinetic call advances the node over the host's step `Δt`. Its cost is the
number of equilibria it solves times the cost of one, and which of the two
dominates is decided by the rate laws, before any choice of integrator:

- **laws that read only the kinetic amounts** (Parrott–Killoh, an imposed
  stoichiometry): the equilibrium partition, when there is one, is solved once
  per *accepted* step and held within it (`speciation = :frozen`, which
  `integrate` takes by itself, [When the partition may be frozen within a step](@ref sec-theory-pe-splitting));
  the problem is mildly stiff, and the integrator that takes the fewest steps
  wins;
- **laws that read the solution** (a saturation ratio, a pH): the partition is
  solved at every evaluation of the right-hand side, and the relaxation to
  equilibrium makes the problem stiff: an explicit integrator is held to steps
  of the order of its relaxation time.

### Integrators and tolerances

Three cases, integrated over their interval by each method that
`using OrdinaryDiffEq` exports, at `reltol` = `10⁻⁴`, `10⁻⁶`, `10⁻⁸` and
`abstol = 10⁻³ reltol` mol; the deviation is the largest one from a Rodas5P run
at `10⁻¹⁰`, relative to the largest magnitude of each component:

- *clinker*: the four Parrott–Killoh laws of
  [Hydration kinetics of a CEM I 52.5 R clinker](@ref sec-clinker-kinetics), products imposed, semi-adiabatic cell, no
  equilibrium partition, 7 days;
- *clinker_eq*: the same clinker, products left to thermodynamics, isothermal, 2
  days: one equilibrium per accepted step;
- *calcite*: calcite dissolving under `r = k(1 − Ω)`, as in
  [Coupling kinetics and equilibrium](@ref sec-coupling), 10⁵ s: one equilibrium
  per evaluation of the right-hand side.

```@example host
rows = bench_rows("integrators.csv")
bench_table(rows, [
    "case" => r -> r["case"],
    "method" => r -> "`" * r["algorithm"] * "`",
    "reltol" => r -> sci(r["reltol"]),
    "time (s)" => r -> fixed(r["time_s"]; d = 3),
    "steps" => r -> r["naccept"],
    "RHS evaluations" => r -> r["nf"],
    "deviation" => r -> sci(r["deviation"]),
])
```

What the table says, case by case:

- without a partition, or with one solved per accepted step, the explicit
  `Tsit5` is the fastest at every tolerance, five to six times faster than
  `Rodas5P` without a partition and 1.5 to 2 times with one: the laws are not
  stiff enough for a stiff method to pay for its Jacobian, and with a partition
  `Tsit5` takes fewer steps, hence fewer equilibria; at `10⁻⁸` the stiff methods
  are the more accurate;
- with the partition read by the law, `Tsit5` takes tens of thousands of steps,
  bounded by stability rather than accuracy, and three orders of magnitude more
  time; `AutoTsit5(Rosenbrock23())`, which starts explicit and switches when it
  detects stiffness, is the fastest there, and `Rodas5P` close behind;
- the deviation follows `reltol`, within a factor of sixty for the low-order
  methods at the loosest tolerance, so the tolerance is set by the accuracy the
  host needs; tightening it by two orders costs `Rodas5P` about twice the time,
  and the second-order `Rosenbrock23` about four times.

So, for a call the host makes at every node and every step:

```julia
using OrdinaryDiffEq

# laws that read only the kinetic amounts
ks = KineticsSolver(; ode_solver = Tsit5(), reltol = 1.0e-4, abstol = 1.0e-7)
# laws that read the solution
ks = KineticsSolver(; ode_solver = AutoTsit5(Rosenbrock23()), reltol = 1.0e-4, abstol = 1.0e-7)
```

`abstol` is in moles: scale it on the amounts of the node, `10⁻³ reltol` times
its largest kinetic amount, so that it means the same on a cubic millimeter of
paste as on a kilogram. `Rodas5P`, the default of `integrate`, is the safe choice
when the stiffness is not known: on these cases it is never more than six times
slower than the fastest.

### A step imposed by the host

A host code does not hand the chemistry an interval to integrate: it imposes
its own steps, one call per step and per node, each from the state the last one
ended at. Two ways to take such a step:

- **an ODE restarted at every step**: a [`KineticsProblem`](@ref) over the
  step, from the state the last step ended at, [`integrate`](@ref), and the
  speciated state at the end of the step, [`speciated_states`](@ref);
- **one implicit step** per host step, [`kinetic_step`](@ref): one Gibbs
  minimization with the extents among its unknowns, first order in the step
  ([The implicit step](@ref sec-theory-implicit-step)).

Two keywords make a restarted run follow the uninterrupted one, at a cost that
does not grow with the length of the history:

Below, the calcite of the nodes above dissolves under a first-order law in its
degree of reaction, `r = k n/n_initial`, over ten host steps of 1000 s, by a run
restarted at each step and by one uninterrupted run:

```@example host
using OrderedCollections, OrdinaryDiffEq

dissolution = Reaction(OrderedDict(cs["Cal"] => 1.0), OrderedDict(cs["Ca+2"] => 1.0, cs["CO3-2"] => 1.0);
                       symbol = "calcite")
dissolution[:rate] = KineticFunc((T, P, t, n, lna, n0) -> 2.0e-6 * n["Cal"] / n0["Cal"], NamedTuple(), u"mol/s")
es = EquilibriumSolver(cs, model, OptimaOptimizer())
ks = KineticsSolver(; ode_solver = Rodas5P(), reltol = 1.0e-6, abstol = 1.0e-12)
problem(st, span; n_initial = nothing) =
    KineticsProblem(cs, [dissolution], st, span; activity_model = model, equilibrium_solver = es, n_initial)

state = first_answer                              # the history starts here
for k in 1:10
    span = (1000.0 * (k - 1), 1000.0 * k)
    kp = problem(state, span; n_initial = first_answer)
    global state = only(speciated_states(integrate(kp, ks), kp; times = [span[2]], start = state))
end
whole = problem(first_answer, (0.0, 1.0e4))
uninterrupted = only(speciated_states(integrate(whole, ks), whole; times = [1.0e4]))
@printf("calcite after 10⁴ s: %.8f mol restarted, %.8f mol uninterrupted\n",
        ustrip(us"mol", moles(state, "Cal")), ustrip(us"mol", moles(uninterrupted, "Cal")))
```

- `n_initial = state0` is what the rate laws read as their `n_initial`, the
  amounts the history started from. Without it, each run counts from the start
  of its own step, and a Parrott–Killoh law, whose degree of hydration is
  `1 − n/n_initial`, restarts at `α = 0` at every step. `kinetic_step` takes the
  same keyword.
- `start = state` hands the replay of the speciated state the answer the last
  step ended with; without it, the replay seeds a chain of certified solves from
  the start of the run.

Measured on the two cases with a partition, over 50 host steps (an hour each
for the clinker, 2000 s for calcite), against the uninterrupted `Rodas5P` run at
`10⁻¹⁰`; the deviation is the largest one over the host steps:

```@example host
bench_table(bench_rows("march.csv"), [
    "case" => r -> r["case"],
    "step" => r -> r["method"],
    "reltol" => r -> sci(r["reltol"]),
    "time per host step (ms)" => r -> fixed(r["time_per_step_ms"]),
    "deviation" => r -> sci(r["deviation"]),
])
```

The implicit step is eight to twenty times cheaper than a restarted ODE, and its
deviation is that of a first-order method over the host's step, of the order of
`10⁻²` on the clinker at one-hour steps, where the restarted ODE meets its tolerance; on
calcite, whose relaxation is far shorter than a host step, both reach the
equilibrium. A restarted ODE pays at every host step what an uninterrupted run
pays once, a start from a small step and the speciated state at its end, so
that a single run over the whole history of the clinker costs less than one of
its host steps; on the clinker the stiff `Rodas5P` is then the cheapest, since
the restart takes away the long steps the explicit methods win with.

So, when the host imposes its steps:

- the **implicit step** when the host step is short against the time over
  which the rates change, or when its accuracy is that of the host's own time
  discretization, first order as well; [`kinetic_step_adaptive`](@ref) when it
  is not, which halves the step until its error estimate is within tolerance
  and the step is certified;
- a **restarted ODE**, with `n_initial` and `start`, when the trajectory within
  a host step must be resolved to a tolerance.


## [Reproducing the figures](@id sec-host-coupling-reproduce)

Each table above is the CSV file a script writes into
`docs/src/assets/bench_coupling/`:

```bash
julia --project=scripts scripts/bench_coupling/equilibrium.jl   # equilibrium.csv
julia --project=scripts scripts/bench_coupling/backends.jl      # backends.csv
julia --project=scripts scripts/bench_coupling/integrators.jl   # integrators.csv, march.csv
```

Times are the minimum of several runs after a warm-up, on one thread with
`OPENBLAS_NUM_THREADS=1`, on an AMD EPYC-Milan processor with Julia 1.13.
The answers do not depend on the machine; the times do.
