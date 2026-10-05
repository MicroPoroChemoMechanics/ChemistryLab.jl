# [Coupling kinetics and equilibrium](@id sec-coupling)

!!! info "Before this page"
    The tutorial [Chemical Kinetics](@ref sec-kinetics). The equations this page
    runs are in [Kinetics under partial equilibrium](@ref sec-theory-pe-kinetics):
    the partition, the right-hand side and its Jacobian, and the implicit step.

One mineral, three ways of advancing it. Calcite dissolves into water under the
law ``r = k(1 - \Omega)``, which reads the saturation ratio of the solution and
therefore the equilibrium partition. The page integrates it with the partition
solved in the right-hand side of the ODE, with the implicit step, and with the
partition frozen within each ODE step, and plots the three trajectories against
the equilibrium they should reach.

## 1. The system and the law

```@example coupling
using ChemistryLab, DynamicQuantities, OptimaSolver, OrderedCollections, OrdinaryDiffEq, Printf, Plots
using Logging # hide

sp = Dict(symbol(x) => x for x in build_species(datapath("slop98-inorganic-thermofun.json"); verbose = false))
cs = ChemicalSystem([sp[x] for x in split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Cal")],
                    ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"])
model = DiluteSolutionModel()

function initial_state()
    st = ChemicalState(cs)
    set_quantity!(st, "Cal", 0.05u"mol")
    set_quantity!(st, "H2O@", 1.0u"kg")
    set_quantity!(st, "H+", 1.0e-7u"mol")
    set_quantity!(st, "OH-", 1.0e-7u"mol")
    return st
end
calcite() = Reaction(OrderedDict(cs["Cal"] => 1.0), OrderedDict(cs["Ca+2"] => 1.0, cs["CO3-2"] => 1.0);
                     symbol = "calcite")

# The saturation ratio of the reaction, from the log activities the law is
# handed and the standard Gibbs energies at 25 °C.
g = [ustrip(us"J/mol", x[:ΔₐG⁰](T = 298.15u"K", P = 1.0e5u"Pa"; unit = true)) / (R_GAS * 298.15) for x in cs.species]
ν = Float64.(KineticReaction(cs, calcite(), KineticFunc((T, P, t, n, lna, n0) -> 0.0, NamedTuple(), u"mol/s")).stoich)
k = 1.0e-4                                       # mol/s
law = KineticFunc((T, P, t, n, lna, n0) -> k * (1 - saturation_ratio(ν, [lna[symbol(x)] for x in cs.species], g)),
                  NamedTuple(), u"mol/s")
dissolution = calcite()
dissolution[:rate] = law

# The equilibrium the dissolution must reach: calcite in water, certified.
eq, cert = equilibrate_certified(initial_state(); model)
cal_eq = ustrip(us"mol", moles(eq, "Cal"))
@printf("at equilibrium: %.6e mol of calcite left of 0.05, certified: %s\n", cal_eq, cert.optimal)
```

## 2. The partition

```@example coupling
kp = KineticsProblem(cs, [dissolution], initial_state(), (0.0, 1.0e5);
                     activity_model = model, equilibrium_solver = EquilibriumSolver(cs, model, OptimaOptimizer()))
println("kinetic:     ", [symbol(kp.system.species[i]) for i in kp.idx_kinetic])
println("equilibrium: ", [symbol(kp.system.species[i]) for i in kp.idx_equilibrium])
```

The state the integrator advances holds the element amounts of the second set
and the amount of calcite, for the reason
[The partition](@ref sec-theory-pe-partition) gives.

## 3. Three routes

The ODE with the partition solved where the right-hand side is evaluated, which
`integrate` chooses by itself for a law that reads it
([When the partition may be frozen within a step](@ref sec-theory-pe-splitting)):

```@example coupling
ks = KineticsSolver(; ode_solver = Rodas5P(), reltol = 1.0e-8, abstol = 1.0e-12)
diagnostics = IOBuffer() # hide
rhs = with_logger(ConsoleLogger(diagnostics)) do # hide
rhs = integrate(kp, ks)
end # hide
occursin("re-speciation failed", String(take!(diagnostics))) && error("a re-speciation failed") # hide
cal_rhs = 0.05 .- vec(reaction_extents(rhs, kp))
@printf("%s in %d steps; calcite at 10⁵ s: %.6e mol\n", rhs.retcode, length(rhs.t), cal_rhs[end])
```

The implicit step ([The implicit step](@ref sec-theory-implicit-step)), marched
by its adaptive controller from one second on:

```@example coupling
kss = KineticStepSolver(cs, model, [KineticReaction(cs, dissolution, law)])
st, t, dt = initial_state(), 0.0, 1.0          # the step lengths are in seconds
t_imp, cal_imp = [0.0], [0.05]
with_logger(ConsoleLogger(diagnostics)) do # hide
while t < 1.0e5
    global st, t, dt
    st, used, dt = kinetic_step_adaptive(kss, st, min(dt, 1.0e5 - t) * u"s")
    t += used
    push!(t_imp, t); push!(cal_imp, ustrip(us"mol", moles(st, "Cal")))
end
end # hide
@printf("%d accepted steps; calcite at 10⁵ s: %.6e mol\n", length(t_imp) - 1, cal_imp[end])
```

And the ODE with the partition frozen within each step, forced:

```@example coupling
frozen = with_logger(NullLogger()) do # hide
frozen = integrate(kp, ks; speciation = :frozen)
end # hide
cal_frozen = 0.05 .- vec(reaction_extents(frozen, kp))
@printf("%s in %d steps; calcite reached %.3g mol, from 0.05\n", frozen.retcode, length(frozen.t), maximum(cal_frozen))
```

## 4. The trajectories

```@example coupling
# The saturation index along the right-hand-side route, on the certified replay;
# the start, in pure water, is left out of the logarithmic time axis.
replay = with_logger(ConsoleLogger(diagnostics)) do # hide
replay = speciated_states(rhs, kp)
end # hide
occursin("could not be certified", String(take!(diagnostics))) && error("an instant was not certified") # hide
logΩ = [saturation_indices(s, model)["Cal"] for s in replay]
on = rhs.t .> 0
p1 = plot(rhs.t[on], cal_rhs[on]; xscale = :log10, label = "ODE, partition in the right-hand side",
          xlabel = "time, s", ylabel = "calcite, mol", lw = 2)
plot!(p1, t_imp[2:end], cal_imp[2:end]; label = "implicit step", marker = :circle, ms = 3, lw = 1)
hline!(p1, [cal_eq]; label = "certified equilibrium", ls = :dash, color = :black)
p2 = plot(rhs.t[on], max.(-logΩ[on], 1.0e-16); xscale = :log10, yscale = :log10, label = false,
          xlabel = "time, s", ylabel = "−log₁₀ Ω of calcite", lw = 2)
fig = plot(p1, p2; layout = (1, 2), size = (900, 330), left_margin = 5Plots.mm, bottom_margin = 6Plots.mm)
savefig(fig, "coupling-calcite.svg"); nothing # hide
```

![](coupling-calcite.svg)

```@example coupling
@printf("calcite at 10⁵ s, mol:  right-hand side %.6e   implicit %.6e   equilibrium %.6e\n",
        cal_rhs[end], cal_imp[end], cal_eq)
@printf("worst log Ω on the replay after 10³ s: %.1e\n", maximum(abs, logΩ[rhs.t .> 1.0e3]))
```

The two routes that treat the partition consistently reach the certified
equilibrium, the saturation index falling to the tolerance of the minimization,
about ``10^{-10}``, where it stays. The frozen one is the route
[When the partition may be frozen within a step](@ref sec-theory-pe-splitting)
explains: with the partition held over a step, the rate is constant over it, the
step overshoots saturation, and the run is returned as a failure rather than a
success.

## See also

- [The silicates of a CEM I clinker, hydrating end to end](@ref sec-coupled-hydration),
  the partition on a clinker, with a law that reads only the kinetic amounts.
- [Validation against Reaktoro](@ref), the equilibrium and the coupling against a
  second code.
- [Writing a kinetic model](@ref sec-kinetics-syntax), the syntax of the routes.
