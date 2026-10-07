# [A slow surface: sorption at a rate](@id sec-tutorial-slow-surface)

!!! info "Before this page"
    [A slow surface](@ref sec-theory-pe-slow-surface) and
    [What a rate law may be](@ref sec-theory-kinetics-admissible) for the theory;
    [Adsorption on a single site family](@ref sec-example-surface-langmuir) for the
    syntax of a surface.

A site family is solved at equilibrium unless a kinetic reaction controls one of
its members. This page builds the smallest slow surface, checks it against its
closed form, sets it beside a fast family, and shows what
[`dissipation`](@ref) says of a law that is not tied to the equilibrium constant.

!!! warning "The rate constant of this page is not a measurement"
    The equilibrium constant below is published. The rate constant is not: the
    papers cited in the theory report equilibrium within milliseconds on oxide
    surfaces, and a slow uptake as an apparent effect of transport. The value
    used here, `k = 10 s⁻¹`, is chosen so that the relaxation takes about a
    quarter of an hour and can be followed. The page shows how the model
    behaves, not how this surface does.

## The surface and the solution

One family of sites, the strong sites of hydrous ferric oxide in the model of
[DzombakMorel1990](@citet), taking calcium:

```math
\equiv\!\mathrm{Hfo_sOH} + \mathrm{Ca^{2+}} \rightleftharpoons \equiv\!\mathrm{Hfo_sOHCa^{2+}} .
```

The constant is read from the copy of `phreeqc.dat` the test oracles use
[ParkhurstAppelo2013](@cite). The family is kept to its free site and this one
complex, and without a surface potential: it is not the full model of Dzombak and
Morel, which also protonates and deprotonates the site and charges the surface.

```@example slow
using ChemistryLab, DynamicQuantities, OptimaSolver, OrderedCollections, OrdinaryDiffEq, Printf, Plots
using Logging # hide
default(framestyle = :box, grid = false)

const RT = R_GAS * 298.15
g0(v) = SymbolicFunc(v * u"J/mol")

dat = read_sorption_model(joinpath(pkgdir(ChemistryLab), "test", "reference", "phreeqc.dat"))
logK = only(reactions_involving(dat, "Hfo_sOHCa+2")).log_K.value

db = Dict(symbol(s) => s for s in
          build_species(datapath("slop98-inorganic-thermofun.json"); verbose = false))
h2o, hp, oh, ca = db["H2O@"], db["H+"], db["OH-"], db["Ca+2"]
G_Ca = ustrip(us"J/mol", ca[:ΔₐG⁰](T = 298.15u"K", P = 1.0e5u"Pa"; unit = true))

function site(sym, g)
    s = Species(sym; aggregate_state = AS_SURFACE, class = SC_SURFCOMPLEX)
    s[:ΔₐG⁰] = g0(g)
    return s
end
# ΔrG° = G(complex) − G(free site) − G(Ca+2) = −RT ln K, with G(free site) = 0
free  = site("XsOH", 0.0)
bound = site("XsOHCa+2", -RT * log(10.0^logK) + G_Ca)
logK
```

The parameters of the run, all defined here:

```@example slow
N_sites = 1.0e-7          # mol of sites, a thousandth of the calcium
n_Ca    = 1.0e-4          # mol of Ca(OH)2 dissolved in one kilogram of water
k       = 10.0            # s⁻¹, illustrative (see the warning above)
t_end   = 2.0e4           # s, some twenty relaxation times
nothing # hide
```

```@example slow
support = SurfaceSupport("hydrous ferric oxide, strong sites", nothing, FixedSurfaceArea(1.0))
family  = SiteFamily("Hfo_s", free, [bound];
                     capacity = TotalSiteAmount(N_sites * u"mol"), support)
cs = ChemicalSystem([h2o, hp, oh, ca, free, bound], [h2o, hp, ca, free];
                    site_families = [family])

idx(s) = findfirst(==(s), symbol.(cs.species))
n0 = Any[fill(0.0u"mol", length(cs.species))...]
n0[idx("H2O@")]  = ustrip(us"mol", 1.0u"kg" / h2o[:M]) * u"mol"
n0[idx("Ca+2")]  = n_Ca * u"mol"
n0[idx("OH-")]   = 2n_Ca * u"mol"
n0[idx("XsOH")]  = N_sites * u"mol"
state0 = ChemicalState(cs, n0)
nothing # hide
```

## A slow family

The reaction is written with the free site first, which makes it the species
the reaction controls. [`sorption_rate`](@ref) gives it the elementary law
``r = k\,n_s\,a_{\mathrm{Ca^{2+}}}(1-\Omega)``, and the whole family goes to the
kinetic side.

```@example slow
rxn = Reaction(OrderedDict(free => 1, ca => 1), OrderedDict(bound => 1);
               symbol = "Ca on Hfo_s", equal_sign = '→')
kr  = KineticReaction(cs, rxn, sorption_rate(k, cs, rxn))

solver = EquilibriumSolver(cs, DiluteSolutionModel(), OptimaOptimizer())
kp  = KineticsProblem(cs, [kr], state0, (0.0, t_end); equilibrium_solver = solver)
ts  = collect(range(0.0, t_end; length = 81))
ks  = KineticsSolver(; ode_solver = Rodas5P(), reltol = 1.0e-10, abstol = 1.0e-20, saveat = ts)
diagnostics = IOBuffer() # hide
sol, states = with_logger(ConsoleLogger(diagnostics)) do # hide
sol = integrate(kp, ks)
states = speciated_states(sol, kp)
sol, states # hide
end # hide
log_text = String(take!(diagnostics)) # hide
occursin("re-speciation failed", log_text) && error("a re-speciation failed") # hide
occursin("could not be certified", log_text) && error("an instant was not certified") # hide
symbol.(cs.species[kp.idx_kinetic])      # the kinetic side: the whole family
```

With a thousand times more calcium than sites, the solution barely moves, and
the occupied amount follows the closed form of the theory page,
``n_c(t) = n_c^{\mathrm{eq}}(1 - e^{-\lambda t})`` with
``\lambda = k\,(a_{\mathrm{Ca^{2+}}} + 1/K)``.

```@example slow
a_Ca = exp(log_activities(states[1], DiluteSolutionModel())["Ca+2"])
λ    = k * (a_Ca + 10.0^-logK)
n_eq = k * N_sites * a_Ca / λ
bound_run    = [ustrip(us"mol", s.n[idx("XsOHCa+2")]) for s in states]
bound_closed = [n_eq * (1 - exp(-λ * t)) for t in ts]
@printf("λ = %.3e s⁻¹, n_eq/N = %.4f, largest gap = %.1e of the site budget\n",
        λ, n_eq / N_sites, maximum(abs.(bound_run .- bound_closed)) / N_sites)
```

```@example slow
plot(ts ./ 60, bound_run ./ N_sites; label = "integrated", lw = 2,
     xlabel = "time (min)", ylabel = "fraction of sites holding Ca")
plot!(ts ./ 60, bound_closed ./ N_sites; label = "closed form", ls = :dash, lw = 2)
```

## Beside a fast family

A second family with the same constant, `Hfo_f`, on a support of its own, that
no reaction controls: it stays on the equilibrium side and follows the solution
at every instant, while the slow one lags.

```@example slow
free_f  = site("XfOH", 0.0)
bound_f = site("XfOHCa+2", -RT * log(10.0^logK) + G_Ca)
fast = SiteFamily("Hfo_f", free_f, [bound_f]; capacity = TotalSiteAmount(N_sites * u"mol"),
                  support = SurfaceSupport("fast sites", nothing, FixedSurfaceArea(1.0)))
cs2 = ChemicalSystem([h2o, hp, oh, ca, free, bound, free_f, bound_f],
                     [h2o, hp, ca, free, free_f]; site_families = [family, fast])
idx2(s) = findfirst(==(s), symbol.(cs2.species))
n2 = Any[fill(0.0u"mol", length(cs2.species))...]
n2[idx2("H2O@")] = n0[idx("H2O@")]
n2[idx2("Ca+2")] = n_Ca * u"mol"
n2[idx2("OH-")]  = 2n_Ca * u"mol"
n2[idx2("XsOH")] = N_sites * u"mol"
n2[idx2("XfOH")] = N_sites * u"mol"
rxn2 = Reaction(OrderedDict(cs2.species[idx2("XsOH")] => 1, cs2.species[idx2("Ca+2")] => 1),
                OrderedDict(cs2.species[idx2("XsOHCa+2")] => 1); symbol = "Ca on Hfo_s", equal_sign = '→')
kp2 = KineticsProblem(cs2, [KineticReaction(cs2, rxn2, sorption_rate(k, cs2, rxn2))],
                      ChemicalState(cs2, n2), (0.0, t_end);
                      equilibrium_solver = EquilibriumSolver(cs2, DiluteSolutionModel(), OptimaOptimizer()))
sol2, states2 = with_logger(ConsoleLogger(diagnostics)) do # hide
sol2 = integrate(kp2, ks)
states2 = speciated_states(sol2, kp2)
sol2, states2 # hide
end # hide
log_text = String(take!(diagnostics)) # hide
occursin("re-speciation failed", log_text) && error("a re-speciation failed") # hide
occursin("could not be certified", log_text) && error("an instant was not certified") # hide
slow_frac = [ustrip(us"mol", s.n[idx2("XsOHCa+2")]) / N_sites for s in states2]
fast_frac = [ustrip(us"mol", s.n[idx2("XfOHCa+2")]) / N_sites for s in states2]
plot(ts ./ 60, fast_frac; label = "fast family", lw = 2,
     xlabel = "time (min)", ylabel = "fraction of sites holding Ca")
plot!(ts ./ 60, slow_frac; label = "slow family", lw = 2)
```

## A law that ignores the database

The same mass action, with a backward constant a tenth of what the equilibrium
constant gives: ``r = k\,n_s a_{\mathrm{Ca^{2+}}} - (k/10K)\,n_c``. It stops
where ``\Omega = 10``, beyond the equilibrium, and from ``\Omega = 1`` on it runs
against its affinity. [`dissipation`](@ref) lists those instants.

```@example slow
K = 10.0^logK
bad = KineticReaction(cs, rxn,
    (T, P, t, n, lna, n0) -> k * n["XsOH"] * exp(lna["Ca+2"]) - k / (10K) * n["XsOHCa+2"])
kp_bad  = KineticsProblem(cs, [bad], state0, (0.0, t_end); equilibrium_solver = solver)
d_good, d_bad = with_logger(ConsoleLogger(diagnostics)) do # hide
sol_bad = integrate(kp_bad, ks)
d_good = dissipation(sol, kp)
d_bad  = dissipation(sol_bad, kp_bad)
d_good, d_bad # hide
end # hide
log_text = String(take!(diagnostics)) # hide
occursin("re-speciation failed", log_text) && error("a re-speciation failed") # hide
occursin("could not be certified", log_text) && error("an instant was not certified") # hide
@printf("consistent law: %d violations; inconsistent law: %d violations\n",
        length(d_good.violations), length(d_bad.violations))
@printf("ln Ω at the end: consistent %.1e, inconsistent %.4f (ln 10 = %.4f)\n",
        -d_good.affinity[end, 1] / RT, -d_bad.affinity[end, 1] / RT, log(10))
```

```@example slow
lnΩ(d) = -d.affinity[:, 1] ./ RT
on = ts .> 0                       # at t = 0 no site is occupied and ln Ω is unbounded
against = [findfirst(==(v.time), ts) for v in d_bad.violations]
plot(ts[on] ./ 60, lnΩ(d_good)[on]; label = "consistent law", lw = 2,
     xlabel = "time (min)", ylabel = "ln Ω", ylims = (-4, 3), legend = :bottomright)
plot!(ts[on] ./ 60, lnΩ(d_bad)[on]; label = "backward constant ÷ 10", lw = 2)
scatter!(ts[against] ./ 60, lnΩ(d_bad)[against]; label = "rate against the affinity", ms = 3)
hline!([0.0, log(10)]; color = :black, ls = :dash, label = nothing)
```

Along the consistent law ``\ln\Omega`` rises to zero and stays there: the rate
has the sign of the affinity throughout and vanishes with it. Along the other it
passes zero, from which instant on the rate runs against the affinity (the
marked points), and settles at ``\ln 10``, where its own rate vanishes: a state
the equilibrium solver would not give.

## What this page does not show

- **That this surface is slow.** The rate constant is illustrative. A rate
  constant fitted on an observed slow uptake is an apparent parameter, which
  usually stands for a transport the zero-dimensional model does not resolve
  ([A slow surface](@ref sec-theory-pe-slow-surface)).
- **A family with fast and slow states.** Every state of a slow family changes
  through a declared reaction; a step much faster than the others is given a
  large constant, which makes the system stiff, rather than treated as
  equilibrated.
