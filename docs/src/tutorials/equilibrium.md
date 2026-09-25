# [Chemical Equilibrium](@id sec-equilibrium)

!!! info "Before this page"
    [Getting started](@ref sec-quickstart), and the Manual pages on
    [Species](@ref sec-species) and on
    [ChemicalSystem and ChemicalState](@ref sec-system-state). Why the
    calculation is the right one is explained in
    [Thermochemistry](@ref sec-theory-thermo).

ChemistryLab computes a thermodynamic equilibrium as the minimum of the Gibbs energy under the element-conservation constraints, and the workflow always follows the same four steps:

1. Build a [`ChemicalSystem`](@ref) (species + stoichiometric matrix).
2. Create an initial [`ChemicalState`](@ref) (temperature, pressure, initial amounts).
3. Call [`equilibrate`](@ref) (or use [`EquilibriumSolver`](@ref) explicitly).
4. Inspect the resulting [`ChemicalState`](@ref).

## Minimal workflow

The convenience function [`equilibrate`](@ref) handles everything with sensible defaults.
The example below computes the equilibrium state of calcite (CaCO₃) dissolving in mildly acidic water — a standard geochemical benchmark.

```@example eq_setup
using Optimization, OptimizationIpopt
using ChemistryLab
using DynamicQuantities

substances = build_species(datapath("slop98-inorganic-thermofun.json"); verbose = false)

# Select the carbonate-system species, calcite and its dissolution product Ca²⁺
dict = Dict(symbol(s) => s for s in substances)
species = [dict[sym] for sym in split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Cal")]

cs = ChemicalSystem(species, ["H2O@", "H+", "Ca+2", "CO3-2", "Zz"])
nothing # hide
```

```@raw html
<details><summary>The chemical system in full</summary>
```

```@example eq_setup
cs
```

```@raw html
</details>
```

```@example eq_setup
state = ChemicalState(cs)

# 1 mmol calcite dissolved in 1 L of acidic water (initial pH ≈ 4)
set_quantity!(state, "Cal",  1e-3u"mol")
set_quantity!(state, "H2O@", 1.0u"kg")

V = volume(state)
set_quantity!(state, "H+",  1e-4u"mol/L" * V.liquid)   # pH = 4
set_quantity!(state, "OH-", 1e-10u"mol/L" * V.liquid)  # charge seed

state_eq = equilibrate(state)
nothing # hide
```

```@raw html
<details><summary>The solved state in full — every species, with its amount</summary>
```

```@example eq_setup
state_eq
```

```@raw html
</details>
```

!!! tip "Quick shortcut"
    Calling `equilibrate(state)` with no extra arguments uses sensible defaults and is usually sufficient for aqueous geochemical problems.

The call above does more than it shows. [`equilibrate`](@ref) starts from every
interior-point back end that is loaded — Ipopt, and the one `OptimaSolver`
provides — refines each answer with the dual Newton method of `OptimaSolver`, and
returns the answer whose optimality it can prove. Which back ends exist, how they compare and how one of
them is imposed are described in [Solving an equilibrium](@ref sec-solving),
together with the constraints other than a fixed temperature and pressure, the
differentiation of an equilibrium with respect to its inputs, the activity models
and the declaration of solid solutions. This tutorial keeps to the calculation
itself.

## Certifying an answer

[`equilibrate`](@ref) takes the certifying route by default, and it returns a
composition together with a **proof** that it is the Gibbs minimum rather than
the point an iteration stopped at. Why such a proof exists — the problem is
convex, so the KKT conditions are sufficient and not merely necessary — what the
three conditions are, and how [`DualEquilibriumSolver`](@ref) aims at them
directly, are in
[Proving that an answer is the answer](@ref sec-theory-certificate).

Driving it explicitly, when the two stages are wanted separately:

```julia
des  = DualEquilibriumSolver(cs, HKFActivityModel())
ipm  = equilibrate(state, OptimaOptimizer())      # into the neighborhood 
dual = solve(des, ipm; b = b)                     # to the KKT conditions
cert = optimality_certificate(des, dual; b = b)
cert.optimal    # true: a proof, for a convex problem
```

## Inspecting the equilibrium state

The returned [`ChemicalState`](@ref) carries all derived thermodynamic quantities:

```@example eq_setup
println("pH      = ", pH(state_eq))
println("pOH     = ", pOH(state_eq))
println("porosity   = ", porosity(state_eq))
println("saturation = ", saturation(state_eq))
```

Phase volumes and mole amounts are accessible via named tuples:

```@example eq_setup
v = volume(state_eq)
println("V liquid = ", v.liquid)
println("V solid  = ", v.solid)
println("V total  = ", v.total)

m = moles(state_eq)
println("n liquid = ", m.liquid)
println("n solid  = ", m.solid)
```

Individual species amounts (in mol):

```@example eq_setup
cs_eq = state_eq.system
for (i, sp) in enumerate(cs_eq.species)
    n_i = state_eq.n[i]
    println(rpad(symbol(sp), 20), ustrip(n_i), " mol")
end
```

## Scaling and normalization

It is often useful to express a composition relative to a reference amount — per mole, per kilogram, or per cubic meter of system. Two mechanisms are provided.

### Scalar multiplication

A [`ChemicalState`](@ref) can be multiplied or divided by a real number. All molar amounts are scaled proportionally; temperature, pressure, and the chemical system are unchanged. The operation is **non-mutating** — a new state is returned:

```julia
state2  = state_eq * 2.0    # double all amounts
state_m = state_eq / 1000   # millimolar scale
```

### `rescale!` — rescale to a target total

[`rescale!`](@ref) scales all molar amounts **in-place** so that the total of the matching physical quantity equals `target`:

| `target` dimension | Quantity brought to `target` |
|:-------------------|:-----------------------------|
| mol                | `moles(state).total`         |
| kg (mass)          | `mass(state).total`          |
| m³ (volume)        | `volume(state).total`        |

All derived quantities (pH, porosity, volume, …) are recomputed automatically after scaling.

```@example eq_setup
# Express the equilibrium composition per kilogram of total system
state_pkg = copy(state_eq)
rescale!(state_pkg, 1.0u"kg")

println("Ca²⁺ = ", moles(state_pkg, "Ca+2"), "  mol/kg")
println("pH   = ", pH(state_pkg))   # intensive quantities are invariant
```

!!! note "Intensive quantities"
    pH, porosity, and saturation are **intensive** — they are invariant under homothety and remain unchanged after `rescale!` or scalar multiplication.

## Using `EquilibriumSolver` explicitly

For batch calculations where many different initial states share the same system and activity model, construct an [`EquilibriumSolver`](@ref) once and reuse it:

```julia
using Optimization, OptimizationIpopt

opt = IpoptOptimizer(
    acceptable_tol        = 1e-12,
    dual_inf_tol          = 1e-12,
    acceptable_iter       = 1000,
    constr_viol_tol       = 1e-12,
    warm_start_init_point = "no",
)

solver = EquilibriumSolver(
    cs,
    DiluteSolutionModel(),
    opt;
    variable_space = Val(:linear),
    abstol  = 1e-10,
    reltol  = 1e-10,
)
```

Once built, `solver` is called with any compatible `ChemicalState`:

```@example eq_setup
using Optimization, OptimizationIpopt #hide

opt = IpoptOptimizer( #hide
    acceptable_tol        = 1e-12, #hide
    dual_inf_tol          = 1e-12, #hide
    acceptable_iter       = 1000, #hide
    constr_viol_tol       = 1e-12, #hide
    warm_start_init_point = "no", #hide
) #hide

solver = EquilibriumSolver( #hide
    cs, #hide
    DiluteSolutionModel(), #hide
    opt; #hide
    variable_space = Val(:linear), #hide
    abstol  = 1e-10, #hide
    reltol  = 1e-10, #hide
) #hide
state_eq2 = solve(solver, state)
```

!!! note "Performance"
    The potential function `μ(n, p)` is compiled once during `EquilibriumSolver` construction. Repeated calls to `solve(solver, ...)` with different states reuse it, avoiding redundant compilation overhead.

## Temperature dependence (10–30 °C)

Calcite solubility varies with temperature. Using the `solver` built above, we sweep from 10 to 30 °C and track pH, dissolved calcium and remaining solid calcite:

```@example eq_setup
using Plots


temperatures = 10:30   # °C

pH_vals   = Float64[]
nCa_vals  = Float64[]  # mmol
nCal_vals = Float64[]  # mmol

i_Ca  = findfirst(sp -> symbol(sp) == "Ca+2", cs.species)
i_Cal = findfirst(sp -> symbol(sp) == "Cal",  cs.species)

# Start from the charged state built above, not from `ChemicalState(cs)`:
# a fresh state holds no matter at all, so every element balance would be zero
# and the sweep would return the same trivial solution at all 21 temperatures.
for θ in temperatures
    s_T = deepcopy(state)
    set_temperature!(s_T, (273.15 + θ) * u"K")
    s_eq = solve(solver, s_T)
    push!(pH_vals,   pH(s_eq))
    push!(nCa_vals,  ustrip(s_eq.n[i_Ca]) * 1e3)
    push!(nCal_vals, ustrip(s_eq.n[i_Cal]) * 1e3)
end
```

The figures can then be drawn.

```@example eq_setup
p1 = plot(collect(temperatures), pH_vals,
    xlabel = "T (°C)", ylabel = "pH", label = "pH",
    marker = :circle, linewidth = 2, title = "pH")
p2 = plot(collect(temperatures), nCa_vals,
    xlabel = "T (°C)", ylabel = "n (mmol)", label = "Ca²⁺",
    marker = :circle, linewidth = 2, title = "Dissolved species")
plot!(p2, collect(temperatures), nCal_vals,
    label = "Cal", marker = :square, linewidth = 2)
plot(p1, p2, layout = (1, 2), left_margin = 8Plots.mm, bottom_margin = 8Plots.mm, size = (900, 400))
```

!!! note "Retrograde Kₛₚ, and yet more dissolved calcium"
    Calcite is **retrograde soluble**: its solubility product falls as temperature
    rises. The sweep above reproduces that — the ionic product ``[\text{Ca}^{2+}]
    [\text{CO}_3^{2-}]`` it settles on goes from ``10^{-8.411}`` at 10 °C to
    ``10^{-8.517}`` at 30 °C, matching the ``K_{sp}`` the database itself gives to
    within 0.003 log units, and ``-8.48`` at 25 °C is the accepted value for
    calcite.

    The dissolved calcium nevertheless **increases**, from 0.146 to 0.158 mmol/L.
    That is not a contradiction: this system is **closed**, with 1 mmol of total
    carbonate and no CO₂ reservoir. As temperature rises the pH drops from 9.83
    to 9.43, which shifts carbonate to bicarbonate and takes free CO₃²⁻ from 26.5
    down to 19.2 µmol/L — more than enough to offset the smaller ``K_{sp}``. The
    familiar statement that retrograde solubility means *less* dissolved calcium
    holds for a system buffered at a fixed CO₂ partial pressure, not for a sealed
    one. Fix the partial pressure instead (add a gas phase) and the sign flips.

## [Reading the aqueous properties back](@id sec-aqueous-properties)

The activity model used so far is the default, [`DiluteSolutionModel`](@ref);
the others, and the reasons for choosing one, are in
[Activity models](@ref sec-activity-models) (syntax) and
[Activity models](@ref sec-theory-activity) (theory). Whatever the model, the
quantities it computes can be read back from a state.

The activity closures compute the molalities, the ionic strength and the
activity coefficients on their way to the log-activities. All of it is readable
off a state:

| call | returns |
|:--|:--|
| [`molalities`](@ref)`(state)` | `mᵢ` of every solute, mol/kg of solvent |
| [`ionic_strength`](@ref)`(state)` | `I = ½ Σ mⱼ zⱼ²`, mol/kg |
| [`activity_coefficients`](@ref)`(state, model)` | `γᵢ` of every aqueous species |
| [`log_activities`](@ref)`(state, model)` | `ln aᵢ` of **every** species |
| [`activities`](@ref)`(state, model)` | `aᵢ` of every species |
| [`pH`](@ref)`(state, model)` | `−log₁₀ a(H⁺)` |
| [`pOH`](@ref)`(state, model)` | `−log₁₀ a(OH⁻)` |

```julia
eq = equilibrate(state; model = HKFActivityModel())

ionic_strength(eq)                      # 0.212 mol/kg
molalities(eq)["K+"]                    # 0.146 mol/kg
activity_coefficients(eq, model)["Ca+2"]
activities(eq, model)["H2O@"]           # water activity
pH(eq, model)                           # activity convention
```

`molalities` and `ionic_strength` need no model: they are properties of the
composition, and every model in the package computes the ionic strength this
way. When comparing against another code, **compare the ionic strength first** —
if it disagrees, the two are not describing the same solution, whatever their
volumes happen to agree on.

### Two conventions of pH, 0.2 units apart

The one-argument [`pH`](@ref)`(state)` and the two-argument
[`pH`](@ref)`(state, model)` are **different quantities**:

- `pH(state)` is `−log₁₀ c(H⁺)`, a **concentration** in mol/L over the computed
  liquid volume, and in an alkaline solution it is reconstructed from OH⁻
  through `pKw`. It is stored on the state and needs no activity model.
- `pH(state, model)` is `−log₁₀ a(H⁺)`, the **activity** on the molality scale.
  This is what GEM-Selektor, PHREEQC and Reaktoro report.

On a Portland cement pore solution at `I ≈ 0.2 mol/kg`, with `γ(H⁺) ≈ 0.61`, the
two differ by about **0.21 units** (13.31 against 13.10). Comparing the wrong
one against another code means chasing a discrepancy that is a convention, not a
result.

### γ comes from the formula, not from a ratio

[`activity_coefficients`](@ref) evaluates the model's own expression, rather than
dividing an activity by a concentration. The ratio agrees for an abundant
solute — the test suite checks that it does — but a species parked at the
solver's `1e-16 mol` lower bound has its log-activity dominated by the closures'
`+ ϵ` regularization, and the ratio then returns values of order `1e300` for a
charge class whose only members are trace. The formula depends on the ionic
strength and the charge alone, so it is exact at any amount.

Relatedly, an activity is a number on a *scale*, and nothing in the number says
which. [`DiluteSolutionModel`](@ref) puts its solutes on the molarity scale and
takes `ρ = 1 kg/L`, so its activities coincide numerically with molalities even
though the scale differs; the other two models are on the molality scale.
[`concentration_scale`](@ref) is the only way to tell them apart.

## Where to go next

The options this tutorial left at their defaults are gathered in
[Solving an equilibrium](@ref sec-solving). The application pages then apply the
same four steps to systems of increasing size: the aqueous cases, beginning with
[CO₂ dissolution and the carbonate system](@ref sec-co2-carbonate), are small
enough to check by hand, and [A CEM I from its clinker phases](@ref sec-cem1-from-clinker)
is the first cement. When the question is how a system evolves in time rather
than where it ends, the next tutorial is [Chemical Kinetics](@ref sec-kinetics).
