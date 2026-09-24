# [Rate laws, and every parameter in them](@id sec-theory-kinetics)

A Gibbs minimization answers *what is stable*. It says nothing about *when*, and
a cement paste is a material whose whole engineering behavior lives in the when.
This page is the catalog of what supplies that: every rate law the package
ships, the physics each one encodes, the parameters each one carries, their
units, and where the numbers come from.

It is written to be read beside [the kinetics tutorial](@ref sec-kinetics)
rather than instead of it: the tutorial says how to call these, this page says
what they mean and what they are worth.

```@example kin
using ChemistryLab
using DynamicQuantities
using Printf
using Plots
default(framestyle = :box, grid = false)
nothing # hide
```

## 0. What a rate law is, here

Every rate in this package is a [`KineticFunc`](@ref) with one signature:

```julia
r = f(T, P, t, n::StateView, lna::StateView, n_initial::StateView)   # mol/s
```

Six arguments, and the choice of six is the point. A rate may depend on the
temperature and pressure, on the elapsed time, on **the current composition**,
on **the activities** of any species in the system, and on the initial
composition — which is what lets a law be written in terms of a degree of
reaction, `α = 1 − n/n₀`, without the caller having to track it.

Two consequences worth stating plainly:

  - **A rate law may read any species, not only the ones in its own reaction.**
    A dissolution rate catalyzed by `H+`, an SCM whose reaction is gated by the
    portlandite left, a mechanism inhibited by the sulfate in solution — all of
    these are written by reading `n` or `lna` at the species concerned. The
    reaction a rate is attached to fixes what the rate *consumes*; it does not
    restrict what the rate may *depend on*. One consequence has its own section
    in the tutorial, because it bites: a rate that reads a species the reaction
    itself consumes to exhaustion — see
    [rate laws that depend on a consumed reactant](@ref kinetics-frozen-species).
  - **The evaluation path is differentiable.** No `Float64` casts, so a
    `ForwardDiff.Dual` propagates through `T`, through `t`, through the
    composition and through the parameters. That is what makes a rate constant
    fittable ([the calibration page](@ref ex-hydration-calibration)) rather than
    only adjustable by hand.

There are two families below, and they answer different questions. The first is
**mechanistic**: a rate proportional to how far the solution is from equilibrium
with the mineral. The second is **empirical**: a degree of reaction as a function
of time, fitted to calorimetry. Clinker and supplementary materials are described
by the second, because nobody can write the first for a multiphase glass.

## 1. Distance from equilibrium: the Palandri–Kharaka form

For a mineral dissolving into or precipitating from a solution, the rate is
driven by the saturation ratio

```math
\Omega = \frac{\mathrm{IAP}}{K}, \qquad
\ln\Omega = \sum_i \nu_i \ln a_i - \ln K,
\qquad \ln K = -\frac{\Delta_r G^\circ}{RT} ,
```

computed by [`saturation_ratio`](@ref) from the stoichiometry, the log-activities
and the standard Gibbs energies the database supplies. `Ω < 1` is
undersaturation and the mineral dissolves; `Ω > 1` and it precipitates; `Ω = 1`
is equilibrium and the rate is exactly zero — which is the property that makes
this form composable with an equilibrium solver rather than in competition with
it.

One mechanism ([`RateMechanism`](@ref)) contributes

```math
r_{\text{mech}} = k(T)\;\Bigl[\prod_j a_j^{\,n_j}\Bigr]\;
                  \operatorname{sign}(1-\Omega)\;\bigl|1-\Omega^{p}\bigr|^{q} ,
```

and a mineral's total rate is the sum over its mechanisms — classically an acid
one, a neutral one and a base one, which is why the products of activities are
called **catalysts** ([`RateModelCatalyst`](@ref)): each contributes `a_j^{n_j}`,
with `n_j = 0.5` on `H+` a typical acid mechanism.

| symbol | meaning | unit |
|:--|:--|:--|
| ``k(T)`` | rate constant, usually [`arrhenius_rate_constant`](@ref) | mol m⁻² s⁻¹ |
| ``a_j, n_j`` | catalyst activity and its exponent | — |
| ``p`` | saturation exponent inside the bracket | — |
| ``q`` | outer exponent | — |
| ``\Omega`` | saturation ratio | — |

The rate constant itself is the Arrhenius law referred to a reference
temperature rather than to an absolute prefactor,

```math
k(T) = k_0 \exp\!\left[-\frac{E_a}{R}\left(\frac 1T - \frac 1{T_{\text{ref}}}\right)\right],
```

so that `k₀` is the measured constant *at* `T_ref` and carries its own unit,
instead of being an extrapolation to infinite temperature. The shipped parameter
sets for this family come from [PalandriKharaka2004](@cite).

## 2. Clinker hydration: Parrot & Killoh

A clinker phase does not dissolve into a solution at a rate set by its own
undersaturation — it is consumed behind a growing layer of hydrate, and the
controlling step changes as the layer thickens. The 1984 description, as reported
by [Lothenbach2008](@cite) and used by [Lavergne2018](@cite), writes three
competing mechanisms in the **degree of hydration** and lets the slowest one
limit:

```math
\dot\alpha_1 = \frac{k_1}{n_1}(1-\xi)\bigl[-\ln(1-\xi)\bigr]^{1-n_1},
\qquad
\dot\alpha_2 = k_2\,\frac{(1-\xi)^{2/3}}{1-(1-\xi)^{1/3}},
\qquad
\dot\alpha_3 = k_3 (1-\xi)^{n_3},
```

```math
\dot\alpha = \min(\dot\alpha_1,\dot\alpha_2,\dot\alpha_3),
\qquad \xi = \alpha/\alpha_{\max} .
```

Each has a physical reading. ``\dot\alpha_1`` is **nucleation and growth** in
Avrami form — the logarithm is the signature of a transformed fraction growing
from nuclei. ``\dot\alpha_2`` is **diffusion** through the hydrate shell, in
Jander's geometry for a shrinking sphere. ``\dot\alpha_3`` is **shell
formation**, a power law in the fraction left. The rate returned is

```math
r = n_{\text{initial}} \times A_T \times \beta_B \times \beta_h \times \dot\alpha
\qquad [\text{mol/s}],
```

the three multipliers being the corrections of section 4.

```@example kin
sets = ["C₃S" => PK84_PARAMS_C3S, "C₂S" => PK84_PARAMS_C2S,
        "C₃A" => PK84_PARAMS_C3A, "C₄AF" => PK84_PARAMS_C4AF]
@printf("%-6s %8s %6s %10s %8s %6s %10s\n",
        "phase", "k1 [1/d]", "n1", "k2 [1/d]", "k3 [1/d]", "n3", "Ea [kJ/mol]")
for (nm, p) in sets
    @printf("%-6s %8.2f %6.2f %10.3f %8.2f %6.1f %10.0f\n",
            nm, ustrip(us"1/d", p.k₁), p.n₁, ustrip(us"1/d", p.k₂),
            ustrip(us"1/d", p.k₃), p.n₃, ustrip(us"J/mol", p.Ea) / 1000)
end
```

Two features of that table are artifacts of the 1984 fit rather than chemistry,
and both are useful as checks on an implementation. With ``n_1 = 1`` the Avrami
branch of **belite** reduces to ``k_1(1-\xi)``, which never limits, so C₂S is
governed by the power law throughout. And **alite** never reaches its
diffusion-controlled stage. The original authors acknowledged both.

!!! warning "Two Parrot–Killoh variants ship, and only one is attributed"
    [`parrot_killoh`](@ref) is a *different*, smoothed variant that predates
    this one in the package (the two are set side by side in
    [the tutorial](@ref pk-variants)) — `min(max(r_{NG}, r_I), r_D)` with a damped
    nucleation term — and its `PK_PARAMS_*` are not transferable to the
    canonical law. It is deprecated and **the attribution to Parrott & Killoh
    was withdrawn** rather than repaired: its nucleation term carries no Avrami
    logarithm, its shell coefficient sits in the diffusion expression, and no
    published set matches its parameters. The primary source is a conference
    proceedings without a DOI that could not be consulted, so the honest action
    was to stop claiming it. Measured consequence: with `PK_PARAMS_*` all four
    phases land on the diffusion branch almost immediately and a CEM I at
    w/c 0.40 reaches α ≈ 0.234 at seven days against the ≈ 0.61 the literature
    reports. Use [`parrot_killoh_avrami`](@ref).

## 3. Supplementary materials: the Waller sigmoid

A slag or a fly ash has no phases, no formula and no single dissolution
mechanism, so its reaction is described by the shape it is observed to have — a
sigmoid in **log** time:

```math
\alpha(t) = \frac{1}{1 + (\tau/t)^{n}} ,
\qquad
\dot\alpha = \frac{n}{\tau}(1-\alpha)^{1+1/n}\,\alpha^{1-1/n} .
```

The second form is what the package evaluates, because writing the rate as a
function of the *current degree* rather than of the clock is what lets the same
temperature, fineness and humidity corrections multiply it as they multiply the
clinker rate. ``\tau`` is the time at which the material is half reacted and
``n`` the sharpness of the transition.

```@example kin
waller_sets = ["fly ash" => WALLER_PARAMS_FLY_ASH, "silica fume" => WALLER_PARAMS_SILICA_FUME,
               "slag" => WALLER_PARAMS_SLAG]
α_waller(p, t_days) = 1 / (1 + (ustrip(us"d", p.τ) / t_days)^p.n)

@printf("%-12s %8s %6s %12s %10s %10s\n",
        "material", "τ [d]", "n", "Ea [kJ/mol]", "α(28 d)", "α(90 d)")
for (nm, p) in waller_sets
    @printf("%-12s %8.0f %6.2f %12.2f %10.3f %10.3f\n",
            nm, ustrip(us"d", p.τ), p.n, ustrip(us"J/mol", p.Ea) / 1000,
            α_waller(p, 28.0), α_waller(p, 90.0))
end
```

Pozzolanic and latent-hydraulic reactions are markedly more **temperature**
sensitive than the clinker's: 83.14 kJ/mol against 21–54 kJ/mol above. That is
why a blended cement gains so much from a warm cure and loses so much in a cold
one, and it is a prediction of the parameter and not an extra rule.

!!! note "These are not the same numbers as a direct measurement of reacted glass"
    The α(28 d) column above is ≈ 0.29 for slag and ≈ 0.32 for fly ash. The
    RILEM TC 238-SCM round robin [Durdzinski2017](@cite), which measured the
    *reacted glass* directly on its own materials rather than inferring it from
    heat, reports 38–49 % for two slags and about 20 % for a siliceous fly ash
    at the same age. The two disagree, and in opposite directions.

    Neither is wrong. Waller's parameters are a fit to particular materials, and
    "a slag" is not a substance — its reactivity depends on its glass content,
    its basicity and its fineness. The round robin's own conclusion is the one to
    keep: the precision of *any* determination of an SCM's degree of reaction is
    "rather low, at best ± 4-5 %". So a reacted fraction is an input to be stated
    and swept, never a constant to be trusted to two digits — which is exactly
    how [the blended binder pages](@ref sec-theory-water-budget) treat it.

Silica fume carries the same τ and n as fly ash: its much higher reactivity is
represented **through the fineness**, at an effective Blaine of 2000 m²/kg
recommended by [Lavergne2018](@cite). Its BET surface, about 20 000 m²/kg, is a
different measurement of a different thing and must not be substituted — a factor
of ten on the rate.

That is now enforced rather than advised. [`BETSurfaceArea`](@ref) and
[`BlaineSurfaceArea`](@ref) are distinct types, and the ratio that produces a
fineness factor is defined only between two of the same kind, so handing a BET
area to [`blaine_factor`](@ref) raises instead of returning a plausible number.
A bare quantity is still read as a Blaine fineness, which is the contract every
existing call relies on; see [Surface areas](@ref sec-manual-surfaces).

## 4. The three multiplicative corrections

Each of the empirical laws above is multiplied by three dimensionless factors.
They are separated because they answer separate questions, and because each can
be switched off by leaving its keyword at `nothing`.

**Temperature**, the same Arrhenius factor as section 1:

```math
A_T = \exp\!\left[-\frac{E_a}{R}\left(\frac 1T - \frac 1{T_{\text{ref}}}\right)\right] .
```

**Fineness** ([`blaine_factor`](@ref)): the rate scales linearly with the
specific surface, relative to the fineness the parameters were fitted at —
385 m²/kg for the clinker phases, 400 m²/kg for the Waller sets.

```math
\beta_B = \frac{B}{B_{\text{ref}}} .
```

**Internal humidity** ([`humidity_factor`](@ref)): hydration stops when the pore
water is no longer available, which is taken to happen below 80 % RH.

```math
\beta_h = \left(\frac{h-0.55}{0.45}\right)^{4} \quad (h > 0.80),
\qquad \beta_h = 0 \quad \text{otherwise.}
```

```@example kin
@printf("β_B at 462 m²/kg (clinker ref) : %.3f\n", blaine_factor(462u"m^2/kg"))
@printf("β_B at 2000 m²/kg (Waller ref) : %.3f\n",
        blaine_factor(2000u"m^2/kg"; blaine_ref = 400u"m^2/kg"))
for h in (0.99, 0.90, 0.81, 0.801, 0.80, 0.75)
    @printf("β_h(%.3f) = %.4f\n", h, humidity_factor(h))
end
```

That cut is a **genuine discontinuity**, and the page says so rather than
smoothing it: the one-sided limit from above is 0.0953 and the value at 0.80 is
exactly 0. It is mild in practice — the rate has already fallen by an order of
magnitude from β_h(0.99) = 0.914 — but a solver stepping across it will feel it,
which is why [`PoreHumidity`](@ref) exists to supply `h` from the current
saturation rather than from a schedule.

### When the area stops being constant

In every law shipped here ``\beta_B`` is a **number**, computed once from the
fineness given at construction and carried through the whole integration. That
is not an oversight: it is the form the Parrot–Killoh and Waller constants were
fitted in, so it is the form in which those constants mean what they say.

A grain that dissolves does not keep its area. [`ShrinkingCoreArea`](@ref) says
so, and passing one where a fineness is expected makes the factor follow the
amount left,

```math
\beta_B(\xi) = \beta_B^0 \, g\!\left(\frac{n}{n_0}\right),
\qquad g(f) \simeq f^{\,p},
\qquad \frac{n}{n_0} = 1-\xi \ \ (\alpha_{\max}=1),
```

with ``p = 2/3`` the geometric value for spheres. The evolving factor is exactly
``\beta_B^0`` at ``n = n_0``, so nothing changes at the first instant; what
changes is everything after it.

**And that freedom is not new everywhere.** Multiply the shell-formation branch
by ``(1-\xi)^p``:

```math
k_3 (1-\xi)^{n_3} \cdot (1-\xi)^{p} = k_3 (1-\xi)^{\,n_3+p} ,
```

so wherever that branch is the active one, ``p`` and ``n_3`` are the same
parameter written twice and no amount of data separates them. The Jander branch
carries no such exponent, so there ``p`` is a genuine new shape. The Waller
sigmoid is the interesting case: its ``n`` sets ``(1-\xi)^{1+1/n}`` and
``\xi^{1-1/n}`` together, in opposite directions, while ``p`` moves only the
first — so ``p`` is distinguishable from ``n``, though correlated with it.

Which of these holds over a given dataset is a measurement and not an argument.
[`identifiability`](@ref) is what makes it, and
[Where the numbers come from](@ref sec-manual-numbers) is where a parameter that
the data turn out not to constrain gets labeled as such instead of quoted.

!!! warning "An evolving area needs its own calibration"
    The published constants were fitted with ``\beta_B`` frozen. Turning the
    area into a function of the state moves the law out of that fit, and the
    honest response is to recalibrate — `scripts/hydration_calibration.jl` is
    the machinery — not to keep the old constants and add a new factor on top.

## 5. How far the reaction can go: `α_max`

Every law above is written in ``\xi = \alpha/\alpha_{\max}``, and
``\alpha_{\max}`` is the one parameter that is not about speed at all. It is the
water-and-space ceiling of [the water budget page](@ref sec-theory-water-budget),
supplied by [`powers_alpha_max`](@ref):

```@example kin
@printf("%6s %12s %12s\n", "w/c", "sealed", "under water")
for wc in (0.25, 0.30, 0.36, 0.40, 0.42, 0.50)
    @printf("%6.2f %12.3f %12.3f\n", wc,
            powers_alpha_max(wc), powers_alpha_max(wc; curing = :saturated))
end
```

Read the last two rows: above 0.42 neither convention does anything, because
water has stopped being what limits the reaction. Read the first: at w/c = 0.25 a
sealed paste cannot pass 60 % hydration however long it is left, and a cured one
cannot pass 69 %. The difference between the columns is the chemical shrinkage,
refilled from the bath in one case and not in the other.

## 6. Which numbers are published and which are fitted

| quantity | where it comes from |
|:--|:--|
| `PK84_PARAMS_*` | Parrott & Killoh (1984) as reported by [Lothenbach2008](@cite), [Lavergne2018](@cite) |
| `WALLER_PARAMS_*` | Waller (1999), as used by [Lavergne2018](@cite) |
| the dissolution rate constants and their exponents | [PalandriKharaka2004](@cite) |
| `blaine_ref` 385 / 400 m²/kg | the finenesses those fits were made at |
| the 0.80 humidity cut and its exponent | Parrot et al., as used by van Breugel |
| `powers_alpha_max` 0.42 / 0.36 | [Powers1948](@cite) |
| `CALIBRATED_THETA` | **fitted here**, on one record — [the calibration page](@ref ex-hydration-calibration) |
| an SCM's reacted fraction at a given age | **an input**, measured or assumed — never a constant of the code |

The last two rows are the ones to keep in view. Everything above them is
somebody's published fit to somebody's materials; the calibration page shows what
happens when those published parameters meet a calorimetry record they were not
fitted to, and how much of the gap a five-parameter fit closes.

## See also

  - [The water budget of a hydrating paste](@ref sec-theory-water-budget) — where
    ``\alpha_{\max}`` comes from, and why it is not thermodynamics
  - [Calibrating a hydration model](@ref ex-hydration-calibration) — the fit, the
    holdout, and the sensitivity
  - [A CEM I from its clinker phases](@ref sec-cem1-from-clinker) — all of this
    running, over ninety days
