# Validation against Reaktoro

!!! info "Before this page"
    The tutorial [Chemical Equilibrium](@ref sec-equilibrium), and [Coupling
    kinetics and equilibrium](@ref sec-coupling) for the last section.

[Reaktoro](https://reaktoro.org) solves the same problem this package does, from
the same literature [Leal2017](@cite), and is mature and widely used. It makes
a good reference — provided the comparison is set up so that a disagreement means
something.

This page reports what agrees, what does not, and by how much. The scripts that
produce the reference values live in `test/reference/` and the assertions in
`test/equilibrium_reference.jl`, so every number below is reproducible and is
checked on every CI run.

## Setting up a comparison that means something

A cross-code comparison has **three** parameters, and each is worth tens of percent:

| parameter | must match | what happens otherwise |
|:--|:--|:--|
| thermodynamic database | same file | different standard-state data, different answer |
| species list | same species | see below — this one is subtle |
| activity model | same convention | `∂Ca²⁺/∂(CO₂)` moves from `+0.152` to `+0.218` |

Reaktoro reads ThermoFun databases, so it can be pointed at **the very file this
package reads**, `cemdata18-thermofun.json` as `datapath` obtains it, with
identical species names.
This package's default is `DiluteSolutionModel()`, matched on the Reaktoro side
by `ActivityModelIdealAqueous()`.

!!! warning "The species list is not a detail"
    Drop the aqueous calcium complexes and free `Ca²⁺` becomes the only aqueous
    home for calcium, so `∂Ca²⁺/∂(CO₂)` and `∂calcite/∂(CO₂)` mirror each other
    exactly. Restore `CaOH⁺`, `Ca(CO₃)@` and `Ca(HCO₃)⁺` and that mirror breaks
    — the difference is what `Ca(HCO₃)⁺` takes up. **Both codes break it the
    same way**, and the element balance closes to `2×10⁻¹⁶` either way. Comparing
    a system that mirrors against one that does not is comparing two different
    chemistries, not validating one.

## The test case

Calcite in water with dissolved CO₂ — small enough to reason about, and it
exercises a solid in equilibrium with an aqueous phase spanning ten orders of
magnitude in amount:

```
H2O@ 55.5 mol   ·   Cal 0.05 mol   ·   CO2@ 0.01 mol   ·   25 °C, 1 bar
species: H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 CaOH+ Ca(CO3)@ Ca(HCO3)+ Cal
```

## What agrees

**Major species** — everything present above about `10⁻⁵` mol — agrees to
`10⁻³` relative or better.

**Sensitivities.** `∂n/∂(CO₂)` obtained by differentiating *through* the solve
(see [Differentiating an equilibrium](@ref)) matches Reaktoro's finite
differences on every species:

| species | ChemistryLab (AD) | Reaktoro (FD) | rel. diff |
|:--|--:|--:|--:|
| H₂O | −0.181336 | −0.181397 | 3.4×10⁻⁴ |
| Ca²⁺ | +0.151920 | +0.151987 | 4.4×10⁻⁴ |
| HCO₃⁻ | +0.333308 | +0.333427 | 3.6×10⁻⁴ |
| CO₂(aq) | +0.818655 | +0.818600 | 6.7×10⁻⁵ |
| Ca(HCO₃)⁺ | +0.029333 | +0.029338 | 1.6×10⁻⁴ |
| calcite | −0.181251 | −0.181324 | 4.0×10⁻⁴ |

Reaktoro's own spread across `h ∈ {10⁻³, 10⁻⁴, 10⁻⁵}` is `7.2×10⁻⁴`, so the
residual difference sits **below the truncation error of the reference itself**. This is as
close as a comparison against a finite-difference reference can get.

**Element balance.** Calcium has no source in this system, so the
sensitivities of all Ca-bearing species must cancel. The implicit-function method
gives `2×10⁻¹⁶`, and a species absent from the solution gets identically zero
(`2×10⁻⁹` in the finite-difference reference, the size of its truncation error).

## The trace species

They agree too, and it took two fixes to get there.

| species | amount | ChemistryLab | Reaktoro | ratio |
|:--|--:|--:|--:|--:|
| H₂O, Cal, Ca²⁺, HCO₃⁻, CO₂, Ca(HCO₃)⁺, Ca(CO₃)@ | > 10⁻⁵ | — | — | 1.0000 ± 0.0002 |
| CO₃²⁻ | 9×10⁻⁷ | 9.3814×10⁻⁷ | 9.3786×10⁻⁷ | 1.0003 |
| H⁺ | 4×10⁻⁷ | 3.6923×10⁻⁷ | 3.6929×10⁻⁷ | 0.9998 |
| OH⁻ | 3×10⁻⁸ | 2.7109×10⁻⁸ | 2.7063×10⁻⁸ | 1.0017 |
| CaOH⁺ | 2×10⁻⁹ | 1.6365×10⁻⁹ | 1.5855×10⁻⁹ | 1.0321 |

`pKw` on this system is **13.9994**, against Reaktoro's **14.0001**. Pure water
gives `[H⁺]/[OH⁻] = 1.0`.

`CaOH⁺`, the smallest amount in the system at 1.6 nmol, is the only species
above 0.2 % out, and it is 3.2 % out. Nothing in the comparison is marked broken.

### How the back-end treats the trace species

Two features of the Newton iteration of `OptimaSolver` matter for the trace species.

The step is computed by the nullspace method. The Schur complement
`S = A H⁻¹ Aᵀ` needs `H` invertible, whereas a pure phase has unit activity,
hence `∂²G/∂nᵢ² = 0` exactly. With `Z` a basis of `null(A)`, the step solves

```math
(\mathbf{Z}^\mathsf{T} \mathbf{H} \mathbf{Z})\, \delta\mathbf{z} = -\mathbf{Z}^\mathsf{T}(\mathbf{e}_x + \mathbf{H}\, \delta\mathbf{n}_p),
```

in which `H` appears only as a product. This is the default of the C++ Optima,
whose `Rangespace` counterpart is documented there as suitable for invertible
diagonal Hessians only. It gives `[H⁺]/[OH⁻] = 1.0` for pure water.

A variable is judged "at its bound" when its slack falls below a threshold
relative to the variable's own bound, not to the largest variable in the
problem. With a solvent at 55 mol, a threshold scaled by the largest slack is
`5.5×10⁻⁷`, so every trace ion below it would be declared to sit on a bound of
`10⁻¹⁶`, nine orders of magnitude away, and its stationarity would never be
enforced.

The convergence test uses the residual evaluated at `μ = 0`, which is the true
KKT error. The residual at the current barrier level `μ` vanishes at the
solution of the barrier subproblem whatever `μ` is, so a test on it can report
success `O(μ)` away from the actual optimum, and the trace species are where
that shows first. The residual at `μ = 0` cannot be met at a loose barrier.

The accuracy obtained on the smallest amounts is:

| species | amount | ratio to Reaktoro |
|:--|--:|--:|
| CaOH⁺ | 3×10⁻⁸ | ×1.032 |
| OH⁻ | 9×10⁻⁸ | ×1.002 |
| H⁺ | 5×10⁻⁷ | ×1.000 |
| CO₃²⁻ | 1×10⁻⁶ | ×1.000 |
| Ca(CO₃)@ | 6×10⁻⁶ | ×1.000 |

Agreement between two back ends run under the same stopping rule is not
evidence of accuracy: two barrier iterations stopped short of stationarity at
`μ = 0` can agree with each other and not with the exact solution. Comparing the
Gibbs energies of two answers, rather than their compositions, tells whether the
solver failed to minimize its own objective.

## The kinetics/equilibrium coupling

**Validated against Reaktoro.** Calcite dissolves at a *constant* rate, so the
kinetic half of the system of [Leal2015](@citet) is analytic —

```math
n_{\text{Cal}}(t) = n_{\text{Cal}}(0) - k t,
\qquad
\mathbf{b}_e(t) = \mathbf{b}_e(0) + k t \, (\text{one CaCO}_3)
```

— which removes every difference of rate-law convention between the codes and
leaves only the question that matters: does the equilibrium partition follow
``\mathbf{n}_e(t) = \varphi(\mathbf{b}_e(t))``? The reference then needs no kinetics integrator at
all, just one equilibrium solve per sample time.

| t (s) | worst relative deviation | on |
|--:|--:|:--|
| 0 | 0.32 % | OH⁻ |
| 600 | 0.44 % | CO₂(aq) |
| 1800 | 4.3 % | CO₂(aq) |
| 3600 | 4.3 % | CO₂(aq) |

The species carrying the largest deviation sits at `2.3×10⁻⁸` mol throughout.
The assertions are in `test/coupling_reference.jl`, the reference solution in
`test/reference/reaktoro_coupling.py`.

Internal consistency is checked separately and independently: over a seven-day
hydration run the partition constraint holds to `‖Aₑnₑ − bₑ‖∞ = 8.7×10⁻⁸` with
no failed re-speciation.

## Reproducing

```sh
conda create -n reaktoro-env -c conda-forge reaktoro thermofun
conda run -n reaktoro-env python test/reference/reaktoro_calcite_co2.py
```

Paste the output into `REAKTORO` in `test/equilibrium_reference.jl`. Every
assertion there is a plain `@test` and all 26 pass — species above `10⁻⁵` mol to
`10⁻³` relative, trace ions to 5 %. There is no `@test_broken` left in the
comparison; if one is ever needed again, the tolerance it hides belongs in the
comment beside it.

## Where to go next

The same comparison on the derivatives of an equilibrium is reported in
[Differentiating an equilibrium](@ref), and the cross-code comparisons of the
surface models are in [Two families of sites, and a metal between them](@ref sec-example-hfo)
and [A published clay sorption model](@ref sec-example-claysor). The
application pages follow, beginning with
[Bogue Calculation](@ref ex-bogue).
