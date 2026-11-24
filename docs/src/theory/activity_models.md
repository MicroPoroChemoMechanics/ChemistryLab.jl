# [Activity models](@id sec-theory-activity)

!!! info "Before this page"
    [Thermochemistry](@ref sec-theory-thermo) and [Standard states](@ref
    sec-theory-standard-states), which define the activity this page models.

Every equilibrium this package computes rests on

```math
\mu_i = \mu_i^\circ + RT \ln a_i ,
```

and on nothing else about the solution. The standard potential ``\mu_i^\circ``
comes from the thermodynamic database; the activity ``a_i`` comes from the
**activity model**, which is therefore the single place where a real solution
stops being ideal. This page is about what those models are, where they come
from, what they need, and where they stop being true.

A model has to answer two questions, not one:

1. the **solute** activity ``a_i = \gamma_i m_i``, which sets every solubility
   and every saturation index;
2. the **solvent** activity ``a_w``, which is what a hydrate reaction consumes,
   and which in a cement paste short of mixing water is the quantity that
   decides how far the reaction can go.

The five built-in models — [`DiluteSolutionModel`](@ref),
[`DaviesActivityModel`](@ref), [`HKFActivityModel`](@ref),
[`SITActivityModel`](@ref) and [`PitzerActivityModel`](@ref) — differ in both;
§5 measures by how much for the first three, and §6a and §6 treat the other two.

!!! note "Two meanings of HKF"
    The initials of Helgeson, Kirkham and Flowers [Helgeson1981](@cite) name two
    different things in this package. [`HKFActivityModel`](@ref) is an
    **activity model**: the extended Debye-Hückel law derived in §1, the term
    ``\dot{B} I`` of §2, the Setschenow form for neutral species and a water
    activity consistent with them; it gives ``\gamma_i``. The **HKF equation of
    state** ([Aqueous solutes](@ref sec-theory-hkf)) gives the standard
    properties of a solute, hence ``\mu_i^\circ``, and knows nothing of the
    other ions. Debye-Hückel is the first term of the former, not a part of the
    latter.

## 1. Where the ``\sqrt{I}`` comes from

An ion in a solution of ions is not in the same state as one alone at the same
molality, and the reason is electrostatic. Around a cation, anions are slightly
more likely to be found than cations — not because of any structure, only
because their energy there is lower — so the ion sits at the center of a diffuse
**ionic cloud** of opposite net charge. That cloud screens the ion's own field
and lowers its energy, which is why activity coefficients of ions are *below* 1.

The classical treatment [Debye1923](@cite) takes four steps, and the
assumptions it makes on the way are exactly the assumptions that later fail. An ion ``i`` of charge
``z_i e`` sits at the origin; ``\varepsilon = \varepsilon_0\varepsilon_r`` is
the permittivity of the solvent, ``\rho_w`` the density of water, and the
molality ``m_j`` of the ion ``j`` far from the origin is a number density
``n_j^0 = N_A\rho_w m_j``.

1. **Poisson's equation** relates the mean potential ``\psi(r)`` around the
   ion to the mean charge density ``\rho_e(r)``, the solvent entering only as a
   **continuum of bulk permittivity**: ``\Delta\psi = -\rho_e/\varepsilon``.
2. **Boltzmann statistics** give the mean densities,
   ``n_j(r) = n_j^0\exp(-z_j e\psi/k_BT)``, treating the ions as *independent*
   in the mean field of all the others.
3. **Linearization**, ``z_j e \psi \ll k_B T``. With electroneutrality,
   ``\sum_j z_j n_j^0 = 0``, the charge density
   ``\rho_e = e\sum_j z_j n_j`` reduces to ``-\varepsilon\kappa^2\psi`` with
   ```math
   \kappa^2 = \frac{e^2}{\varepsilon k_B T}\sum_j z_j^2 n_j^0
            = \frac{2 e^2 N_A \rho_w}{\varepsilon_0\varepsilon_r k_B T}\, I ,
   \qquad I = \tfrac12\sum_j m_j z_j^2 ,
   ```
   and Poisson's equation becomes ``\Delta\psi = \kappa^2\psi``. The length
   ``\kappa^{-1}`` is the **Debye length**; ``\kappa \propto \sqrt{I}``.
4. **Exclusion**: no other ion comes closer to the center than a distance
   ``\mathring{a}``, so that for ``r < \mathring{a}`` the only charge is the
   central one and ``\Delta\psi = 0``.

### The potential around the ion

The radial solutions, the inner one with the central charge and the outer one
bounded at infinity, are

```math
\psi(r) = \frac{z_i e}{4\pi\varepsilon\, r} + C \quad (r < \mathring{a}),
\qquad
\psi(r) = D\,\frac{\mathrm{e}^{-\kappa r}}{r} \quad (r > \mathring{a}) .
```

The potential and the radial field are continuous at ``r = \mathring{a}``,
where there is no charge. The field gives
``D\,\mathrm{e}^{-\kappa\mathring{a}}(1+\kappa\mathring{a}) = z_i e/(4\pi\varepsilon)``,
then the potential gives ``C``:

```math
D = \frac{z_i e}{4\pi\varepsilon}\,\frac{\mathrm{e}^{\kappa\mathring{a}}}{1+\kappa\mathring{a}} ,
\qquad
C = \frac{z_i e}{4\pi\varepsilon\,\mathring{a}}\Bigl(\frac{1}{1+\kappa\mathring{a}} - 1\Bigr)
  = -\frac{z_i e}{4\pi\varepsilon}\,\frac{\kappa}{1+\kappa\mathring{a}} .
```

``C`` is the potential that the cloud creates at the center of the ion. It is
negative for a cation, positive for an anion: the cloud screens.

### Charging the ion

The excess chemical potential of the ion is the reversible work needed to
bring its charge from ``0`` to ``z_i e`` with the cloud present, less the work
needed to charge it alone, which belongs to the standard state. Only the
central ion is charged, the others keeping their charges: this is the
assumption of the charging process, under which ``\kappa``, set by the other
ions, stays fixed while ``C`` is proportional to the charge ``q`` being built
up. Per ion,

```math
k_B T\ln\gamma_i = \int_0^{z_i e} C(q)\,\mathrm{d}q
  = -\frac{\kappa}{4\pi\varepsilon(1+\kappa\mathring{a})}\int_0^{z_i e} q\,\mathrm{d}q
  = -\frac{z_i^2 e^2}{8\pi\varepsilon}\,\frac{\kappa}{1+\kappa\mathring{a}} .
```

Writing ``\kappa = B\sqrt{I}`` and passing to decimal logarithms gives the
**extended law**, which [`HKFActivityModel`](@ref) implements,

```math
\log_{10}\gamma_i = -\frac{A z_i^2\sqrt{I}}{1 + B\mathring{a}_i\sqrt{I}} ,
\qquad
B = \sqrt{\frac{2 e^2 N_A\rho_w}{\varepsilon_0 k_B}}\,\frac{1}{\sqrt{\varepsilon_r T}} ,
\qquad
A = \frac{e^2}{8\pi\varepsilon_0 k_B\ln 10}\,\frac{B}{\varepsilon_r T} ,
```

and, as ``\mathring{a} \to 0``, the **limiting law**

```math
\log_{10}\gamma_i = -A z_i^2 \sqrt{I} .
```

As ``I \to 0`` the cloud spreads over a distance ``\kappa^{-1}`` much larger
than the size of the ion, where the potential is small, so that the
linearization of step 3 and the neglect of the size become accurate together:
the limiting law is the dilute end of the theory. The ``\sqrt{I}`` is not
fitted: it is ``\kappa``.

The three relations can be checked symbolically: the two potentials solve their
equations, they match at ``r = \mathring{a}``, and the charging integral gives
the law above.

```@example debye_huckel
using Symbolics
@variables r a κ q ε            # a stands for å, ε for ε₀εᵣ
ψ_in = q / (4π * ε * r) - q * κ / (4π * ε * (1 + κ * a))
ψ_out = q / (4π * ε) * exp(κ * a) / (1 + κ * a) * exp(-κ * r) / r
∂r = Differential(r)
laplacian(ψ) = expand_derivatives(∂r(r^2 * expand_derivatives(∂r(ψ)))) / r^2   # radial
at_a(f) = substitute(f, Dict(r => a))
point = Dict(r => 1.7, a => 0.4, κ => 2.3, q => 1.1, ε => 0.9)   # any values
residuals = [
    laplacian(ψ_out) - κ^2 * ψ_out,                                # Δψ = κ²ψ outside
    laplacian(ψ_in),                                               # Δψ = 0 inside
    at_a(ψ_in - ψ_out),                                            # continuous potential
    at_a(expand_derivatives(∂r(ψ_in - ψ_out))),                    # continuous field
]
value_at(f) = Symbolics.value(substitute(f, point; fold = Val(true)))
maximum(abs ∘ value_at, residuals)
```

The work up to a charge ``Q`` is the primitive of ``C`` that vanishes at
``Q = 0``:

```@example debye_huckel
C(q) = -q * κ / (4π * ε * (1 + κ * a))                 # the cloud at the center
W(Q) = -Q^2 * κ / (8π * ε * (1 + κ * a))               # the work up to the charge Q
@variables Q
simplify(expand_derivatives(Differential(Q)(W(Q))) - C(Q)), W(0)
```

### ``A`` and ``B`` are properties of water, not fitting constants

``A`` and ``B`` hold four fundamental constants and two properties of water,
``\rho_w`` and ``\varepsilon_r``, at the temperature and the pressure of the
solution. With ``\rho_w`` in g/cm³ and ``B`` in Å⁻¹(kg/mol)^½, the constants
combine into two numbers, computed by the package from the CODATA values:

```@example debye_huckel
using ChemistryLab, DynamicQuantities
κ₁ = sqrt(2 * ELEMENTARY_CHARGE_Q^2 * AVOGADRO_Q * 1u"g/cm^3" * 1u"mol/kg" /
          (VACUUM_PERMITTIVITY_Q * BOLTZMANN_Q * 1u"K"))   # κ for unit ρ_w, I and ε_r T
B₁ = ustrip(us"1/Å", κ₁)
A₁ = ustrip(ELEMENTARY_CHARGE_Q^2 / (8π * VACUUM_PERMITTIVITY_Q * BOLTZMANN_Q * 1u"K") * κ₁) / log(10)
B₁, A₁
```

[`hkf_debye_huckel_params`](@ref) evaluates ``A`` and ``B`` from these two
numbers and the density and permittivity given by this package's own equation of
state for water. So the ``A = 0.5114`` and ``B = 0.3288`` that the
models carry as defaults at 25 °C are **derived** as well as tabulated: they
agree to 0.01 % with the values of the LLNL aqueous model that
[ParkhurstAppelo2013](@citet) tabulate, from 0 to 100 °C — compared at three
temperatures in [What the choice of activity model costs](@ref sec-app-activity-models).
[Helgeson1981](@citet), Table 1, computed from the water properties of the time,
gives 0.5091 and 0.3283 at 25 °C.

Both rise with temperature, because water's dielectric constant falls faster than
``T`` rises: hot water screens worse, so the same ionic strength costs more.

### The screening length is the size of a gel pore

The Debye length is ``\kappa^{-1} = 1/(B\sqrt{I})`` in ångström when ``B`` is in
Å⁻¹(kg/mol)^½. At 25 °C it runs from
9.6 nm at ``I = 10^{-3}`` mol/kg to 0.30 nm at 1 mol/kg, passing **0.55 nm at
``I = 0.3`` mol/kg** — the values are tabulated in
[What the choice of activity model costs](@ref sec-app-activity-models).

A cement pore solution sits around ``I \approx 0.1``–``0.5`` mol/kg, so its
screening length is a **few ångström** — the thickness of two or three water
molecules. That is the same scale as the water films in the gel pores of C-S-H,
and both assumptions of step 1 and step 2 above are strained there: a continuum
of bulk permittivity, and ions independent in a mean field. Nothing in the
formulas announces it. It is the physical reason to treat an extended
Debye-Hückel model as a correlation valid in bulk solution rather than as a
theory of confined water.

## 2. The ``\dot{B} I`` term is a deviation function, not a physical term

Measured activity coefficients turn back *upwards* at high ionic strength, which
no screening argument produces. The B-dot model adds a linear term for it,

```math
\log_{10}\gamma_i = -\frac{A z_i^2\sqrt{I}}{1 + B\mathring{a}_i\sqrt{I}}
                    + \dot{B} I .
```

Its status is that of an empirical correlation. [AndersonCrerar1993](@citet) (§17.7.1,
pp. 445–446) record that [Helgeson1969](@citet) defined ``\dot{B}`` as a **deviation
function**: the difference between the *observed* activity coefficient of an
electrolyte — NaCl — and what the extended Debye-Hückel expression predicts for
it. So it carries short-range ion-solvent and ion-ion interaction *and* whatever
the first two terms failed to capture, together, in one number fitted to one
salt. They add that [Helgeson1981](@cite) later split it into a hydration term
from the Born equation and a residual short-range term.

That is why this model has a *ceiling* rather than an asymptote, and why the
ceiling is quoted vaguely as "about a molal": the term is not wrong so much as
it is standing in for physics it does not contain.

Away from 25 °C the model keeps ``\dot{B}`` at its 25 °C value. With
`temperature_dependent = true`, ``A`` and ``B`` are recomputed at the
temperature and pressure of the solution from the water properties of §1,
which the equation of state of water carries to 1000 °C; ``\dot{B}``, fitted to
NaCl at 25 °C, has no temperature function here, and a calculation far from
25 °C in a solution concentrated enough for the term to matter should supply
its own value of it, a keyword away.

The Davies equation [Davies1962](@cite) makes the same move with one more
simplification. It drops the ion size, replacing ``1 + B\mathring{a}_i\sqrt{I}``
by ``1 + \sqrt{I}``, and fixes the linear term once for every ion:

```math
\log_{10}\gamma_i = -A z_i^2\left(\frac{\sqrt{I}}{1 + \sqrt{I}} - b\,I\right) ,
\qquad b = 0.3 .
```

No per-species datum is needed, and every ion of a given charge has the same
activity coefficient. The ``bI`` term plays the part of ``\dot{B} I``, with a
weight proportional to ``z_i^2``, and its status is the same: a correction
fitted to data, not derived from screening.

Neutral species get the **Setschenow** form, ``\log_{10}\gamma_i = K_n I``:
water engaged in the solvation shells of ions is water unavailable to solvate a
neutral molecule, so its activity rises with ionic strength and its solubility
falls. This is salting out, and `CO₂(aq)` is the case that matters for
carbonation.

## 3. The water activity, and why it cannot be assumed separately

The solvent is not a solute and its activity is not obtained by the same
formula. At constant ``T`` and ``P``,

```math
\sum_i n_i \,\mathrm{d}\mu_i = 0 ,
```

which is not an optional refinement: it is the statement that the solvent and
the solutes are parts of one thermodynamic system. Once the activities of the
solutes are chosen, it leaves the solvent none to choose. Integrating it over
the solutes' terms gives the osmotic coefficient ``\varphi``,

```math
\ln a_w = -M_w \varphi \sum_j m_j ,
```

and that is how the solvent row of the ideal model and of the two extended
Debye-Hückel models is built. The ideal model has ``\varphi = 1``, the exact
partner of ``\ln m_i``. Davies, whose ions all carry one function of ``I``
times ``z_i^2``, has a closed form for ``\varphi`` (see
[`DaviesActivityModel`](@ref)). [`HKFActivityModel`](@ref) integrates its
``A``, ``B`` and ``\dot{B}``. Raoult's mole fraction ``x_w`` agrees with the
ideal ``-M_w \sum_j m_j`` to first order, the two differing by
``(M_w \sum_j m_j)^2/2``, but breaks the relation by ``1 - x_w``; a model that
corrected its solutes and left its solvent at ``x_w`` would set the two halves of
one model against each other. [`TruesdellJonesActivityModel`](@ref) and
[`SITActivityModel`](@ref) (§6a) keep it.

The one approximation in the B-dot method is that ``\varphi`` uses a single
charge-weighted mean radius,
``\mathring{a}_{\text{eff}} = \sum_i m_i z_i^2\mathring{a}_i / \sum_i m_i z_i^2``,
where the ``\gamma_i`` use per-ion radii. §5 measures that too.

## 4. Inputs and outputs, model by model

What each model **needs** and what it **returns** — the full parameter tables,
with the provenance of every default, are in the docstrings
([`HKFActivityModel`](@ref), [`DaviesActivityModel`](@ref)); this is the summary
that lets you choose.

| | [`DiluteSolutionModel`](@ref) | [`DaviesActivityModel`](@ref) | [`HKFActivityModel`](@ref) | [`SITActivityModel`](@ref) | [`PitzerActivityModel`](@ref) |
|:--|:--|:--|:--|:--|:--|
| solute scale | molarity | molality | molality | molality | molality |
| ``\gamma_i`` | ``\equiv 1`` | Davies | extended D-H + ``\dot{B} I`` | D-H + ``\sum_k \varepsilon_{ik} m_k`` | virial expansion |
| ``a_w`` | ``-M_w \sum m`` | osmotic coefficient, closed form | osmotic coefficient | Raoult (§6a) | osmotic coefficient |
| per-species data | none | none | ion radii ``\mathring{a}_i`` (tabulated, overridable) | **one parameter per ion pair — caller input** | **a parameter per ion pair and per triplet — caller input** |
| scalar inputs | none | ``A``, ``b``, ``b_n`` | ``A``, ``B``, ``\dot{B}``, ``K_n``, ``\mathring{a}_{\text{default}}`` | ``A`` and the convention ``b = 1.5`` | the shape constants ``\alpha_1``, ``\alpha_2``, ``b`` |
| ``T``, ``P`` dependence | none | ``A(T,P)`` on request | ``A(T,P)``, ``B(T,P)`` on request | ``A(T,P)`` on request; the ``\varepsilon`` set is fitted at one temperature | ``A_\varphi(T,P)`` only; the ``\beta`` set is fitted at one temperature |
| returns | ``\ln a_i`` for every species | same | same | same | same |
| ``\gamma`` useful to | ``I \lesssim 10^{-4}`` (1 % on a monovalent ion) | ``I \lesssim 0.5`` | ``I \lesssim 1`` | ``I \lesssim 3`` to ``4`` | the range its set was fitted over, a few mol/kg |
| Gibbs-Duhem consistent | **exactly** | **exactly** on ions; not for a neutral with ``b_n \neq 0`` (§5b) | to ``10^{-5}`` (§5) | see §6a | **exactly, by construction** (§6) |

All of them return the same object — a vector of ``\ln a_i`` indexed like
`cs.species`, covering solutes, solvent, pure crystals (``0``), gases and
solid-solution end-members — so they are interchangeable at every call site, and
`concentration_scale` tells the accessors which convention was used.
[`activity_coefficients`](@ref) evaluates the per-ion formula of the first three
directly. SIT and Pitzer couple the solutes and have none; their coefficient is
read back from the closure, which writes ``\ln a_i = \ln\gamma_i + \ln m_i`` with
``m_i`` built on its own floored amounts, so that
``\ln\gamma_i = \ln a_i - \ln m_i`` is exact for a trace solute as for a major one.

One limit is common to all of them and is not one of accuracy. A cation cannot
be added to a solution without an anion, so a measurement gives the mean
coefficient of a neutral combination,
``\gamma_\pm = (\gamma_+^{\nu_+}\gamma_-^{\nu_-})^{1/\nu}``, and the ratios that
appear in balanced reactions, never the coefficient of one ion
[Richet2001; Secs. 11.6b–c, pp. 273–275](@cite). How the mean is shared among
the ions is the model's, and so is any quantity that rests on one ion alone,
the pH, ``-\log_{10} a_{\mathrm{H^+}}``, among them
[Richet2001; Sec. 12.1b, p. 287](@cite).

## 5. What the difference is worth, measured

The models are compared, on an imposed NaCl composition and without solving
anything, in
[What the choice of activity model costs](@ref sec-app-activity-models). Three
results from that page belong here, because they are about the theory rather
than about the numbers:

**The ideal model is already several percent off at a millimolal.** Ideality is
not a safe default that degrades gracefully; it is exact only in a limit.

**Davies and the B-dot model part company above a tenth molal**, and by
3 mol/kg Davies returns ``\gamma > 1`` while the B-dot model is still below 1.
That is the ``bI`` term of §2 taking over from the screening term — the ceiling
of a deviation function arriving, visible in the numbers.

**The values of ``a_w`` barely separate below a molal**, and in a single salt
the Gibbs-Duhem relation holds for every model to rounding: each builds its
solvent row from its own solutes' terms, and two ions of one charge and one size
leave nothing to break it. What breaks it shows in a mixed solution, through the
symmetry measured in §5b: ions of different sizes or charges under a linear
term, and neutral species carrying a salting-out coefficient the ions'
coefficients do not return. The cost of such a break is not a large error in
``\mu_w``: an equilibrium is set by the values of the potentials, and two water
activities within one percent of each other shift ``\mu_w`` by less than
``0.01\,RT``. It is that the solvent and the solutes no longer draw their
potentials from one Gibbs energy, so that there is no function whose minimum
the equilibrium is; §5b tests this exactly, and the certificate depends on it.

The B-dot model's own approximation — a single charge-weighted mean ion size in
its osmotic coefficient — costs nothing in NaCl, whose two ions carry the ion
size of the salt, and appears only in a solution whose ions differ in size.

## [5b. Does a model come from one excess Gibbs energy at all?](@id sec-theory-potential)

§5 measures a Gibbs-Duhem residual **along one path**. A small residual there is
necessary and not sufficient: a model can cancel along a dissolution and fail
badly in another direction. The question underneath is whether the activities are
the gradient of *any* function of ``\mathbf{n}``, and that has an exact test.

If ``\mu_i = \partial G/\partial n_i`` for some ``G``, then second derivatives
commute, so

```math
\frac{\partial \ln a_i}{\partial n_j} = \frac{\partial \ln a_j}{\partial n_i}
\qquad\text{(Maxwell)},
```

and, ``\mu_i`` being homogeneous of degree zero,

```math
\sum_i n_i \frac{\partial \ln a_i}{\partial n_j} = 0
\qquad\text{(Gibbs-Duhem)} .
```

The first is the sharper: it fails exactly when **no** ``G`` has these activities
for its gradient, whatever ``G`` might be. Both are read off the Jacobian of the
log-activities, the closure [`activity_model`](@ref) returns, by automatic
differentiation — no ``G`` is constructed, no solve
runs, nothing is differenced against a second solve.

### What is measured

On `H2O@ / Na⁺ / Cl⁻ / Ca²⁺ / SO₄²⁻` at 0.1, 0.12, 0.02 and 0.01 mol, worst
relative asymmetry over pairs, ion/ion and solvent/ion separated:

| model | ion/ion | solvent/ion |
|:--|--:|--:|
| Pitzer (on NaCl) | ``1.3\times10^{-15}`` | ``1.3\times10^{-15}`` |
| ideal dilute | ``< 10^{-12}`` | ``< 10^{-12}`` |
| Davies (no neutral species) | ``< 10^{-12}`` | ``< 10^{-10}`` |
| B-dot, ``\mathring{a}`` per ion, ``\dot{B}\neq0`` (default) | ``1.8\times10^{-1}`` | ``2.9\times10^{-2}`` |
| B-dot, common ``\mathring{a}``, ``\dot{B}\neq0`` | ``1.2\times10^{-1}`` | ``2.6\times10^{-2}`` |
| B-dot, ``\mathring{a}`` per ion, ``\dot{B}=0`` | ``2.4\times10^{-1}`` | ``6.0\times10^{-2}`` |
| **B-dot, common ``\mathring{a}`` and ``\dot{B}=0``** | ``0`` | ``1.8\times10^{-14}`` |
| **Debye-Hückel limiting law** (``\mathring{a}=0``, ``\dot{B}=0``) | ``0`` | ``3.4\times10^{-13}`` |
| ``\mathring{a}=0``, ``\dot{B}=0.0976`` (one GEM-Selektor run [Kulik2013](@cite), not the prescription of Cemdata18 [Lothenbach2019](@cite)) | ``1.2\times10^{-1}`` | ``3.6\times10^{-1}`` |

Two readings, and the second is the one that was not expected.

**The limiting law is exact.** It is the classical result, and it is the anchor:
that model *does* come from one excess energy, to machine precision.

**It takes both corrections to break it, and removing either alone does not
mend it.** Dropping ``\dot{B}`` while keeping per-ion radii makes the ion/ion
asymmetry *worse*.

### Why, in two lines

Every ``\gamma_i`` of §1 depends on composition only through ``I``, and, at a
fixed amount of water, ``\partial I/\partial n_j = z_j^2/(2\,n_w M_w)``. So

```math
\frac{\partial \ln\gamma_i}{\partial n_j}
   = \ln 10 \; f'(I;z_i,\mathring{a}_i)\; \frac{z_j^2}{2\,n_w M_w} ,
```

and Maxwell demands that ``f'(I;z_i,\mathring{a}_i)/z_i^2`` not depend on ``i``.

- With a **common** ``\mathring{a}`` and no ``\dot{B}``, ``f'`` is proportional
  to ``z_i^2`` exactly, and the ratio is the same for every ion. Hence the zeros
  in the table.
- A **per-ion** ``\mathring{a}_i`` puts ``i`` inside the denominator
  ``1+B\mathring{a}_i\sqrt{I}``, where it does not factor out.
- The **``\dot{B} I`` term is added with the same coefficient to every ion**, so
  it contributes ``\dot{B}`` and not ``\dot{B}z_i^2``. An ion of charge 2 and an
  ion of charge 1 then disagree by construction.

Neither is a defect of this implementation: both are the published extended form,
which every geochemical code uses. §2 already says what ``\dot{B}`` is — a
deviation function fitted to one salt — and this is the same fact seen from the
thermodynamics rather than from the fit. The model is an approximation with a
stated domain; outside a fairly dilute solution it does not have a potential
behind it, and saying so is not a criticism of the model or of the codes built
on it.

### What it means for a certificate

[`equilibrate_certified`](@ref) proves a KKT point of the problem it is handed.
When the activities are not the gradient of a single ``G``, that problem is not
the minimization of a Gibbs energy, and `optimal = true` says less than its name
suggests: a stationary point of the residuals, not a minimum of a potential.

In practice this is the regime every geochemical calculation on extended
Debye-Hückel has always worked in, and the table says how far from a potential
one is standing. Where the water activity carries the answer — a sealed paste,
a concentrated brine — [`PitzerActivityModel`](@ref) is the model that has one,
and §6 is about what that costs.

`test/activities.jl` pins the zeros and the orders of magnitude above, so a change
to any ``\gamma`` shows up here rather than in a result.

The two lines of the derivation are checked as such rather than only through
their consequence. [`log10_gamma_expression`](@ref) returns the kernel as a
symbolic expression, so ``f'(I;z_i,\mathring{a}_i)/z_i^2`` can be differentiated
and compared between ions: with a common ``\mathring{a}`` and no ``\dot{B}`` the
two agree to ``10^{-14}``, and with either correction in place they do not. The
same expression is asserted to agree with the compiled kernel, which is what
makes the formula written in §1 the formula that runs.

## [6a. SIT: one parameter per pair, and no radius](@id sec-theory-sit)

Between the corrected Debye-Hückel laws above and the virial expansion below
sits the **Specific ion Interaction Theory** (SIT). It keeps the screening term and replaces the empirical deviation term
with a sum over *pairs*:

```math
\log_{10}\gamma_i \;=\; -z_i^2\,\frac{A\sqrt{I}}{1 + 1.5\sqrt{I}}
   \;+\; \sum_k \varepsilon(i,k)\, m_k
```

the sum running over ions ``k`` of charge **opposite** to ``i``.

Read that against §2. Davies and the ``\dot{B} I`` form both write the deviation
as *one number times the ionic strength*, which says the correction depends on
how much salt there is and not on which salt. SIT says it does depend on which:
``\varepsilon(\mathrm{Na}^+,\mathrm{Cl}^-)`` and
``\varepsilon(\mathrm{Na}^+,\mathrm{SO}_4^{2-})`` are ``+0.03`` and ``-0.12``,
of opposite sign. That is the whole content of the model, and it is why a
sodium sulfate solution and a sodium chloride solution of the same ionic
strength do not have the same activity coefficients.

The ``1.5`` is a **convention**, not a fit. Every published ``\varepsilon`` was
obtained with it, so changing it makes a compilation inconsistent with the model
that produced it; [`SITActivityModel`](@ref) exposes it as a field only to make
that visible.

### What it gives, and where it stops

``\gamma`` under SIT does not fall forever: the pairwise term is positive for
many pairs and eventually overcomes the screening term, so the curve **turns
around**. In sodium chloride at 25 °C, ``\log_{10}\gamma_{\mathrm{Na}^+}``
reaches ``-0.174`` near ``I = 1`` and has come back to ``-0.156`` at ``I = 3``.
A one-parameter deviation term cannot do that and match the dilute end too,
which is why the range is the test rather than a single point.

The model is a **truncation**: the pairwise term is the first order of the same
virial expansion §6 carries to second order, so at high ionic strength the terms
it drops stop being small. The range usually quoted, ``I \lesssim 3``–``4``
mol/kg (the table of §4), is a rule of thumb from practice, not a bound derived
from the model. It also expects genuine ion pairs to appear in
the *speciation* and not inside ``\varepsilon``.

### Why it matters here, beyond being one more option

The thermodynamic databases of the NEA reviews and of ANDRA's ThermoChimie are
**calibrated in SIT**. A ``\log K`` taken from those and used under Davies is
not the constant that was fitted, and the difference is not a rounding: it is
the same category of error as using a surface constant fitted with a diffuse
layer in a model without one. Reproducing published work from those compilations
needs the model they were written in.

### What it costs

The water activity here is Raoult's mole fraction, as in
[`TruesdellJonesActivityModel`](@ref), not an osmotic coefficient: a departure
from a full SIT treatment, in which the solvent follows from the terms of the
solutes through ``\ln a_w = -M_w\varphi\sum_j m_j``
[Richet2001; Sec. 11.6d, Eq. (11.116), p. 276](@cite), and one that does not
enter a comparison made at prescribed proton activity.

And an unlisted ``\varepsilon`` is taken as **zero**, which is the literature's
convention and never a measurement. `missing_epsilon_pairs` reports which pairs
of a given system are resting on it, because a convention silently applied is
indistinguishable from a coefficient somebody determined.

## 6. The ion-interaction model: a different kind of object

Everything above is a **corrected Debye-Hückel law**: one screening term derived
from electrostatics, one size correction, and one empirical term standing in for
the rest. [`PitzerActivityModel`](@ref) is not that. It starts from the excess
Gibbs energy and expands it as a **virial series in the molalities**, exactly as
the pressure of a non-ideal gas is expanded in its density —
[AndersonCrerar1993](@cite) (§17.8) derive it as a cluster expansion with
osmotic pressure in place of pressure, which is where the analogy is exact:

> the first term in a virial equation always represents ideal behavior; in the
> second term ``B_2`` represents the non-ideal contribution from pairwise
> interactions of molecules; ``B_3`` gives the interactions of triples.

Hence the shape of the parameter set: one coefficient per **ion pair**, one per
**triplet**, and no per-species radius at all.

```math
\frac{G^{\text{ex}}}{RT\, w_w} =
  f(I)
  + \sum_{c}\sum_{a} m_c m_a\, B_{ca}(I)
  + \sum_{c<c'} m_c m_{c'}\,\theta_{cc'}
  + \sum_{a<a'} m_a m_{a'}\,\theta_{aa'}
  + \text{triplets } (\psi)
  + \text{neutrals } (\lambda)
```

where ``w_w = n_w M_w`` is the mass of water in kilograms, the basis of the
molalities.

The long-range term keeps the Debye-Hückel physics — it must, since the limiting
law is exact — in the form

```math
f^{\gamma} = -A_\varphi\left[\frac{\sqrt{I}}{1 + b\sqrt{I}}
             + \frac{2}{b}\ln\left(1 + b\sqrt{I}\right)\right] ,
```

and the pair term carries the ionic-strength dependence that a single constant
cannot:

```math
B_{ca}(I) = \beta^{(0)}_{ca}
          + \beta^{(1)}_{ca}\,g\!\left(\alpha_1\sqrt{I}\right)
          + \beta^{(2)}_{ca}\,g\!\left(\alpha_2\sqrt{I}\right) ,
\qquad
g(x) = \frac{2\left[1 - (1+x)e^{-x}\right]}{x^{2}} .
```

``A_\varphi`` is the same Debye-Hückel slope as §1, on the osmotic basis:
``A_\varphi = A\ln 10/3``, so it comes from the water model rather than from a
second implementation. ``\alpha_1``, ``\alpha_2`` and ``b`` are not fitted —
they fix the functional form, and a published ``\beta`` table is only valid with
the values it was fitted against.

### Why the whole construction matters here

``\gamma_i`` and ``\varphi`` are both partial derivatives of **one** function.
So the Gibbs-Duhem relation between the solutes and the solvent is an identity
of the algebra, not an approximation — and that is the defect §3 measured in the
other two models. It is verified rather than asserted: with the derivative taken
analytically, ``\sum_i n_i\,\mathrm{d}\mu_i`` comes out **exactly zero** for
this model along three composition directions at 0.1, 1 and 3 mol/kg, where the
B-dot model's residual is not small.

### What it costs

The parameters are **caller input**. The ThermoFun databases read by this
package carry none; PHREEQC distributes a set in its `pitzer.dat`, and the
package ships the cement set of [Reardon1990](@citet) in
`data/pitzer-reardon1990.toml`, read by [`build_pitzer_parameters`](@ref). What
the package refuses is to invent a missing one. A
[`PitzerParameters`](@ref) takes every table as a keyword without a default, and
completeness is checked against the species list rather than against the set
alone — every cation-anion pair present must have a ``\beta^{(0)}``, and the
error names those that do not.

Two further terms belong to the model, and both are carried:

  - the **higher-order electrostatic terms** ``{}^E\theta(I)`` and
    ``{}^E\theta'(I)`` of [Pitzer1975](@citet), for two ions of like sign and
    unlike charge, Na⁺ with Ca²⁺ in a cement pore solution. They vanish for a
    symmetrical pair, enter ``\gamma`` and the osmotic coefficient from the same
    excess energy, and use Pitzer's function ``J(x)`` by the Chebyshev
    approximation PHREEQC uses [Plummer1988](@cite). A set whose ``\theta`` were
    fitted without them is used with `etheta = false`;
  - **temperature terms** of the coefficients, in the form PHREEQC's databases
    give them, used when the model is built with `temperature_dependent = true`.
    A published set fitted at one temperature, as the shipped one is, carries
    none, and is used at its own values whatever the temperature.

A set in PHREEQC's format is read by [`build_pitzer_parameters`](@ref) with
`format = :phreeqc`, from the file the caller has.

### The speciation a set assumes is part of the set

This one is easy to miss and changes results. A Pitzer set absorbs ion
association into its ``\beta`` coefficients: [HarvieMollerWeare1984](@citet)
fit the Ca–SO₄ interaction rather than postulating a `CaSO₄⁰` complex. A species list
that carries both the free ions **and** the ion pairs therefore counts the same
association twice, and CEMDATA18 does carry them — `Ca(SO4)@`, `CaOH+`,
`Na(SO4)-`, `NaOH@` among others. It also names the silica species differently
(`HSiO3-`, `SiO2@`) from the set's `H3SiO4-` and `H2SiO4-2`.

So the shipped set of [Reardon1990](@citet) does not drop into a CEMDATA18 cement calculation,
and the completeness check refuses such a system rather than returning a number.
That refusal is the correct outcome and not a limitation of the implementation:
combining a dissociated parameterization with an associated speciation is not a
defensible calculation in any code.

## 7. Outside the domain

None of this survives arbitrary concentration, and the failure is silent: the
formulas go on returning finite, plausible numbers. Two guards exist rather than
one warning.

[`solvent_fraction`](@ref) reports the mole fraction of water *within* the
aqueous phase, and [`SOLVENT_FRACTION_FLOOR`](@ref) is where
[`equilibrate_certified`](@ref) refuses to call the answer a solution at all.
The regime it catches is real: a cement paste below ``w/c \approx 0.30``, on the
species lists used here, drives the free water to the solver's floor and the
ionic strength to hundreds of mol/kg — a composition at which every molality is
per kilogram of a solvent that is no longer there. See
[the w/c example](@ref sec-wc-ratio) and
[the water budget](@ref sec-theory-water-budget).

A Pitzer model does **not** rescue that regime. Its fitted range is a few
mol/kg, not hundreds, and no activity model repairs a composition that has left
the physical picture of a solution.

## Where to go next

The chapter continues with [Solid solutions](@ref sec-theory-solid-solutions),
the other place where a mixture departs from ideality. The expressions of this
page are evaluated side by side in
[What the choice of activity model costs](@ref sec-app-activity-models), the
ion-interaction model is loaded and checked in
[The Pitzer model, and what it can be used with](@ref sec-app-pitzer), and the
syntax for choosing a model is in [Activity models](@ref sec-activity-models).
