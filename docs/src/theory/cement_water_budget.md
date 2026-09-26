# [The water budget of a hydrating paste](@id sec-theory-water-budget)

!!! info "Before this page"
    [Thermochemistry](@ref sec-theory-thermo) §6, where volumes and porosity are
    defined, and [Rate laws, and every parameter in them](@ref
    sec-theory-kinetics).

Mix one gram of cement with ``w/c`` grams of water, seal it, and wait. Some of
the water ends up inside hydrate formulas, some stays in pores the hydrates
leave behind, and the reaction stops well before the cement is gone. Powers'
rule of thumb summarizes the arrest as

```math
\alpha_{\max} = \frac{w/c}{0.42} ,
```

and [`powers_alpha_max`](@ref) supplies it to the rate laws. This page is about
what that 0.42 is made of, which part of it a Gibbs energy minimization can
predict, and why the other part is **not a thermodynamic quantity**.

## 1. Three destinations for the mixing water

At a degree of hydration ``\alpha``, the mixing water has gone to three places:

| destination | symbol | what it is |
|:--|:--|:--|
| **chemically bound** | ``b\,\alpha`` | water written into hydrate formulas — C-S-H, portlandite, AFm, AFt. It is no longer water |
| **the volume the reaction lost** | ``s\,\alpha`` | hydrates occupy less than the water and clinker they came from, so voids open: **chemical shrinkage** |
| **still liquid** | the remainder | pore solution, held more or less tightly depending on the pore it sits in |

Both coefficients are **computed**, not assumed. From certified equilibria on a
CEM I at imposed ``\alpha`` — the calculation is executed on the
[self-desiccation page](@ref sec-self-desiccation) — this species list gives

```math
b = 0.3095 \ \mathrm{g/g} , \qquad s = 0.0639 \ \mathrm{cm^3/g} ,
```

each **constant in ``\alpha`` to four digits**, which is what makes a budget
linear in ``\alpha`` legitimate for this system rather than merely convenient.
``s \approx 0.064`` cm³ per gram of cement is the chemical shrinkage that cement
chemistry reports for an ordinary Portland cement, and here it is a consequence
of the standard molar volumes in the database.

## 2. What the minimization predicts: the stoichiometric threshold

A Gibbs minimization under a conservation budget answers one question: what is
the lowest-energy composition reachable *with the atoms available*. Water
appears in it as hydrogen and oxygen, and nothing more. So the threshold it
predicts is the point at which the hydrogen runs out — that is, ``w/c \approx b``.

Measured on the [w/c example](@ref sec-wc-ratio)'s species list, at equilibrium:

| ``w/c`` | clinker left at equilibrium |
|:--|:--|
| 0.15 | 54.6 % |
| 0.28 | 2.0 % |
| 0.30 and above | 0 % |

That is a real prediction and it is a **stoichiometric** one. Below ``w/c \approx
0.30`` clinker survives because there is not enough water to convert it, and the
number 0.30 is nothing but ``b`` plus the little the solution retains.

!!! warning "Below that threshold the calculation also leaves its own domain"
    At ``w/c = 0.28`` the free water goes to the solver's floor, the solvent mole
    fraction falls to 0.21 and the ionic strength reaches hundreds of mol/kg. The
    composition is arithmetically a solution and physically not one — every
    molality is per kilogram of a solvent that is no longer there. This is what
    [`solvent_fraction`](@ref) and [`SOLVENT_FRACTION_FLOOR`](@ref) exist to
    name, and no activity model repairs it.

## 3. What Powers' 0.42 actually contains

Powers splits his coefficient himself: about **0.23 g/g** of non-evaporable water
— what survives D-drying, his operational definition of chemically bound — plus
about **0.19 g/g** of **gel water**, which is water *present in the paste and
unusable by the reaction*.

So the gap between the two thresholds,

```math
0.30 \;\;(\text{stoichiometric, predicted})
\qquad\longleftrightarrow\qquad
0.42 \;\;(\text{Powers, measured}) ,
```

**is exactly the gel water**, and it is the whole of the question. A paste at
``w/c = 0.35`` has enough hydrogen to consume all its clinker — the minimization
says so, correctly — and does not do it.

Written as a budget, the arrest condition is

```math
\frac{w}{c} = \underbrace{b\,\alpha_{\max}}_{\text{bound}}
            + \underbrace{s\,\alpha_{\max}\frac{S^\ast}{1-S^\ast}}_{\text{the water left at arrest}} ,
\qquad\text{so}\qquad
\alpha_{\max} = \frac{w/c}{k} ,
\quad k = b + s\,\frac{S^\ast}{1-S^\ast} ,
```

where ``S^\ast`` is the degree of saturation at which hydration stops. Two of the
three ingredients are thermodynamic and computed; ``S^\ast`` is neither.

## 4. Why the gel water is not a thermodynamic effect

The natural hypothesis is that confined water is thermodynamically unavailable:
water held in a fine pore has a reduced activity, so hydrates that consume it
become less stable, and the arrest should fall out of the minimization.
**The magnitude is wrong by two orders of magnitude**, and the arithmetic is
short enough to do here:

```@example wb
using ChemistryLab, Printf          # `R_GAS` is the package's, not a literal
R, T = R_GAS, 298.15
a_w = 0.80                     # the internal humidity a sealed paste arrests at
ΔG_water = R * T * log(a_w)                    # J per mole of water
n_water_per_alite = 3.3                        # mol H₂O per mol C3S → C-S-H + CH
@printf("RT ln a_w                        = %8.1f J/mol of water\n", ΔG_water)
@printf("× %.1f mol water per mol alite    = %8.2f kJ per mol of alite\n",
        n_water_per_alite, n_water_per_alite * ΔG_water / 1000)
@printf("hydration ΔG, order of magnitude = %8.0f kJ per mol of alite\n", -100.0)
@printf("\nactivity that would null it      = %.2e\n", exp(-100_000 / (n_water_per_alite * R * T)))
```

An activity of that size corresponds, through Kelvin, to a meniscus radius
**smaller than a water molecule**. And the measurement agrees with the
arithmetic: imposing a water activity from saturation down to 0.80 through
[`CapillaryWater`](@ref) leaves the equilibrium assemblage of a CEM I paste
unchanged to six digits, with a certificate on every answer — that is the first
negative control on the [self-desiccation page](@ref sec-self-desiccation).

So the capillary lowering of the water activity is real, representable, and
energetically negligible. What stops a real paste is elsewhere.

### It is a question of scale, which a 0D calculation does not have

The reaction needs water *at the surface of a clinker grain*, and dissolved ions
carried *away* from it. As hydration proceeds:

  - the coarse capillary pores empty first, because they hold water most loosely;
  - what remains sits in gel pores — nanometers wide, inside the C-S-H — as films
    one or two molecules thick;
  - the permeability of the paste falls by orders of magnitude, and the liquid
    path between the remaining water and an unreacted grain, tens of microns
    away, is broken.

The reactants stop meeting. That is a statement about **geometry and transport
across four orders of magnitude of length**, and a 0D model has no length at
all: it assumes one well-mixed solution in contact with every solid, so a water
molecule in a 1 nm gel pore is the same object as one in a 10 μm capillary pore.
No activity model, however good, can distinguish them — the distinction is not
thermodynamic.

There is a numerical hint of that boundary inside the theory itself. The
[Debye screening length](@ref sec-theory-activity) of a cement pore solution at
``I \approx 0.3`` mol/kg is **0.55 nm** — the thickness of two water molecules,
and the width of the pores in question. Where the screening length and the pore
are the same size, the continuum-dielectric and mean-field assumptions that
produce every activity coefficient on that page are both strained. The formulas
keep returning numbers; they have simply left the picture they were derived in.

## 5. So what *can* be predicted, and how

Three routes, and they differ in what has to be assumed.

**Impose ``\alpha``.** Give the minimization a degree of hydration and let it
compute the assemblage — the construction of [LothenbachWinnefeld2006](@cite),
used by the [w/c example](@ref sec-wc-ratio) below its stoichiometric bound and
by the self-desiccation page throughout. Honest, and ``\alpha`` is an input:
either measured, or taken from `powers_alpha_max`. This is also what GEM-Selektor
and Reaktoro offer, since they are 0D equilibrium codes too.

**Close the budget with a measured isotherm.** Take ``b`` and ``s`` from the
minimization and ``S^\ast`` from a published desorption isotherm at an assumed
arrest humidity. Then ``k`` is predicted up to that one empirical threshold — and
inverting it is more informative than quoting it: with Powers' own
``w_n = 0.23``, ``k = 0.42`` corresponds to an arrest at **77.5 % relative
humidity**, the window sealed pastes are independently reported to stop in. The
[self-desiccation page](@ref sec-self-desiccation) does this, and shows why the
*proportionality* ``\alpha_{\max}\propto w/c`` is structural and therefore no
evidence at all.

**Integrate a rate law that reads the humidity.** [`humidity_factor`](@ref)
implements the empirical cut that stops hydration below about 80 % RH, and
[`PoreHumidity`](@ref) computes that humidity from the current saturation, so a
kinetic run arrests on its own. The arrest becomes a result of the integration
rather than a criterion applied afterwards — while the *criterion* stays
empirical. That is the honest ceiling of a 0D framework.

**Beyond 0D.** Making the arrest a genuine prediction needs a model with a length
in it: a microstructural one, where the arrest emerges because the liquid path
percolates no longer, or a transport calculation on the pore network. That is a
different package.

## 6. The same argument for a blended binder

Nothing in sections 3 and 4 mentions alite. The coefficient

```math
k = b + s\,\frac{S^\ast}{1-S^\ast}
```

is built from the water the hydrates bind and the water the gel holds at arrest;
it is a property of the **pore space**, not of the grain that made it. A slag
particle or a fly ash sphere sitting in the same paste is in the same predicament
as an unhydrated clinker core: the water that remains is in gel pores two
molecules wide, and it is not going to reach a grain a thousand times larger.

So the ceiling transposes, and [`powers_alpha_max`](@ref) is applied to every
constituent of a blended binder on the pages for
[CEM II](@ref ex-cem2-blended), [CEM III](@ref ex-cem3-slag),
[CEM IV](@ref ex-cem4-pozzolanic) and [CEM V](@ref ex-cem5-composite). Two
qualifications go with it, and both matter:

**It is the water/*binder* ratio that is used.** Powers measured ``b`` and ``s``
on Portland cement pastes, and a slag binder's C-A-S-H is not his C-S-H: it has a
lower Ca/Si, it binds a different amount of water, and nobody has published the
equivalent coefficients per supplementary material. Using ``w/b`` treats every
constituent as drawing on the same water, which is the conservative reading and
is stated as an approximation rather than presented as a measurement.

**For a slag or an ash the water is usually not what binds.** The ceiling says
how far the reaction *can* go. What decides how far it *has* gone at an ordinary
age is the glass's own dissolution rate, and that is much the lower of the two:
the RILEM TC 238-SCM round robin [Durdzinski2017](@cite) reports, at 28 days and
by SEM image analysis, 38–49 % for two ground granulated slags and about 20 % for
a siliceous fly ash — against a water ceiling of 0.95 at ``w/b = 0.40``. So the
reacted fraction is the **smaller** of the two ceilings, and for the glasses it
is the kinetic one. The same study is worth quoting on how well any of this is
known: the precision of determination is "rather low, at best ± 4-5 %".

!!! danger "This is the assumption an equilibrium code hides best"
    A Gibbs minimization reacts whatever budget it is handed, without comment. A
    page that hands it the whole binder has asked what the paste becomes after
    every grain has dissolved — a question about geological time — and will get a
    confident, certified, completely unphysical answer: on a CEM V at 48 %
    replacement, a pH of 14.4 and an element balance that cannot close because no
    assemblage in the database will hold the alkalis and aluminum released. The
    failure is in the question. [The CEM V page](@ref cem5-dor) shows both.

## 7. Curing is a boundary condition, and there are two of them

Everything above is a **sealed** specimen: it exchanges nothing with its
surroundings, the volume that chemical shrinkage empties becomes gas-filled
porosity, the saturation falls, and the paste desiccates itself. That is one
boundary condition. The other is a specimen kept under water after setting,
which draws in what the shrinkage empties and never desiccates.

Powers gives both, and the difference between them is exactly the chemical
shrinkage:

```math
k_{\text{sealed}} = 0.42 ,
\qquad
k_{\text{saturated}} = 0.36 ,
\qquad
k_{\text{sealed}} - k_{\text{saturated}} \approx 0.06 \;\text{g/g} ,
```

the 0.06 g of water per gram of cement being the volume the reaction loses
because the hydrates are denser than the reagents that made them. Sealed, that
volume has to come out of the paste's own water; immersed, it comes from the
bath. So the same mix reaches full hydration from a lower mixing water content
when it is cured under water — which is why curing is specified, and why a
strength result quoted without its curing regime is incomplete.

```julia
powers_alpha_max(0.32)                        # 0.762 -- sealed
powers_alpha_max(0.32; curing = :saturated)   # 0.889 -- under water
```

`curing = :saturated` moves the **ceiling**, which is what a kinetic run needs.
It does not by itself open the specimen: a coupled run is a closed system, so the
water the bath supplies is not in its balance.

For an equilibrium the specimen *can* be opened, and that is
[`SaturatedCuring`](@ref) — the mirror image of [`CapillaryWater`](@ref). Where
the sealed constraint lets the saturation fall and lowers the water activity by
the Kelvin term, this one holds the specimen's total volume at the fresh paste's
and draws water in to make up what chemical shrinkage empties:

```math
\sum_i \bar V_i\, n_i = V_{\text{ref}} ,
```

one equation, linear in the composition, with the amount of water imbibed as its
unknown. That amount is not a numerical device: it **is** the chemical shrinkage,
the quantity a chemical-shrinkage test measures by watching a specimen drink.

!!! tip "And the shrinkage is computed, which makes 0.06 a prediction rather than a constant"
    ``\Delta V`` is a difference of standard molar volumes — the same data that
    fixed the assemblage — so the package can be asked a question it was never
    fitted to answer: does it reproduce the coefficient Powers measured? On a
    w/c = 0.40 paste, per gram of reacted cement:

    | | |
    |:--|--:|
    | shrinkage from the molar volumes | 0.0606 cm³ |
    | water the cured specimen drew in | 0.0604 g |
    | Powers, as ``k_{\text{sealed}} - k_{\text{saturated}}`` | 0.0600 g |

    The two internal routes agree with each other to 0.3 %, and both land within
    1 % of a number measured on pastes in 1948. They are not obliged to: the
    assemblage is a declared species list and not a real paste's, the molar
    volumes are ideal, and Powers' coefficient is an average over the cements he
    had. Which is what makes the agreement a check on the volume data rather
    than a restatement of it — and what makes an empirical coefficient
    *intelligible* rather than merely used. [The w/c example](@ref sec-wc-ratio)
    runs it.


!!! warning "A cure is not `FixedActivity("H2O@", 1.0)`"
    The tempting way to write "kept under water" is to prescribe unit water
    activity, and it is wrong. A cement pore solution sits near ``a_w = 0.98``
    from its dissolved salts alone, so prescribing 1 would draw water in until
    the solution was dilute enough to reach it — which never happens, and the
    constraint would imbibe without bound. A bath does not fix the activity
    inside the specimen. It fixes the **availability**: the pore space stays
    full. That is a volume statement, and it is why the constraint is written on
    the volume.

## 8. Measuring it: what thermogravimetry gives, and what it needs

Bound water is the quantity a thermogram integrates to, which makes
thermogravimetry the natural second observable beside calorimetry — and the one
[the calibration example](@ref ex-hydration-calibration) asks for by
name, because heat constrains three combinations of six kinetic parameters and
a measurement that sees the *phases* breaks correlations heat cannot.

[`ignition_loss`](@ref) computes the total from the formulas alone:

```math
m_{\mathrm{H_2O}} = M_{\mathrm{H_2O}}\sum_{\text{solids}} \frac{n_i H_i}{2},
\qquad
m_{\mathrm{CO_2}} = M_{\mathrm{CO_2}}\sum_{\text{solids}} n_i C_i
```

It counts **hydrogen**, not formula water, and the difference is not pedantry:
portlandite is `Ca(OH)₂`, has no `H₂O` written in it, and loses one water per
formula unit on ignition. A rule that searched for `H₂O` would report zero for
the second most abundant hydrate in a paste.

The aqueous phase is excluded, which is the distinction this whole page is
about: pore solution is water and is not bound water.

### From the total to a curve

A total is not a thermogram, and the difference is what identifies phases:
C-S-H, AFt and AFm all release below 200 °C and are told apart by the *shape* of
the release, not by its size. That needs, per phase, a temperature window —

```math
f_i(T) = \frac{1}{1 + \exp\!\left(-\dfrac{T - T_{1/2,i}}{w_i}\right)}
```

— and a window is **not** a consequence of a formula. It comes from one of two
places, and [`DecompositionWindow`](@ref) makes a curve say which:

 1. **A publication.** Then the two numbers are `PROV_PUBLISHED` and carry their
    source.
 2. **A measured thermogram.** The windows are *identifiable* from one, which is
    the whole reason [`thermogram`](@ref) is written as a smooth function of
    them: `window_parameters` hands them to an optimizer, and
    [`identifiability`](@ref) says afterwards which of them the curve actually
    determined.

Neither is invented. A window given as a bare number is `PROV_UNSTATED` — the
weakest claim there is — and one written to get a picture on the screen should
say `PROV_PLACEHOLDER` and keep saying it until a measurement replaces it.

**A phase can need more than one window.** Gypsum loses its two waters in two
steps, `CaSO₄·2H₂O → CaSO₄·½H₂O → CaSO₄`, and C-S-H does not leave in one piece
either; a window per stage with `fraction` splitting the release is how that is
written. The fractions are checked to sum to one, because two windows each
accounting for all of a phase would release its mass **twice** and the only
symptom would be a curve integrating to more than `ignition_loss` — a silent
doubling rather than an error.

!!! warning "Overlapping peaks are where this earns its keep"
    Two phases releasing in the same window is the ordinary case in a paste, and
    a fit that reported four numbers there would be reporting two. Running
    [`identifiability`](@ref) on the windows is not a formality: for two peaks
    5 K apart it returns a rank of two out of four, which is the measurement
    saying so ([the executed case](@ref sec-example-tga)).

[`phases_without_windows`](@ref) is the other half of the honesty: a phase with
no window contributes to the starting mass and never leaves, so a curve computed
without noticing integrates to less than `ignition_loss` and says nothing about
it. It reports `phase => product` pairs rather than phases, because **a phase
can need two windows** — a carbonated hydrate carries hydrogen and carbon,
releases water and carbon dioxide, and does so at different temperatures.
Counting coverage per phase would call such a phase done when half of it is, and
a hemicarboaluminate is not an exotic case in a cement.

[`windows_without_phases`](@ref) is the mirror, and it is the one that catches a
typo: a window on a phase that releases nothing contributes nothing and raises
nothing, so its only symptom is a peak that is not there.

## See also

  - [Self-desiccation](@ref sec-self-desiccation) — the budget closed, with the
    numbers computed and the negative controls executed
  - [The w/c ratio](@ref sec-wc-ratio) — the stoichiometric threshold, scanned
  - [Activity models](@ref sec-theory-activity) — why the screening length
    matters here
  - [`powers_alpha_max`](@ref), [`PoreHumidity`](@ref), [`CapillaryWater`](@ref)
