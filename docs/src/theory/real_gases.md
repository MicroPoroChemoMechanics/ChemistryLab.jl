# [Real gases and pressure](@id sec-theory-real-gases)

!!! info "Before this page"
    [Standard states](@ref sec-theory-standard-states), §2, for the reference
    state of a gas, the pure ideal gas at ``P^\circ``.

A gas keeps the pure ideal gas at ``P^\circ = 1`` bar as its standard state
whatever its pressure; what departs from the ideal gas is its activity, through
a fugacity coefficient given by an equation of state. Pressure also moves the
standard state of the condensed species, solids, solutes and solvent, through
their volumes. This page states both, the equation of [PengRobinson1976](@citet)
for the gas and the equation of [Haar1984](@citet) for the solvent, and what the
equilibrium solver sees of them.

## 1. The fugacity coefficient

The chemical potential of a gas in a mixture at ``T``, ``P`` and mole fractions
``\mathbf{y}`` is written

```math
\mu_i = \mu_i^\circ(T) + RT \ln\frac{f_i}{P^\circ}, \qquad f_i = \varphi_i\,y_i\,P ,
```

the fugacity ``f_i`` replacing the partial pressure ``y_i P`` of the ideal gas, so
that the activity is ``a_i = \varphi_i y_i P/P^\circ``. The fugacity coefficient
measures what the real mixture adds to the ideal one at the same ``T``, ``P``
and composition. With ``G^{\mathrm{res}}`` the residual Gibbs energy of the
phase, its Gibbs energy minus that of the ideal gas mixture,

```math
\ln\varphi_i = \frac{1}{RT}\left(\frac{\partial G^{\mathrm{res}}}{\partial n_i}\right)_{T, P, n_{j \ne i}} .
```

Three consequences follow from this one definition and are what the tests of
the package check, by automatic differentiation, for any constants:

  - **Gibbs–Duhem**: ``\sum_i n_i\,\partial\ln\varphi_i/\partial n_j = 0``, the
    fugacity coefficients being the partial molar quantities of one function;
  - **the partial molar volume**: ``RT\,\partial\ln\varphi_i/\partial P = \bar V_i - RT/P``,
    since ``\partial G/\partial P = V``;
  - **the ideal limit**: ``\varphi_i \to 1`` as ``P \to 0``, as
    ``\ln\varphi = B_2 P/RT + O(P^2)`` for a pure gas, ``B_2`` its second virial
    coefficient.

## 2. The equation of Peng and Robinson

The equation of state relates pressure, temperature and molar volume through two
parameters,

```math
P = \frac{RT}{V - b} - \frac{a\,\alpha(T)}{V(V + b) + b(V - b)} ,
```

the attraction ``a`` and the covolume ``b``, which the critical temperature
``T_c`` and pressure ``P_c`` set, and a temperature function ``\alpha`` which the
acentric factor ``\omega`` shapes:

```math
a = \Omega_a\frac{R^2 T_c^2}{P_c}, \qquad b = \Omega_b\frac{R T_c}{P_c}, \qquad
\alpha(T) = \left[1 + \kappa\left(1 - \sqrt{T/T_c}\right)\right]^2, \qquad
\kappa = 0.37464 + 1.54226\,\omega - 0.26992\,\omega^2 .
```

In the compressibility factor ``Z = PV/RT`` and the reduced parameters
``A = a\alpha P/(RT)^2`` and ``B = bP/RT``, it is a cubic:

```math
Z^3 - (1 - B)\,Z^2 + (A - 3B^2 - 2B)\,Z - (AB - B^2 - B^3) = 0 .
```

**The two constants.** At the critical point the vapor and the liquid merge,
and the cubic has a triple root ``Z_c``: ``(Z - Z_c)^3`` matched to it
coefficient by coefficient gives ``Z_c = (1 - B)/3``,
``A = 3Z_c^2 + 3B^2 + 2B``, and for ``B`` the cubic
``64B^3 + 6B^2 + 12B - 1 = 0``. Its root is ``\Omega_b = 0.0777960739``, and
then ``\Omega_a = 0.4572355289`` and ``Z_c = 0.3074013087``. The paper prints
``0.45724`` and ``0.07780``, the same numbers rounded; the package takes the
roots themselves. Rounded, they move the triple root by the cube root of the
rounding, about ``0.01`` in ``Z``, so the critical point of the equation would
no longer be the one it is given.

**The residual energy.** Integrated along the isotherm from the ideal gas, the
equation gives, per mole of gas,

```math
\frac{G^{\mathrm{res}}}{NRT} = Z - 1 - \ln(Z - B)
    - \frac{A}{2\sqrt 2\,B}\ln\frac{Z + (1 + \sqrt 2)B}{Z + (1 - \sqrt 2)B} ,
```

which is ``\ln\varphi`` for a pure gas.

## 3. Which root

Below the critical temperature and between the spinodals the cubic has three
real roots above ``B``: the vapor, the liquid and, between them, a root that
describes no stable state. The stable phase is the root of least Gibbs energy,
and that is the one taken. Where the two energies are equal, the fugacities of
the vapor and the liquid are, which is the saturation: for carbon dioxide at
25 °C, with the constants of `phreeqc.dat` [ParkhurstAppelo2013](@cite), the
equation places it at 64.4 bar, against 64.35 bar for the reference equation of
pure carbon dioxide that [Spycher2003](@citet) quote. Across it, the root jumps
from the vapor to the liquid and the fugacity stays continuous.

The root is found on plain numbers, in closed form, then lifted into the dual
numbers of whatever is being differentiated by two Newton steps on the cubic
``F(Z) = 0``. The value does not move, and the derivatives come out as the
implicit-function theorem gives them, ``\partial Z/\partial\theta =
-(\partial F/\partial\theta)/(\partial F/\partial Z)``, exact to second order. At
the critical point itself ``\partial F/\partial Z`` vanishes and the derivatives
are infinite, which is the physics: the compressibility of a fluid diverges
there.

## 4. Mixtures

A mixture is described as one fluid whose parameters are those of the van der
Waals mixing rules,

```math
a\alpha = \sum_i \sum_j y_i y_j (1 - k_{ij})\sqrt{a_i\alpha_i\,a_j\alpha_j}, \qquad
b = \sum_i y_i b_i ,
```

the binary parameters ``k_{ij}`` correcting the geometric mean of the
attractions; they are zero unless given. Differentiating ``N G^{\mathrm{res}}``
with respect to the amount of each gas gives

```math
\ln\varphi_k = \frac{b_k}{b}(Z - 1) - \ln(Z - B)
    - \frac{A}{2\sqrt 2\,B}\left(\frac{2\sum_j y_j (a\alpha)_{jk}}{a\alpha} - \frac{b_k}{b}\right)
      \ln\frac{Z + (1 + \sqrt 2)B}{Z + (1 - \sqrt 2)B} ,
```

with ``(a\alpha)_{jk} = (1 - k_{jk})\sqrt{a_j\alpha_j\,a_k\alpha_k}``. A gas phase
follows the equation as a whole: when some of its gases carry critical
constants and others do not, there is no mixing rule between the equation and
its absence, and the package refuses the system rather than mix them.

## 5. What the equilibrium sees

The activity of a gas becomes

```math
\ln a_i = \ln y_i + \ln\varphi_i(T, P, \mathbf{y}) + \ln\frac{P}{P^\circ},
```

in every activity model, and the gas phase occupies ``V = Z N R T/P`` instead of
``NRT/P``. Two consequences for the solver are worth stating.

  - **A gas alone** is a pure phase whose activity, ``\varphi(T, P)\,P/P^\circ``,
    does not depend on its amount: its Gibbs energy is linear in that amount
    whatever the equation of state, the choice of root having already settled
    whether it is a vapor or a liquid. The problem keeps the convexity it had,
    and a certificate can still prove a global minimum
    ([Proving that an answer is the answer](@ref sec-theory-certificate)).
  - **A mixture of real gases** can split into two fluids of different
    compositions, a liquid rich in one gas and a vapor rich in the other. Its
    energy is not convex in the composition, and a certificate proves no more
    than a local minimum.

A constraint that counts the volume of the system as ``\sum_i n_i V_i^\circ``,
linear in the amounts, has no place for a phase whose molar volume depends on
its composition, and refuses a real gas.

## 6. Pressure in the condensed phases

The standard state of a condensed species is the pure substance, or the
hypothetical one-molal solution, at the pressure of the system, so its standard
Gibbs energy moves with pressure by

```math
\mu_i^\circ(T, P) - \mu_i^\circ(T, P^\circ) = \int_{P^\circ}^{P} V_i^\circ(T, P')\,\mathrm{d}P' .
```

For a solute, the HKF equation of state gives ``V_i^\circ(T, P)``; for a solid of
constant molar volume, the integral is ``V_i^\circ (P - P^\circ)``. The solvent
of the ThermoFun databases follows the equation of state of water of
[Haar1984](@citet), its volume at ``P^\circ`` being the one its record tabulates:

```math
V_w(T, P) = V_w^\circ\,\frac{\rho(T, P^\circ)}{\rho(T, P)}, \qquad
\mu_w^\circ(T, P) - \mu_w^\circ(T, P^\circ) = V_w^\circ\,\rho(T, P^\circ)\,\left[g(T, P) - g(T, P^\circ)\right],
```

with ``\rho`` the density and ``g = f + P/\rho`` the specific Gibbs energy of the
equation, ``f`` its specific Helmholtz energy; the second follows from the first
since ``\partial g/\partial P = 1/\rho``. The enthalpy, the entropy and the heat
capacity gain the derivatives of the same term, so that the four functions remain
those of one Gibbs energy. Against a constant volume, the compressibility of
water, about ``4.5 \times 10^{-10}`` per pascal near 25 °C, takes 1.1 % from the
volume at 250 bar and 10 J/mol from the Gibbs energy at 500 bar. Every term
vanishes at ``P^\circ``, where nothing changes.

## Where to go next

[Carbon dioxide in water under pressure](@ref sec-validation-co2-solubility)
computes the measured solubility of carbon dioxide from 25 to 500 atm with the
gas ideal and real, and separates what the fugacity coefficient and the volume
of the solute each contribute. The syntax is in
[Gases under pressure](@ref sec-manual-real-gases).
