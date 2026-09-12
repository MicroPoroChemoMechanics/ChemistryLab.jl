# [Oxidation state, and the potential conjugate to it](@id theory-redox)

A Portland cement can be modeled without ever mentioning oxidation states. Its
sulfur is sulfate, its iron is ferric, and nothing in the paste changes either.
A cement made with blast-furnace slag cannot: the slag brings sulfur as
**sulfide**, S(-II), into a pore solution whose sulfur is **sulfate**, S(+VI),
and a calculation that cannot hold both has already decided the answer.

This page is about what has to be added to a Gibbs minimization for it to hold
both, and — as importantly — about what that addition does *not* buy.

## Charge is conserved separately from the elements, but not always

The constraint of an equilibrium calculation is `A n = b`: the amounts must
reproduce the component totals. The question is what the components are.

Take a system whose every element has a single valence — sodium always Na(I),
chlorine always Cl(-I). Then the charge of a species is a fixed linear function
of its composition: one sodium contributes `+1`, one chloride `−1`, always. The
charge row of the stoichiometric matrix is a **combination of the element rows**,
it adds no information, and carrying it would impose the same constraint twice.

Now add sulfur at two valences, say `SO4-2` and `HS-`. No fixed multiplier turns
a sulfur atom into a charge any more: it depends on which sulfur. The charge row
is then **independent**, and dropping it would let the calculation move electrons
around for free — turning sulfide into sulfate at no cost, which is precisely
the wrong answer.

So the test is a rank test, and it is the one the package performs:

```math
\text{redox} \iff
\operatorname{rank}\begin{pmatrix} M_{\text{elements}} \\ M_{\text{charge}}\end{pmatrix}
> \operatorname{rank} M_{\text{elements}}
```

When it holds, `Zz` — the unit-charge pseudo-component — is kept as a row of the
conservation matrix, and the **oxidation state of the system is conserved**
alongside its elements. When it does not, the row is dropped as redundant.

!!! note "This happens more often than a reader expects"
    Species lists are usually built with [`speciation`](@ref), which derives
    species from the **atoms** of what it is given. Ask for `SO4-2` and it
    brings `HS-`, `SO3-2` and `S2O3-2` along, because they are made of the same
    atoms. Any system containing sulfate therefore has redox freedom whether or
    not its author wanted any — which is correct, and worth knowing before
    reading a result.

## The potential conjugate to it

Conserving a quantity is half of a thermodynamic description; the other half is
the intensive variable conjugate to it. For the elements those are the chemical
potentials. For charge it is the **electron activity**, written

```math
pe = -\log_{10} a_{e^-}, \qquad
E_h = \frac{RT\ln 10}{F}\, pe
```

with ``F`` the Faraday constant — 0.05916 V per pe unit at 25 °C. `pe` and
``E_h`` are the same statement in different units, as `pH` and the proton
chemical potential are.

### There is no electron to read it from

Here the analogy with `pH` breaks, and the break is structural rather than an
inconvenience of implementation. `pH` is read off a species: `H+` is in the
system, it has an activity, take its logarithm. **No database tabulates a free
electron in solution**, because there is no such thing to tabulate — the
electron in a redox reaction is an accounting device, not a dissolved substance.

So the electron activity is *inferred*, from a couple that is in the system. Any
redox couple obeys a half-reaction, which the package balances from the element
and charge balance alone:

```math
\ce{SO4^2- + 9H+ + 8e- = HS- + 4H2O}
```

Its equilibrium constant relates the activities of the members to the activity of
the electron, so with the members' activities known the electron's follows:

```math
pe = \frac{1}{n}\left(\log_{10}K
     - \sum_{\text{products}}\nu_i\log_{10}a_i
     + \sum_{\text{reactants}\neq e^-}\nu_i\log_{10}a_i\right)
```

The ``\log_{10}K`` is computed from the same standard Gibbs energies as
everything else, with the electron at the conventional standard state
``\Delta_f G^0 = 0`` — the same convention that puts `H+` at zero, and equally a
convention rather than a measurement. That it is the *usual* convention is what
makes the numbers comparable with published half-reaction constants:

| half-reaction | computed here | published |
|:--|--:|--:|
| ``\ce{SO4^2- + 9H+ + 8e- = HS- + 4H2O}`` | 33.69 | 33.66 |
| ``\ce{Fe^3+ + e- = Fe^2+}`` | 13.02 | 13.03 |

Computed from CEMDATA18's own data, so the agreement also checks that CEMDATA18
and the sources of those published constants share a reference state.

## What this does not buy: one potential per system

A single `pe` describes a system only if **every** couple in it is at mutual
equilibrium. That is a strong assumption, and in a cement it is false.

The demonstration costs nothing. On one solution carrying sulfate, sulfide and
both irons at comparable amounts, the two couples report

| couple | `pe` | ``E_h`` |
|:--|--:|--:|
| ``\ce{Fe^3+}/\ce{Fe^2+}`` | +13.0 | +0.77 V |
| ``\ce{SO4^2-}/\ce{HS-}`` | −3.7 | −0.22 V |

almost a volt apart. Neither is wrong. They are the potentials **of those
couples**, and they differ because the solution is not at redox equilibrium with
itself: sulfate reduction is kinetically frozen on any time scale a cement cares
about, so the sulfur couple retains whatever state the slag gave it while the
iron couple relaxes.

This is why [`pe`](@ref) takes the couple as an argument and reports which one it
used. Computing it from two couples and comparing them is a measurement of how
far the single-potential assumption is from holding, and it is worth making
before trusting either number.

## Two ways to pose a redox calculation

**Conserve it.** The default, and what a sealed paste does: the charge row is in
`b`, the oxidation state is whatever the reactants brought, and the potential
comes out as a result. Nothing needs to be said; it happens because the
component is there.

**Prescribe it.** [`FixedpE`](@ref) and [`FixedEh`](@ref) impose the potential
and let the oxidation state follow, for a system genuinely open to an outside
buffer — a controlled atmosphere, an electrode, a measured ``E_h`` one wants to
reproduce. The mechanism is the implicit titrant of [`FixedpH`](@ref), with the
difference forced by the absence of an electron species: what is prescribed is
not one species' activity but the **linear combination** the half-reaction
gives,

```math
\sum_i \nu_i \log_{10} a_i = \log_{10}K - n\,pe
```

one equation, one unknown (the titrant amount), one extra column in the
conservation rows.

Prescribing a potential on a **closed** paste is a modeling error rather than a
modeling choice: it says the system exchanges electrons with something, and a
sealed specimen does not. The first form is the one a cement wants.

## What a slag cement actually needs

Putting the pieces together, the reason this chapter exists:

1. the slag brings S(-II), the clinker brings S(+VI), so the species list spans
   two valences and the charge row survives — the conservation is automatic;
2. the equilibrium then **distributes** sulfur between the two, and how far it
   goes is a thermodynamic result rather than an input;
3. the resulting ``E_h`` is strongly reducing, which is why a slag cement
   protects embedded steel differently from a Portland one, and why its pore
   solution can carry sulfide;
4. but the answer is an *equilibrium* answer, and sulfate reduction is slow. A
   real slag paste is somewhere between the sulfur its slag brought and the
   sulfur equilibrium would give it, and nothing in a Gibbs minimization knows
   where.

Point 4 is the honest limit of this chapter. The calculation is now able to pose
the question; a kinetic description of sulfate reduction, which this package does
not have, is what would answer it.
