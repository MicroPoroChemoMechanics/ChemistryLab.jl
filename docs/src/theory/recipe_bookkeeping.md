# [What a recipe puts into the equilibrium](@id sec-theory-recipes)

!!! info "Before this page"
    [Proving that an answer is the answer](@ref sec-theory-certificate), where the budget
    ``\mathbf{b}`` and the conservation matrix ``\mathbf{A}`` enter the
    minimization, and the manual page [Recipes](@ref man-recipes), which shows
    the same objects at work.

A cement paste is described by the materials it was mixed from, their masses and
the fraction of each that has reacted, whereas an equilibrium calculation
expects the amounts of the components it must conserve. This page writes the
conversion from the first description to the second as [`budget`](@ref)
performs it, together with what is kept out of the equilibrium, and the
identities that follow and that the tests check.

## 1. Masses and extents

A recipe holds ``m_\mathrm{b}`` grams of binder, made of materials in the mass
fractions ``f_\mathsf{m}``, and the water ``(w/b)\,m_\mathrm{b}``. A material
``\mathsf{m}`` is a list of constituents ``c``, each a mass fraction ``x_c`` of
it, so that the mass of a constituent in the recipe is

```math
m_c = f_\mathsf{m}\, m_\mathrm{b}\, x_c .
```

An addition outside the binder (a salt, an admixture) enters with its own mass
in place of ``f_\mathsf{m} m_\mathrm{b}``. The degree of reaction of a
constituent at the time ``t`` is the product of that of its material and its
own,

```math
\alpha_c(t) = \alpha_\mathsf{m}(t)\,\alpha'_c(t) ,
```

which lets a degree measured on a material as a whole, that of a slag by image
analysis for instance, apply to the constituents that react, while the crystals
of the same material are held at ``\alpha'_c = 0``
([`effective_extent`](@ref)).

## 2. The budget

A constituent is described either by a database species ``s`` of molar mass
``M_s``, or by its oxides alone, the mass fraction ``\omega_{c,k}`` of each
oxide ``k`` of molar mass ``M_k``. The reacted part of the first kind enters the
equilibrium as that species,

```math
n^0_s = \sum_{c \to s} \frac{\alpha_c\, m_c}{M_s} ,
```

the water as ``\mathrm{H_2O}``, and the reacted part of the second kind as the
elements of its oxides, each oxide decomposed over the primaries of the system
by the vector ``\mathbf{d}_k`` of [`primary_decomposition`](@ref). The budget is
then

```math
\mathbf{b} = \mathbf{A}\,\mathbf{n}^0
  + \sum_{c} \alpha_c\, m_c \sum_{k \in K_c} \frac{\omega_{c,k}}{M_k}\,\mathbf{d}_k ,
```

where ``K_c`` is the set of the oxides of ``c`` whose element the system has a
primary for. The decomposition of a basic oxide takes protons from the water,
``\mathrm{CaO} = \mathrm{Ca^{2+}} + \mathrm{H_2O} - 2\,\mathrm{H^+}``, and that
of an acidic one gives them back, so that ``\mathbf{d}_k`` holds negative
entries; only ``\mathbf{b}`` as a whole has to be a budget the species can
reach. The analysis is not renormalized: ``\sum_k \omega_{c,k}`` falls short of
one by what the analysis does not report as an oxide, a loss on ignition for
instance, and scaling it up would put into the paste matter that the analysis
does not attest.

## 3. What is kept aside

The unreacted part of each constituent, ``(1 - \alpha_c)\,m_c``, stays in the
paste as a solid that takes no part in the chemistry. Its volume is that of its
species, ``(1 - \alpha_c)\,m_c V^\circ_s / M_s``, or its mass over the density
of the constituent, and its enthalpy is the enthalpy of formation of its
species. A constituent known by its oxides has neither a molar volume nor an
enthalpy in any database, and both are left unknown unless its source gives
them. The reacted mass of an oxide outside ``K_c``, ``\alpha_c m_c
\omega_{c,k}``, is kept aside in the same way, with the reason
`:not_in_system`: the titanium of a cement in a system without titanium has
reacted, but nothing in the system can hold it.

## 4. The mass balance

Since ``\mathbf{d}_k`` decomposes one formula unit of the oxide without loss,
the equilibrium holds the mass of every species put into it and of every oxide
of ``K_c``. It follows that

```math
m_\mathrm{b}\left(1 + w/b\right) + m_\mathrm{add}
  = m_\mathrm{eq} + m_\mathrm{res}
  + \sum_c \alpha_c\, m_c \Bigl(1 - \sum_k \omega_{c,k}\Bigr) ,
```

with ``m_\mathrm{eq}`` the mass of the equilibrium, ``m_\mathrm{res}`` that of
what is kept aside and ``m_\mathrm{add}`` that of the additions. The last term,
the part of the analyses not reported as oxides, vanishes for a constituent
described by a species and for an analysis that sums to one.

## 5. Volume and heat

The initial volume ``V_0`` is that of the reacted species, of the water and of
what is kept aside; the porosity ``\phi`` and the volume fractions of
[`volume_fractions`](@ref) are referred to it, the void ``V_0 - V_\mathrm{eq} -
V_\mathrm{res}`` being the chemical shrinkage that the products do not refill
([Recipes](@ref man-recipes)). The enthalpy of the paste is that of the
equilibrium plus that of what is kept aside, ``H = H_\mathrm{eq} +
H_\mathrm{res}``, and the heat released at constant temperature and pressure
between two states of one paste is ``Q = H_1 - H_2``. A part kept aside whose
mass is the same in both states contributes the same enthalpy to both and
cancels, known or not; a part whose mass changes and whose enthalpy is unknown
makes ``Q`` unknown, and [`heat_release`](@ref) returns `NaN` rather than a
heat that leaves it out.

## 6. A constituent given a rate

When a constituent is given a rate law ([`KineticsProblem`](@ref) of a recipe),
the rate decides its degree of reaction, and the constituent enters the problem
whole and unreacted, ``n_s(t_0) = m_c / M_s``, whatever ``\alpha_c(t_0)``.
A constituent known by its oxides has no formula to dissolve, so it is first
given one: [`glass_species`](@ref) writes the elements of the oxides of
``K_c`` that ``M`` grams of it carry into one formula unit, ``M`` being the mass
of material the unit stands for, and [`with_species`](@ref) makes it a
constituent of that species. Dissolving ``\alpha_c\, m_c / M`` formula units
then adds to the budget

```math
\frac{\alpha_c\, m_c}{M}\; M \sum_{k \in K_c} \frac{\omega_{c,k}}{M_k}\,\mathbf{d}_k ,
```

the same contribution as the reacted part of the oxide constituent in Section 2,
so that the kinetic route and the route of imposed extents put the same elements
into the equilibrium at the same degree of reaction. The two masses differ,
however: the solid loses ``\alpha_c m_c`` while the equilibrium gains
``\alpha_c m_c \sum_{k \in K_c} \omega_{c,k}``, the rest being what Section 3
keeps aside on the route of imposed extents and what the formula unit does not
model on the kinetic one (`modeled_mass_fraction`).
