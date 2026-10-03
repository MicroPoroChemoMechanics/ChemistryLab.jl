# [Solid solutions](@id sec-theory-solid-solutions)

!!! info "Before this page"
    [Thermochemistry](@ref sec-theory-thermo) and [Standard states](@ref
    sec-theory-standard-states) §2, where the standard state of an end-member is
    set.

A pure solid has ``a = 1``: it is its own standard state, its chemical potential
does not depend on how much of it there is, and it either precipitates or it does
not. A **solid solution** is a single crystalline phase whose composition varies
continuously between end-members, and that changes the thermodynamics
qualitatively: the activity of an end-member depends on the composition of the
phase, so the phase can absorb a component gradually instead of appearing all at
once.

This matters in cement more than anywhere else. C-S-H is not a compound with a
formula but a phase of variable Ca/Si; the AFm phases exchange sulfate, carbonate
and hydroxide on the same interlayer site; alkalis partition into hydrates rather
than staying in solution. A model that only knows pure phases predicts sharp
appearances and disappearances that are not observed.

## 1. Ideal mixing: the activity *is* the mole fraction

Put ``n_k`` moles of each end-member into one phase, with mole fractions
``x_k = n_k/\sum_j n_j``. Even with no interaction at all, mixing lowers the
Gibbs energy, because there are more ways to arrange a mixture than a pure
solid. Per mole of phase,

```math
G^{\text{mix}} = RT\sum_k x_k \ln x_k \;\; (<0) ,
```

which is entropy alone — the ``-T\Delta S`` of mixing. Differentiating for one
end-member gives its chemical potential and hence its activity:

```math
\mu_k = \mu_k^\circ + RT\ln x_k
\qquad\Longrightarrow\qquad
\boxed{\;\ln a_k = \ln x_k\;}
```

That is [`IdealSolidSolutionModel`](@ref), and it is the default. Note what it
already buys: an end-member at ``x_k = 10^{-3}`` has ``\ln a_k = -6.9``, so its
saturation index is shifted by ``-3`` decades relative to the pure phase. A trace
component is *stabilized* by being diluted in a host, which is why a solid
solution can take up an ion that would never precipitate as its own phase.

## 2. Non-ideal mixing: the excess Gibbs energy

Real end-members interact. Everything beyond the ideal term is collected into
the **excess** Gibbs energy ``G^{\text{ex}}``, and the activity coefficient is
its partial molar derivative:

```math
G = \underbrace{\sum_k x_k\mu_k^\circ}_{\text{mechanical mixture}}
  + \underbrace{RT\sum_k x_k\ln x_k}_{\text{ideal mixing}}
  + \underbrace{G^{\text{ex}}(x,T)}_{\text{interaction}} ,
\qquad
\ln\gamma_k = \frac{\partial}{\partial n_k}\!\left(\frac{n\,G^{\text{ex}}}{RT}\right) ,
```

and then ``\ln a_k = \ln x_k + \ln\gamma_k``. This is exactly what the code
does: for these models `_ss_log_activities!` writes
`log(x[k] + ϵ) + _excess_ln_gamma(model, k, x, T)` into the log-activity vector,
so the two terms are visibly separate and a model supplies only the second. The
model of section 7 is written differently, and has its own method.

## 3. `RegularSolutionModel` — one interaction energy per pair

The simplest non-ideal form keeps only pairwise contacts, with one energy per
pair of end-members:

```math
\frac{G^{\text{ex}}}{RT} = \sum_{i<j}\frac{W_{ij}}{RT}\,x_i x_j ,
\qquad
\ln\gamma_k = \sum_{j\neq k}\frac{W_{kj}}{RT}x_j \;-\; \sum_{i<j}\frac{W_{ij}}{RT}x_i x_j ,
```

which is the symmetric multicomponent Margules expression, and it is
`_excess_ln_gamma(::RegularSolutionModel, …)` line for line. For a **binary** it
collapses to the form usually quoted:

```math
\ln\gamma_1 = \frac{W_{12}}{RT}x_2^2 ,
\qquad
\ln\gamma_2 = \frac{W_{12}}{RT}x_1^2 .
```

``W_{ij}`` (J/mol) is the energy of an ``i``–``j`` contact *relative to the
average* of ``i``–``i`` and ``j``–``j``, so its sign is physical:

  - ``W_{ij} > 0`` — unlike neighbors are unfavorable. ``\gamma > 1``, the phase
    resists mixing, and above a threshold it unmixes;
  - ``W_{ij} < 0`` — unlike neighbors are favorable. ``\gamma < 1``, mixing is
    stabilized beyond ideal, and the phase may order.

`W` must be square and symmetric, which a **validating inner constructor**
enforces; the reason it is inner rather than outer is recorded in the source, and
it is a mistake worth not repeating.

### The threshold is ``W = 2RT``, and it is checkable

For a symmetric binary, the molar Gibbs energy of mixing is
``G^{\text{mix}}/RT = x\ln x + (1-x)\ln(1-x) + (W/RT)x(1-x)``. A phase is
unstable to unmixing where that function is concave, and

```math
\frac{\partial^2}{\partial x^2}\frac{G^{\text{mix}}}{RT}
 = \frac{1}{x} + \frac{1}{1-x} - \frac{2W}{RT} ,
```

which at ``x = 1/2`` is ``4 - 2W/RT``: negative as soon as ``W > 2RT``. So
``W/RT = 2`` — about **5 kJ/mol at 25 °C** — is the critical point of a symmetric
regular solution: below it one homogeneous phase, above it a miscibility gap
that widens as ``W`` grows. The threshold is evaluated, and the sign of the
second derivative tabulated on either side of it, in
[Solid solution models, in numbers](@ref sec-app-solid-solutions).

This matters in practice because **nothing detects it**. A solid solution is
entered as one phase, the activity expression goes on returning values past the
threshold, and those values describe a metastable single phase. Checking
``W/RT`` against 2 is the caller's job.

## 4. `RedlichKisterModel` — a binary that need not be symmetric

A regular solution is symmetric by construction: swapping the end-members leaves
``G^{\text{ex}}`` unchanged. Real binaries often are not, so the excess is
expanded in a polynomial of the composition difference. As implemented, for two
end-members,

```math
\ln\gamma_1 = \frac{x_2^2}{RT}\Big[a_0 + a_1(3x_1 - x_2) + a_2(x_1 - x_2)(5x_1 - x_2)\Big] ,
```
```math
\ln\gamma_2 = \frac{x_1^2}{RT}\Big[a_0 - a_1(3x_2 - x_1) + a_2(x_2 - x_1)(5x_2 - x_1)\Big] ,
```

with ``a_0, a_1, a_2`` in J/mol. Read the roles off the formulas: ``a_0`` is the
symmetric term and is exactly a regular solution's ``W_{12}``; ``a_1`` and
``a_2`` are the asymmetric corrections, so setting them to zero must reproduce
the regular model. That is an identity between two independently written
methods, and it is checked rather than trusted in
[Solid solution models, in numbers](@ref sec-app-solid-solutions), where the two
agree to the last digit.

`RedlichKisterModel` **requires exactly two end-members**, and
[`SolidSolutionPhase`](@ref) refuses the combination at construction rather than
letting a ternary reach an expression written for a binary. For three or more
end-members the choices are the ideal model or `RegularSolutionModel` with a
full ``\mathbf{W}`` matrix.

## 5. What the three models do to an activity

Evaluated over the whole composition range in
[Solid solution models, in numbers](@ref sec-app-solid-solutions), with
``W = \pm 4`` kJ/mol — inside the stable range of §3, so that the numbers
describe a phase that stays homogeneous. Two things to read off it.

A positive ``W`` pushes the activity of a dilute end-member **above** its mole
fraction: the host is rejecting it. A negative ``W`` pulls it below: the host is
stabilizing it. And at ``x_k \to 1`` all three models converge, because
``\gamma_k \to 1`` as the phase becomes pure — the standard state of an
end-member is the pure end-member, so they agree there by construction.

The ideal case is the one to keep in mind for cement. An end-member at
``x_k = 10^{-3}`` has ``a_k = 10^{-3}``, so its saturation index sits three
decades below the pure phase: **a trace component is stabilized simply by being
diluted in a host.** That is how a solid solution takes up an ion which would
never precipitate as a phase of its own, and it is the whole reason C-S-H and
the AFm phases can absorb alkalis, sulfate and carbonate continuously instead of
in jumps.

## 6. The spinodal, the common tangent, and why one amount is not enough

Sections 3 and 4 gave models whose excess term can be strong enough to make the
molar Gibbs energy of mixing **concave** over an interval. This section is about
what that means for the equilibrium, because it is not a detail of the model: it
changes what the answer *is*.

### Convexity is the hypothesis under everything else

For a binary write the molar Gibbs energy of mixing, in units of ``RT``, as

```math
\frac{g(x)}{RT} = x\ln x + (1-x)\ln(1-x) + \frac{g^{\mathrm{ex}}(x)}{RT} .
```

The ideal part is convex everywhere — its second derivative is
``1/x + 1/(1-x) > 0`` — so unmixing is always the excess term's doing. Where

```math
\frac{\mathrm{d}^2 g}{\mathrm{d}x^2} < 0
```

the phase is **inside its spinodal**, and there the straight line joining two
compositions lies *below* the curve between them. A system at an intermediate
overall composition therefore lowers its energy by separating into those two, and
the minimum of ``G`` is not a point but a **pair**.

For the symmetric regular solution of section 3 the criterion is exactly
``\mathrm{d}^2g/\mathrm{d}x^2 = 4 - 2W/RT`` at ``x = \tfrac12``, so unmixing
begins at ``W = 2RT`` — the threshold section 3 already quoted, here derived from
the same inequality. [`spinodal_interval`](@ref) evaluates the second derivative
on a grid rather than in closed form, which is why the same routine covers the
asymmetric Redlich-Kister case without a separate derivation.

### The answer inside a gap is the common tangent

Which pair? The two compositions ``x_\alpha < x_\beta`` at which the chemical
potentials of *both* end-members agree:

```math
\left.\frac{\mathrm{d}g}{\mathrm{d}x}\right|_{x_\alpha}
  = \left.\frac{\mathrm{d}g}{\mathrm{d}x}\right|_{x_\beta}
  = \frac{g(x_\beta) - g(x_\alpha)}{x_\beta - x_\alpha} ,
```

which says geometrically that one straight line is tangent to the curve at both
points — the **common tangent**. Between ``x_\alpha`` and ``x_\beta`` the
equilibrium state is a mixture of the two, in the proportions the lever rule
gives, and the energy follows the tangent line rather than the curve.

The common tangent construction is wider than the spinodal: the *binodal*
``[x_\alpha, x_\beta]`` contains the spinodal, and between the two the phase is
metastable rather than unstable. A minimization sees only the tangent.

### Computing the pair: what PHREEQC does, and what it costs

The pair is not found by minimizing. It is **computed from the model alone**, by
[`common_tangent`](@ref): two equations in two unknowns, solved by Newton with
`ForwardDiff` supplying the derivatives. Nothing about the rest of the chemical
system enters — not the aqueous solution, not the other phases, not the element
budget — which is why it costs microseconds and can be used as the starting point
of a full equilibrium rather than as its result.

This is the construction PHREEQC uses for binary solid solutions, after
[GlynnReardon1990](@cite). It was validated here against a case with a closed
form: for a symmetric model the pair must be symmetric about ``x = 1/2``, and it
comes out at ``(0.070720,\ 0.929280)`` with a residual of ``2.3\times10^{-13}``
and the symmetry exact to the last bit.

[`miscibility_split`](@ref) then applies the lever rule. Inside the gap the two
**compositions are fixed** and only their proportions move with the overall
composition — which is what makes a miscibility gap flat in a phase diagram, and
what a single-composition answer cannot reproduce at any resolution.

!!! note "Where the three codes agree, and where this one adds a step"
    The criterion is the **same object** in GEM-Selektor and here, arrived at
    independently from the same KKT conditions — its phase stability index
    ``\Lambda_k = \log_{10}\Omega_k`` is term for term what `phase_split_measure`
    computes [Kulik2013](@cite). The *construction* of the pair is PHREEQC's,
    after Glynn & Reardon. Neither this package nor the others invented either.

    What differs is only where the duplication comes from. GEM-Selektor's users
    get the right answer inside a gap because CEMDATA18 ships the AFm and AFt
    binaries under two names each, so the declaration is already doubled in the
    database; here it is asked for by a keyword. Both are sound, and the database
    route has the advantage of being the published one.

    The step this package adds is the **refusal**: `SolidSolutionPhase` evaluates
    the second derivative at construction and declines a model that unmixes,
    naming the interval — and `OptimaSolver`'s certificate applies the
    tangent-plane test to phases that are **present**, not only to absent ones.
    Assuming convexity and leaving the duplication to the caller is a reasonable
    design for a general-purpose code; checking it is worth the few lines here,
    because a phase sitting inside its own spinodal otherwise certifies on the
    stationarity of its members alone and says nothing.

### Why a formulation with one amount per species cannot hold it

The composition vector of a `ChemicalSystem` has one entry per species, so a
species belongs to exactly one phase and a phase has exactly one composition. A
common-tangent pair is *two compositions of one substance*, and there is no way
to write that down with one amount each.

So the formulation has to be given the substance twice. That is what
`SolidSolutionPhase(...; instances = 2)` does: `ChemicalSystem` builds a second
copy of each end-member under a derived symbol (`monosulphate12#2`) sharing the
same thermodynamic record, and each copy belongs to its own phase. The groups
stay disjoint, the duplicated columns of the conservation matrix are copies of
columns already there — so the row rank, and with it conservation, is untouched —
and the minimization is free to populate either lobe or both.

It is refused for a convex model, and the reason is worth stating because it is
not conservatism. Two instances of a convex phase are **degenerate**: the energy
is the same however the amount is split between them, so the minimum becomes a
flat manifold and the optimizer is asked to pick a point on it arbitrarily.
Inside a spinodal the common-tangent pair is unique and no such direction exists.

### Michelsen's tangent-plane distance, which is what detects it

Representation is one problem; knowing that it is needed is another. A phase
sitting inside its own spinodal satisfies every first-order condition — its
end-members are stationary, the element balance closes, nothing absent is
supersaturated — so stationarity alone certifies a non-minimum.

The test that separates them is Michelsen's: for a trial composition
``\hat{\mathbf{x}}`` of the phase, the **tangent-plane distance**

```math
D(\hat{\mathbf{x}}) = \sum_k \hat{x}_k\,
   \bigl[\,\mu_k(\hat{\mathbf{x}}) - \mu_k(\mathbf{x}^\star)\,\bigr]
```

measures how far the Gibbs surface at ``\hat{\mathbf{x}}`` lies above the tangent
plane drawn at the reported composition ``\mathbf{x}^\star``. If ``D`` is
negative anywhere, some other composition lies *below* that plane and the
reported state is not a minimum. At ``\hat{\mathbf{x}} = \mathbf{x}^\star`` the
distance is exactly zero and is a stationary point of the search, which is why
the measure has to be started from the **end-member corners** — the other lobe is
where the negative value lives.

`OptimaSolver`'s `phase_split_measure` computes this for every present
mole-fraction phase and folds the result into the certificate. Concretely:
`check_convexity = false` stops being a silent loss of the proof — a concave
declaration that does unmix now fails to certify, with the phase named.

### The same criterion, in another code

GEM-Selektor's `PhaseSelection` computes, for every phase, a stability index

```math
\Lambda_k = \log_{10}\Omega_k , \qquad
\Omega_k = \sum_{j\in l_k}\pi_j , \qquad
\pi_j = \frac{\omega_j(\hat{\mathbf{u}})}{\gamma_j(\hat{\mathbf{n}})} ,
```

the activity from the **dual** solution divided by the activity coefficient from
the **primal** one — an estimate of the mole fraction. That is term for term what
this package computes, in `_repair_start`'s ``\Omega = \sum_i 10^{SI_i}`` (the
ideal case ``\gamma = 1``) and in `phase_split_measure`'s log-sum-exp of
``d_j = u_i - g_i - \ln\gamma_i``. The GEMS3K paper calls ``\Omega_k``
*"a generalization of the saturation index"* and notes that it derives from the
KKT conditions [Kulik2013](@cite).

So the criterion is one object, arrived at independently. What differs is the
declaration. CEMDATA18 ships the AFm and AFt binaries under two names each, so a
GEMS user represents a gap by declaring the binary twice in the database;
`instances = 2` asks for the same thing with a keyword. `instances = :auto` asks
for less: the phase is solved with one composition, and a second instance is added
only when the certificate finds it wanting to split. On a binary whose overall
composition the element budget holds inside the gap, calcite and magnesite in
equal amounts, that finds the common-tangent pair:

```@example split
using ChemistryLab, DynamicQuantities, Printf
sp = Dict(symbol(s) => s for s in build_species(datapath("slop98-inorganic-thermofun.json"); verbose = false))
gap = RedlichKisterModel(a0 = 14_000.0)
carbonate = SolidSolutionPhase("carbonate", [sp["Cal"], sp["Mgs"]]; model = gap, instances = :auto)
cs = ChemicalSystem([sp[s] for s in split("H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 Mg+2 Cal Mgs")],
                    ["H2O@", "H+", "Ca+2", "Mg+2", "CO3-2", "Zz"]; solid_solutions = [carbonate])
st = ChemicalState(cs)
set_quantity!(st, "H2O@", 1.0u"kg"); set_quantity!(st, "Cal", 0.025u"mol"); set_quantity!(st, "Mgs", 0.025u"mol")
eq, cert = equilibrate_certified(st; b = Float64.(cs.SM.A) * ustrip.(us"mol", st.n))
n = ustrip.(us"mol", eq.n)
@printf("certified %s, instances added to %s\n", cert.optimal, cert.instances_added)
for (g, ph) in zip(eq.system.ss_groups, eq.system.solid_solutions)
    @printf("  %-12s %.5f mol, x(calcite) = %.4f\n", name(ph), sum(n[g]), n[g[1]] / sum(n[g]))
end
@printf("common tangent of the model: %.4f and %.4f\n", common_tangent(gap)...)
```

It does not make a nonconvex problem easy: the certificate names a phase to split
only once the one-composition solve has converged.

The executed counterpart of this section is
[the miscibility-gap page](@ref ex-miscibility-gap), which runs one cement three
ways and reports the certificate each time.

## 7. Mixing on several sites: `SublatticeModel`

Sections 1 to 4 mix **end-members**. In a C-S-H the entities that actually mix
are smaller: a few structural positions of the silicate chain and of the
interlayer, each held by one species or another. An end-member is then one
particular filling of those positions, and mixing is random on each position
separately. That is the model of Kulik (2011) for CSH3T and of Myers et al.
(2014) for the CNASH gel, and it is what [`SublatticeModel`](@ref) implements.

**The model.** A formula unit has sites ``s``, site ``s`` counted ``m_s`` times,
and end-member ``k`` puts the species ``\sigma_s(k)`` on site ``s``. The fraction
of site ``s`` held by species ``i`` is the total mole fraction of the end-members
that put it there, ``y_{s,i} = \sum_k [\sigma_s(k) = i]\,x_k``, and the
configurational Gibbs energy is the ideal one of each site:

```math
\frac{G^{\text{mix}}}{RT} = \sum_s m_s \sum_i y_{s,i}\ln y_{s,i}
\qquad\Longrightarrow\qquad
\boxed{\;\ln a_k = \sum_s m_s \ln y_{s,\sigma_s(k)}\;}
```

The implication is ``\partial(n\,G^{\text{mix}}/RT)/\partial n_k``, exact because
every end-member fills every site once. One site counted once is section 1 again.
The papers write the same model as a *fictive activity coefficient*
``\lambda_k = a_k/x_k``, which is what [`excess_ln_gamma_expression`](@ref)
returns for this model.

**CSH3T.** Two bridging-tetrahedral sites hold either Si or Ca. TobH puts Si on
both, T2C puts Ca on both, and the ordered T5C puts Ca on the first and Si on the
second (Kulik 2011, Eq. 19). The Cemdata18 records are half Kulik's formula
units, so on them each site counts one half:

```@example sublattice
using ChemistryLab
subs = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-thermofun.json"); verbose = false))
csh3t = sublattice_model("Kulik2011:csh3t", [subs[n] for n in ("CSH3T-TobH", "CSH3T-T5C", "CSH3T-T2C")])
csh3t.multiplicity
```

`sublattice_model` finds that factor itself, by comparing each record's formula
with the one the paper prints, and refuses a record that is not such a multiple.
At ``\mathbf{x} = (0.2, 0.5, 0.3)`` the first site holds Si only through TobH, the second
through TobH and T5C:

```@example sublattice
site_fractions(csh3t, [0.2, 0.5, 0.3])
```

and the fictive activity coefficient of T5C is Kulik's Eq. (20), halved:

```math
\ln\lambda_{\mathrm{T5C}} = \tfrac12\ln(x_{\mathrm{TobH}} + x_{\mathrm{T5C}})
  + \tfrac12\ln(x_{\mathrm{T2C}} + x_{\mathrm{T5C}}) - \ln x_{\mathrm{T5C}}
```

The block below builds that expression from the model and subtracts the formula
above; the difference simplifies to zero:

```@example sublattice
using Symbolics
x = [Symbolics.variable(:x, i) for i in 1:3]   # x₁ = TobH, x₂ = T5C, x₃ = T2C
λ = excess_ln_gamma_expression(csh3t, 2, 3)
simplify(expand(λ - (0.5log(x[1] + x[2]) + 0.5log(x[3] + x[2]) - log(x[2]))))
```

**Two properties the solver relies on.**

  - *Convexity.* ``G^{\text{mix}}`` is a sum of convex functions of the site
    fractions, which are linear in ``\mathbf{x}``, so it is convex, and strictly so when
    the occupancy matrix (one row per site and species, one column per
    end-member) has full column rank, stored as `rank`. The activities are the
    gradient of that energy, so a certificate keeps its `:global_minimum` scope.
  - *The dilute limit.* As ``x_k \to 0`` only the site fractions of the species
    ``k`` alone carries vanish with it: ``\ln a_k \simeq e_k \ln x_k``, with
    ``e_k`` the multiplicity of those sites, stored as `exponents`. T5C owns no
    species of its own, ``e_k = 0``: its activity stays finite as it disappears,
    and it can be absent from a CSH3T that is present.

```@example sublattice
csh3t.exponents
```

**What the solver does with it.** The composition of a mixing phase is found by
successive substitution, and on a sublattice model that iteration diverges: its
rate is set by the sum of the multiplicities, nine for the CNASH gel. Such a phase
is therefore inverted by Newton's method, and the members with ``e_k = 0`` are
declared to the certificate as able to leave the phase, which then tests them by
the inequality a pure phase obeys (OptimaSolver's `SolutionPhase`, keywords
`newton` and `bounded_members`). A system with a sublattice phase is first solved
with that phase mixing its end-members ideally, and the answer is the start of
the sublattice solve.

`data/solid_solutions.toml` ships `CNASH_ss` as Myers' model and `CSH3T` as
Cemdata18 ships it, ideal between end-members; `sublattice_model("Kulik2011:csh3t",
members)` gives the site form.

## 8. Sites with the energies of the compounds: `CompoundEnergyModel`

Section 7 mixes ideally on the sites and keeps the end-members' own standard
energies, ``\sum_k x_k G^\circ_k``. That is exact only when the end-members are
independent in the site fractions. The CASH+ model of C-S-H [Kulik2022](@cite)
has two mixing sites: the bridging tetrahedron holds a silicate, a vacancy or a
calcium (S, v, C), and the interlayer a vacancy or a calcium (v, C). Its six
end-members are the six ways of filling them, TSvh to TCCh, but three site
fractions fix the composition, and the six members are not independent. The
reciprocal reaction

```math
\mathrm{TSvh} + \mathrm{TCCh} = \mathrm{TSCh} + \mathrm{TCvh}
```

leaves every site fraction unchanged, yet its standard Gibbs energy is not zero.
The compound energy formalism, in which Kulik et al. write the model, gives such
a reaction its energy. [`CompoundEnergyModel`](@ref) implements it.

**The Gibbs energy.** The end-members must be every *compound* of the sites, each
once. Compound ``j`` puts species ``j_s`` on site ``s`` (the ``\sigma_s(j)`` of
section 7). Per formula unit,

```math
G = \underbrace{\sum_j G^\circ_j \prod_s y_{s,j_s}}_{G_\text{ref}}
  + RT \sum_s m_s \sum_i y_{s,i}\ln y_{s,i}
  + \sum_s \sum_{i<l} W_{s,il}\, y_{s,i}\, y_{s,l} .
```

The first term, the *reference surface*, equals ``G^\circ_j`` at compound ``j``
and is linear in the fractions of each site, so a reaction between compounds that
changes no site fraction keeps its energy. The last term is a regular interaction
between two species of one site (Berman's symmetric form, with the parameters
``W`` of the paper). ``G`` depends on the site fractions alone.

**The activities.** With ``y_{s,i} = n_{s,i}/n``, the chemical potential of
member ``k`` is ``\mu_k = \partial(nG)/\partial n_k``. For a function of the site
fractions, ``\partial y_{s,i}/\partial n_k = ([k_s = i] - y_{s,i})/n``, so

```math
\mu_k = G + \sum_s\Big(\frac{\partial G}{\partial y_{s,k_s}} - \sum_i y_{s,i}\frac{\partial G}{\partial y_{s,i}}\Big).
```

On the reference surface ``\sum_i y_{s,i}\,\partial G_\text{ref}/\partial y_{s,i} = G_\text{ref}``
for every site, because ``G_\text{ref}`` is linear in each site's fractions. With
``n_\text{site}`` mixing sites, and writing ``\ln a_k = (\mu_k - G^\circ_k)/RT``,

```math
\ln a_k = \sum_s m_s \ln y_{s,k_s}
  + \frac{1}{RT}\Big(\sum_s \frac{\partial G_\text{ref}}{\partial y_{s,k_s}} - (n_\text{site}-1)\,G_\text{ref} - G^\circ_k\Big)
  + \frac{1}{RT}\sum_s\Big(\sum_l W_{s,k_s l}\, y_{s,l} - \sum_{i<l} W_{s,il}\, y_{s,i}\, y_{s,l}\Big).
```

The middle term vanishes at every compound. It also does not change when each
``G^\circ_j`` is shifted by the energies of its elements, since those are a sum
over the sites of the species' shares, which the reference surface interpolates
exactly. So the convention of the database does not matter, but the energies of
the members do: the solver hands them to the model at the temperature of the
solve.

**Three consequences.**

  - *The amounts are not unique.* Every split of the same site fractions between
    the members has the same ``G`` and the same element content: six members for
    four independent site fractions in CASH+. Left to choose, the solver returns one
    of them, which one depending on where the search started. It is therefore given
    one split, the product of the site fractions, ``x_j = \prod_s y_{s,j_s}``, by adding to the energy it
    minimizes ``RT\,D(\mathbf{x})``, with

    ```math
    D(\mathbf{x}) = \sum_j x_j \ln x_j - \sum_s \sum_i y_{s,i}\ln y_{s,i} .
    ```

    ``D`` is the divergence of ``\mathbf{x}`` from the product of its own site fractions:
    it is never negative, and it vanishes only at that product, where its gradient
    ``\ln x_k - \sum_s \ln y_{s,k_s}`` vanishes too. The minimum over the splits
    of one set of site fractions is therefore the product, with the energy and
    the chemical potentials of the model. The activity of a vanishing member now
    goes as its mole fraction, as in a mixture of end-members. Measured on twelve
    pastes of CASH+NK, the equilibria certify with or without ``D``, at the same
    cost; what ``D`` adds is an answer that does not depend on the start.
  - *Convexity must be decided.* The reference surface is multilinear, and a
    reciprocal energy can make ``G`` concave somewhere. In the site fractions, each
    site contributes a curvature of at least ``2m_s + \mu_s``, the configurational
    bound of section 9 plus the smallest tangent eigenvalue ``\mu_s`` of its
    interactions ``\mathbf{W}/RT``, and with two sites the reference surface couples them
    through a constant matrix ``\mathbf{C}``, the ``G^\circ_j/RT`` of the compounds projected
    on the two tangent spaces. So ``G`` is convex when
    ``(2m_1 + \mu_1)(2m_2 + \mu_2) > \lVert \mathbf{C} \rVert^2``. That proves the CASH+
    core convex (34.3 against 13.2 at 25 °C), and a certificate on it keeps its
    `:global_minimum` scope as far as the gel goes; the bound does not decide
    CASH+NK, whose sampled Hessian shows no concave point, and its scope is then
    `:local_minimum` or `:kkt_point`, according to the curvature at the answer
    ([`mixing_convexity`](@ref)).
  - *Euler's relation holds.* ``\sum_k x_k\,\mu_k = G`` for any amounts with the
    site fractions ``\mathbf{y}``, as it must for a Gibbs energy.

The block below checks the gradient and the last point on the core model, at a
composition of the six members. The activities are compared with the gradient of
``n(G/RT + D)`` coded from the formulas above and differentiated by ForwardDiff;
at the product split, ``D`` and its gradient vanish:

```@example cef
using ChemistryLab, DynamicQuantities, ForwardDiff, LinearAlgebra
subs = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-cashplus.json"); verbose = false))
members = ["TSvh", "TSCh", "Tvvh", "TCvh", "TvCh", "TCCh"]
cash = compound_energy_model("Kulik2022:cashplus", [subs[m] for m in members])
T = 298.15
g = [ustrip(us"J/mol", subs[m][:ΔₐG⁰](T = T * u"K", P = 1.0e5u"Pa"; unit = true)) for m in members] ./ (ChemistryLab.R_GAS * T)
o = cash.lattice.occupancy
function nG(n; split = true)
    x = n ./ sum(n)
    y = site_fractions(cash, x)
    G = sum(g[j] * y[1][o[1, j]] * y[2][o[2, j]] for j in 1:6) + sum(v * log(v) for ys in y for v in ys)
    G += sum(W / (ChemistryLab.R_GAS * T) * y[s][i] * y[s][l] for (s, i, l, W) in cash.interactions)
    split && (G += sum(v * log(v) for v in x) - sum(v * log(v) for ys in y for v in ys))   # D(x)
    return sum(n) * G
end
n = [0.10, 0.25, 0.05, 0.20, 0.15, 0.25]
lna = ChemistryLab._ss_log_activities!(zeros(6), 1:6, n ./ sum(n), cash, T, 0.0, g)
μ = ForwardDiff.gradient(nG, n)
y = site_fractions(cash, n ./ sum(n))
xp = [y[1][o[1, j]] * y[2][o[2, j]] for j in 1:6]          # the product split
(gradient = maximum(abs, lna .+ g .- μ), euler = dot(n, μ) - nG(n),
 D_at_product = nG(xp) - nG(xp; split = false))
```

All three are at the rounding of the arithmetic. The model and its data are
in `data/literature/Kulik2022.json`, the sodium and potassium of Miron et al.
[Miron2022a, Miron2022b](@citet) in `Miron2022a.json` and `Miron2022b.json`, and the
twelve end-members in the database `cemdata18-cashplus.json`.
[The CASH+ page](@ref ex-cashplus-csh) computes the C-S-H in water and in alkali
solutions with it.

## 9. More than two end-members

### Convexity with more than two end-members

A one-dimensional scan decides a binary. With ``n`` end-members the question is
whether the Hessian of ``g/RT`` is positive on the tangent space of the simplex,
``\{\mathbf{d} : \sum_i d_i = 0\}``, at every composition, and
[`mixing_convexity`](@ref) answers it by model.

For a regular model, ``g/RT = \sum_i x_i \ln x_i + \sum_{i<j} w_{ij} x_i x_j`` with
``w_{ij} = W_{ij}/RT``, and the Hessian is ``\operatorname{diag}(1/\mathbf{x}) + \mathbf{W}/RT``. The first term
is at least 2 on the tangent space: for a unit ``\mathbf{d}`` with ``\sum_i d_i = 0``, the
Cauchy–Schwarz inequality gives
``\sum_i d_i^2/x_i \ge (\sum_i |d_i|)^2 / \sum_i x_i = (\sum_i |d_i|)^2``, and the
ℓ1 norm of such a ``\mathbf{d}`` is at least ``\sqrt 2`` (its positive and negative parts
have equal sums). Hence

```math
\lambda_{\min}\!\left(\mathbf{Q}^\mathsf{T}\,\frac{\mathbf{W}}{RT}\,\mathbf{Q}\right) \ge -2
\quad\Longrightarrow\quad \text{convex},
```

``\mathbf{Q}`` an orthonormal basis of the tangent space. For two end-members this is
``W \le 2RT``, the threshold of section 3; for more it is a sufficient
condition. In the other direction, a pair with ``W_{ij} > 2RT`` is a witness:
at the middle of that edge the second derivative along it is ``4 - 2w_{ij} < 0``.
Between the two, the smallest eigenvalue of the projected Hessian is searched on
a lattice of compositions, and a negative one is a witness too; finding none
proves nothing, and the verdict is then `:undecided`.

```@example ss_convexity
using ChemistryLab, DynamicQuantities, Printf
RT = ChemistryLab.R_GAS * 298.15
W(a, b, c) = [0.0 a b; a 0.0 c; b c 0.0] .* RT
[mixing_convexity(RegularSolutionModel(W(1.9, 1.9, 1.9)), 3).verdict,
 mixing_convexity(RegularSolutionModel(W(3.0, 0.0, 0.0)), 3).verdict]
```

A phase whose verdict is `:nonconvex` is refused at construction, as a concave
binary is, and admitted with two instances. The certificate of an answer is not
scoped `:global_minimum` unless every mixing phase is proved convex; it is then
`:local_minimum` where the Hessian over the directions that conserve matter is
positive definite at the answer, and `:kkt_point` otherwise.

### A member that is a mixture of two others

Two different things can make one end-member the mixture of two others in
composition. In an **ordered** member the difference of Gibbs energy is the point:
the pentameric T5C of CSH3T is the average of TobH and T2C less an ordering
energy (Kulik 2011, Eq. 18), and the siliceous hydrogarnet
C3AS0.41H5.18 lies between C3AH6 and C3AS0.84H4.32. When the Gibbs energy is
the average as well, the phase holds one substance twice, once as a member and
once as a mixture of two, and ideal mixing counts its configurations twice.
[`SolidSolutionPhase`](@ref) warns when that holds within ``0.1\,RT``; the
MgAl-OH-LDH ternary of Cemdata18 is the case, its M6A member the average of M4A
and M8A:

```@example ss_convexity
subs = Dict(symbol(s) => s for s in build_species(datapath("cemdata18-thermofun.json"); verbose = false))
g(s) = ustrip(us"kJ/mol", subs[s][:ΔₐG⁰](T = 298.15u"K", P = 1.0e5u"Pa"; unit = true))
@printf("MgAl-OH-LDH  M6A - (M4A + M8A)/2  = %+.4f kJ/mol\n",
        g("M6A-OH-LDH") - (g("M4A-OH-LDH") + g("M8A-OH-LDH")) / 2)
@printf("CSH3T        T5C - (TobH + T2C)/2 = %+.2f kJ/mol\n",
        g("CSH3T-T5C") - (g("CSH3T-TobH") + g("CSH3T-T2C")) / 2)
```

The LDH difference, 0.3 J/mol, is about 10⁻⁴ RT, where the −4.35 kJ/mol of T5C
is an ordering energy. The shipped LDH entry keeps the published model with
`acknowledge_degenerate = true`.

## 10. How a solid solution is declared

A [`SolidSolutionPhase`](@ref) names its end-members and carries a model:

```julia
phase = SolidSolutionPhase("CSHQ", [sp["CSHQ-JenD"], sp["CSHQ-JenH"],
                                    sp["CSHQ-TobD"], sp["CSHQ-TobH"]];
                           model = IdealSolidSolutionModel())
```

Every end-member must be `AS_CRYSTAL`, and the phase is checked against its
model's arity **and against the convexity of its mixing energy** at construction;
`instances` is how a declaration asks for more than one coexisting composition,
which section 6 is about. Sets of end-members with their interaction
parameters are read from a TOML through
[`build_solid_solutions`](@ref), which is also the only interaction-parameter
file format in the package — `data/solid_solutions.toml` is the shipped example.

Two consequences to keep in mind when reading results:

  - a `ChemicalSystem` groups end-members into `ss_groups`, and **every** aqueous
    activity model must call the solid-solution branch or the end-members come
    back with ``\ln a = 0``, i.e. treated as pure phases;
  - an end-member's [`saturation_indices`](@ref) entry is relative to its
    **current** mole fraction, since its activity is ``x_k\gamma_k``. For an
    end-member sitting at the solver's lower bound that is a statement about a
    vanishing phase, not about whether the solid solution would form.

## See also

  - [Thermochemistry](@ref sec-theory-thermo) — where ``\ln a`` enters ``\mu``
  - [Activity models](@ref sec-theory-activity) — the aqueous half
  - [Solving an equilibrium](@ref sec-solving) — declaring and solving with
    solid solutions
