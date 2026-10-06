# [Kinetics under partial equilibrium: the right-hand side and the calorimeters](@id sec-theory-pe-kinetics)

!!! info "Before this page"
    [Rate laws, and every parameter in them](@ref sec-theory-kinetics), which
    writes the rates ``\mathbf{r}``; [Proving that an answer is the
    answer](@ref sec-theory-certificate), for the minimization and its
    certificate; and [The two laws, and what the Gibbs energy
    measures](@ref sec-theory-laws), where the enthalpy a calorimeter measures
    is defined.

A hydrating paste holds reactions on two time scales. The aqueous speciation,
protonation, complexation and the autoprotolysis of water, reaches equilibrium
in microseconds; alite dissolves over days. Integrating every reaction with a
rate law would need kinetic constants nobody measures for the fast ones and
would bring the integrator down to their time scale; equilibrating every phase
would dissolve the clinker at once. Partial equilibrium integrates the slow
reactions and solves the fast ones by a minimization of the Gibbs energy at each
instant [Leal2015, Leal2017](@cite).

The equilibria computed this way are conditional on the kinetic amounts of the
instant. They are not jumps over activation barriers, and they do not say that
the paste has reached its stable equilibrium: they are minima of ``G`` on the
compositions the slow reactions have made accessible, as [Metastable does not
mean a local minimum of G](@ref sec-theory-metastable) explains.

Section 1 sets out the matrix every later equation is written with, the
stoichiometric matrix, over the elements and over the primary species. The
package then advances such a system in two ways, which answer two questions.
Sections 2 to 4 write the differential system of [Leal2015](@citet), whose
right-hand side holds the minimization and whose Jacobian is exact; Section 5
the implicit step of [Leal2017](@citet), where the step itself is one
minimization; Section 6 extends the first to the energy balance of a
calorimeter.

## [1. The stoichiometric matrix](@id sec-theory-pe-matrix)

Everything below rests on one matrix: what each species is made of. The
partition, the state of the integration and the reactions it advances are all
written with it, and the two ways of writing it, over the elements or over a
few chosen species, are the source of the one error a kinetic scheme makes most
easily, a budget kept in one basis and spent in the other.

```@raw html
<figure style="margin: 1.2em 0; text-align: center;">
<svg viewBox="0 0 860 250" width="100%" style="max-width: 860px; font-family: inherit;" role="img" aria-label="From the formulas of the species to their reactions">
  <defs>
    <marker id="pe-arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">
      <path d="M 0 0 L 10 5 L 0 10 z" fill="currentColor"/>
    </marker>
  </defs>
  <g fill="none" stroke="currentColor" stroke-width="1.4">
    <rect x="10" y="40" width="185" height="150" rx="10" fill="rgba(70,130,180,0.12)"/>
    <rect x="230" y="40" width="185" height="150" rx="10" fill="rgba(70,130,180,0.12)"/>
    <rect x="450" y="40" width="185" height="150" rx="10" fill="rgba(46,139,87,0.14)"/>
    <rect x="670" y="40" width="180" height="150" rx="10" fill="rgba(205,133,63,0.16)"/>
    <line x1="195" y1="115" x2="228" y2="115" marker-end="url(#pe-arrow)"/>
    <line x1="415" y1="115" x2="448" y2="115" marker-end="url(#pe-arrow)"/>
    <line x1="635" y1="115" x2="668" y2="115" marker-end="url(#pe-arrow)"/>
  </g>
  <g fill="currentColor" font-size="14" text-anchor="middle">
    <text x="102" y="30" font-weight="bold">formulas</text>
    <text x="322" y="30" font-weight="bold">over the elements</text>
    <text x="542" y="30" font-weight="bold">over the primaries</text>
    <text x="760" y="30" font-weight="bold">reactions</text>
    <text x="102" y="75">N species,</text>
    <text x="102" y="97">each a column of</text>
    <text x="102" y="119">atom counts and</text>
    <text x="102" y="141">its charge</text>
    <text x="322" y="75" font-style="italic">W</text>
    <text x="322" y="97">E + 1 rows, N columns</text>
    <text x="322" y="119">b = W n</text>
    <text x="322" y="141">rank C ≤ E + 1</text>
    <text x="322" y="168" font-size="12">a dependent row is redundant</text>
    <text x="542" y="75">C independent columns:</text>
    <text x="542" y="97">the primary species</text>
    <text x="542" y="125" font-style="italic">A = [ I  |  ν̂ᵀ ]</text>
    <text x="542" y="150">W = W<tspan baseline-shift="sub" font-size="10">P</tspan> A</text>
    <text x="542" y="172" font-size="12">coordinates in the basis</text>
    <text x="760" y="75">the null space of A,</text>
    <text x="760" y="97">dimension M = N − C</text>
    <text x="760" y="125" font-style="italic">ν = [ −ν̂  |  I ]</text>
    <text x="760" y="150">A νᵀ = 0</text>
    <text x="760" y="172" font-size="12">one per secondary species</text>
    <text x="430" y="228" font-size="13">The columns are ordered primaries first; ν̂ holds the coordinates of each secondary species on the primaries.</text>
  </g>
</svg>
</figure>
```

### Over the elements

A species is a column of numbers, the atoms of each element in its formula and
its charge. Placed side by side, the columns of the ``N`` species of a system
form its **formula matrix** ``\mathbf{W}``, with one row per element and one for
the charge, `Zz` in the package, and the amounts of the elements in a
composition ``\mathbf{n}`` are

```math
\mathbf{b} = \mathbf{W}\,\mathbf{n} ,
```

the balance every equilibrium and every reaction conserves
[Leal2017](@cite). The carbonate system, water with dissolved carbon dioxide and
its ions, is the smallest one that shows everything that follows:

```@example pe_matrix
using ChemistryLab, LinearAlgebra

aq(f, c = SC_AQSOLUTE) = Species(f; aggregate_state = AS_AQUEOUS, class = c)
H2O = aq("H2O", SC_AQSOLVENT)
Hp, OHm, CO2, HCO3m, CO3mm = aq("H+"), aq("OH-"), aq("CO2"), aq("HCO3-"), aq("CO3-2")
carbonate = [H2O, Hp, OHm, CO2, HCO3m, CO3mm]
W = CanonicalStoichMatrix(carbonate)
pprint(W; label = :symbol)
```

Four rows, six columns, and yet only three of the rows are independent:

```@example pe_matrix
println("rank W = ", rank(Float64.(W.A)))
```

The fourth row is a combination of the other three. Weighting the hydrogen by
``+1``, the oxygen by ``-2`` and the carbon by ``+4``, their oxidation numbers,
gives the charge of every species, so the rows satisfy
``\mathrm{H} - 2\,\mathrm{O} + 4\,\mathrm{C} - \mathrm{Z} = 0``, and conserving
three of the four quantities conserves the fourth. That is the statement that no
electron changes hands between these species, which
[Oxidation state](@ref theory-redox) turns into a test. A minimization over all four
rows would carry a constraint that adds nothing and a matrix that is not of full
rank; the rank, not the number of elements, counts the independent balances.

### Primary species: a basis of the columns

Linear algebra says more. The ``N`` columns of ``\mathbf{W}`` live in a space of
dimension ``C = \operatorname{rank}\mathbf{W}``, the space of the compositions
the system can express, and any ``C`` of them that are linearly independent are a
**basis** of it. The species of such a basis are the **primary species**; the
``N - C`` others are the **secondary species**. Being a basis, the primary
columns decompose every other column in exactly one way,

```math
\mathbf{w}_s = \sum_{p} \hat\nu_{sp}\,\mathbf{w}_p ,
```

and each decomposition is a balanced reaction, the formation of the secondary
species ``s`` from the primaries, ``\sum_p \hat\nu_{sp}\,\mathrm{P}_p \rightleftharpoons
\mathrm{S}_s``. The coefficients ``\hat\nu_{sp}`` form the **canonical
stoichiometric matrix** ``\hat{\boldsymbol\nu}``, one row per secondary species.

The same decomposition rewrites the whole formula matrix. With ``\mathbf{W}_P``
the columns of the primaries, ``\mathbf{W} = \mathbf{W}_P\,\mathbf{A}``, where the
column ``i`` of ``\mathbf{A}`` holds the coordinates of species ``i`` in the
basis: a unit vector for a primary, ``\hat{\boldsymbol\nu}`` for a secondary.
Ordering the primaries first,

```math
\mathbf{A} = \begin{bmatrix} \mathbf{I}_C & \hat{\boldsymbol\nu}^\mathsf{T} \end{bmatrix},
\qquad
\mathbf{W}\,\mathbf{n} = \mathbf{W}_P\,(\mathbf{A}\,\mathbf{n}) .
```

``\mathbf{A}`` is the **conservation matrix over the primaries**, `SM.A` in the
package, of ``C`` rows, all independent. Since ``\mathbf{W}_P`` has independent
columns, ``\mathbf{W}\mathbf{n} = \mathbf{W}\mathbf{n}'`` holds exactly when
``\mathbf{A}\mathbf{n} = \mathbf{A}\mathbf{n}'``: the two matrices state the same
conservation, ``\mathbf{A}`` with no redundant row, and the budget over the
primaries ``\mathbf{b}_A = \mathbf{A}\mathbf{n}`` is the budget over the elements
written in another basis, ``\mathbf{b} = \mathbf{W}_P\,\mathbf{b}_A``. With water,
the proton and dissolved carbon dioxide as primaries:

```@example pe_matrix
SM = StoichMatrix(carbonate, [H2O, Hp, CO2])
pprint(SM; label = :symbol)
```

The column of the carbonate ion reads ``\mathrm{CO_3^{2-}} = \mathrm{CO_2} +
\mathrm{H_2O} - 2\,\mathrm{H^+}``: a negative coordinate is no difficulty, the
decomposition being one of vectors, not of matter. The basis is a choice, and
another one gives other coordinates for the same species,

```@example pe_matrix
pprint(StoichMatrix(carbonate, [H2O, Hp, HCO3m]); label = :symbol)
```

while the budgets it describes are the same compositions:

```@example pe_matrix
# mol, in the order of `carbonate`; the proton balances the anions
n = [55.5, 1e-4 + 1e-10 + 2e-9, 1e-10, 1e-2, 1e-4, 1e-9]
A₁ = Float64.(SM.A)
A₂ = Float64.(StoichMatrix(carbonate, [H2O, Hp, HCO3m]).A)
println("over H₂O, H⁺, CO₂:   b = ", round.(A₁ * n; digits = 6))
println("over H₂O, H⁺, HCO₃⁻: b = ", round.(A₂ * n; digits = 6))
println("over the elements:  b = ", round.(Float64.(W.A) * n; digits = 6))
```

Three budgets, one composition. In the first basis the row of the proton is the
row of the charge, so its coordinate is the charge of the solution, zero here,
and it would be negative for a solution holding a base: a coordinate, not an
amount. Which basis is chosen changes the numbers, not the constraints; what it
may not do is change in the middle of a calculation. The package uses the primaries of
the system throughout, a choice made once ([Formation from primary
species](@ref sec-theory-primaries) for the energies that go with it).

### Reactions: the null space

A reaction is a vector of coefficients ``\boldsymbol\nu``, negative for what it
consumes and positive for what it produces, and it is balanced when it conserves
every component, ``\mathbf{A}\boldsymbol\nu = \mathbf{0}``. The balanced reactions
are therefore the **null space** of ``\mathbf{A}``, of dimension
``M = N - C``: there are as many independent reactions as secondary species,
and the canonical reactions are a basis of them. Stacked as rows, they form the
**stoichiometric matrix** ``\boldsymbol\nu`` of ``M`` rows and ``N`` columns,

```math
\boldsymbol\nu = \begin{bmatrix} -\hat{\boldsymbol\nu} & \mathbf{I}_M \end{bmatrix},
\qquad
\mathbf{A}\,\boldsymbol\nu^\mathsf{T}
= \begin{bmatrix} \mathbf{I}_C & \hat{\boldsymbol\nu}^\mathsf{T} \end{bmatrix}
  \begin{bmatrix} -\hat{\boldsymbol\nu}^\mathsf{T} \\ \mathbf{I}_M \end{bmatrix}
= \mathbf{0} ,
```

and any other balanced reaction is a combination of these. A minimization of the
Gibbs energy needs no list of reactions at all: it moves ``\mathbf{n}`` within
this null space, which ``\mathbf{A}\mathbf{n} = \mathbf{b}`` defines
([Proving that an answer is the answer](@ref sec-theory-certificate)). A kinetic
scheme needs some, the slow ones, and writes them in the same way.

```@example pe_matrix
for r in reactions(SM)
    println(r.equation)
end
```

```@example pe_matrix
ν = [r[s] for r in reactions(SM), s in carbonate]     # one row per reaction, r[s] < 0 for a reactant
println("A νᵀ = 0: ", iszero(SM.A * transpose(ν)))
```

## [2. The partition, and why the state carries element amounts](@id sec-theory-pe-partition)

The species are split into a **kinetic partition**, of amounts
``\mathbf{n}_k``, the phases whose transformation a rate law controls, and an
**equilibrium partition**, of amounts ``\mathbf{n}_e``, everything taken as
equilibrated at every instant: the aqueous species and the phases free to
precipitate or dissolve. The kinetic reactions, of rates ``\mathbf{r}``, have the
stoichiometric matrix ``\boldsymbol{\nu} = [\boldsymbol{\nu}_e\;\boldsymbol{\nu}_k]``
split column-wise between the two. Each is oriented so that its controlling
mineral carries the coefficient ``-1``, a positive rate being a dissolution: a
reaction generated from the nullspace of the conservation matrix comes out with
an arbitrary sign, and taken as it is it would grow the clinker.

The obvious state, every amount ``(\mathbf{n}_e, \mathbf{n}_k)``, does not work.
Advancing ``\mathbf{n}_e`` along ``\boldsymbol{\nu}_e^\mathsf{T}\mathbf{r}`` moves the
equilibrium species along the kinetic reactions without re-equilibrating them,
so that after one step the hypothesis that defined the partition no longer
holds; and re-equilibrating afterwards does not repair it, because what the
fast reactions conserve is not a composition. A dissolution reaction written
from the primary species is written in ``\mathrm{H^+}``, which a cement paste
does not contain: along the way the amount of ``\mathrm{H^+}`` "consumed" would
go negative, whereas the amount of hydrogen in the partition never does. The
state therefore carries the element amounts of the equilibrium partition,

```math
\mathbf{b}_e = \mathbf{A}_e\,\mathbf{n}_e ,
```

with ``\mathbf{A}_e`` the conservation matrix restricted to the partition. Every
fast reaction conserves ``\mathbf{b}_e`` by construction, only the kinetic
reactions move it, and the minimizer, not the caller, distributes it over a
feasible composition.

```@raw html
<figure style="margin: 1.2em 0; text-align: center;">
<svg viewBox="0 0 860 300" width="100%" style="max-width: 860px; font-family: inherit;" role="img" aria-label="The kinetic and the equilibrium partitions exchanging elements">
  <defs>
    <marker id="pe-arrow2" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">
      <path d="M 0 0 L 10 5 L 0 10 z" fill="currentColor"/>
    </marker>
  </defs>
  <g fill="none" stroke="currentColor" stroke-width="1.4">
    <rect x="15" y="40" width="300" height="190" rx="12" fill="rgba(205,133,63,0.16)"/>
    <rect x="545" y="40" width="300" height="190" rx="12" fill="rgba(70,130,180,0.12)"/>
    <line x1="318" y1="110" x2="542" y2="110" marker-end="url(#pe-arrow2)"/>
    <line x1="542" y1="160" x2="318" y2="160" marker-end="url(#pe-arrow2)" stroke-dasharray="5 4"/>
  </g>
  <g fill="currentColor" font-size="14" text-anchor="middle">
    <text x="165" y="30" font-weight="bold">kinetic partition</text>
    <text x="695" y="30" font-weight="bold">equilibrium partition</text>
    <text x="165" y="75" font-style="italic">n<tspan baseline-shift="sub" font-size="10">k</tspan></text>
    <text x="165" y="100">the minerals a rate law holds:</text>
    <text x="165" y="122">clinker phases, glasses, calcite</text>
    <text x="165" y="158">dn<tspan baseline-shift="sub" font-size="10">k</tspan>/dt = ν<tspan baseline-shift="sub" font-size="10">k</tspan>ᵀ r</text>
    <text x="165" y="195" font-size="12">integrated by the ODE solver</text>
    <text x="695" y="75" font-style="italic">n<tspan baseline-shift="sub" font-size="10">e</tspan> = φ(b<tspan baseline-shift="sub" font-size="10">e</tspan>, T)</text>
    <text x="695" y="100">the solution, the hydrates,</text>
    <text x="695" y="122">solid solutions, surface sites</text>
    <text x="695" y="158">min G  subject to  A<tspan baseline-shift="sub" font-size="10">e</tspan> n = b<tspan baseline-shift="sub" font-size="10">e</tspan></text>
    <text x="695" y="195" font-size="12">solved at every evaluation</text>
    <text x="430" y="100" font-size="13">elements released by the rates r</text>
    <text x="430" y="182" font-size="13">composition read by the rates</text>
    <text x="430" y="268" font-size="14">state  u = (b<tspan baseline-shift="sub" font-size="10">e</tspan>, n<tspan baseline-shift="sub" font-size="10">k</tspan>, ξ),   db<tspan baseline-shift="sub" font-size="10">e</tspan>/dt = A<tspan baseline-shift="sub" font-size="10">e</tspan> ν<tspan baseline-shift="sub" font-size="10">e</tspan>ᵀ r = −A<tspan baseline-shift="sub" font-size="10">k</tspan> ν<tspan baseline-shift="sub" font-size="10">k</tspan>ᵀ r</text>
  </g>
</svg>
</figure>
```

On calcite dissolving in the carbonate solution, calcite alone is kinetic. Its
canonical reaction on the primaries is its formation; turned so that calcite
carries ``-1``, it is its dissolution, which consumes two protons, so that the
equilibrium partition receives the calcium, the carbon and the oxygen of the
calcite and gives up two protons per mole dissolved:

```@example pe_matrix
Cal = Species("CaCO3"; aggregate_state = AS_CRYSTAL, class = SC_COMPONENT)
Cap = aq("Ca+2")
species = [carbonate; Cap; Cal]
SMc = StoichMatrix(species, [H2O, Hp, CO2, Cap])
formation = only(r for r in reactions(SMc) if r[Cal] != 0)
println(formation.equation)
ν = [formation[s] for s in species]
ν = ν ./ -ν[end]                                        # the dissolution: calcite carries −1
e, k = 1:7, 8:8                                         # the two partitions
show(b) = join(("$(symbol(p)) $(x)" for (p, x) in zip(SMc.primaries, b)), ",  ")
println("A_e ν_eᵀ  = ", show(SMc.A[:, e] * ν[e]))      # per mole dissolved
println("−A_k ν_kᵀ = ", show(-SMc.A[:, k] * ν[k]))
```

The budget moves by the calcite's column, the two expressions are one, and no
amount of any equilibrium species is touched: the minimization decides where the
released calcium and carbon go. Advancing the amounts themselves instead would
remove two protons per mole from a solution that holds almost none. In a cement
paste at pH 13 a liter holds ``10^{-13}`` mol of them:

```@example pe_matrix
n_H = 1e-13                                             # mol of H⁺ in the solution
ξ = 1e-3                                                # mol of calcite dissolved
println("H⁺ advanced along the reaction: ", round(n_H + ν[2] * ξ; sigdigits = 3), " mol")
```

The coordinate of the budget on the proton may well be negative, it is a
coordinate; the amounts the minimization returns for it are not, the protons
"consumed" being taken from the hydroxide and the water, as the equilibrium
says.

Three properties of the construction follow, each of which an implementation can
get wrong while conserving matter:

- **``\mathbf{A}_e`` is the matrix the minimization is posed on.** The formula
  matrix over the elements and the matrix over the primary species of the system
  are different matrices, the same conservation in two bases (§1); a budget
  built in one basis and a minimization posed in the other is infeasible at
  every step. The equilibrium sub-system inherits
  the primary species of the parent system, so that ``\mathbf{b}_e``, its
  derivative and the minimization share one basis.
- **The minimization runs over the equilibrium partition only.** Posed on the
  whole system it would equilibrate the kinetic minerals instantaneously, which
  is what a kinetic description exists to prevent. The sub-system keeps the
  phases of the parent whose members are all in the partition, its solid
  solutions and its site families, on an all-or-nothing rule: a family split
  between the two partitions is refused, its members sharing one budget.
- **Surface sites belong to the equilibrium partition.** A site family's budget
  is a conservation row of the minimization; moving one member across would take
  the row with it. The states of a site are taken to redistribute as fast as the
  aqueous speciation; an adsorption slow enough to need a rate law is another
  model, declared by naming its species as kinetic.

## [3. The state and the right-hand side](@id sec-theory-pe-rhs)

The state is

```math
\mathbf{u} = (\mathbf{b}_e,\ \mathbf{n}_k,\ \boldsymbol{\xi}),
```

``\boldsymbol{\xi}`` being the extents of the kinetic reactions, and it evolves
according to [Leal2015; Eqs. 2.25–2.26](@cite)

```math
\frac{\mathrm{d}\mathbf{n}_k}{\mathrm{d}t} = \boldsymbol{\nu}_k^\mathsf{T}\mathbf{r},
\qquad
\frac{\mathrm{d}\mathbf{b}_e}{\mathrm{d}t} = \mathbf{A}_e\,\boldsymbol{\nu}_e^\mathsf{T}\mathbf{r},
\qquad
\frac{\mathrm{d}\boldsymbol{\xi}}{\mathrm{d}t} = \mathbf{r},
\qquad
\mathbf{r} = \mathbf{r}(\mathbf{n}_e, \mathbf{n}_k, T, t),
```

the amounts of the equilibrium partition being those of the minimum of the Gibbs
energy at the budget and the temperature of the state,

```math
\mathbf{n}_e = \varphi(\mathbf{b}_e, T)
= \arg\min_{\mathbf{n}\,\ge\,0}\ G(\mathbf{n}; T, P)
\quad\text{subject to}\quad \mathbf{A}_e\,\mathbf{n} = \mathbf{b}_e .
```

The right-hand side ``\mathbf{f}(\mathbf{u}, t)`` is thereby a function of the
state, since ``\varphi`` is: the same state always gives the same derivative,
whichever point of the trajectory the integrator evaluates it at, a stage, a
rejected step or a Jacobian probe. The initial state is equilibrated before the
first step, so that the trajectory starts on the constraint
``\mathbf{n}_e = \varphi(\mathbf{b}_e, T)`` rather than drifting onto it.

### [When the partition may be frozen within a step](@id sec-theory-pe-splitting)

When no rate law reads the equilibrium partition, as for a Parrott–Killoh or a
Waller law, ``\mathbf{f}`` does not depend on ``\varphi`` at all. The
minimization then changes nothing in the trajectory, and solving it once per
accepted step, only to report it, is exact: the splitting is exact
(`speciation = :frozen`). [`speciated_states`](@ref) certifies the partition at
any instant afterwards.

When a rate law reads it, through an activity, a saturation ratio or the amount
of an equilibrium species, the partition has to be solved where the integrator
evaluates ``\mathbf{f}`` (`speciation = :rhs`). Frozen within a step instead, the
rate is constant over the step, and the stiff method integrates the extent
explicitly however implicit it is. A rate ``r = k(1-\Omega)`` relaxes ``\Omega``
to one over a time of order ``1/k`` divided by the amount it acts on; a step
longer than that time overshoots the equilibrium, ``\Omega`` then exceeds one by
orders of magnitude, the rate reverses, and the run precipitates the mineral
from nothing while the integrator, whose error estimate sees a constant rate,
reports no error. Bounding the step does not cure it, since the bound would have
to be the relaxation time of the fastest such law. `integrate` tells the two
cases apart from the rate laws themselves (`speciation = :auto`), and a
trajectory that reaches amounts the system cannot hold is returned as a failure
rather than a success.

A semi-adiabatic calorimeter takes the second route whatever its laws read: its
temperature is solved with the partition (Section 6).

## [4. The Jacobian](@id sec-theory-pe-jacobian)

A stiff method needs ``\mathbf{J} = \partial\mathbf{f}/\partial\mathbf{u}``.
Applying the chain rule to the rates gives [Leal2015; Eqs. 2.34–2.35](@cite)

```math
\frac{\partial\mathbf{r}}{\partial\mathbf{u}}
= \frac{\partial\mathbf{r}}{\partial\mathbf{n}_e}\,\frac{\partial\varphi}{\partial\mathbf{b}_e}\,\frac{\partial\mathbf{b}_e}{\partial\mathbf{u}}
+ \frac{\partial\mathbf{r}}{\partial\mathbf{n}_k}\,\frac{\partial\mathbf{n}_k}{\partial\mathbf{u}},
```

in which only ``\partial\varphi/\partial\mathbf{b}_e`` is not explicit. It
follows from the optimality conditions of the minimization. Writing them
``\mathbf{F}(\mathbf{n}_e, \mathbf{y}; \mathbf{b}_e, T) = \mathbf{0}``, the
stationarity of the Lagrangian of the species present and the conservation
``\mathbf{A}_e\mathbf{n}_e = \mathbf{b}_e``, with ``\mathbf{y}`` the multipliers
of the conservation rows, and differentiating them at the answer with its active
set held yields

```math
\frac{\partial\mathbf{F}}{\partial(\mathbf{n}_e, \mathbf{y})}
\begin{bmatrix} \partial\mathbf{n}_e/\partial\mathbf{b}_e \\ \partial\mathbf{y}/\partial\mathbf{b}_e \end{bmatrix}
= -\frac{\partial\mathbf{F}}{\partial\mathbf{b}_e},
```

and the same with ``T`` in place of ``\mathbf{b}_e``. The package does not form
these matrices column by column. The minimization is solved on plain numbers,
and when the state carries the dual numbers of automatic differentiation the
answer is lifted by one Newton correction with the Jacobian of ``\mathbf{F}``
at the answer,

```math
\mathbf{n}_e^\star \;\leftarrow\; \mathbf{n}_e^\star
- \left[\frac{\partial\mathbf{F}}{\partial(\mathbf{n}_e,\mathbf{y})}\right]^{-1}
\mathbf{F}(\mathbf{n}_e^\star, \mathbf{y}^\star; \mathbf{b}_e, T),
```

whose value part vanishes and whose dual part is exactly the derivative above,
in every direction the dual numbers carry. The Jacobian of ``\mathbf{f}`` is thus
exact to first order, which is all a Rosenbrock method or the Newton iteration of
a BDF method requires, and it holds the derivative of the partition rather than
the zero a partition frozen within a step would contribute.

## [5. The implicit step](@id sec-theory-implicit-step)

The differential system of Sections 3 and 4 lets the integrator choose the step.
The implicit step of [Leal2017](@citet) makes the step itself one minimization,
in which the extents are unknowns beside the amounts:

```math
\min_{\mathbf{n} \ge 0,\ \Delta\boldsymbol{\xi}}\ G(\mathbf{n})
\quad\text{subject to}\quad
\mathbf{A}\,\mathbf{n} = \mathbf{b}_0,
\qquad
\mathbf{K}^\mathsf{T}\mathbf{n} - \Delta\boldsymbol{\xi} = \boldsymbol{\xi}_0,
\qquad
\Delta\boldsymbol{\xi} - \Delta t\,\mathbf{M}\,\mathbf{r}(\mathbf{n}) = \mathbf{0},
\qquad
\mathbf{M} = \mathbf{K}^\mathsf{T}\mathbf{K}.
```

``\mathbf{K}`` holds the stoichiometry of each reaction restricted to its
non-aqueous participants: an aqueous product re-speciates, and holding it would
stop holding the mineral. On calcite ``\mathbf{K} = [-1]``, so that
``\mathbf{M} = 1`` and ``\Delta\xi = \Delta t\, r``. The rate is evaluated at the
composition of the end of the step, which makes the step a backward Euler step of
the extents, with no frozen speciation and no lag between the kinetic and the
equilibrium species. The reactivity rows are linear in the unknowns and join the
conservation block, so that the algebraic cost of the kinetics is the number of
reactions, not of species.

Backward Euler makes the step unconditionally stable. Under ``r = k(1-\Omega)``
the end-of-step rate vanishes as the solution reaches saturation, and the step
approaches ``\Omega = 1`` from below at any ``\Delta t``, where an explicit step of
``\Delta t\,k`` would dissolve far more mineral than is present. Its accuracy is
that of a first-order method: a step ten times too long is ten times less
accurate.

**Two constraints on the products.** One reactivity row per reaction
(`coupling = :reactions`) leaves the solids free to rearrange, the assemblage
being whatever minimizes ``G``. One row per kinetic species
(`coupling = :species`), ``\mathbf{n}_i = \mathbf{n}_i(0) + \sum_j \nu_{ij}\Delta\xi_j``,
imposes the products of a stoichiometric scheme while the solution still
minimizes ``G`` under element conservation. These species are eliminated rather
than constrained: their amounts are an affine function of the extents, so they
are removed from the unknowns, their element content subtracted from the budget,
and an equilibrium over what remains is solved inside a Newton iteration on the
few extents. Holding them by a linear row instead keeps their stationarity in the
system, satisfied only by a multiplier that must reach the mineral's own
chemical potential, of order ``10^2``–``10^3`` in ``RT``, and the Newton
iteration stalls at a fixed point that its line search cannot leave.

**Choosing the step, and why the certificate decides.** For a rate that
vanishes at equilibrium the step has two solutions: the equilibrium root, and
the composition with the mineral wholly dissolved, which satisfies the element
balance and the reactivity row ``\mathbf{K}^\mathsf{T}\mathbf{n} - \Delta\boldsymbol{\xi}``
while violating ``\Delta\boldsymbol{\xi} - \Delta t\,\mathbf{M}\,\mathbf{r}(\mathbf{n})``
by the whole extent. Which root a Newton iteration reaches from far away depends
on its path. [`kinetic_step_adaptive`](@ref) estimates the error of a step of
``\Delta t`` from two steps of ``\Delta t/2`` (Richardson's estimate, which for a
first-order method is the error of the coarse step) and accepts a step on the
estimate **and** on the certificate of the augmented problem. The estimate alone
is blind to an error both resolutions share: the coarse step and both half-steps
can each dissolve the whole mineral and agree to every digit. The tolerance is
relative to the amount each reaction acts on, not to the extent: since
``\Delta\xi \propto \Delta t``, a tolerance on the extent vanishes with the step
while the noise of the minimization does not, and the controller would halve the
step without end.

**What the certificate of a step proves.** It is taken on the augmented problem,
whose reactivity rows carry the multipliers that make the stationarity of a
kinetically held mineral satisfiable. The same composition tested against the
unconstrained equilibrium fails, because a mineral held back by a rate law is
supersaturated by construction, which is what being held back means. Two choices
of the step are left to the certificate rather than to the user because neither
works on every problem: whether the starting guess of a system with solid
solutions is equilibrated first (a mixing phase is admitted by a tangent-plane
test, which a composition without the phase cannot pass), and whether the
kinetic minerals are held in the active set (a mineral exhausted at equilibrium
is otherwise removed by the active-set rule, after which nothing enforces its
reactivity row). The problem being convex, the certificate decides between the
candidates exactly.

## [6. The calorimeters](@id sec-theory-pe-calorimeters)

The enthalpy of the paste is the sum of the standard enthalpies of its species,

```math
H(\mathbf{n}, T) = \sum_i n_i\, h_i^\circ(T),
```

the excess enthalpies of the activity models being left out, as
[What the package counts as heat](@ref sec-theory-heat-output) states. Under
partial equilibrium the kinetic reactions only release ions, the hydrates being
precipitated by the minimization, so that the heat of those reactions is not the
heat of the paste; the heat is the change of ``H`` over the whole composition,
``\mathbf{n}_e = \varphi(\mathbf{b}_e, T)`` included. The reference is the
partition at the start, ``H_0 = H(\varphi(\mathbf{b}_e(0), T_0), \mathbf{n}_k(0), T_0)``,
so that the first equilibrium releases nothing.

### The isothermal cell

The bath holds the paste at ``T_0`` and takes the heat ``Q`` the paste releases,
so that the enthalpy of the paste and of what the bath has received does not
change, ``H + Q = H_0``. The heat is therefore a function of the state,

```math
Q(\mathbf{u}) = H_0 - H\bigl(\varphi(\mathbf{b}_e, T_0), \mathbf{n}_k, T_0\bigr),
```

and requires no integration. Its rate along the trajectory is the derivative of
the same function,

```math
\dot q = -\frac{\partial H}{\partial\mathbf{n}_e}\,\frac{\partial\varphi}{\partial\mathbf{b}_e}\,\frac{\mathrm{d}\mathbf{b}_e}{\mathrm{d}t}
- \frac{\partial H}{\partial\mathbf{n}_k}\,\frac{\mathrm{d}\mathbf{n}_k}{\mathrm{d}t},
```

evaluated by lifting the partition along ``\mathrm{d}\mathbf{u}/\mathrm{d}t``
([`heat_flow`](@ref)), and the certified replay
([`heat_release`](@ref)) computes both from the same partition.

### The semi-adiabatic cell

The vessel and what it holds besides the paste, of heat capacity ``C_v``, take
the temperature of the paste, and heat leaves through the walls at the rate
``\mathcal{L}(T - T_{\rm env})`` ([`SemiAdiabaticCalorimeter`](@ref)). The
enthalpy of the cell, ``H + C_v(T - T_0)``, changes by the losses alone, so that
the state is extended by its change since the start,

```math
\mathbf{u} = (\mathbf{b}_e,\ \mathbf{n}_k,\ \boldsymbol{\xi},\ \Delta H),
\qquad
\frac{\mathrm{d}\,\Delta H}{\mathrm{d}t} = -\mathcal{L}(T - T_{\rm env}),
```

and the temperature is the root of the energy balance

```math
\mathcal{R}(\mathbf{u}, T) = H\bigl(\varphi(\mathbf{b}_e, T), \mathbf{n}_k, T\bigr) - H_0 + C_v\,(T - T_0) - \Delta H = 0 .
```

The temperature is solved for, not integrated. Its derivative with respect to
``T`` at fixed state,

```math
C_{\rm eq} = \frac{\partial\mathcal{R}}{\partial T}
= C_v + \sum_i n_i\,\frac{\mathrm{d}h_i^\circ}{\mathrm{d}T}
+ \sum_{i \in e} h_i^\circ\,\frac{\partial\varphi_i}{\partial T},
```

is the heat capacity of the cell at equilibrium, that of the vessel and of the
paste at fixed composition, plus the heat the partition takes up as it shifts
with the temperature. It is positive wherever the balance has a unique root, and
the temperature is found by Newton's method,
``T \leftarrow T - \mathcal{R}(\mathbf{u}, T)/C_{\rm eq}``, from the temperature
of the last accepted step, each iteration solving the minimization at the
current ``T`` and lifting it in temperature to obtain ``C_{\rm eq}``; a balance
whose ``C_{\rm eq}`` is not positive makes the evaluation fail, and the
integrator rejects the step.

The implicit-function theorem applied to ``\mathcal{R}`` then gives the
derivatives of the temperature the Jacobian needs,

```math
\frac{\partial T}{\partial\mathbf{u}} = -\frac{1}{C_{\rm eq}}\,\frac{\partial\mathcal{R}}{\partial\mathbf{u}},
\qquad
\frac{\partial\mathcal{R}}{\partial\mathbf{b}_e} = \sum_{i \in e} h_i^\circ\,\frac{\partial\varphi_i}{\partial\mathbf{b}_e},
\quad
\frac{\partial\mathcal{R}}{\partial\mathbf{n}_k} = \mathbf{h}_k^\circ,
\quad
\frac{\partial\mathcal{R}}{\partial\,\Delta H} = -1,
```

and the package obtains them, as in Section 4, by a single correction lifted
into the dual numbers of the state, ``T \leftarrow T^\star - \bigl(\mathcal{R}(\mathbf{u},
T^\star) - \mathcal{R}^\star\bigr)/C_{\rm eq}``, which keeps the value ``T^\star``
at which the partition was solved. The partition at the state is then lifted in
``\mathbf{b}_e`` and in ``T`` together, so that every entry of the Jacobian holds
the first derivatives of ``\varphi`` and of the root, and nothing else.

That is the reason for carrying ``\Delta H`` rather than ``T``. Written on the
temperature, the balance would read
``C_{\rm eq}\,\mathrm{d}T/\mathrm{d}t = -\partial\mathcal{R}/\partial\mathbf{b}_e \cdot
\mathrm{d}\mathbf{b}_e/\mathrm{d}t - \dots - \mathcal{L}``, a right-hand side that
holds ``\partial\varphi/\partial\mathbf{b}_e`` and ``\partial\varphi/\partial T``
themselves, so that its Jacobian would need the second derivatives of the
minimization. Following the derivatives of a partition frozen at the last
accepted step instead, and correcting the temperature by the heat of the next
re-speciation, makes the right-hand side depend on the history of the
integration and puts part of the energy balance outside its error control.

The heat and the temperature the accessors report are functions of the state:
[`temperature_profile`](@ref) solves the root at each instant,
[`cumulative_heat`](@ref) returns ``Q = C_v(T - T_0) - \Delta H``, the enthalpy
the paste has lost, and [`heat_flow`](@ref) its rate,
``\dot q = C_v\,\mathrm{d}T/\mathrm{d}t + \mathcal{L}(T - T_{\rm env})``.

### Without an equilibrium partition

When the kinetic reactions produce the hydrates themselves and no minimization is
solved, every amount follows from the extents,
``\mathbf{n} = \mathbf{n}(0) + \boldsymbol{\nu}^\mathsf{T}\boldsymbol{\xi}``, independently of the temperature, and
the same balance reduces to an equation the state can carry directly,

```math
\Bigl(C_v + \sum_i n_i\,c_{p,i}^\circ(T)\Bigr)\frac{\mathrm{d}T}{\mathrm{d}t}
= \sum_j r_j\,\bigl(-\Delta_r H_j^\circ(T)\bigr) - \mathcal{L}(T - T_{\rm env}),
```

whose right-hand side holds no derivative of a minimization; this formulation
integrates the temperature, and an isothermal cell the heat
``\mathrm{d}Q/\mathrm{d}t = \sum_j r_j(-\Delta_r H_j^\circ)``. The reaction
enthalpy ``\Delta_r H_j^\circ`` is computed from the enthalpies of formation of the
participants, so that a reaction must be balanced for its heat to mean anything.

## [7. What is checked](@id sec-theory-pe-checks)

The formulation is held to identities that a route sharing nothing with it
computes:

- **The partition carries the integrated budget.** At every accepted step
  ``\|\mathbf{A}_e\mathbf{n}_e - \mathbf{b}_e\|_\infty`` is recorded, and so is
  the number of re-speciations that failed; `integrate` reports both. A run in
  which the minimization never succeeded would otherwise be indistinguishable
  from a healthy one.
- **The coupling reaches the equilibrium it claims.** A solid releasing calcium
  at a constant rate onto a sorbent (`test/kinetics/test_surface_coupling.jl`)
  holds the site budget at every instant, accounts for every mole released, and
  satisfies the mass action of the binding reaction on the integrated state; the
  equilibrium and the coupling along a constant-rate dissolution agree with
  Reaktoro ([Validation against Reaktoro](@ref)).
- **The Jacobian holds the derivative of the partition.** Its columns in
  ``\mathbf{b}_e`` equal ``(\partial\mathbf{r}/\partial\mathbf{n}_e)
  (\partial\varphi/\partial\mathbf{b}_e)`` computed through the certified
  solver (`test/kinetics/test_rhs_speciation.jl`).
- **The calorimeters conserve energy** (`test/kinetics/test_calorimetry.jl`):
  the enthalpy of an adiabatic cell, ``H + C_v(T - T_0)`` evaluated on the
  certified replay at the temperature the run solved, stays ``H_0`` to within
  ``10^{-7}`` of the heat released; ``1/(\partial T/\partial\,\Delta H)``, taken
  through the lifted root, equals ``C_{\rm eq}`` computed with the certified
  equilibrium differentiated with respect to its temperature, and
  ``\partial T/\partial\mathbf{b}_e`` the same with respect to the budget; the
  heat of an isothermal cell equals the enthalpy difference of the certified
  replay, and its rate the rate of the replay; a Rosenbrock, a BDF and an
  explicit Runge–Kutta method integrate the same temperature, to ``10^{-5}`` K,
  with losses through the walls equal to ``-\Delta H``.

## [8. Limits](@id sec-theory-pe-limits)

- The excess enthalpies of the activity models are not in ``H``, whereas the
  minimization carries the temperature dependence of the activity coefficients;
  ``C_{\rm eq}`` then holds a term that the excess enthalpies would have made a
  quadratic form, and its positivity is not guaranteed by the stability of the
  equilibrium alone. A balance with a non-positive ``C_{\rm eq}`` is refused at
  the evaluation rather than corrected.
- ``\varphi`` is smooth only while the assemblage does not change. Where a phase
  appears or vanishes it has a kink, which the integrator meets as a loss of
  accuracy and answers by shortening its step; no event is located there.
- The partition is the answer of the certified solver to its tolerance, and the
  temperature is solved to ``10^{-9}`` K.
- The implicit step is first order in time; its accuracy is controlled by the
  adaptive march, not by the step alone.

## See also

- [Coupling kinetics and equilibrium](@ref sec-coupling), the three routes run
  on one mineral, with their trajectories.
- [The silicates of a CEM I clinker, hydrating end to end](@ref sec-coupled-hydration),
  the partition on a clinker.
- [A CEM I 52.5 N mortar in a semi-adiabatic calorimeter, inside the kinetics](@ref ex-semiadiabatic),
  the semi-adiabatic cell on a cement.
- [Writing a kinetic model](@ref sec-kinetics-syntax), the syntax of the
  routes and of the calorimeters.
