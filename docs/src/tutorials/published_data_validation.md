# Validation against published data

[Validation against Reaktoro](@ref) compares this package against another code
reading the same database. That comparison cannot see an error in the database
itself: both codes would make it together, agree perfectly, and be wrong.

This page is the other half. It checks the **vendored data and the package
against the papers the data came from** — first the Cemdata18 article
([Lothenbach2019](@cite)), which publishes its solubility products *and* the
standard properties they were derived from, so the two can be required to agree;
then a set of measured solution compositions ([Atkins1992](@cite)), which can
disagree with the database and does.

Every number below is pinned by an assertion in `test/cemdata18_reference.jl`,
`test/lothenbach2010_reference.jl`, `test/lothenbach2008_reference.jl`,
`test/atkins1992_reference.jl`,
`test/limestone_blending_reference.jl`, `test/chloride_binding_reference.jl`,
`test/hong_glasser1999_reference.jl` or `test/duan2016_reference.jl`, so it is
checked on every CI run — **pinned, not merely bounded**.

The distinction is what makes a page like this worth reading. A threshold that
says a disagreement is *under* `0.05` leaves the `0.041` printed beside it free
to drift to `0.047` while the suite stays green and the page goes quietly
stale. Every displayed value is therefore asserted at the precision it is
printed to, which is half of its last digit.

!!! note "The one class of number that is not pinned"
    **The timings.** The seconds quoted for a cold and a warm solve measure one
    machine. Asserting them would make the suite fail on a faster one without
    anything being wrong, which is the opposite of what a test is for. Every
    rung of every sweep in the tables below is computed and pinned, including
    those that only illustrate a trend.

## What has been checked

Nine sources, more than eight hundred assertions, and the coverage is uneven on purpose:
the database is checked exhaustively because it is cheap to check exhaustively,
while the equilibrium cases are checked one composition at a time because each
one costs seconds to minutes.

| what | against | how closely | where it stops |
|:--|:--|:--|:--|
| **every solubility product** in the shipped CEMDATA18 — 52 phases | [Lothenbach2019](@cite) Tables 2-3 | 50 close to `0.041`, 48 to `0.005` | two M-S-H end members are `0.48` and `0.40` out; the source disagrees with itself |
| **standard properties and HKF coefficients**, 19 aqueous species and 7 gases | [Lothenbach2019](@cite) Tables D.1-D.2 | exact, all seven coefficients each | nitrite-AFm and Fe-Friedel's salt are not in the file |
| **`ΔₐG⁰` rebuilt from `ΔfH°` and `S°`** | the file's own `ΔfG°` | exact for 220 of 228 substances | the eight exceptions are all phases whose `S°` Cemdata18 estimated |
| **HKF away from 298.15 K, 1 bar** | [Duan2016](@cite) Table 4, HKF column | `0.03 %` at the reference point; `0.3 %` across a factor 2.9 in pressure | their two constants in one row sit at two different pressures |
| **calcite `log Ksp`** at 25 °C | the accepted value, [PlummerBusenberg1982](@cite) | `−8.480` against `−8.48` | diverges from Duan's non-HKF method by 1.4 log units at 478 K |
| **a measured solution** over a two-phase assemblage | [Atkins1992](@cite) Table 2 | Al and pH agree, robustly | Si is `×5` and Ca has a floor the model cannot leave; 9 of their 10 mixtures are not usable at all |
| **the carboaluminate sequence** under limestone | [Kulik2021](@cite) Fig. 7A, [Lothenbach2019](@cite) Figs. 13-14 | order and thresholds reproduce; iron partition exact | only appears when Al exceeds Fe — see below |
| **chloride binding** and the AFm → Friedel transition | [Guo2018](@cite) Fig. 1(b) | plateau to `1 %`, both conservation laws close | pH not reproducible (their alkalis are unpublished); above 2 % NaCl the activity model is out of range |
| **the previous generation**, Cemdata07's 29 solubility products | [Lothenbach2010](@cite) Table 1 | 11 unchanged to the printed digit, from the package's own energies | 9 revised, by up to 1.6 log units; 9 have no CEMDATA18 solid of the same composition |
| **the molar volumes of 28 solids** | [Lothenbach2008](@cite) Table 4, cemdata2007 | 24 within half the printed digit | the other 4 by up to `1.7 cm³/mol`, three of them iron phases |
| **the formation energies of the same solids** | [Lothenbach2008](@cite) Table 4 | 21 solids and water the same to `0.01 kJ/mol` | 7 revised by up to `10 kJ/mol`; another hydration state costs the water's energy, within 3 % |
| **an aged pore solution** against portlandite and ettringite | [Lothenbach2010](@cite) Table 2 | within a quarter of a log unit of saturation | which side depends on the activity model at `I ≈ 0.5 mol/kg` |
| **alkali uptake by C-S-H**, 4 Ca/Si × 6 concentrations × Na and K | [HongGlasser1999](@cite) Tables 1-2 | pH to `0.072` from 15 to 100 mM below Ca/Si 1.8 | alkali over-bound at 46 of 48 points; the end members were fitted to these data |

Two things this chapter deliberately does **not** do. It does not check the
kinetics — [Validation against Reaktoro](@ref) covers the coupling, and no
published hydration curve is reproduced here. And it does not check anything
above about `1 mol/kg` ionic strength, because `HKFActivityModel`'s B-dot term
was not fitted there.

## Traps

Seven things that cost time to find, collected because each is reusable well
beyond the case that established it.

**A solid entered as oxides loses its water.** Feeding a C-S-H as lime and
silica gives the solver its calcium and its silicon and none of the 2.1 H₂O per
formula unit it carries. On [Guo2018](@cite)'s inventory that is 44.5 g per liter
of concrete: the solver takes it back out of the pore solution, the pore volume
falls from 146 mL to 84, and every concentration is wrong by nearly a factor of
two. See [Chloride binding](@ref).

**Walk a sweep downhill.** A warm start survives only while the phase
*assemblage* holds; at a point where a phase appears or vanishes the certificate
refuses it and the full multi-start cascade runs again. Starting from the
composition-rich end and walking toward the simple one put one limestone sweep
from 152 s to 52 s — and removing a single intermediate rung, so that one step
crossed the boundary in a jump, put it back to two and a half minutes. See
[Blending a CEM I with limestone](@ref).

**Suppress gibbsite below 60 °C.** Cemdata18 §2.1 says so, and it is not a
harmless over-specification: left in the phase list it takes over from
microcrystalline `AlOHmic` in the sulfate-poor mixtures and moves the aluminum.
See [Atkins et al. (1992): a measurement, not a calculation](@ref).

**The iron decides whether carboaluminates form at all.** `C3AFS0.84H4.32` takes
one aluminum per iron, so the siliceous hydrogarnet is capped by the iron
available. On a clinker with `Fe₂O₃ = 4.49 %` there is enough iron to pair with
every aluminum, no monocarbonate forms at any limestone content, and the
published sequence simply does not appear. Below about 3.5 % it does. See
[The sequence, and the condition nobody states](@ref).

**Temperature belongs to the state, not to the solve.** The solves take `T` and
`P` from the `ChemicalState`, and `equilibrate`, `equilibrate_certified` and
`EquilibriumSolver` refuse a `T` or a `P` passed to them, with an error that
names `set_temperature!`. They refuse it because the mistake is otherwise silent
and not small: a keyword the solve does not know goes to the optimizer's options
and is dropped there, and at a fixed hydroxide the pH of a pore solution moves
by the change in `pKw` between the two temperatures, `0.17` from 20 to 25 °C —
larger than most of the agreements this chapter reports. A measurement made at
20 °C is compared at `ChemicalState(cs; T = 293.15u"K")`. See
[The HKF model away from 25 °C and 1 bar](@ref).

**A surface species' standard energy is on the aqueous gauge, not its own.**
The site occupancies are driven by `ΔₐG⁰` like everything else, so a complex
formed from an aqueous ion has to carry that ion's formation energy:

```
≡SOH + M^z ⇌ ≡SOM^(z−1) + H⁺     ⟹     G(≡SOM) = G(M^z) − RT ln K
```

with `G(≡SOH) = G(H⁺) = 0` fixing the gauge. Setting `G(≡SOM) = −RT ln K` and
dropping `G(M^z)` is off by 553 kJ/mol for calcium — and it does not produce a
wrong answer, it produces `0.00000` mol of every cation complex while the solve
still certifies. Deprotonation keeps working throughout, because `≡SO⁻` is the
one complex whose reaction involves no aqueous species but H⁺, so the omission
is invisible on the row most likely to be checked first. `test/surface_complexation.jl`
carries the same warning at its zinc case.

`G(H⁺) = 0` is the database's own convention. `G(≡SOH) = 0` is a gauge **only
while the site budget is fixed**: both sides of every surface reaction carry a
site, so the choice cancels, and `test/surface_complexation.jl` checks that the
answer does not move over a 20 kJ/mol shift of the whole family. Under
`SITES_FOLLOW_HOST` it stops canceling — the host carries `−ν` of the site
component, so the choice reaches the host's own solubility — and the free site's
reference is then not a convention but the energy of the matter it carries,
`μ°(H₂O) − μ°(H⁺) = −237.2 kJ/mol` for an oxide. [`host_coupling_bias`](@ref)
computes it. See [The surface half of Guo](@ref).

**A rebuilt Gibbs energy is not always the tabulated one.** `ΔₐG⁰(T)` is formed
from `ΔfH°` and `S°`, so at 298.15 K it must reproduce `ΔfG°` — and does, for
220 of the 228 substances. The eight that do not are exactly those whose entropy
and heat capacity Cemdata18 estimated rather than measured. If a result turns on
one of them, check which number you are standing on. See
[The same two phases disagree with themselves](@ref).

## Cemdata18 Tables 2 and 3: the solubility products

Table 2 gives `log Ks0` for each phase **and** writes out the dissolution
reaction it refers to; Table 1 gives the `ΔfG°` those constants were computed
from. Recomputing one from the other through Cemdata18's own aqueous Gibbs
energies closes a loop that can only close if the transcription is exact and the
two halves of the paper share a reference state.

### Transcribing a reaction without transcribing it twice

Table 2 writes its reactions over `Al(OH)₄⁻`, `Fe(OH)₄⁻` and `SiO(OH)₃⁻`.
Cemdata18's GEMS primaries are `AlO2-`, `FeO2-` and `HSiO3-`, which differ from
those by 2, 2 and 1 H₂O. Transcribing both the substitution and the paper's
water coefficients would be two chances to get it wrong, and the second would
silently absorb the first.

So only the **non-water** products are transcribed, and the water is recovered
from the hydrogen balance. Oxygen and charge then become free checks — they can
only close if the transcribed products are right. They close to `10⁻⁹` on all 52
rows.

### Result

**50 of the 52 rows close.** The largest disagreement among them is `0.041` on
`M8A-OH-LDH`, whose published value is quoted to one decimal; 48 are inside
`0.005`.

Two do not:

| phase | from tabulated `ΔfG°` | published `log Ks0` | Δ | in energy |
|:--|--:|--:|--:|--:|
| `M075SH` (M₁.₅S₂H₂.₅) | −28.317 | −28.80 | **+0.483** | 2.76 kJ/mol |
| `M15SH` (M₁.₅SH₂.₅) | −23.175 | −23.57 | **+0.395** | 2.26 kJ/mol |

Both are M-S-H end members. With every other phase in the table closing to
better than `0.05`, this is a property of the source and not of the
transcription.

### The same two phases disagree with themselves

`ΔₐG⁰(T)` is formed from `ΔfH°` and `S°` rather than read off the file, so at
`T = 298.15 K` it must reproduce the tabulated `ΔfG°`. For **220 of the 228**
substances it does, to the last bit. For eight it does not:

| substances | gap on `ΔfG°` | implied inconsistency in `S°` |
|:--|--:|--:|
| `ECSH1-KSH`, `ECSH2-KSH` | −243.7 J/mol | 0.82 J/K/mol |
| `ECSH1-NaSH`, `ECSH2-NaSH` | −207.6 J/mol | 0.70 J/K/mol |
| `KSiOH`, `NaSiOH` | −208.6 J/mol | 0.70 J/K/mol |
| `M15SH` | −1088.6 J/mol | 3.65 J/K/mol |
| `M075SH` | −1364.8 J/mol | 4.58 J/K/mol |

The eight are not a random selection, and — this was the first guess and it is
wrong — they are **not** the ones missing a heat-capacity block: six of them
have one. What they share is that Cemdata18 says their entropy and heat capacity
were **estimated rather than measured**: the six alkali C-S-H end members by the
linear Ca/Si relations of Table 4 (the paper's Eqs 2a and 2b), the two M-S-H
end members from talc, chrysotile and water (Table 1, footnote r, after
[Nied2016](@cite)).

An estimated `S°` that was never reconciled with the tabulated `ΔfG°` leaves the
triplet `(ΔfG°, ΔfH°, S°)` inconsistent, and rebuilding the third from the other
two is what makes it visible. The control is the five zeolites: they carry no
heat-capacity block at all and still land on their tabulated value exactly. So
this is the data, not the code path.

For `M075SH` that makes **three** mutually inconsistent numbers for one phase in
one paper — the tabulated `ΔfG°`, the value its own `ΔfH°` and `S°` imply, and
the value its published `log Ks0` implies — spread over 2.8 kJ/mol. Nothing in a
calculation announces this; the phase simply sits half a log unit off wherever
M-S-H matters.

## Appendix D: what the solubility products cannot test

Table 2 exercises the `ΔfG°` of the aqueous primaries hard — a drift in any one
of them breaks dozens of reactions at once — but says nothing about `S°`, `Cp°`,
`V°` or the HKF equation-of-state coefficients, which are what carry the
database away from 25 °C and 1 bar. Those are checked directly against
Tables D.1 and D.2, on 19 aqueous species and all seven gases.

Two conventions of the printed tables have to be undone first.

**The HKF columns are scaled.** Table D.1 prints `a₁·10`, `a₂·10⁻²`, `a₄·10⁻⁴`,
`c₂·10⁻⁴` and `ω⁰·10⁻⁵`, in the calorimetric units the HKF papers use. For
Ca²⁺:

| | `a₁` | `a₂` | `a₃` | `a₄` | `c₁` | `c₂` | `ω⁰` |
|:--|--:|--:|--:|--:|--:|--:|--:|
| printed | −0.1947 | −7.2520 | 5.2966 | −2.4792 | 9.0000 | −2.5220 | 1.2366 |
| in the file | −0.01947 | −725.20 | 5.2966 | −24792 | 9.0 | −25220 | 123660 |

All seven match on all 19 species, which is what makes the reading of the
scaling more than a guess.

**Table D.1's volume column is mislabelled.** It is headed `V⁰ (J/bar)` and
carries cm³/mol: Ca²⁺ is listed at −18.44, and −18.44 J/bar would be
−184.4 cm³/mol. The file stores −1.8439 J/bar, which is the same −18.44 cm³/mol.

Table D.2's gas column, under the same header, really is J/bar — 2479 J/bar is
24.79 L/mol, the ideal-gas molar volume at 298.15 K and 1 bar. The two tables
share a header and not a unit.

## What the shipped file does not carry

Two rows of Table 2 cannot be checked, because the phase is in the printed table
and not in the vendored file. The test asserts their absence, so that a database
update which adds them turns this note red rather than leaving it stale.

| row | why not |
|:--|:--|
| nitrite-AFm | the solid `mononitrite` is present; `NO2-` is not, so the reaction cannot be written over the file's own primaries |
| Fe-Friedel's salt (`C4FCl2H10`, `log Ks0 = −28.62`) | absent altogether, as are amorphous and microcrystalline `Fe(OH)₃` |

## Atkins et al. (1992): a measurement, not a calculation

Their Table 2 reports **analyzed** solution compositions for slurries of
synthetic hydrates in CO₂-free water at 25 °C, two to four solids at a time.
Ten mixtures are tabulated. **One is usable quantitatively**, and the reasons the
others are not are worth knowing before anyone spends an afternoon on them.

### The one that is: ettringite + C-S-H

Nominal C-S-H at Ca/Si = 0.9, measured at six months (mmol/L):

| | measured | calculated | |
|:--|--:|--:|:--|
| Al | 0.136 | 0.149 | agreement |
| pH | 11.0 | 11.33 | agreement |
| SO₄ | 1.08 | 1.18 | close, but see below |
| Si | 0.076 | 0.396 | **×5.2** |
| Ca | 1.95 | 2.72 | **+40 %** |

**Aluminum is the strong result.** Across the sweep of the C-S-H Ca/Si below,
from 0.75 to 1.0, it stays between 0.142 and 0.157 mmol/L. That is the database
answering, not a fit — and it is where this package does better than Atkins' own
1992 model, which gave 0.398 on the same mixture.

**Silicon is over-predicted for a reason the paper gives.** Atkins' own model
gave 0.448 mmol/L on the same mixture, and they explain it: electron microscopy
showed the ettringite had taken up silicon, *"substituting on average for 40 %
of the available SO₄ sites"*, which *"would tend to lower aqueous Si
concentrations, and increase SO₄ levels as observed"*. Cemdata18 carries no
Si-bearing AFt end member either, so the same silicon has nowhere to go and
stays in solution.

**Calcium is a disagreement the model cannot absorb**, and the obvious escape
does not work. Atkins report that their C-S-H dissolved incongruently over the
test, which would have lowered its Ca/Si below the nominal 0.9 — and a lower
Ca/Si ought to mean less calcium in solution. It does not:

| C-S-H Ca/Si | Ca | Al | Si | SO₄ | pH |
|--:|--:|--:|--:|--:|--:|
| 0.75 | 3.46 | 0.157 | 0.504 | 2.51 | 11.00 |
| 0.80 | 3.06 | 0.154 | 0.465 | 1.92 | 11.13 |
| 0.85 | 2.83 | 0.152 | 0.429 | 1.50 | 11.24 |
| 0.90 | **2.72** | 0.149 | 0.396 | 1.18 | 11.33 |
| 1.00 | 2.72 | 0.142 | 0.339 | 0.76 | 11.47 |
| measured | 1.95 | 0.136 | 0.076 | 1.08 | 11.0 |

The sweep holds the lime and varies the silica, and every row is pinned to the
digits above. Calculated calcium is flat between Ca/Si 0.9 and 1.0 and rises
below it, the direction incongruent dissolution would take it; it never
approaches 1.95. The floor is a property of the CSHQ model, not of the mixture:
Atkins' own 1992 model, on the same nominal C-S-H, gave 2.04.

Sulfate lands close, but it is the one number here that moves freely with the
Ca/Si — 0.76 at 1.0, 2.51 at 0.75 — so the agreement is not evidence of much.
Atkins' own 1992 model gave 0.597; Cemdata18 is nearer, and that is the whole
claim.

!!! note "Gibbsite has to be suppressed, and it changes the answer"
    Cemdata18 §2.1 says its precipitation *"should be suppressed for
    calculations at ambient temperatures, where microcrystalline Al(OH)₃ will
    form instead"*. Leaving `Gbs` in the phase list is not a harmless
    over-specification: it takes over from AH₃ in the sulfate-poor mixtures and
    moves the aluminum.

### Why the other nine are not usable

| reason | evidence |
|:--|:--|
| the solid-to-water ratio is not reported, and these assemblages are **not invariant** | wherever a C-S-H buffers the calcium, a fourfold change in loading moves the answer by a factor of two: mixture 1 (AFt + CSH 1.7) gives Ca = 8.0 mmol/L at one loading and 16.0 at four times it, against 2.67 measured |
| four mixtures grew a phase Atkins identifies as **metastable**, so the measured solution is not an equilibrium one | AFm in mixture 3, C₄AH₁₃ in 8 and 21, siliceous hydrogarnet in 30 — the paper says so itself and uses them to argue the point |
| two are in **0.4 M NaOH** | past the range the B-dot term in [`HKFActivityModel`](@ref) was fitted for (`I ≲ 1 mol/kg`, and less for the water activity); a Pitzer model would be the honest instrument |

Their Table 3 — pore fluids of five-year-old OPC, 30 % BFS and 30 % FA pastes —
is a better target for a future test: alkali and pH are reported and the binder
compositions are given. It needs the alkali uptake of the CSHQ Na/K end members
and belongs with the CEM I/II examples.

## Blending a CEM I with limestone

[Kulik2021](@cite) Fig. 7A and [Lothenbach2019](@cite) Figs. 13-14 report the
same thing: what happens to the aluminate phases as limestone replaces clinker.
CemGEMS is the closest oracle this package has — same CEMDATA18, same
Gibbs-energy minimization, and the paper states that it and GEM-Selektor agree
exactly. What it does not give is the cement it ran on, so what is checked is
the part that does not depend on it: the order the phases appear in, where added
limestone stops reacting, and where the iron goes.

### The iron

Fig. 14D is the sharpest quantitative claim in either paper: with Cemdata18
*"close to 100 % of the iron is bound by the siliceous hydrogarnet solid
solution"*, against about 80 % sitting in hemi-/monocarbonate under Cemdata07.

Measured here: **100.0 %**, at every limestone content and every iron content
tried, from 0.5 to 4.49 % `Fe₂O₃`, to within `10⁻³`. That one reproduces without
qualification.

### The sequence, and the condition nobody states

Both papers report

```
monosulfate  →  hemicarbonate  →  monocarbonate + ettringite  →  free calcite
```

It does not appear on every CEM I. On the oxide analysis this repository uses
elsewhere — `Fe₂O₃ = 4.49 %` — **no carbonate AFm forms at any limestone
content**, the added calcite stays inert, and every aluminum ends in the mixed
Al-Fe siliceous hydrogarnet. That is precisely the behavior [Lothenbach2019](@cite)
§3.1 attributes to *Cemdata07* and says Cemdata18 cures.

It is not a defect, and the mechanism is stoichiometric. `C3AFS0.84H4.32` takes
**one aluminum per iron**, so the phase is capped by the iron available. Sweeping
`Fe₂O₃` at 4 g of limestone, the balance going to CaO:

| `Fe₂O₃` % | 4.49 | 3.50 | 2.50 | 1.50 | 1.00 | 0.50 |
|:--|--:|--:|--:|--:|--:|--:|
| siliceous hydrogarnet | 0.0487 | 0.0396 | 0.0283 | 0.0170 | 0.0113 | 0.0057 |
| monocarbonate | **0** | 0.0031 | 0.0083 | 0.0135 | 0.0161 | 0.0187 |

The test walks the six iron contents at 4 g of limestone, downhill from the
ferriferous clinker, and every column is pinned — the hydrogarnet amount, which
this section turns on, as well as the monocarbonate.

The hydrogarnet tracks the iron one for one. Where there is enough iron to pair
with every aluminum, nothing is left to make a carboaluminate; below about
3.5 % there is, and the published sequence returns.

At `Fe₂O₃ = 2.5 %` the figure comes back step by step (mol per 100 g of binder):

| g limestone | monosulfate | hemicarbonate | monocarbonate | ettringite | free calcite |
|--:|--:|--:|--:|--:|--:|
| 0.0 | 0.0090 | 0 | 0 | 0.0067 | 0 |
| 0.5 | **0** | 0.0028 | 0.0036 | 0.0097 | 0 |
| 1.0 | 0 | **0** | 0.0085 | 0.0096 | 0.0015 |
| 4.0 | 0 | 0 | 0.0083 | 0.0093 | 0.0317 |

Half a gram destroys the monosulfate outright and stabilizes the ettringite by
more than 40 % — the mechanism the whole figure is about. By one gram the
hemicarbonate is spent and calcite begins to survive undissolved; past that the
limestone is a filler.

!!! tip "Walk the sweep downhill"
    A cold solve of this system costs about 26 s; a warm one, started from a
    neighboring answer, 0.2-0.7 s. But a warm start only helps while the phase
    *assemblage* holds: where a phase appears or vanishes the certificate
    refuses the warm answer and the full multi-start cascade runs again.

    Ascending from 0 g, the monosulfate-to-carboaluminate switch at 0.5 g cost
    **125 s on its own**. Descending from 4 g, where the assemblage is simple and
    stable, the same sweep costs **52 s in total**. Same answers, three times
    faster — and dropping a single intermediate rung, so that one step crosses
    the boundary in one jump, puts it back to two and a half minutes.

## Chloride binding

[Guo2018](@cite) Fig. 1(b) is unusually easy to reproduce, because the paper
states its inventory outright instead of leaving it to be reconstructed: per
liter of concrete, C-S-H 225 g, CH 90 g, AFm 9 g, AFt 22.5 g, porosity 14.6 %.

**Their constants are not CEMDATA18's.** Their §2 says so — *"Cemdata2007 gives
`Kp` and `Δ_r G_T^0` for nearly all phases in cement hydrate"* — and the note
under their dissolution table points at Lothenbach, Matschei, Möschner and
Glasser (2008), which is Cemdata07. The Cemdata07 table of [Lothenbach2010](@cite)
confirms it: its ettringite, monosulfate and jennite-type C-S-H are the `−44.9`,
`−29.26` and `−13.17` that Guo print as `−44.9085`, `−29.2628` and `−13.1659`
(see [Cemdata07, and an aged pore solution](@ref)). So this is a comparison
**across database versions**, and the interesting question is where the two
still agree.

The phases do map one to one, and the molar masses of their Table 3 settle it:
AFm 622.5, AFt 1255.1, Friedel's salt 561.3, CH 74.1 and C-S-H 191.4 all
reproduce, to better than 0.03 g/mol, the molar masses the package computes from
the CEMDATA18 formulas and from their own C-S-H formula. (Their CaCO₃, printed at 100.9, is
a typo for 100.09; it touches nothing in this figure.)

!!! warning "The C-S-H carries 44.5 g of water, and it is easy to drop"
    Guo write the C-S-H dissolution for `(CaO)₅(SiO₂)₃(H₂O)₆.₃`, `M = 574.1`,
    while their Table 3's `191.4` is **one third** of that. So 225 g is
    1.1755 mol of `(CaO)₁.₆₆₇(SiO₂)(H₂O)₂.₁` and carries **2.469 mol — 44.5 g —
    of water inside the solid**.

    Entering the C-S-H as lime and silica alone leaves the system 44.5 g short,
    and the solver then takes that water back out of the pore solution to
    hydrate the C-S-H: the pore volume comes out at **84 mL** against the
    146 mL the stated porosity implies, every concentration is wrong by nearly a
    factor of two, and Friedel's salt appears at 0.1 % NaCl where the paper says
    it cannot. Counting the water moves the pore volume to 129 mL and the onset
    back above 0.1 %.

### On the phase list that drew the figure

Guo's dissolution table carries five solids: C-S-H, CH, AFm, AFt and Friedel's salt.
Restricted to those, the figure comes back (mol per liter of concrete):

| | calculated | Fig. 1(b) |
|:--|--:|--:|
| AFm at 0 % NaCl | 0.014455 | 9 g / 622.5 = 0.014458 |
| AFt at 0 % | 0.017928 | 22.5 g / 1255.1 = 0.017927 |
| AFm at 1 % | **0** | consumed by ≈ 1 % |
| AFt plateau | 0.02274 | 0.023 |
| Friedel's salt plateau | 0.00965 | 0.010 |

What makes the agreement more than two curves happening to sit on top of each
other is that both conservation statements behind it close as well. Writing
`Δ` for what the AFm loses:

- **sulfate** — the AFm gives up one, ettringite takes three, so
  `ΔAFt = ΔAFm / 3`: 0.00482 against 0.014455/3 = 0.004818;
- **aluminum** — what the AFm held is split between Friedel's salt and the
  ettringite that grew, both carrying two aluminums per formula, so
  `FS + ΔAFt = ΔAFm`: 0.00965 + 0.00482 = 0.01447, against 0.014455.

### What a fuller phase list adds

Two things Guo's five solids cannot express, and one place where it stops
mattering.

**The hydration state.** CEMDATA18 carries monosulfate at 9, 10.5, 12, 14 and 16
waters, and at these conditions the stable one is the **14-hydrate**, not the 12
Guo wrote. That is a difference between database versions, not a slip of theirs:
their `log K = −29.2628` is Cemdata07's value for the 12-hydrate, and Cemdata18
*recalculated* it to −29.23 — its Table 2 marks that entry `***`, "recalculated
in this paper from ΔfG° values" — while carrying −29.26 for the 14-hydrate.
Their number therefore lands within 0.003 of CEMDATA18's 14-hydrate by
arithmetic coincidence. The amount is the same either way.

**Kuzel's salt**, half a chloride and half a sulfate per AFm layer, is the phase
the transition actually passes through. It holds the entire low-chloride range
and peaks at 0.5 % NaCl — where Guo has Friedel's salt barely starting:

| % NaCl | 0.1 | 0.25 | 0.5 | 0.75 | 1.0 | 1.5 |
|:--|--:|--:|--:|--:|--:|--:|
| monosulfate (14-hydrate) | 0.0131 | 0.0085 | 0.0010 | 0 | 0 | 0 |
| Kuzel's salt | 0.0011 | 0.0047 | 0.0107 | **0.0116** | 0.0058 | 0 |
| Friedel's salt | 0 | 0 | 0 | 0 | 0.0048 | 0.0096 |

The sweep solves `0`, `0.1`, `0.25`, `0.5`, `0.75`, `1`, `1.5` and `2 %`, and
the six columns above are pinned. The largest value is at `0.75 %`, but between
the rungs the true maximum is not located more finely than that.

**And where it stops mattering.** Once the chloride is high enough to take the
last sulfate out of the AFm layer, Kuzel's salt is gone and the two phase lists
agree to three digits. Guo's plateau is right even though the path to it is not
theirs.

### Not checked, and why

- **The pore volume**, which comes out at 129 mL against the 146 mL their
  14.6 % porosity states. The 17 mL is the difference between Guo's C-S-H,
  carrying 2.1 H₂O per `(CaO)₁.₆₆₇(SiO₂)`, and CEMDATA18's CSHQ end members,
  which carry more and take the rest out of the pore. Nothing closes that gap
  without replacing the C-S-H model.
- **pH.** Guo reports 13.213 falling to 13.128; this system sits at 12.5, which
  is portlandite in alkali-free water. The difference is their pore solution's
  alkalis, whose concentrations the paper takes from a reference it does not
  reproduce.
- **The 5 % end of their abscissa.** At that loading the ionic strength reaches
  1.4 mol/L, past the range the B-dot term of [`HKFActivityModel`](@ref) was
  fitted for; the solve stops certifying and returns pH 15.3. The sweep here
  stops at 2 %, where `I ≈ 0.4`. Pitzer is the instrument for the rest.
- **The surface complexation half of the paper** — the five `≡SiOH` reactions,
  the diffuse layer, the Ca²⁺ > Cl⁻ > Na⁺ > K⁺ ordering. Four of the five rows
  can be entered as printed; the surface is not built in this chapter. See
  [The surface half of Guo](@ref).

## The surface half of Guo

[Guo2018](@cite) closes their chloride model with a `≡SiOH` surface on the
C-S-H — 500 m²/g, 4·10⁻³ mol of sites per gram, a diffuse layer, and five
reactions in their Table 1. Four of the five rows can be entered as printed; one
cannot.

| row | reaction as printed | `log K` | status |
|:--|:--|:--|:--|
| 1 | `≡SiOH + OH⁻ → ≡SiO⁻ + H₂O` | `−12.7` | constant belongs to a different reaction |
| 2 | `≡SiOH + Ca²⁺ → ≡SiOCa⁺ + H⁺` | `−9.4` | usable |
| 3 | `≡SiOH + Cl⁻ → ≡SiOHCl⁻` | `−0.35` | usable |
| 4 | `≡SiOH + Na⁺ → ≡SiONa + H⁺` | `−13.6` | usable |
| 5 | `≡SiOH + K⁺ → ≡SiOK + H⁺` | `−13.6` | usable |

The table is [Elakneswaran2010](@cite)'s: its Eq. (11) is row 1 with the same
`−12.7`, its Eq. (13) row 3 with the same `−0.35`, and its Eqs. (15) and (16)
rows 4 and 5. So row 1 is not a slip of Guo's; the source writes it that way.

As written, row 1 cannot be what the authors used. Against `OH⁻`, a constant of
`10⁻¹²·⁷` leaves `≡SiO⁻/≡SiOH = 10⁻¹³·⁷` at pH 13 — an uncharged surface — while
the source reports a negative ζ at that pH and explains it by this very
dissociation. In the proton form, `≡SiOH ⇌ ≡SiO⁻ + H⁺`, the same constant gives a
ratio of 2 at pH 13, two thirds of the sites deprotonated before electrostatics,
which is the surface they describe; written against `OH⁻` it would carry
`log K ≈ +1.3`. The evidence points one way, but repairing the row is still
choosing which half of it to believe, which is a modeling decision rather than a
transcription.

Row 3, the chloride row, is balanced: the product carries the charge of the
chloride it took up, `≡SiOHCl⁻`, as both papers print it. It is the row the
surface exists to provide.

Two things of the source are not in Guo's table. Its Eq. (14),
`≡SiOH + Ca²⁺ + Cl⁻ ⇌ ≡SiOCaCl + H⁺`, holds chloride on a site that has already
taken up calcium. And its double layer is written in the Donnan approximation,
the ions of the layer related to those of the free solution by a Boltzmann
factor and the layer kept electrically neutral (its Eqs. (9)-(10)). The surface
is not built here: it belongs with a model of chloride binding by C-S-H, and
this chapter checks the dissolution half only.

### The machinery is not what is missing

Worth separating, because the two failures look alike from the outside. The
site-family route reproduces the chemistry it is given. On a 1 kg pore solution
holding 20 mmol Ca, 100 mmol NaCl and 10 mmol KCl, with rows 1, 2, 4 and 5 and
`DiffuseLayer`, every point certifies and the occupancies behave:

| sites, site model | `≡SiOH` | `≡SiO⁻` | `≡SiOCa⁺` | `≡SiONa` | `≡SiOK` |
|:--|--:|--:|--:|--:|--:|
| Guo's 0.90 mol, `DiffuseLayer` | `0.939` | `0.038` | `0.022` | `0.0007` | `0.0001` |
| 1/100 of that, ideal mixing | `0.175` | `0.120` | `0.704` | `0.0011` | `0.0001` |
| 1/100 of that, `DiffuseLayer` | `0.217` | `0.376` | `0.406` | `0.0013` | `0.0001` |

Three checks on those numbers, none of which needs the paper.

**The inventory bounds the calcium.** At Guo's loading there are 0.9 mol of
sites and 0.02 mol of calcium in the entire system, so the calcium fraction
cannot exceed `0.02/0.9 = 0.0222` — and it sits at `0.0219`. The surface has
taken essentially all of it, and the ceiling is the inventory, not the affinity.

**The alkali competition is exact.** Rows 4 and 5 share a `log K`, and both
cations are monovalent, so the competitive closed form collapses to

```
≡SiONa / ≡SiOK = a(Na⁺) / a(K⁺)
```

with the site capacity, the shared denominator and the proton activity all
canceling. At one hundredth of the loading the surface removes `9.5·10⁻⁶` mol
of sodium out of `0.1`, so the free concentrations still stand at `10.000` — and
the occupancy ratio comes out `9.864`. The residual `1.4 %` is `γ(Na⁺)/γ(K⁺)`,
the one term the cancellation leaves behind. That is a check on the site
mixing itself, not on any constant.

**The diffuse layer pushes the right way.** At 70 % `≡SiOCa⁺` against 12 %
`≡SiO⁻` the surface carries `+0.58` charge per site, so the electrostatic term
must repel Ca²⁺ and favor deprotonation — which is the `0.704 → 0.406` shift in
the last two rows. The sign is worth checking explicitly, because the intuition
that a silanol surface is negative is wrong here: calcium has reversed it.

### What a budget that follows the C-S-H would cost

The numbers above hold a **fixed** site budget, which is the wrong idealization
for this system: the C-S-H is a reaction product, so its amount moves and the
sites should move with it. `SITES_FOLLOW_HOST` expresses that, and it is worth
knowing before reaching for it here, because Guo's site density is high enough
to make the reference energy of the free site matter a great deal.

At `4·10⁻³ mol/g` on a `(CaO)₁.₆₆₇(SiO₂)(H₂O)₂.₁` of `191.4 g/mol`, the coupling
ratio is `ν = 0.766` mol of sites per mol of host — nearly four times Dzombak
and Morel's `ν = 0.2` for hydrous ferric oxide, where
[`host_coupling_bias`](@ref) records `8.3` log units. The bias is linear in `ν`,
so here it is

```
ν · |μ°(H₂O) − μ°(H⁺)| / (RT ln 10) = 0.766 × 237.2 / 5.708 = 31.8 log units
```

which would not perturb the C-S-H, it would annihilate it. With the reference
set to `−237.2 kJ/mol` the bias is zero and the coupling is free; left at zero,
the family is refused and the number printed. So the fixed budget used above is
not merely a simplification here — it is the reason the free site's energy
could be left at zero without consequence.

For a surface case that *is* checked against an external oracle, see
`test/surface_complexation.jl`, which runs against PHREEQC at matched proton
activity, and `test/diffuse_layer.jl`.

## The HKF model away from 25 °C and 1 bar

The Cemdata18 section above pins every HKF equation-of-state coefficient in the
shipped file. Nothing there pins what the model *does* with them once the
temperature and pressure leave the reference point — and that is the half that
carries a burial, an autoclave or a steam-cured calculation.

[Duan2016](@cite) is a weak source in most respects: their carbonate data for
ferrocalcite and ankerite are estimated rather than measured, their tables
disagree with their own text about which is which, and their conclusion prints
"−2.24 mmol·L⁻¹mmol·L⁻¹". But their Table 4 carries one column this package can
be held to — two constants computed by **HKF itself**, taken from Yu, Dong and
Ruan (2008), at two temperatures. That is an HKF oracle, and it is the only one
in this chapter.

| | ChemistryLab | Duan Table 4 | difference |
|:--|--:|--:|--:|
| `K₃` (HCO₃⁻ ⇌ H⁺ + CO₃²⁻), 298.15 K, 1 bar | 4.6885×10⁻¹¹ | 4.69×10⁻¹¹ | −0.03 % |
| `K₄` (H₂O ⇌ H⁺ + OH⁻), 298.15 K, 1 bar | 9.9971×10⁻¹⁵ | 1.00×10⁻¹⁴ | −0.03 % |

### Their pressure column is in bar, and `K₃` proves it

Table 4 is headed `P(Pa)` and reads `1.0` and `1000.0`. Neither can be pascals —
1 Pa is a hard vacuum, and no aqueous constant is tabulated there. Bar fits, and
the second row settles it. At 373.15 K the published `K₃` is 2.42×10⁻¹⁰, against

| pressure | `K₃` | difference |
|:--|--:|--:|
| 1 bar | 8.240×10⁻¹¹ | −66 % |
| **1000 bar** | **2.412×10⁻¹⁰** | **−0.3 %** |

Reproducing a factor of 2.9 in pressure to three parts in a thousand is a real
check of the volume terms, and it is the strongest single piece of evidence in
this chapter that the HKF implementation is right away from the reference point.

!!! warning "The same row mixes two pressures"
    `K₄` does not agree at 1000 bar. It matches at **1 bar** — 5.504×10⁻¹³
    against their 5.38×10⁻¹³, 2.3 % — and is 14.7 % out at 1000 bar. Water
    autoprotolysis at 100 °C and 1 bar is `pKw = 12.26`, so the 1 bar value is
    the physical one and the row carries its two constants at different
    pressures.

### Calcite, and a method that is not HKF

Their Table 6 gives calcite over 298–478 K and 0.1–70 MPa, computed with their
own Gibbs-energy integration and SRK volumes rather than with HKF. So a
disagreement is a difference of method, not a defect — but at the reference
point there is an independent arbiter, and it favors this package: the accepted
`log Ksp` of calcite at 25 °C and 1 bar is **−8.48**, which is what comes out
here. Duan print −8.53.

| T (K) | P (MPa) | ChemistryLab | Duan | Δ |
|--:|--:|--:|--:|--:|
| 298.15 | 0.1 | **−8.480** | −8.53 | +0.05 |
| 301.15 | 15 | −8.440 | −8.53 | +0.09 |
| 301.15 | 70 | −8.261 | −7.82 | −0.44 |
| 343.15 | 15 | −8.855 | −8.69 | −0.17 |
| 418.15 | 40 | −9.827 | −9.06 | −0.77 |
| 478.15 | 15 | −11.055 | −9.69 | −1.37 |

Both agree on the two signs a burial calculation turns on — heating dissolves
less, compressing dissolves more — and diverge steadily with temperature, this
package giving the lower solubility. The divergence is worth recording because
Duan's own Table 4 puts their method within 3 % of HKF on `K₃` and `K₄`, which
is 0.01 log units; 1.4 log units on calcite is far outside what they claim for
it. Their estimated carbonate data are the likeliest reason — the ferrocalcite
row of their Table 3 carries a heat-capacity coefficient of `+2.09×10⁶` where
every other carbonate in the table has zero or a large negative.

## Cemdata07, and an aged pore solution

[Lothenbach2010](@cite) reviews what equilibrium calculations do for cements,
and its Table 1 lists the solubility products it computed with: those of
Cemdata07, the generation before the one this package ships, for 29 solids. It
prints no Gibbs energy and no molar volume, only the constant and its reaction.
Read through the package's own energies, the same reactions separate what
Cemdata18 kept from what it revised.

| | phases |
|:--|:--|
| the same to the printed digit (11) | ettringite, tricarboaluminate, Fe-ettringite, monocarbonate, hemicarbonate, strätlingite, hydrotalcite (`M₄AH₁₀`), tobermorite-type C-S-H, amorphous silica, syngenite, amorphous `Al(OH)₃` |
| revised (9), ours minus Cemdata07 | Fe-monosulfate `+1.63`, C₃AH₆ `+0.34`, C₄AH₁₃ `+0.15`, monosulfate `+0.026`, jennite-type C-S-H `+0.007`, thaumasite `−0.10`, CAH₁₀ `−0.10`, C₃FH₆ `−1.14`, C₄FH₁₃ `−1.35` |
| no CEMDATA18 solid of the same composition (9) | siliceous hydrogarnet `C₃AS₀.₈H₄.₄`, C₂AH₈, C₂FH₈, Fe-monocarbonate, Fe-hemicarbonate, Fe-strätlingite, `M₄AC̄H₉`, `M₄FH₁₀`, `Fe(OH)₃` |

Seven of the eight values the paper marks as tentative are among the revised or
the absent; hydrotalcite is the one that survived. The phases are paired on
their composition, never on their constant, which takes two precautions:
amorphous silica shares its formula with quartz, and amorphous `Al(OH)₃` with
gibbsite, so those two are paired by name; and CEMDATA18 carries several solids
also at half their formula unit, for its solid solutions, so the whole unit is
preferred and the constant scaled when only a half exists.

### Cemdata07's energies and molar volumes

The same generation is printed in full by [Lothenbach2008](@cite), Table 4:
the Gibbs energy, enthalpy, entropy and molar volume of 37 solids and water,
with the solubility products of [Lothenbach2010](@cite) phase for phase. Its
energies can be set against the package's without any reaction in between, and
its volumes, most of them calculated from unit cells or densities, have no
reason to move with a revision of the solubility products.

| | energies | molar volumes |
|:--|:--|:--|
| the same to the printed digit | 21 solids and water, to `0.01 kJ/mol` | 24 of 28 solids, to `0.5 cm³/mol` |
| different | 7 revised: Fe-monosulfate `+9.30`, C₃AH₆ `+1.94`, C₄AH₁₃ `+0.87`, monosulfate `+0.14`, Fe-monocarbonate `+5.19`, C₃FH₆ `−6.53`, C₄FH₁₃ `−7.71` kJ/mol | Fe-monocarbonate `+1.67`, monosulfate `+1.10`, Fe-monosulfate `−0.86`, Fe-ettringite `+0.56` cm³/mol |
| no CEMDATA18 solid of the same composition | 9 | 9 |

The revisions concern the iron phases and the hydrogarnets, which Cemdata18
revisited, and the monosulfate, recalculated from its Gibbs energy.

Three of the nine unpaired rows do have a CEMDATA18 counterpart at another
hydration state, and the difference in Gibbs energy between the two is the
energy of the water that separates them:

| Cemdata07 | CEMDATA18 | water | kJ/mol per H₂O |
|:--|:--|--:|--:|
| C₂AH₈ | C₂AH₇.₅ | 0.5 | `−234.4` |
| C₄FC̄₀.₅H₁₂ | Fe-hemicarbonate | 2 | `−243.7` |
| Fe(OH)₃ (mic.) | FeOOH (mic.) | 1 | `−231.5` |

Each is within 3 % of liquid water's `−237.18 kJ/mol`, which identifies the three
as compositional differences rather than revisions. The water is counted on the
hydrogen: the CEMDATA18 formula of Fe-hemicarbonate ends in `(H₂O)9.5`, and it
carries ten waters, not 9.5.

### An aged pore solution

Table 2 of the same paper is the pore solution of an ordinary Portland cement
paste (w/c 0.4) after 69 days, expelled at five pressures, and the same page
reports that the pore solutions of old pastes are often found a little above
saturation for portlandite and for ettringite. At an ionic strength near 0.5 mol/kg, that
statement depends on the activity model. With the analyzed totals at 25 °C (the
table gives no temperature), the hydroxide from the charge balance, and lithium
and strontium left out, the effective saturation indices — the index divided by
the number of ions, 3 and 15, as the paper defines them — are:

| activity model | portlandite | ettringite |
|:--|--:|--:|
| Davies | `+0.14` to `+0.16` | `+0.18` to `+0.23` |
| B-dot, as identified from GEM-Selektor (`å = 0`) | `−0.11` to `−0.10` | `−0.18` to `−0.15` |

Both put the solution within a quarter of a log unit of equilibrium, and
neither can say on which side of it. With the second model the free hydroxide
the paper reports comes back within 6 % from the charge balance of the totals.

## Alkali uptake by C-S-H

[HongGlasser1999](@cite) equilibrated synthetic C-S-H of Ca/Si 0.85, 1.2, 1.5 and
1.8 with NaOH and KOH solutions of 1 to 300 mM, 0.6 g of gel in 9 mL at 20 °C,
and analyzed 48 solutions for alkali, calcium, hydroxide and pH. Nothing else
enters: no clinker, no kinetics, no unpublished input.

!!! warning "The alkali was fitted to these data"
    The Cemdata18 authors fine-tuned the Gibbs energies of their Na and K C-S-H
    end members on these isotherms ([Lothenbach2019](@cite), §2.7 and Fig. 10).
    Agreement on the alkali checks that the model is applied as it was fitted.
    The calcium and the pH were not fitted, and they are where it can disagree.

The gel enters as lime, silica and its own water: 14 % of the soft-dried mass at
Ca/Si 0.85 and 1.2, as measured, and 18 % above, the value the authors adopt from
Taylor. The 9 mL of solution are taken as 9 g of water. The temperature is set
on the state, at 20 °C, and the activity model is the B-dot form identified from
a GEM-Selektor run of CEMDATA18, as for [the Atkins case](@ref "Atkins et al. (1992): a measurement, not a calculation").
Each Ca/Si ratio is walked down from 300 mM, and every point certifies.

At the 100 mM rung, measured / computed:

| Ca/Si | alkali | alkali, mmol/L | Ca, mmol/L | pH | Rd, mL/g |
|:--|:--|--:|--:|--:|--:|
| 0.85 | Na | 79.0 / 69.7 | 0.10 / 0.15 | 12.85 / 12.805 | 3.82 / 6.41 |
| 1.2 | Na | 87.0 / 83.1 | 0.63 / 0.43 | 12.96 / 12.940 | 2.08 / 2.99 |
| 1.5 | Na | 95.3 / 91.5 | 5.20 / 3.81 | 13.01 / 12.999 | 0.60 / 1.30 |
| 1.8 | Na | 96.9 / 93.7 | 4.88 / 8.33 | 13.00 / 13.033 | 0.34 / 0.93 |
| 0.85 | K | 78.2 / 76.1 | 0.10 / 0.15 | 12.87 / 12.835 | 4.33 / 4.94 |
| 1.2 | K | 89.5 / 87.6 | 1.02 / 0.41 | 12.98 / 12.964 | 1.89 / 2.36 |
| 1.5 | K | 94.0 / 94.6 | 4.80 / 3.62 | 13.00 / 13.015 | 1.09 / 1.03 |
| 1.8 | K | 98.2 / 96.4 | 3.90 / 8.06 | 13.04 / 13.046 | 0.40 / 0.74 |

Over the 48 points:

  - **The pH**, which nothing was fitted to, is within `0.072` of the electrode
    from 15 to 100 mM on the three gels below Ca/Si 1.8. At 300 mM it is low at
    every Ca/Si, by `0.08` to `0.21`.
  - **The alkali is over-bound.** Less stays in solution than was measured at 46
    of the 48 points, and the computed Rd at 100 mM is 0.95 to 2.7 times the
    measured one. The measured Rd hardly depends on the concentration, which is
    the paper's main finding; the computed one falls as it rises.
  - **The calcium** is low at Ca/Si 1.2 and 1.5 and high at 1.8. At 1.8 the CSHQ
    model cannot hold all the calcium of the gel, and portlandite precipitates,
    3.7 to 6.8 % of the solid, where the preparation carried 0.2 %.
  - Aluminum is outside the model: [HongGlasser2002](@cite) replaced 6 to 7 % of
    the silicon by aluminum and found a markedly higher Rd, and the CSHQ solid
    solution carries no aluminum.

What the measurements themselves can carry is read from the paper's own columns.
It prints the cation sum `Na⁺ + 2Ca²⁺` beside the measured hydroxide, and the gap
between them is the silicate that was not measured. The authors call ±5 %
satisfactory at high ionic strength; two of the sixteen solutions at 100 and
300 mM exceed it, the Ca/Si 0.85 gel in NaOH at 100 mM (9.1 %) and in KOH at
300 mM (5.7 %). In the sixteen dilute solutions, at 1 and 5 mM, the gap reaches
55 % and is 16 % on average, so nothing compared there carries much. One printed
sum disagrees with its own row: 45.0 for `14.6 + 2 × 15.5 = 45.6` (KOH, 15 mM,
Ca/Si 1.8).

## Reproducing

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

The Cemdata18 checks are pure arithmetic over the shipped file and need no
solver. The Atkins case runs one `equilibrate_certified` and takes about a
second once compiled.
