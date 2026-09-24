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

Every number below is pinned by an assertion in `test/cemdata18_reference.jl`
(584 of them), `test/atkins1992_reference.jl` (12), `test/limestone_blending_reference.jl` (25)
`test/chloride_binding_reference.jl` (25) or
`test/duan2016_reference.jl` (11), so it is checked on every CI run.

## What has been checked

Five sources, seven hundred assertions, and the coverage is uneven on purpose:
the database is checked exhaustively because it is cheap to check exhaustively,
while the equilibrium cases are checked one composition at a time because each
one costs seconds to minutes.

| what | against | how closely | where it stops |
|:--|:--|:--|:--|
| **every solubility product** in the shipped CEMDATA18 — 52 phases | [Lothenbach2019](@cite) Tables 2-3 | 50 close to `0.041`, 48 to `0.005` | two M-S-H end members are `0.48` and `0.40` out; the source disagrees with itself |
| **standard properties and HKF coefficients**, 19 aqueous species and 7 gases | [Lothenbach2019](@cite) Tables D.1-D.2 | exact, all seven coefficients each | nitrite-AFm and Fe-Friedel's salt are not in the file |
| **`ΔₐG⁰` rebuilt from `ΔfH°` and `S°`** | the file's own `ΔfG°` | exact for 220 of 228 substances | the eight exceptions are all phases whose `S°` Cemdata18 estimated |
| **HKF away from 298.15 K, 1 bar** | [Duan2016](@cite) Table 4, HKF column | `0.03 %` at the reference point; `0.3 %` across a factor 2.9 in pressure | their two constants in one row sit at two different pressures |
| **calcite `log Ksp`** at 25 °C | the accepted value | `−8.480` against `−8.48` | diverges from Duan's non-HKF method by 1.4 log units at 478 K |
| **a measured solution** over a two-phase assemblage | [Atkins1992](@cite) Table 2 | Al and pH agree, robustly | Si is `×5` and Ca has a floor the model cannot leave; 9 of their 10 mixtures are not usable at all |
| **the carboaluminate sequence** under limestone | [Kulik2021](@cite) Fig. 7A, [Lothenbach2019](@cite) Figs. 13-14 | order and thresholds reproduce; iron partition exact | only appears when Al exceeds Fe — see below |
| **chloride binding** and the AFm → Friedel transition | [Guo2018](@cite) Fig. 1(b) | plateau to `1 %`, both conservation laws close | pH not reproducible (their alkalis are unpublished); above 2 % NaCl the activity model is out of range |

Two things this chapter deliberately does **not** do. It does not check the
kinetics — [Validation against Reaktoro](@ref) covers the coupling, and no
published hydration curve is reproduced here. And it does not check anything
above about `1 mol/kg` ionic strength, because `HKFActivityModel`'s B-dot term
was not fitted there.

## Traps

Five things that cost time to find, collected because each is reusable well
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

**Aluminum is the strong result.** Across a sweep of the C-S-H Ca/Si from 0.75
to 1.0 and a fourfold change in solid loading it never leaves 0.138–0.159
mmol/L. That is the database answering, not a fit.

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
| 0.75 | 3.39 | 0.156 | 0.498 | 2.41 | 11.02 |
| 0.80 | 3.04 | 0.154 | 0.463 | 1.90 | 11.14 |
| 0.85 | 2.83 | 0.152 | 0.429 | 1.50 | 11.24 |
| 0.90 | **2.72** | 0.149 | 0.396 | 1.18 | 11.33 |
| 1.00 | 2.73 | 0.142 | 0.335 | 0.74 | 11.47 |
| *measured* | *1.95* | *0.136* | *0.076* | *1.08* | *11.0* |

Calculated calcium goes through a **minimum near Ca/Si 0.85–0.90** and rises
again below it. It never approaches 1.95. The floor is a property of the CSHQ
model, not a free parameter, so the gap is real.

Sulfate lands close, but it is the one number here that moves freely with the
Ca/Si — 0.74 at 1.0, 2.41 at 0.75 — so the agreement is not evidence of much.
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

Measured here: **100.0 %**, at every limestone content and both iron contents
tried, to within `10⁻³`. That one reproduces without qualification.

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
Glasser (2008), which is Cemdata07. So this is a comparison **across database
versions**, and the interesting question is where the two still agree.

The phases do map one to one, and the molar masses of their Table 3 settle it:
AFm 622.5, AFt 1255.1, Friedel's salt 561.3, CH 74.1 all reproduce from the
CEMDATA18 formulas to better than 0.03 g/mol. (Their CaCO₃, printed at 100.9, is
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
| Friedel's salt plateau | 0.00964 | 0.010 |

What makes the agreement more than two curves happening to sit on top of each
other is that both conservation statements behind it close as well. Writing
`Δ` for what the AFm loses:

- **sulfate** — the AFm gives up one, ettringite takes three, so
  `ΔAFt = ΔAFm / 3`: 0.00482 against 0.014455/3 = 0.004818;
- **aluminum** — what the AFm held is split between Friedel's salt and the
  ettringite that grew, both carrying two aluminums per formula, so
  `FS + ΔAFt = ΔAFm`: 0.00964 + 0.00482 = 0.01446.

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
  the diffuse layer, the Ca²⁺ > Cl⁻ > Na⁺ > K⁺ ordering. `SC_SURFCOMPLEX` and
  the site families arrived in v0.20.0, so it is now expressible; it is simply
  not attempted yet.

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

## Reproducing

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

The Cemdata18 checks are pure arithmetic over the shipped file and need no
solver. The Atkins case runs one `equilibrate_certified` and takes about a
second once compiled.
