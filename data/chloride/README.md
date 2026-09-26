# A chloride end member for CSHQ

`data/cemdata18-chloride.json` is CEMDATA18 with one end member appended,
`CSHQ-Cl` = (CaCl2)0.5, and `CSHQ_Cl` in `data/solid_solutions.toml` is CSHQ
with that end member. The file is **generated**, and `regenerate.jl` in this
directory is what generates it:

    julia --project=docs data/chloride/regenerate.jl

It is deterministic; `test/cshq_chloride.jl` checks that the shipped file is
what it writes.

## What it is for, and what it is not

CSHQ describes the calcium of the C-S-H and its alkalis, not chloride. The
published surface model of the C-S-H binds chloride on calcium adsorbed at its
silanol sites, and cannot be put on CSHQ without counting that calcium twice.
It can be put on a C-S-H frozen after a first equilibrium, which is what
`freeze_solid_solution` does, but only while portlandite fixes the composition
of the gel. An end member of CSHQ holds chloride at any Ca/Si, portlandite or
not.

It is an **effective** description. Plusquellec and Nonat (2016) found that
chloride does not adsorb specifically on C-S-H: the chloride a depletion measurement
counts as bound accompanies the calcium the surface adsorbs, and sits in the
diffuse layer that screens it. An end member carrying CaCl2 lumps both, and says
nothing about where in the gel the chloride is.

## The data

Hirao, Yamada, Takahashi and Zibara (2005), *Chloride binding of cement estimated
by binding isotherms of hydrates*, Journal of Advanced Concrete Technology 3,
77-84, [10.3151/jact.3.77](https://doi.org/10.3151/jact.3.77), Fig. 5: C-S-H from
alite hydrated 56 days, with the portlandite it formed, 1 g in 10 cm³ of NaCl
solution at 20 °C, the bound chloride from what the solution lost. The figure is
a vector drawing, and the points were read from its coordinates
(`data/literature/Hirao2005.json`). Only the three points up to 1 mol/L enter the
fit: above, the solution leaves the range of the B-dot activity model the
equilibria use.

The criteria written before the data were read asked for synthetic C-S-H at
several Ca/Si, with the solution analyzed. No set in hand meets them in full.
Hirao et al. report one gel, at portlandite saturation, and the chloride of the
solution only; portlandite fixes its calcium and its pH, so the test can be
recomputed, and it is the one used. The Ca/Si dependence is therefore not
fitted: it is what the model of an ideal solid solution gives.

Examined and not used:

| source | why not |
|:--|:--|
| Beaudoin, Ramachandran & Feldman (1990), Cem. Concr. Res. 20, 875, [10.1016/0008-8846(90)90049-4](https://doi.org/10.1016/0008-8846(90)90049-4) | CaCl2 solutions up to 3.8 %, an ionic strength of about 1 mol/kg, with no analysis of the solution but its chloride; figures raster only |
| Zibara, Hooton, Thomas & Stanish (2008), Cem. Concr. Res. 38, 422, [10.1016/j.cemconres.2007.08.024](https://doi.org/10.1016/j.cemconres.2007.08.024) | raster figures; its trend with the Ca/Si is used as a criterion, below |
| Yoshida, Elakneswaran & Nawa (2021), Cem. Concr. Compos. 121, 104109, [10.1016/j.cemconcomp.2021.104109](https://doi.org/10.1016/j.cemconcomp.2021.104109) | CaCl2 solutions from 0.5 mol/L, an ionic strength of 1.5 mol/kg and more, above the range of the activity model |
| Plusquellec & Nonat (2016), Cem. Concr. Res. 90, 89, [10.1016/j.cemconres.2016.08.002](https://doi.org/10.1016/j.cemconres.2016.08.002) | a finding about the mechanism, stated above, rather than a set to fit |
| Hirao et al. (2005), Fig. 12 | pastes whose test conditions are those of an earlier paper, not in hand |

All DOIs were resolved against the Crossref REST API.

## The fit

One number is fitted: δ, the Gibbs energy of forming the end member from
½ Ca²⁺ + Cl⁻ at 293.15 K, NaSiOH's reference temperature and that of the tests.
Its entropy and volume are those of the two ions, so that the dissolution has
neither, and its heat capacity is zero. These are estimates, and they only
enter between 20 °C and the temperature a calculation runs at.

The model of the test reproduces the measurement, not the binding alone: 1 g of
hydrate split between C-S-H and portlandite as the two capacities Hirao et al.
print imply (0.47001/0.61602 of C-S-H), the gel at the Ca/Si the model gives it
at portlandite saturation with 18 % of water (the value Hong and Glasser adopt
for dried gels of high Ca/Si), 10 cm³ of solution, and the bound chloride
computed from the concentration the model leaves in its liquid, as the depletion
computes it. The gel takes up water as it rehydrates, which concentrates the
solution: at 1 mol/L this hides 0.12 mmol/g of the 0.41 the end member holds.

The result is recorded in the end member's `chloride_provenance`: δ, its
uncertainty (from the residuals through `identifiability`, and from the water of
the gel, refitted at 14 %), the points and the residuals. One parameter cannot
follow the shape of the measured points: the fit is within 0.03 mmol/g at
0.5 mol/L and 0.05 low at 1 mol/L, and four times the measurement at 0.1 mol/L,
where Hirao's own Langmuir curve is 2.7 times it.

## Why this end member

Two were written and fitted. The first follows the sodium end member of CSHQ:
two formula units of NaSiOH with their NaOH exchanged for half as much CaCl2,
((CaCl2)1.25SiO2H2O)0.4. It fits the tests as well, and binds **less** chloride
on a gel of Ca/Si 1.4 than on one of 1.0, where Zibara et al. (2008) and Beaudoin
et al. (1990) measured more: the silica it carries is favored by the silica
activity of a low-Ca/Si gel. The one shipped carries calcium chloride alone, and
binds more at the higher Ca/Si. Both fits, and the trend of each, are recorded
by the generator; the trend is a criterion the end member was chosen on, not a
validation.

What neither the fit nor the choice saw: at equal chloride, a CaCl2 solution
makes the gel bind twice as much as a NaCl one, the ordering Tran et al. (2018,
[10.1016/j.conbuildmat.2018.10.058](https://doi.org/10.1016/j.conbuildmat.2018.10.058))
report.
