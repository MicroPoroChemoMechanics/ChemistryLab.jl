# Experimental data

Measured data, vendored so that the test suite and the documentation build
without network access. Unlike everything else under `data/`, these files are
**not** thermodynamic databases and **not** covered by the package's license —
see [`LICENSE`](LICENSE) in this directory.

| record | cement | w/b | Blaine (m²/kg) | duration | final *Q* (J/g) | role |
|--- |--- |--- |--- |--- |--- |--- |
| `122-cemI-52.5R-cizkovice` | CEM I 52.5 R Cizkovice | 0.50 | 397 | 262 h | **376** | calibration target |
| `116-cemI-52.5R-ladce` | CEM I 52.5 R Ladce | 0.45 | 415 | 617 h | **355** | holdout |
| `165-cemII-A-LL-42.5R-hranice` | CEM II-A-LL 42.5 R Hranice | 0.45 | 423 | 384 h | **329** | limestone-blended Portland |
| `149-cemII-B-S-32.5R-mokra` | CEM II-B-S 32.5 R Mokra | 0.40 | 380 | 504 h | **296** | slag-blended Portland |
| `184-cemIII-A-42.5N-hranice` | CEM III-A 42.5 N Hranice | 0.40 | 408 | 359 h | **261** | blastfurnace, 36-65 % slag |
| `200-cemV-A-S-V-32.5R-prachovice` | CEM V-A (S-V) 32.5 R Prachovice | 0.40 | 444 | 501 h | **259** | composite, slag + fly ash |
| `121-cemIII-B-32.5N-mokra` | CEM III-B-32.5 N LH-SR Mokra | 0.45 | 495 | 663 h | **234** | blastfurnace, 66-80 % slag |

!!! warning "The deposit's internal `Cement name` field is not reliable — the file name is"
    Each source record carries a `Cement name:` line, and it disagrees with the
    file name often enough that it cannot be used. Record 122 gives
    `Cizkovice 52.5R`, with no EN 197-1 designation at all. Record 116 gives
    `CEM I 42.5R Ladce` where its own file name says **52.5 R** — the two
    contradict each other, and nothing in the deposit settles which is right.

    The designations in the table above are therefore taken from the **file
    name**, which is the deposit's own index and is consistent across all 65
    records. The internal field is left in each vendored file, unedited, so a
    reader can see the discrepancy rather than take our word for it.

Sorted by the heat released, which is the physical ordering: **every joule comes
from the clinker**, so replacing clinker with slag, fly ash or limestone lowers
it. The 376 J/g of a CEM I and the 234 J/g of a CEM III/B are the two ends of the
EN 197-1 range, measured on the same instrument at the same temperature.

Both come from the CC-BY-4.0 Zenodo deposit
[10.5281/zenodo.15212785](https://doi.org/10.5281/zenodo.15212785) of Šmilauer
and Reiterman, and are used by
[`scripts/hydration_calibration.jl`](../../scripts/hydration_calibration.jl).

## Format

A `#`-commented provenance block, then a plain comma-separated table with the
header `time_h,heat_flow_W_per_g,heat_J_per_g`. `read_calorimetry` in the script
reads it, and reads the metadata lines it needs out of the comment block.

Any file in this shape can be substituted: that is the seam for calibrating
against your own measurements. The four comment lines the reader looks for are
`blaine`, `wb`, `temperature` and `Released heat up to`.

## Why these two records

Of the 65 cements in the deposit, 14 are plain CEM I. Record `122` was chosen as
the calibration target because its w/b of 0.50 and Blaine of 397 m²/kg sit within
a few percent of the CEM I 52.5 N formulation that
[`scripts/ionic_hydration.jl`](../../scripts/ionic_hydration.jl) already runs
after Lavergne et al. (2018) — so the one input the deposit does not report, the
clinker phase composition, is the least of an extrapolation available.

Record `116` is the holdout: same nominal strength class, but a different w/b
**and** a different Blaine. Predicting it with parameters fitted on `122`, without
refitting, is a direct test of the two corrections that claim to carry a
calibration across mixes — `powers_alpha_max(w/b)` and `blaine_factor(blaine)`.

## Caveats, all of them

These matter for interpreting a fit and are stated in the documentation page too.

- **The records do not start at *t* = 0.** Both are truncated at the calorimeter's
  thermal equilibration — 1.02 h for `122`, 1.21 h for `116`. The header line
  `Released heat up to 45 minutes (J/g of binder)` gives the heat already gone
  before the record begins (12.0 J/g, "estimated", for both). A model compared
  against these curves must be offset by that amount, not zeroed at the first
  sample.
- **No phase composition, no oxide analysis.** The deposit reports the cement
  type, the Blaine fineness, the w/b ratio and the temperature — not the clinker
  mineralogy. The composition therefore has to be assumed, and only the products
  of a fitted rate multiplier and its assumed phase fraction are identifiable.
- **One temperature.** Everything was measured at 20 °C, so activation energies
  are not identifiable from these data and must be held at published values.
- **Filename and internal header disagree for the Ladce records.** The file named
  `116-CEM I 52.5 R Ladce-415.csv` in the deposit carries
  `Cement name: CEM I 42.5R Ladce-415` internally; the same is true of `115`.
  Which of the two is right is not resolvable from the deposit, which is why
  `116` is used as a plausibility check rather than as a validation.
- **The depositors' own reference fit is quoted twice, inconsistently.** The
  provenance block reports the affinity-model parameters from the deposit's
  `Fit-affinity.zip`, which give `Ea 38300 J/mol` for every cement. The header of
  the source data file for `122` carries a different line, with `AEn 40000
  J/mol`. The value from the fit archive is the one kept, since it is the one that
  goes with the fitted `B1`, `B2` and `eta`.

## Regenerating

```
julia data/experimental/regenerate.jl            # download and rebuild
julia data/experimental/regenerate.jl --check    # rebuild and compare only
```

Needs network access and `unzip`; neither is needed to *use* the files. The
subset is deterministic — nearest row to each of 500 log-spaced times — so
`--check` is an exact byte comparison, and it is how you tell that nothing here
was hand-edited.
