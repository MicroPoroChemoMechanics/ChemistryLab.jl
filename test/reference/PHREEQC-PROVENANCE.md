# `phreeqc.dat` — where it came from, and the terms it comes under

`phreeqc.dat` in this directory is the thermodynamic database of **PHREEQC**,
redistributed here **unmodified**.

| | |
|:--|:--|
| file | `database/phreeqc.dat` |
| source | <https://github.com/phreeqc-dev/phreeqc3>, tag **`v3.7.3`** |
| md5 | `4d3f4378c46498730d569cacc37c4fa8` |
| size | 48 364 bytes |
| note | byte-identical to tag `v3.7.1`; the two releases ship the same database |
| authors | D.L. Parkhurst and C.A.J. Appelo |
| publisher | U.S. Geological Survey |
| terms | see `PHREEQC-NOTICE.txt` beside this file |

`PHREEQC-NOTICE.txt` is the **User Rights Notice** shipped as `doc/NOTICE.TXT`
in the same distribution, reproduced verbatim. It is here because that notice
requires it: *"If you distribute copies or modifications of the software and
related material, make sure the recipients receive a copy of this notice and
receive or can get a copy of the original distribution."* It covers data as
well as code — *"This software and related material (data and (or)
documentation), contained in or furnished in connection with PHREEQC"* — which
is what makes this database file redistributable at all.

**No modification has been made.** The notice requires a prominent statement of
any change, its author and its date; there is none to state, and the md5 above
is the check.

## Why this version, and not the newest

The oracle generators run against the IPhreeqc engine bundled with
`phreeqpython`, which is **3.7.3-15968**. A database and an engine are a matched
pair: `phreeqc.dat` from PHREEQC 3.8 or 3.9 **fails to load** in this engine —

```
ERROR: Equation has no equal sign.  MEAN_GAMMAS
ERROR: Equation has no equal sign.  H2O(g)CO2(g)0.19
```

— because later editions use keywords it does not know. Tag `v3.7.3` is the
database of exactly that engine's release, so the pair is coherent rather than
merely working.

## Why not the copy `phreeqpython` bundles

That was the obvious choice and it is the wrong one. Its md5 is
`fca384eb269fda632fcefa6f74c8f90e` and **it matches no upstream tag** — not
v3.6.1, v3.6.2, v3.7.0, v3.7.1 or v3.7.3, all of which were checked. It is an
earlier edition, possibly modified, and nothing beside it says which. The USGS
notice requires a redistributor either to ship an unmodified file or to state
precisely what was changed, by whom and when. Neither is possible for a file
whose provenance cannot be established, so it is not the one that is shipped
here.

## Citation

Attribution is a condition of the notice, not a courtesy:

> Parkhurst, D.L., and Appelo, C.A.J., 2013, *Description of input and examples
> for PHREEQC version 3 — A computer program for speciation, batch-reaction,
> one-dimensional transport, and inverse geochemical calculations*: U.S.
> Geological Survey Techniques and Methods, book 6, chap. A43, 497 p.
> <https://doi.org/10.3133/tm6A43>

The pressure- and temperature-dependence data this particular database carries
are from Appelo, Parkhurst and Post (2014), *Geochimica et Cosmochimica Acta*
**125**, 49-67, as its own header states.
