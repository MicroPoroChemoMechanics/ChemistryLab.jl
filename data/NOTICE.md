# Data files and databases

## What this directory holds

Everything here is covered by the LGPL-2.1-or-later license of ChemistryLab.jl,
with one exception: the CSV files of `experimental/` come under their own terms,
stated in `experimental/LICENSE`.

- `solid_solutions.toml`, `gel_models.toml`, `pitzer-reardon1990.toml`: models
  and parameters, each entry naming its source.
- `literature/`: values transcribed from published papers, one file per source
  (see `literature/README.md`).
- `chloride/` and `zeolites/`: the data ChemistryLab adds to Cemdata18. That is
  the fitted chloride end member of CSHQ (`chloride/cshq_cl.json`) and the
  scripts that produce or check each extension.

## Thermodynamic databases published by others

ChemistryLab reads the thermodynamic databases of their publishers. It obtains
them the first time they are needed, checks the SHA-256 of each file against
the version it was validated with, and keeps them in a cache in the Julia
depot. `datapath("<name>")` returns the path of any of them. `database_info()`
says where each one currently resolves, and `fetch_databases()` obtains them all
in advance. The manual page *Databases* explains the procedure and every
message it can print.

| file | database | obtained from | license stated by the publisher | cite |
|:--|:--|:--|:--|:--|
| `cemdata18-thermofun.json` | Cemdata18, ThermoFun format | ThermoHub, release v1.1.1 (github.com/thermohub/thermohub, doi:10.5281/zenodo.7385311) | GPL-3.0 (repository) | Lothenbach et al. (2019), Cem. Concr. Res. 115, 472-506, doi:10.1016/j.cemconres.2018.04.018 |
| `psinagra-12-07-thermofun.json` | PSI/Nagra 12/07 | ThermoHub, release v1.1.1 | GPL-3.0 (repository) | Thoenen, Hummel, Berner & Curti (2014), PSI Bericht Nr. 14-04 |
| `aq17-thermofun.json` | aq17 | ThermoHub, release v1.1.1 | GPL-3.0 (repository) | Miron, Wagner, Kulik & Lothenbach (2017), Am. J. Sci. 317, 755-806, doi:10.2475/07.2017.01 |
| `slop98-inorganic-thermofun.json`, `slop98-organic-thermofun.json` | SUPCRT slop98 | ThermoHub, release v1.1.1 | GPL-3.0 (repository) | Johnson, Oelkers & Helgeson (1992), Comput. Geosci. 18, 899-947, doi:10.1016/0098-3004(92)90029-Q |
| `CEMDATA18-31-03-2022-phaseVol.dat` | Cemdata18, PHREEQC format | Empa, https://www.empa.ch/web/s308/thermodynamic-data, downloaded by hand and installed with `install_database` | the terms of Empa's download page | Lothenbach et al. (2019) |

## Databases ChemistryLab builds

These three are assembled on first use from the downloaded Cemdata18 file and
from ChemistryLab's own data in this directory. The entries of the Cemdata18
file are copied through unchanged, except `CaSiO3@` in `cemdata18-cashplus.json`,
whose standard properties are those refitted with the CASH+ model; the file lists
it under `replaced_substances`. They come under the license of the Cemdata18 file
they extend.

| file | what is added | from |
|:--|:--|:--|
| `cemdata18-zeolites.json` | 28 zeolites | Ma & Lothenbach (2020), doi:10.1016/j.cemconres.2020.106111; Ma & Lothenbach (2021), doi:10.1016/j.cemconres.2021.106537 |
| `cemdata18-chloride.json` | the chloride end member of CSHQ | fitted on Hirao et al. (2005), doi:10.3151/jact.3.77 |
| `cemdata18-cashplus.json` | the end members of the CASH+ and CASH+NK C-S-H, the aqueous species of the extended model, and `CaSiO3@` refitted | Kulik, Miron & Lothenbach (2022), doi:10.1016/j.cemconres.2021.106585; Miron, Kulik, Yan, Tits & Lothenbach (2022), doi:10.1016/j.cemconres.2021.106667; Miron, Kulik & Lothenbach (2022), doi:10.1617/s11527-022-02045-0 |
