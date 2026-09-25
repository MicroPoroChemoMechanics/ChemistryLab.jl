# Databases

```@index
Pages = ["databases.md"]
```

```@autodocs
Modules = [ChemistryLab]
Pages = ["databases/paths.jl", "databases/thermofun_json.jl", "databases/phreeqc_dat.jl", "databases/merge_dat_json.jl"]
```

## Published values

The values taken from articles, one file per source under `data/literature/`;
see [Published values live in data files](@ref man-literature-data).

```@autodocs
Modules = [ChemistryLab]
Pages   = ["databases/literature.jl"]
```

## Published sorption models

Reading a sorption model written in the PHREEQC format — its site families, its
exchangers, and the `log K` of each reaction with the reference and the
uncertainty the compilation states beside it.

None ships here. ClaySor 2023, the model these were written against, is
CC-BY-4.0 and freely available from its Zenodo deposit; a published model is
also written against a particular **aqueous** database, and its constants are
that database's, so importing the reactions is not the same as being able to
reproduce the model.

```@autodocs
Modules = [ChemistryLab]
Pages   = ["databases/phreeqc_sorption.jl"]
```
