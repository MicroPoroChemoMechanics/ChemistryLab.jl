# Published values

One JSON file per published source, named after the source's key in
[`docs/src/refs.bib`](../../docs/src/refs.bib). Each file holds the numbers the
package takes from that source (fitted parameters, measured values, tables,
compositions), with the place in the source they were read from, how they were
transcribed and checked, their unit and, when the source gives one, their
uncertainty.

Values from articles are kept here and nowhere else. They are not typed into
`src/`, `test/`, the documentation or `scripts/`: those read them with
`literature(key)`, `literature_value(key, name)` and `literature_table(key,
name)`. Coefficients that define an equation of state or a published model
(the water equation of state, the HKF constants, the Pitzer α and b) are part of
the model and stay in the code, with their source in a comment. Species and
phases with thermodynamic properties go into a dedicated database in the
ThermoFun format instead, as `cemdata18-zeolites.json` does.

## Format (`"schema": "chemistrylab-literature/1"`)

| field | content |
|:--|:--|
| `source.key` | the `refs.bib` key; equal to the file name |
| `source.doi` | the DOI, as in `refs.bib`, or `null` for a source without one |
| `source.citation`, `source.location` | a readable citation, and the table or page the values come from |
| `transcription` | where the values were transcribed from (`from`), whether they were checked against the source itself (`checked_against_source`) and how (`note`) |
| `quantities.<name>` | `value`, `unit` (unit arithmetic, `"1"` for a pure number), `kind` (`published`, `measured`, `fitted`, `estimated`, `placeholder` or `unstated`), `uncertainty` in the same unit or `null`, optional `location` and `description` |
| `tables.<name>` | `columns`, `units` (one per column, `null` for a text column), `kind`, `location`, `rows` |
| `notes` | free remarks |

The test suite (`test/literature.jl`) reads every file, checks it against this
format, and checks that its key is an entry of `refs.bib` whose DOI it repeats.
