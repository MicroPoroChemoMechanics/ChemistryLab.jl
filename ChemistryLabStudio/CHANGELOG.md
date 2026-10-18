# Changelog of ChemistryLabStudio

## v0.1.0 — ChemistryLab Studio, a graphical interface for cement pastes

The first version of ChemistryLab Studio, a local web interface to ChemistryLab
for cementitious materials, started with `clstudio()`. It is written for a
student who has never used Julia: the forms build the calculation, the results
are drawn and tabulated, and the Julia script of each calculation is written as
the forms change, so that the studio is also a way of learning the package.

### What it computes

- **One equilibrium** of a paste, certified, its degrees of reaction read at
  the age given.
- **Paths of equilibria**, each point started from the previous one: hydration
  in time, the replacement of the binder by a supplementary material,
  carbonation, a salt or any species added, leaching by renewals of pure water
  or of an attacking solution, ranges of temperatures and of water-to-binder
  ratios.
- **Kinetic runs**: the constituents given a rate law react by it (Parrott and
  Killoh for the clinker phases, Waller for a fly ash or a silica fume, a
  first-order law, or a Julia expression of the user's), the others are at
  equilibrium from the start; isothermal, or in a semi-adiabatic cell.
- **Constraints** on an equilibrium: a pH held, the activity of a species held
  (a partial pressure of CO₂), saturated curing, adiabatic conditions.

### What it is built from

- **Databases**: those ChemistryLab ships (Cemdata18 by default) or a file of
  the user's in any format ChemistryLab imports.
- **Phases**: a published phase list, with phases left out or added from a
  browser of the species of the database (search by name, formula or cement
  notation, filter by elements and state), or a list of the user's with the
  solid solutions it chooses; the model of the C-S-H gel.
- **Materials**: a published template, an oxide analysis (by Bogue's
  calculation, or kept as oxides) or a mineral composition, entered in a table
  that reads a block pasted from a spreadsheet or a CSV file; a law of the
  degree of reaction per constituent (constant, logistic, tabulated, Parrott
  and Killoh, capped).
- **Activity models**: the extended Debye-Hückel of Cemdata18, ideal dilute,
  Davies, Truesdell-Jones, B-dot, SIT and Pitzer.

### What it shows

- Phase masses and volume fractions (stacked along a path or in time), the
  pore solution (pH, ionic strength, dissolved elements), porosity and bound
  water, saturation indices, the certificate of the solver, the degrees of
  reaction and the heat of a kinetic run, a table to copy.
- **A diagnosis** before computing, from ChemistryLab's `diagnose_system`: a
  constituent the system cannot hold, a phase counted twice, data read outside
  their temperature range, missing volumes or enthalpies; and, in the page, a
  material whose minerals the system lacks, with a button that adds them.
- **The script**: seven commented sections, the studio running the first five
  itself so that the page shows what the script computes; saved, and reopened
  with every form restored from the model embedded in it. A script edited by
  hand is read back from its model and never overwritten.

### Around it

- A library of formulations published in papers, loaded in one click.
- English or French, a light or a dark theme, the expert fields hidden by
  default.
- The calculations run in a separate Julia process, in an environment the
  first start prepares (`@chemistrylab-studio`). A run is refused when the
  machine cannot spare the memory, and stopped when it outgrows what was left;
  *Stop* interrupts it, or restarts the process. The compute process ends with
  the studio, however the studio ends.
- The page answers only to the token of the session, from the studio's own
  address; the studio listens on the loopback interface, prints its address
  for a remote session, and takes the next free port when 8766 is taken.
