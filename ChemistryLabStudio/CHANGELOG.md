# Changelog of ChemistryLabStudio

## 0.1.0 (unreleased)

The first version of ChemistryLab Studio, a local web interface to ChemistryLab
for cementitious materials, started with `clstudio()`.

- One certified equilibrium of a mix of published materials, at chosen degrees
  of reaction, water-to-binder ratio and temperature, with the phase list and
  the activity model to choose.
- Results: phase assemblage, volume fractions, porosity, pore solution,
  saturation indices and the certificate of the solver.
- The Julia script of the calculation, written as the forms change, saved and
  reopened; it runs on its own in the studio's environment.
- A library of published formulations, loaded in one click.
- Interface in English or French, light or dark.
- The calculations run in a separate process that can be stopped, under a
  memory guard; the page is reachable only with the token of the session.
