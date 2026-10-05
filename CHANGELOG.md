# Changelog

## Unreleased

### Added

- **The two other C-S-H models in a coupled blended run**
  (`dw11k_setup(; gel)` in `scripts/deweerdt2011_kinetics.jl`, Section 7 of the
  page *CEM II/B-V and CEM II/B-M (V-LL), with their CEM I, integrated in
  time*): the fly-ash pastes of De Weerdt et al. (2011) integrated over six
  months with `CNASH_ss` and with `CASH+NK` in place of `CSHQ`, against the
  portlandite and ettringite of their Table 7 and the Ca/Si and Al/Si of their
  SEM-EDX analyses. In time as at the measured extents, neither fixes the
  pastes: `CASH+NK` gives the gel and the portlandite of `CSHQ` within a point,
  with no aluminum; `CNASH_ss` takes aluminum (Al/Si 0.09 to 0.11, 0.06 to 0.13
  measured) at a Ca/Si of 1.12 to 1.16 where the paste's gel falls from 1.7 to
  1.4, and leaves 1.7 to 3.9 points too much portlandite. `CSHQ` stays the gel
  of these pages, and the page says why.
- **Alkali uptake by C-A-S-H against L'Hôpital et al. (2016)** (the page
  *Alkali uptake by C-A-S-H*, `scripts/lhopital2016_alkali.jl`): their 49
  batches of C-S-H and C-A-S-H, Ca/Si 0.6 to 1.6 in water or in KOH or NaOH
  from 0.01 to 0.5 mol/L, transcribed whole from their Appendices A to C into
  `data/literature/LHopital2016b.json` (the key `LHopital2016` reads on, with a
  deprecation), and computed with CSHQ and with CASH+NK.
  The comparison is a prediction for CSHQ, whose alkali end members were fitted
  on other isotherms (Hong and Glasser 1999), and a check for CASH+NK, which
  was fitted on these. CSHQ takes up too little alkali where the gel is poor in
  calcium, as little as an eighth to a third of what was measured at Ca/Si 0.6
  to 1.0, reaching it only in 0.5 mol/L KOH: the opposite of what it does on
  the isotherms it was fitted on. It also dissolves too much silicon at a high
  pH, 2.5 to 3.6 times the measurement at the median, as the authors found with
  the same model. CASH+NK reproduces its fitting data: the alkali of the gel at
  0.62 to 1.55 times the measured one up to Ca/Si 1.2, the dissolved alkali and
  silicon at the median within 4 % below Ca/Si 1.1. In a gel poor in calcium,
  CASH+NK is the model whose alkali uptake was fitted there.
- **Aluminum uptake by C-S-H against L'Hôpital et al. (2016a)** (the page
  *Aluminum uptake by C-S-H*, `scripts/lhopital2016_aluminum.jl`): their 34
  syntheses without alkali, Ca/Si 0.6 to 1.6 and Al/Si 0 to 0.33, transcribed
  from their Appendices A, B and D into `data/literature/LHopital2016a.json`,
  and computed with CNASH_ss, the gel that takes aluminum, and with CSHQ.
  CNASH_ss reproduces the uptake: all the aluminum up to Al/Si 0.05, then an
  Al/Si held between 0.10 and 0.12 whatever the Ca/Si, the paper's finding
  (0.15 ± 0.05). It forms no strätlingite, where the syntheses hold 3 to 13 wt.%
  of it, putting the rest in gibbsite below Ca/Si 1.2 and in katoite above, and
  its gel stops at a Ca/Si of 1.19, the calcium it does not take precipitating
  as portlandite. CSHQ holds the calcium of the gel within 0.052 and puts all the
  aluminum in the hydrates from the smallest Al/Si. All 34 certify with both
  models; with CNASH_ss they need OptimaSolver 0.8.1.

### Changed

- **OptimaSolver 0.8.1 is required** (`[compat] OptimaSolver = "0.8.1"`). Its
  dual Newton left a member of a sublattice phase that owns no species on any
  site, and is barely unstable, at 1e-16 to 1e-28 mol instead of the floor, and
  the certificate then refused the answer: 14 of the 34 syntheses above did not
  certify with CNASH_ss, whose member 5CA is such a member. 0.8.1 holds it at
  the floor; every answer that certified before is the same.

### Fixed

- **A kinetic run whose partition holds a C-S-H under the compound energy
  formalism stopped before its first step.** The activities of the members of a
  `CompoundEnergyModel` phase, the CASH+ and CASH+NK gels of Kulik et al. (2022)
  and Miron et al. (2022a, b), depend on their standard Gibbs energies, which the
  equilibrium solvers hand to the activity model and the parameters of a kinetic
  run did not carry: the first evaluation of the activities raised "none were
  passed". The parameters of a run now hold those energies over RT when a phase
  of the system reads them, at the temperature of the run, or at the cell's own
  under a semi-adiabatic calorimeter. Runs on any other system are unchanged.

## v0.32.1 — Blended cements and their hydrates against temperature

Blended cements are now checked against what was measured on them at their
curing temperature, and their hydrates against temperature. The clinker and
slag laws meet the degrees of reaction of a slag-limestone cement cured at 5,
20 and 40 °C (Snellings et al. 2022): the slag's Waller law needs a ceiling
that rises with the water, and its activation energy, fitted, is the authors'
67 kJ/mol; the same cement integrated at each temperature gives its bound water
and portlandite within two points from four weeks on. The hydrates of two
cements from 0 to 60 °C (Lothenbach et al. 2008), the chloride AFm phases from
0 to 99 °C (Balonis 2019) and twelve equilibrium constants from 5 to 90 °C
(against the fits PHREEQC ships) are compared with published calculations and
measurements; the ternary cements of De Weerdt et al. (2011) are integrated in
time and their alite calibrated. Nothing the package computed before changes:
the release adds material templates, Rietveld phase names, published data and
pages.

### Added

- **The ternary cements of De Weerdt et al. (2011) integrated in time**
  (`scripts/deweerdt2011_kinetics.jl`, the page *CEM II/B-V and CEM II/B-M
  (V-LL), with their CEM I, integrated in time*): a clinker and gypsum with
  limestone powder, siliceous fly ash or both, the four clinker phases under
  Parrott–Killoh, the fly-ash glass at the degree of reaction the authors
  measured (their fit, written on the state as the Waller law is), against the
  clinker phases, portlandite and ettringite of their Table 7 over six months.
  The clinker as a whole reacts at about the measured rate, but the published
  parameters dissolve the alite too slowly and the belite too fast; the
  portlandite of the pastes without fly ash is within two points of the
  measurement; with fly ash the glass consumes it, and without limestone the
  ettringite is lost, CSHQ taking no aluminum. The fineness of Table 1 is
  transcribed. The page then fits the law of the alite on the plain cement,
  with the identifiability of its constants (exact sensitivity): two
  combinations are determined, those of the diffusion term (seven times the
  published constant); the interaction constant and the critical degree are
  only bounded from below. The limestone cement, not fitted on, follows within
  1.2 points; with fly ash the measured alite is faster still, the filler effect
  the law does not carry; and the belite, whose rate under this law is greatest
  at the mixing, cannot be held by any value of its constants.
- **The clinker and slag laws at 5, 20 and 40 °C** (`scripts/snellings2022_kinetics.jl`,
  the page *CEM I 52.5 R with slag and limestone at 5, 20 and 40 °C*), against
  the degrees of reaction Snellings et al. (2022) measured by X-ray diffraction
  on a 50:40:10 cement, slag and limestone blend at w/b 0.4 to 0.6 for six
  months (`data/literature/Snellings2022.json`: Tables 1 and 2 transcribed,
  Fig. 6 digitized from its raster image to 0.3 point). With its published
  constants and activation energies, Parrott–Killoh brings the clinker within
  four points of the measurement at three and six months, but it is 20 points
  short at one day at 20 and 40 °C. No Waller constant is published for a slag:
  the fly-ash shape, which runs to complete reaction, misses the slag by 15 to
  18 points, while one time, one exponent and one activation energy with a
  ceiling that rises with the water describe the 54 degrees to 2.6 points, all
  six constants determined. The activation energy, 67 kJ/mol, is the authors'
  own (their Table 2, by another law), below the 83 kJ/mol of the fly-ash set.
  On the two slags of the RILEM round robin (Durdziński et al. 2017), at an
  assumed 20 °C, the same constants hold within 7.2 points in three series of
  four. They are fitted on one slag and are not shipped.
- **The hydrates of two cements from 0 to 60 °C** (`scripts/lothenbach2008_temperature.jl`,
  the page *CEM I 52.5 N HTS and CEM II/A-L 42.5 R from 0 to 60 °C*), against
  the calculation of Lothenbach et al. (2008) with cemdata2007, whose Tables 1
  to 3 are now transcribed and whose Figs. 5 and 6 are read at 5 and 58 °C
  (`data/literature/Lothenbach2008.json`). With Cemdata18 and the package's
  Portland phase list, ettringite and monocarbonate give way to monosulfate at
  53.1 °C in the CEM I and 52.6 °C in the CEM II, where the article finds about
  48 °C and allows 42 to 54 °C for an uncertainty of 0.1 log units on its
  solubility products. Portlandite and calcite agree within 1.1 cm³ per 100 g of
  cement; the monocarbonate is a third to a half of the article's, whose AFm
  and AFt phases hold the iron that this list puts in an iron hydroxide.
- **The slag-limestone cement of Snellings et al. (2022) integrated in time at
  5, 20 and 40 °C** (`scripts/snellings2022_pastes.jl`, the page *CEM I 52.5 R
  with slag and limestone, integrated in time at 5, 20 and 40 °C*): its three
  materials as templates built from Table 1, the clinker under Parrott–Killoh,
  the slag glass under the Waller law fitted on its degree of reaction, the
  rest at equilibrium, against the bound water and portlandite of Fig. 8 and
  the hydrates of Fig. 10, digitized for w/b 0.5. The bound water is within two
  points from 28 days on and the portlandite within 1.8 at every age, though
  its fall with the temperature at six months is not reproduced. The ettringite
  is bounded by the sulfate of the cement, 9.9 g per 100 g of binder, where the
  diffraction finds 12 at 20 °C and 14 at 5 °C; the carboaluminate is the
  stable monocarbonate where the paste holds hemicarbonate; the hydrotalcite
  forms at once from magnesia released at the mixing.
- **Friedel's and Kuzel's salts from 0 to 99 °C** (`scripts/balonis2019_temperature.jl`,
  the page *Friedel's and Kuzel's salts from 0 to 99 °C*), against the four
  model mixtures of Balonis (2019): C3A, portlandite and water with sulfate,
  chloride and calcite. Her two salts are the Cemdata18 records (her Table 1,
  checked to the last digit). With chloride plentiful the calculation follows
  hers: Friedel's salt goes between 90 and 95 °C (80 and about 90 °C), the
  solids lose 22.8 and 26.0 % of their volume (23 and 25 %), monocarbonate
  goes between 50 and 55 °C (about 50). With chloride scarce, Kuzel's salt
  stays up to 99 °C where hers goes above 28 °C: of the same records, Kuzel's
  salt is stable against half monosulfate and half Friedel's salt by 1.6 kJ/mol
  at 25 °C and 0.3 at 99 °C, a margin her ideal solid solutions of Friedel's
  salt, which the package does not declare, close.
- **Equilibrium constants from 5 to 90 °C** (`scripts/logk_temperature_check.jl`,
  the page *Equilibrium constants from 5 to 90 °C*): twelve reactions of
  Cemdata18 (water, the carbonate system, bisulfate, CO2(g), calcite,
  aragonite, dolomite, gypsum, anhydrite, quartz, amorphous silica) against
  the analytical expressions of `phreeqc.dat`. Ten share their constant at
  25 °C within 0.01 (anhydrite 0.08 apart, quartz 0.23); their temperature
  dependences agree within 0.07 up to 50 °C, within 0.24 up to 90 °C, the
  largest difference being calcite's.
- The XRD-Rietveld names of a material template now include `Alite`, `Belite`,
  `C3A ortho` and `Aphthitalite` (by its oxides), as Table 1 of Snellings et al.
  (2022) prints them.

### Documentation

- **Which C-S-H for a fly-ash blend**, on the validation page of De Weerdt et
  al. (2011): the four pastes at measured extents with each of the three gels
  the package ships. `CASH+NK` behaves as `CSHQ` (no aluminum); `CNASH_ss` takes
  aluminum (Al/Si 0.11 against 0.13 measured) but at a Ca/Si of 1.16 in every
  paste, its end members reaching 1.5 at most where a Portland cement's gel is
  at 1.8, so that it puts the calcium in portlandite (31.7 against 21.8 wt.% in
  the CEM I). With fly ash, `CSHQ` and `CNASH_ss` bracket the measured
  portlandite, and neither keeps the ettringite of the paste without limestone:
  what is missing is a gel that takes aluminum at a Portland cement's Ca/Si.
- **The Guggenheim parameters and the temperature**, in the theory of solid
  solutions: the dimensionless parameters Cemdata18 prints for its AFt and AFm
  binaries become interaction energies in J/mol at 298.15 K and are held there
  at every temperature, which the source does not specify. Said, and what it
  does shown: the gap of the AFm SO₄/OH binary narrows from 0.48–0.98 at 5 °C to
  0.57–0.94 at 80 °C, where holding the dimensionless parameters would keep it
  at the printed 0.50–0.97.
- Two more pages built a dictionary of states by a comprehension, which Julia
  specializes on the type of a state: the validation page of De Weerdt et al.
  (2011) and the carbonation page.

## v0.32.0 — The calorimeter's energy balance inside the ODE, and blended cements in time

The semi-adiabatic calorimeter under partial equilibrium now integrates as Leal
et al. (2015) write the amounts: the right-hand side is a function of the
state, the temperature is the root of the cell's energy balance solved with
the partition at every evaluation, and the Jacobian is exact, so any
integrator applies. Blended cements are integrated in time for the first time,
slag (Gruyaert et al. 2010) and slag, fly ash and limestone together (Schöler
et al. 2015), and pore solutions are checked from 7 to 80 °C (Deschner et al.
2013). On the 32 cement pastes of a thesis the answers are those of 0.31.1, in
three quarters of the time.

### Breaking changes

- **The compatibility bound.** Below 1.0 a minor release is breaking for the
  registry: a package bounding ChemistryLab at `"0.31"` does not accept 0.32 and
  has to widen its bound.
- **Under partial equilibrium, the last entry of the state of a calorimeter's run
  is the change of the enthalpy of the cell**, not the temperature of a
  semi-adiabatic cell nor the heat of an isothermal one. Read the temperature
  with `temperature_profile(sol, cal; times)` and the heat with
  `cumulative_heat` or `heat_release`; `sol(t)[end]` no longer gives either. The
  stoichiometric formulation, without an equilibrium solver, is unchanged.
- **`WALLER_PARAMS_SLAG` is removed.** Its characteristic time of 100 days was
  attributed to Waller (1999) since the law was added, and recorded as unstated
  in 0.29. Read on 2026-10-04, the thesis does not contain it: it names a slag
  once, p. 219, among the additions an adiabatic test may contain, and gives no
  kinetic parameter for one. A slag's time is now the caller's, from a source of
  their own: `waller(merge(WALLER_PARAMS_FLY_ASH, (τ = τ_slag,)), "GGBS")`.

### Added

- **The dissolution mechanisms of Palandri and Kharaka (2004) for 34 minerals**:
  `palandri_kharaka(mineral; mechanisms, pco2, assume_Ea)` builds the acid,
  neutral, base and carbonate mechanisms the report tabulates, ready for
  `transition_state`, and `palandri_kharaka_minerals()` lists them. Tables 4,
  6, 26, 31, 32 and 34 of the report (quartz and amorphous silica,
  pyroxenoids, oxides, hydroxides, sulfates) are transcribed beside the
  carbonates of Table 33, read from the USGS PDF. A mechanism the report does not
  give is absent, and an activation energy it does not give (the neutral
  mechanism of gypsum, printed 0 where its text says the data could not
  determine it) is refused unless the caller states one. The report's base
  mechanism of quartz prints a pre-exponential factor of 10 and a log k that
  follows from the 491 it was adjusted from; both are stored as printed and the
  note says which.
- **A page of theory on recipes**, *What a recipe puts into the equilibrium*:
  the budget, what is kept aside, the mass balance with what an analysis does not
  report, the heat, and what a kinetic constituent changes in them.
- `glass_species(constituent, system; symbol, M)` gives a glass known by its
  oxides the formula a rate law can dissolve, from the oxides the system can
  hold, and `with_species(material, Dict(name => species))` puts it in its place:
  the glass of a template can now be integrated in time.
- A cell left empty by a source (`n.d.`, a dash) is `null` in a literature file
  and reads as `missing`.
- **Fe-Friedel's salt**, Ca₄Fe₂Cl₂(OH)₁₂·4H₂O, in the chloride extension
  `cemdata18-chloride.json`, with `Friedel_AlFe`, its ideal solid solution with
  Friedel's salt: the Cemdata18 paper tabulates it (Tables 1 and 2), its
  ThermoFun export does not carry it, and the iron of a slag or a fly ash could
  not bind chloride. Its record is the row of Table 1, transcribed in
  `data/literature/Lothenbach2019.json`, and the build refuses to write it unless
  its log Ks0 recomputed through the aqueous species of Cemdata18 is the one of
  Table 2 (it is within 0.011). The nitrite AFm, whose solid ships and whose
  NO₂⁻ does not, is checked through the NO₂⁻ of slop98 (within 0.001): the two
  rows of Table 2 the shipped file could not be checked on are checked now.
- **A blended cement integrated in time**, for the first time in the package: the
  CEM I 52.5 N and slag pastes of Gruyaert et al. (2010), from the mixing, the
  clinker under Parrott–Killoh and the slag glass under the Waller law, through
  `KineticsProblem(recipe, …)`, with the certified replay
  (`scripts/gruyaert2010_kinetics.jl`, the page *CEM I 52.5 N and slag pastes,
  integrated in time*). One number is fitted, the slag's characteristic time,
  on one measurement. Unfitted, the heat of the plain paste at 2 days is the
  measured one to 0.1 %, and the cement of the blends reaches the 94 % the image
  analysis found, more than in the plain paste, through the water/cement factor
  below. Two measurements are missed, and the page says what they need: the
  slag of a paste that is 85 % slag, which the portlandite stops activating,
  and the slow start of its cement.
- **Quaternary cements integrated in time**: the ten pastes of Schöler et al.
  (2015), a CEM I 52.5 R with blast-furnace slag, siliceous fly ash and
  limestone, from the mixing to six months (`scripts/scholer2015_kinetics.jl`,
  the page *CEM I 52.5 R with slag, fly ash and limestone, integrated in time*):
  the clinker under Parrott–Killoh, each glass under the Waller law with its
  time set on the degree the authors assume after a year, the limestone and the
  sulfates at equilibrium, against their thermogravimetry (Table 8), which
  nothing was fitted to. The bound water follows it within one to three points
  of the dry mass. The portlandite does not: the glasses consume it at
  equilibrium where the pastes keep 12 %, as in the authors' own calculation
  (Table 7, transcribed with the rest of the paper).
- The apparent activation energies of the cement and the slag that Gruyaert et
  al. (2010) fit against the cement-to-binder ratio (Eqs. 2 and 3), in
  `data/literature/Gruyaert2010.json`; the slag's rate of the page carries the
  latter.
- **The water/cement factor of Parrott and Killoh in the rate law**:
  `parrott_killoh_avrami(…; w_c, H)` and `pk_wc_factor`, the slowdown of a
  clinker phase once its degree passes `1.333 w/c`, or `H w/c` with the critical
  degree per phase of Lothenbach et al. (2008). It existed only in
  `ParrottKillohExtent`, an imposed extent, with its constants written in the
  code; both now read them from the papers' records in `data/literature/`.

### Fixed

- **A species without a rate law and without a standard Gibbs energy was put in
  the equilibrium partition**, where the minimization cannot price it: a glass
  known only by its formula, in a mix that did not hold it, made every
  re-speciation of the run fail, each step kept a frozen composition,
  and the run ended `Success` with one warning. `KineticsProblem` now refuses
  such a species by name when it is given an equilibrium solver.
- **Every kinetic run under partial equilibrium printed a warning of SciMLBase**,
  that parameters held in arrays of different types "can hurt performance": the
  system of the partition, a vector of species, was among the parameters of the
  ODE problem. It is held in a reference, read only to build states, and the
  warning is gone from every page that integrates.
- **`transition_state` read a catalyst the system did not hold as an activity of
  one**, without a word: a mechanism catalyzed by a species absent from the
  system ran at its rate constant. Such a catalyst is now refused at
  construction, with its name.
- **The templates of a Rietveld analysis gave two constituents to one species**
  where the analysis tells polymorphs apart (α′ and β C₂S, cubic and
  orthorhombic C₃A): two rate laws would have dissolved the species twice. They
  are summed into one constituent named by the database symbol, and
  `KineticsProblem` of a recipe refuses two kinetic constituents of one species.

### Fixed: the energy balance of a calorimeter under partial equilibrium

- **The temperature of a semi-adiabatic cell under partial equilibrium depended
  on the history of the integration, not on its state.** Its rate followed the
  equilibrium partition through derivatives taken at the last accepted step,
  frozen until the next one, and the heat of what they had not predicted was
  added to the temperature as a jump at each accepted step, outside the
  integrator's error control; the heat capacity of the shift of the partition
  was clamped at zero where it came out negative. The right-hand side is now a
  function of the state, as Leal et al. (2015) write it for the amounts: the
  state carries the change of the enthalpy of the cell, which only the losses
  through its walls move, and the temperature is the root of the cell's energy
  balance, solved with the partition at every evaluation. Its derivatives are
  first derivatives of the minimization, which the implicit-function theorem
  gives exactly, so the Jacobian is exact and any integrator applies: a
  Rosenbrock method, a BDF method and an explicit Runge–Kutta method now give
  the same temperature at two days to a few microkelvin. On the alite paste of the tests, the
  enthalpy of an adiabatic cell is conserved to 6e-9 of the heat released
  (1e-3 was allowed before), and the heat capacity the Jacobian implies equals
  the one of the certified equilibrium differentiated in temperature to every
  printed digit. The theory page *Kinetics under partial equilibrium* writes
  the formulation.
- **The heat of an isothermal cell under partial equilibrium** is the enthalpy
  the paste has lost, `H₀ − H`, at the partition of each state, and no longer
  the integral of a linearized rate corrected by jumps: it equals the heat of the
  certified replay (`heat_release`) to the third decimal of a joule, where it was
  10 J off out of 2 kJ near the first hour.
- In a semi-adiabatic run, the log-activities the rate laws read were evaluated
  at the initial temperature; they are now evaluated at the cell's.
- **The partition solved in the right-hand side started from the last accepted
  one lifted to 1e-10 mol**, the floor of the interior point, which moves the
  potentials of the traces: on the C100 mortar at 1.6 h the dual Newton failed
  from there and the run stopped, where the partition itself, floored at 1e-16
  mol, certified at once. It now starts there, the lifted start second. With the
  root of the energy balance kept for the passes of a Jacobian, the five days
  of that mortar in its semi-adiabatic cell take 220 s.

### Documentation

- **Validation on pore solutions from 7 to 80 °C**, the pastes of Deschner et al.
  (2013), a CEM I with half quartz powder or siliceous fly ash, cured at 7, 23,
  40, 50 and 80 °C (`data/literature/Deschner2013.json`, Tables A.1, A.2 and
  B.1): every solution speciated at its hydroxide and its temperature, the
  effective saturation indices of portlandite, ettringite, monosulfate and
  strätlingite against the authors'. From 7 to 50 °C they differ by 0.08 at
  most; at 80 °C the fly-ash paste, the richest in sulfate, differs by up to
  0.41, which the page states without settling the cause.
- **Six pages compiled a dictionary for minutes.** A dictionary built by a
  comprehension over states or runs is specialized on their types, long enough
  that compiling it took 402 of the 413 s of a block that solved 54 pore
  solutions. The pages build them with `Any` values and say why.

### Changed: the inversion of the aqueous phase three times lighter

- **The inversion of the aqueous phase allocated and dispatched at run time in
  its innermost loop.** The function whose root is the ionic strength, evaluated
  some twenty times per inversion and an inversion per trial step of the dual
  Newton, assigned the ionic strength to a name the enclosing function also
  assigned, and Julia boxed that variable: every evaluation went through a
  `Core.Box`. Measured on the C100 mortar of Lavergne et al. (2018), from the
  same start, a certified solve that needs many inversions took 150 ms and
  328 MB and takes 54 ms and 99 MB; one that converges at once allocates 4.9 MB
  instead of 12.8, in the same time. The answers are the same to the last bit.
  The same defect is removed from the Newton inversion of the SIT and Pitzer
  models, from the Donnan potential of a diffuse layer and from the implicit
  kinetic step.

### Changed: kinetic runs twice as fast

- **Every evaluation of an activity model in a kinetic run computed a number
  type by walking the run's parameters at run time.** The helper that decides
  whether the parameters carry dual numbers was written to be resolved by the
  compiler, but it mapped over the values of the tuple it was given; a kinetic
  run hands the model its own parameters, some fifty fields, a `map` over a
  tuple that long is not unrolled, and each field was dispatched at every
  evaluation. On the pore-humidity test that was 44 % of the integration. The
  type is now computed from the types alone, folded into a constant: that
  integration takes 76 s instead of 148 s (130 s with 0.30.0), the trajectory
  unchanged. Present since 0.29.0, in every activity model.
- **The certified replay of a kinetic run returned compositions it had not
  certified, on its first instant above all.** `speciated_states`, and
  `heat_release` through it, certify each instant from three starts and, from
  the second instant on, by a continuation from the last certified one; the
  first has no neighbor to walk from. On a slag paste at 2 and 7 days (Gruyaert
  et al. 2010) nothing certified, and the interior-point composition the
  replay fell back to had a pH of 15.3. The full certified search, with its
  start from the linear program, is now the last resort: it certifies those
  instants at once, at pH 12.82.
- The searches that ask a back end only for a start pass it `polish = false`,
  a new keyword of `solve(::EquilibriumSolver, state)`, instead of running the
  solve under a scoped value: another 1.4 s of compilation on the first
  equilibrium of a session.

## v0.31.1 — The first equilibrium of a session compiles as fast as with 0.30

A patch release: the answers are those of 0.31.0, the same on the 32 cement
pastes of a thesis, and the first certified equilibrium of a session compiles
again within 2 s of the time 0.30 took. Measured on the same machine, the 32
pastes take 396 s with 0.30.0, 421 s with 0.31.0 and 398 s with 0.31.1.

### Fixed

- **The first equilibrium of a session compiled 13 s longer than with 0.30.**
  0.31.0 ran an explicit back end whose answer was to be polished under a
  scoped value that suspended the strict convergence check, and inferring the
  solve through it cost 13 s of compilation on the first certified equilibrium
  of a CEM I paste: 125 s against 112 s, where 0.30.0 took 109 s. The check is
  now simply not made when the polish decides on the answer.
- With it goes a warning that misled: a back end stopping short, at `MaxIters`
  for instance, was reported as not converged and counted in `NONCONVERGED`,
  although the answer returned was the one the polish then certified. A polish
  that fails is reported, as before.

## v0.31.0 — One equilibrium whatever the back end, a scope that says what is proved, and kinetics that read the speciation

An external audit of 0.29.0 reported seven defects, and the code confirmed
each: back ends that computed another composition than the equilibrium, and
derivatives of a map they did not return; a saddle certified a global minimum;
an ODE trajectory that created matter and reported success; a gas whose state
could not be built; pressure absent from the gases and from the condensed
phases; and a validation page that read a reference temperature as an
inconsistent entropy. This release corrects them, makes the ideal model and
Davies derive from one Gibbs energy, and completes the Pitzer model and the
heat capacities the audit listed as limits.

### Breaking changes

- **The compatibility bound.** Below 1.0 a minor release is breaking for the
  registry: a package bounding ChemistryLab at `"0.30"` does not accept 0.31 and
  has to widen its bound.
- **A back end's answer is polished.** `equilibrate(state, solver)`,
  `equilibrate(state; certify = false)` and `solve(::EquilibriumSolver, state)`
  return the composition the dual Newton certifies from the back end's answer
  when OptimaSolver is loaded and the system has an aqueous phase with `H2O@`.
  Where that answer was not at the equilibrium, the numbers move to it.
- **Without OptimaSolver, a back end minimizing `n⋅μ(n)` refuses** an activity
  model that breaks the Gibbs–Duhem relation (the B-dot model, Davies with a
  neutral solute and `bₙ ≠ 0`, Truesdell–Jones, SIT), with an `ArgumentError`.
- **`certify` defaults to `nothing`**, which behaves as `true` did; an explicit
  `certify = true` now raises where the certified search does not apply.
- **The solvent rows of `DiluteSolutionModel` and `DaviesActivityModel`** are
  the Gibbs–Duhem partners of their solutes' terms: `ln a_w` moves at second
  order for the first, and by the osmotic correction for the second.
- **Scopes.** `:local_minimum` is a new value of `certificate.scope`, and
  `:global_minimum` needs a convexity proved over the whole domain: an ideal
  answer now has it, a Pitzer or SIT answer never does, and one under the
  B-dot or Davies model only with one ion size, no linear or salting-out term,
  and under the Debye–Hückel bound. The certificate has a field
  `reduced_curvature`.
- **Pressure.** A gas's activity carries `ln(P/P°)` and its molar volume is
  `RT/P`; a species declared at constant volume carries `V⁰ (P − P°)` in its
  `ΔₐG⁰` and `ΔₐH⁰`, which are then `NumericFunc`s of `T` and `P`. Nothing moves
  at 1 bar.
- **Heat capacity on several intervals.** The functions of such a species are
  piecewise `NumericFunc`s; they do not move inside the interval of the
  reference temperature.
- **Pitzer.** The higher-order electrostatic terms are on by default
  (`etheta = false` restores the previous model); a mixture of unlike charges
  moves. `PitzerParameters` has a field `temperature` and
  `PitzerActivityModel` a field `etheta`; code calling their positional inner
  constructors has to pass them.
- **Kinetics.** A rate law that reads the speciation is integrated with the
  partition solved in the right-hand side; its trajectory changes, to the right
  one. A trajectory that reaches kinetic amounts the system cannot hold is
  returned with `retcode = Unstable` instead of `Success`.

### Documentation

- The README described the state of 0.18. It called Ipopt the default backend
  and OptimaSolver optional, when OptimaSolver takes over whenever it is loaded
  and is what certifies; it called the certificate a proof of the global
  minimum, when the certificate states its scope; it said Cemdata18 was shipped
  with the package, when `datapath` obtains it on first use; and it listed 12 of
  the 21 phases of `data/solid_solutions.toml`. The feature list now names the
  activity models, solid-solution models, surfaces, recipes and derivatives that
  came since.
- The home page and Getting started repeated that the database is distributed
  with the package, and *Database Interoperability* counted two built databases
  where there are three, `cemdata18-cashplus.json` included.
- `CITATION.cff` and `.zenodo.json` carry the same abstract again, brought up to
  date, with the keyword "surface complexation".
- The Pitzer pages said that no thermodynamic database ships Pitzer parameters.
  PHREEQC distributes a set in its `pitzer.dat`, and the package itself ships
  the cement set of Reardon (1990); what the ThermoFun databases lack is said
  of them alone.
- A recipe constituent known by its oxides and given a rate is refused, as
  before; the message now names the way round it, a pseudo-species built by
  `glass_species` and declared as a mineral constituent.
- **The published-data page blamed an estimated entropy for what is a reference
  temperature.** Eight CEMDATA18 records, the alkali C-S-H and M-S-H end
  members, are tabulated at 293.15 K, and `ΔₐG⁰` is anchored at that
  temperature, as GEMS anchors them. The page and its test compared their 20 °C
  tabulated energy with the package's 25 °C one and read the 5 K step as an
  inconsistency in `S°`, and the solubility-product check put the same 20 °C
  energy into a 25 °C constant, which doubled the M-S-H offsets. Each record is
  now compared at its own temperature; the M-S-H end members miss Table 2 by
  `0.244` and `0.205` rather than `0.48` and `0.40`. The genuine data check,
  `ΔfG° = ΔfH° − Tst (S° − Σ S°el)`, is added: the alkali C-S-H records satisfy
  it exactly at 293.15 K, the M-S-H records at neither temperature. The page's
  counts (220 of 228, 52 phases) were those of an older file, and its test now
  checks the numbers the page prints.
- The derivative example of *Solving an equilibrium* printed `0.15193` for a
  system the page did not define. The system is written out, and the value is
  the one it gives, `0.163095`, through an explicit back end and the certified
  route alike. The explicit row of the calcite table of *Writing a kinetic
  model* was measured again (84 586 steps, 250 s).
- The w/c page said an Ipopt answer comes without a certificate and leaves the
  absent phases at its lower bound; with OptimaSolver loaded it is polished and
  certified, and the page says so.

### Changed: one equilibrium, whatever the back end

- **Ipopt computed another composition than the equilibrium.** It minimizes
  `n⋅μ(n)`, whose gradient is `μ + Jᵀn`, and `Jᵀn = 0` is the Gibbs–Duhem
  relation, which the B-dot model, Davies with a neutral solute and SIT do not
  satisfy: the minimum of `n⋅μ(n)` is then not where `μ(n) = −Aᵀy`, the
  conditions the dual Newton solves and the certificate audits. On calcite and
  carbon dioxide in a sodium chloride solution under Davies, Ipopt's dissolved
  calcium was 7e-4 away from the equilibrium's, against 2e-6 under the ideal
  model, while the element balance held to 1e-15 mol either way. The
  logarithmic route of OptimaSolver differentiated the same scalar.
- **Every back end's answer is now polished by the dual Newton** when
  OptimaSolver is loaded and the system has an aqueous phase with `H2O@`: the
  dual Newton is started from it, and the composition it certifies is returned,
  by `equilibrate(state, solver)`, by `equilibrate(state; certify = false)` and by
  `solve(::EquilibriumSolver, state)` alike. The logarithmic route of
  OptimaSolver is handed the gradient `n ∘ μ`, as the linear one is handed `μ`.
  Where a back end only supplies a start, to the certified search, the homotopy,
  a kinetic run or its replay, its answer is not polished twice.
- **The derivatives through a back end were those of another map.** They were
  lifted with the conditions of the dual Newton at an answer that did not
  satisfy them; they are now lifted at the polished answer, which does, and equal
  those of the certified route.
- **Without OptimaSolver**, nothing can polish, and a back end that minimizes
  `n⋅μ(n)` refuses an activity model that breaks the Gibbs–Duhem relation, with
  an error that names the species where it fails.
- **`certify` has three values.** `nothing`, the default, certifies where the
  certified search applies and takes the single back end elsewhere, as `true`
  did; an explicit `true` now refuses a system where it does not apply rather
  than return an answer that was never certified; `false` takes the single back
  end, polished. `equilibrate(…; certificate = Ref{Any}())` receives the
  certificate of the answer returned, or `nothing` when none was computed.

### Changed: Gibbs–Duhem exact for the ideal model and for Davies

- **The ideal model broke the Gibbs–Duhem relation on its solvent.** Its solutes
  are `ln mᵢ`, and their partner is `ln a_w = −M_w Σ m`, where it took Raoult's
  mole fraction `ln x_w`, off by `1 − x_w`: its activities were not the
  gradient of a Gibbs energy, the certificate of an ideal equilibrium was
  `:self_consistent`, and Ipopt alone minimized another function than the
  energy. The solvent row is now `−M_w Σ m`; `ln a_w` moves by `(M_w Σ m)²/2`,
  2e-8 on the calcite reference, where the comparison with Reaktoro is unchanged.
- **So did Davies.** Its ions carry one function of `I` times `zᵢ²`, whose
  Gibbs–Duhem partner for the solvent has a closed form, now used in place of
  Raoult's: with ions alone, or with `bₙ = 0`, Davies derives from one Gibbs
  energy. A neutral species with `bₙ ≠ 0` still breaks the symmetry, and says
  so. The page comparing the activity models measured the Gibbs–Duhem residual
  by finite differences and concluded that Davies was less consistent than the
  ideal model; it is measured by automatic differentiation, and in a single
  salt every model now satisfies the relation to rounding.

### Changed: what a certificate proves

- **A saddle could be certified a global minimum.** The scope rested on the
  symmetry of the Jacobian alone, which says that an energy exists, not that it
  is convex: a synthetic Pitzer set with `β⁽⁰⁾ = −1`, under which dissolving
  halite into a sodium chloride solution lowers the Gibbs energy, was scoped
  `:global_minimum`, and so was HKF with an ion size of 1 Å.
- `:global_minimum` now rests on a convexity proved over the whole domain: ideal
  mixing, convex solid solutions, ideal site mixing, the ideal dilute model, and
  a Debye–Hückel form with one function of `I` under the bound
  `ln(10)·A·z_max²/(8·B·å) ≤ 1` (0.45 for divalent ions of 4 Å at 25 °C), whose
  derivation is on the page of the certificate. Elsewhere the Hessian of the
  energy over the directions that conserve matter decides: positive definite, the
  new scope `:local_minimum`; otherwise `:kkt_point`, with the reason, a saddle
  or a flat direction. The certificate reports its smallest eigenvalue as
  `reduced_curvature`, and documents `stationarity_abs`.

### Changed: a rate law that reads the speciation

- **The ODE route returned an impossible trajectory as a success.** With the
  equilibrium partition frozen within a step, a rate law reading it (a
  saturation ratio, an activity) is constant over the step, so the stiff method
  integrated the extent explicitly: on calcite under `r = k(1 − Ω)`, `Rodas5P`
  stepped past the second or so over which `Ω` relaxes, `Ω` then exceeded one by
  orders of magnitude, and the run ended on hundreds of moles of calcite from
  0.05 with `retcode = Success` and a warning. With the partition solved in the
  right-hand side, the same run takes 82 steps and ends at the equilibrium.
  The source said in one place that a missing Jacobian was the cause and in
  another that it was not.
- `integrate` now reads the rate laws (`_rates_read_speciation`): when one reads
  the partition, the right-hand side solves it at the state it is evaluated at,
  by the certified solver warm-started from the last accepted step, and lifts its
  derivative with respect to `bₑ` into the Jacobian by the implicit-function
  theorem; the calcite case reaches the equilibrium. A law that reads only the
  kinetic amounts keeps the split route, which is exact for it and unchanged.
  `integrate(…; speciation = :rhs | :frozen)` forces either.
- **Feasibility is checked on the whole trajectory**, against the element totals
  of the system, each kinetic amount at its own scale, rather than on the final
  state against twice the total amount of matter, which the water dominated. In
  `:rhs` mode a step leaving those bounds is rejected; in any mode a trajectory
  that reaches them is returned with `retcode = Unstable`, or raises under
  `STRICT_CONVERGENCE`.
- The step callback declares the integrator's derivative information stale
  whenever the re-speciation changes what the right-hand side reads (the heat of
  a calorimeter, a law run frozen), where it said nothing had changed.

### Added: the rest of the Pitzer model, and heat capacities past a transition

- **The higher-order electrostatic terms of Pitzer (1975)** for two ions of like
  sign and unlike charge, Na⁺ with Ca²⁺, Cl⁻ with SO₄²⁻, which the model left out
  and its docstring said so. `J(x)` is evaluated by the Chebyshev approximation
  of Harvie that PHREEQC uses, checked against the values Reaktoro tabulates
  independently (4e-8); the terms enter `γ` and the osmotic coefficient from one
  excess energy, so the model keeps satisfying the Gibbs–Duhem relation exactly.
  `PitzerActivityModel(…; etheta = false)` leaves them out, for a set fitted
  without them. A single salt is unchanged to the last bit; a mixture of unlike
  charges moves.
- **Temperature terms of the Pitzer coefficients**, in PHREEQC's six-term form,
  used under `temperature_dependent = true`, zero at 298.15 K exactly;
  `PitzerParameters(…; temperature)` and the TOML reader take them.
- **A reader of the `PITZER` block of a PHREEQC database**,
  `build_pitzer_parameters(path; format = :phreeqc)`, for the file the caller
  has: none is shipped. Neutral species take the `@` this package names them
  with; the identifiers the model has no counterpart for are skipped and named.
- **A heat capacity given on several intervals is followed past the first.**
  ThermoFun records give one polynomial per interval and the transitions between
  them; only the interval holding the reference temperature was kept, so that
  hematite past its transition at 950 K extrapolated the wrong polynomial. Each
  interval is now anchored on the one before, the transition adding its
  enthalpy and entropy; inside the reference interval nothing moves.

### Changed: pressure enters the gases and the condensed phases

- **A gas's activity ignored the pressure.** Every activity model gave a gas
  `ln a = ln xᵢ`, the mole fraction, which is its activity at 1 bar only: a gas
  over water dissolved the same amount at 1 and at 10 bar, and its chemical
  potential did not grow with pressure. It is now `ln xᵢ + ln(P/P°)`, its
  fugacity over the standard pressure, in Dilute, HKF, Davies, Truesdell–Jones,
  SIT and Pitzer alike, so that `∂μᵢ/∂P = RT/P`. `P_STANDARD` and
  `P_STANDARD_Q` (1 bar) are exported.
- **A gas has the ideal gas's volume, `RT/P`.** The databases' gases carried a
  constant 24.79 L/mol, their volume at 298.15 K and 1 bar, at every
  temperature and pressure; their records declare the ideal gas (`mv_pvnrt`),
  and that is what they now get. With `ln(P/P°)` in the activity, a gas's
  chemical potential and its volume satisfy the Maxwell relation, which a
  volume constraint (`FixedVolume`, `SealedVolume`) needs.
- **A condensed species declared at constant volume had no pressure in its
  standard energy.** A record declaring `mv_constant` (the solids, the solutes
  outside HKF) now carries `V⁰ (P − P°)` in its `ΔₐG⁰` and `ΔₐH⁰`, so that
  `∂G⁰/∂P = V⁰`; so does the solvent, whose equation of state is not
  implemented, with its compressibility neglected. Three crystals CEMDATA18
  marks as ideal gases (`CA`, `CA2`, `C12A7`) are treated as the crystals they
  are. At 1 bar nothing changes: the term is an exact zero there. Calcite
  under pressure is the visible case: its reaction volume, ions minus crystal,
  is now about −56 cm³/mol between 15 and 70 MPa at 301 K, and its `log K` at
  70 MPa agrees with Duan et al. (2016) to 0.01, where it was 0.44 short.

### Fixed

- **A gas built without a molar volume made its state impossible to build.**
  The ideal-gas fallback converted to `u"m^3"`, which DynamicQuantities refuses
  as a target, so `ChemicalState` threw for any system holding such a gas. A gas
  without `V⁰` now takes `RT/P` wherever a volume is read.

## v0.30.0 — A trace held to its own amount, and SIT and Pitzer solved by Newton's method

The certified equilibrium judged its element balance in moles, against `1e-10`:
a component of a nanomole could be 10 % wrong and certified. Each balance row is
now judged against what it holds, as PHREEQC and GEMS judge a mass balance
against its element total (OptimaSolver 0.8). And the aqueous solutes under SIT
and Pitzer, recovered one by one until now, are solved by Newton's method on the
model's own Jacobian; the first equilibria computed with these models are in the
test suite, the solubility of halite under Pitzer within 0.27 % of its
measurement.

### Breaking changes

- **The compatibility bound.** Below 1.0 a minor release is breaking for the
  registry: a package bounding ChemistryLab at `"0.29"` does not accept 0.30 and
  has to widen its bound.
- **A certificate is stricter.** `optimal` asks each balance row that has a
  budget to be met relative to what it holds, below one mole, and in moles above
  it as before. An answer 0.29 certified with a trace off by more than `1e-10` of
  itself is now solved further, or refused. `balance` stays in moles;
  `balance_relative` is new.
- **SIT at the activity floor.** The SIT closure took the logarithm of a molality
  regularized twice, `log(mᵢ + ϵ)` on amounts already floored, which put a
  species at the floor at `ln 2` above it; HKF, Davies and Pitzer had lost that
  offset in 0.28. The log-activities of floored species under SIT move by `ln 2`,
  and nothing else does.
- **`DecompositionWindow` has a new field, `shape`** (`:logistic` or
  `:interval`); code calling its positional inner constructor has to pass it.
  Its keyword constructors are unchanged.
- **`scripts/ionic_hydration.jl`** runs its report as `report_ionic_hydration()`,
  no longer `main()`: `scripts/hydration_calibration.jl` includes that file and
  has a `main` of its own, which replaced it.

### Changed: the balance judged against what each row holds

- `optimality_certificate` reports the balance twice: `balance`, the worst row in
  moles, and `balance_relative`, the worst row with a budget relative to what it
  holds. `optimal` follows OptimaSolver's `kkt_certificate`, which judges the
  larger of the two below one mole. A row whose budget is zero within rounding,
  such as the electron row of a redox pair, and a component nobody supplies (a
  budget below `1e-12` of the largest) are judged in moles.
- The ranking of candidates (`_kkt_error`), that of a kinetic step, and the rule
  that offers a phase its second instance read the balance as the certificate
  judges it, and the refusal messages print both figures.
- Measured on the same machine as 0.29.0, with OptimaSolver 0.8.0: the 32 cement
  pastes of a thesis in 388 s instead of 349 s, every printed value unchanged and
  their balance residuals down from `1e-11` to `1e-15`, the time added on the one
  paste whose nanomole of carbon was left `8e-6` of itself off. The three-hour
  hydration of `scripts/ionic_hydration.jl` takes 5.1 s instead of 5.8 s, on an
  identical trajectory, and its 28-day certified replay 2.6 s instead of 2.5 s.

### Changed: SIT and Pitzer, solved by Newton's method

- **The aqueous solutes under SIT and Pitzer are recovered by Newton's method**
  on their log-amounts, with the Jacobian of the model, exact by forward
  differentiation: the specific ion interaction terms `ε(i,k) mₖ` and the Pitzer
  pair and triplet sums depend on the molalities themselves, and no single
  equation gives the solutes back as it does for the Debye–Hückel family. The
  iteration starts from the better of the composition the solve holds and the
  one the model's Debye–Hückel part gives through the ionic strength, and its
  steps are capped at 30 in any log-amount and halved until the squared residual
  falls. A solute below the activity floor is placed from its activity
  coefficient at the composition found. No composition found is reported as
  such, which OptimaSolver handles as it handles the limiting law past its range.
- **The first equilibria with these models.** Until now no test, page or script
  solved one; their closures had only been evaluated at given compositions. The
  suite certifies a sodium chloride brine at 3 mol/kg under SIT, and halite in
  water under Pitzer: the saturated molality comes out at 6.1605 mol/kg against
  the 6.144 Hamer & Wu (1972) measured, from the standard Gibbs energy of
  halite in slop98 and the Na–Cl parameters of the Reardon set, neither fitted to
  it. The Pitzer page shows that equilibrium. On these systems the sweeps
  converged too, to the same compositions; the inversion is exact where they
  were not guaranteed to be.

### Added: thermograms as the cement literature reports them

- **A decomposition window given as a temperature interval**,
  `DecompositionWindow(phase; between = (T₁, T₂))`: the phase releases all of its
  water or carbon dioxide between the two temperatures and none outside, along a
  smooth step whose rate vanishes at both ends. That is how the papers attribute a
  loss to a phase (portlandite between 350 and 500 °C, De Weerdt et al. 2011;
  the carbonate from about 300 to 850 °C, Shi et al. 2016; the water of ettringite
  between 30 and 150 °C, Möschner et al. 2009), and the loss between the two ends
  is then exactly the phase's content, where a logistic window loses part of it
  outside. Both forms mix in one set, fit to a curve (`window_parameters` gives an
  interval's ends) and serve `bound_water` over a range. `window_interval`
  returns an interval's ends, or the 1 % and 99 % points of a logistic.
- **`thermogram(...; relative_to)`** returns the curve in percent of a reference
  mass, `mass_percent` and `loss_percent` (the loss counted from the first
  temperature of the grid), with `reference_mass`. The reference is the sample at
  the start (`:initial`, the default), at a temperature (a dry mass: 500 °C for
  Schöler et al. 2015, 800 °C for Shi et al. 2016), or ignited (`:ignited`; Shi &
  Lothenbach 2020 give bound water in percent of the sample ignited at 980 °C),
  computed at its own temperature. The thermogravimetry page sets these
  conventions out, each with its source.
- The intervals and reference temperatures are transcribed in `data/literature`
  (`DeWeerdt2011`, `Shi2016`, and the new `Moschner2009`, `Scholer2015`,
  `LHopital2016`, `ShiLothenbach2020`), with the page and the sentence they come
  from.

### Fixed

- The bibliography gave the first authors of Nied et al. (2016) as "David" Nied
  and "Eléonore" L'Hôpital; Crossref and the article give Dominik and Emilie.
- `SITActivityModel` declares its `concentration_scale`, molality, which the
  interface asks of every activity model and which it alone lacked.
- The last method redefinitions of the test suite are gone:
  `scripts/validation_common.jl`, included by four validation scripts, defines
  its helpers once per session, and the page extractor is included by its test in
  a module of its own. With `--warn-overwrite=yes`, as `Pkg.test` runs, 0.29.0
  still printed eight.

### Changed: requirements

- ChemistryLab requires **OptimaSolver 0.8** (`OptimaSolver = "0.8"`), whose
  certificate and dual Newton judge each balance row against what it holds.

## v0.29.0 — Derivatives through every forward model, and no difference quotient left

This release takes every derivative of the package by forward-mode
differentiation, and lets one flow through everything a forward model is built
from: an equilibrium, a kinetic run, a recipe, a surface. It is what an inversion
needs whose forward model is an equilibrium or a hydration, differentiated with
`ForwardDiff` and differentiated again for its Hessian. It also showed that the
difference quotients had been producing results of their own: a degeneracy half
hidden, a rate constant credited with an influence it does not have, a spectrum
flattened. And the certified solver now recovers the aqueous solutes through
the ionic strength, as one equation, where it had swept them one by one without
converging: the same answers, several times faster, and a model with no solution
at given potentials recognized as such.

### Breaking changes

- **The compatibility bound.** Below 1.0 a minor release is breaking for the
  registry: a package bounding ChemistryLab at `"0.28"` does not accept 0.29 and
  has to widen its bound.
- **Six exported types gained a type parameter**, so that they can hold dual
  numbers: `Recipe{R}`, `RecipeState{S, C, M, B, I}`,
  `MineralConstituent{S, E, F}`, `OxideConstituent{E, F}`,
  `ParrottKillohExtent{V}` and `DonnanLayer{T}`. Their constructors are
  unchanged; code that names their type parameters (an inner constructor called
  as `MineralConstituent{S, E}(…)`, a method signature listing them all) has to
  be updated.
- **`identifiability` reads its rank differently with an observation**: from the
  standard errors of the singular directions rather than from the largest gap
  of the spectrum (see below). The same inputs can return another `rank`, and
  `as_traced` can then mark another set of parameters as placeholders. Its
  `stderr` are larger by `√(n/(n − p))`, the noise level now being estimated on
  the degrees of freedom the parameters leave, and `Identifiability` has a
  `noise` field. The correlation matrix is formed from the singular value
  decomposition, which changes its sign at an exact degeneracy, where the
  inversion of `JᵀJ` had none to give.

### Changed: the aqueous solutes are solved through the ionic strength

- **One equation where the solver swept.** For the Debye–Hückel (`HKFActivityModel`,
  with or without an ion size), Davies, Truesdell–Jones and dilute models, the
  inner level of the certified solver recovers the solutes from their potentials
  through the ionic strength: every solute is explicit at a given `I`, and `I`
  solves one equation, as PHREEQC carries it, an unknown of its own. OptimaSolver
  had recovered them one by one, which cycles where multivalent ions couple
  through `I`: on cement pastes most inversions had ended unconverged, and every
  one under the limiting law past its range. The root taken is the first above
  the dilute limit, the branch connected to it, so that the answer does not
  depend on the path; a dip of the equation below zero narrower than a step is
  found from the sign change of its derivative. A dip that stays above zero ends
  the dilute branch, and so does 1e4 mol/kg: under the limiting law the root past
  such a dip lay at 6500 mol/kg for a paste loaded with sodium chloride, and a
  solve that started from it never recovered. No root says that the potentials
  hold no composition: a trial step from an iterate that held one is then
  passed over, and an iterate that holds none is recovered by the sweeps, as are
  its trials, until its potentials hold one (OptimaSolver 0.7.8). Pitzer and SIT
  keep the sweeps.
- Measured with the same answers. The 32 cement pastes of a thesis, solved
  with the certified search and every printed value unchanged, take 349 s
  instead of 1238 s with OptimaSolver 0.7.3; the attempt under a limiting law
  past its range, before the fall to an ion size per ion, takes seconds where it
  took minutes. The three-hour hydration of `scripts/ionic_hydration.jl` takes
  5.8 s instead of 7.8 s, on a trajectory identical to the last digit, and its
  28-day certified replay (`speciated_states`) 2.5 s instead of 10.3 s.
- A constraint whose `hq` changes the parameters of the activity model (the
  temperature of an adiabatic solve, the pressure of a fixed-volume one) declares
  them with `pq`, so that the inversion sees what the model sees; one that does
  not declare it keeps the sweeps.

### Changed: derivatives through a certified equilibrium are exact at every level

- **The derivative of a certified equilibrium is that of the problem solved.**
  An equilibrium whose amounts, temperature, pressure or budget carry
  `ForwardDiff` dual numbers is solved on their values and lifted by the
  implicit-function theorem at the certified answer (OptimaSolver's
  `dual_newton_tangent`), on the problem `solve` poses: the constraint's
  unknowns and the surface potentials included. The derivative of the titrant
  a prescribed pH needs, or of the temperature an adiabatic solve reaches, is
  now returned with the answer. Until 0.28.2 the derivative came from the
  optimality conditions of the unconstrained problem, whatever the constraint.
- **Nested differentiations are exact.** Each level of duals is stripped in
  turn, so a derivative taken inside another is solved on the outer one's
  duals. Stripping every level at once, as 0.28.2 did, gave a second derivative
  of exactly zero. This is what an inversion needs whose forward model is an
  equilibrium: a gradient-based fit differentiated again for its Hessian, or a
  sensitivity of a fitted parameter.
- `SciMLBase.solve(::DualEquilibriumSolver, state; b)` and `solve_certified`
  accept a state or a budget carrying dual numbers, by the same route.
- **The implicit kinetic step has no difference quotient left.** Its Newton over
  the reaction extents (`coupling = :species`) took its Jacobian by differences,
  one certified equilibrium per extent. It now takes the residual and its
  Jacobian from one evaluation on dual numbers, through the certified
  equilibrium of the free species. The rate closures of the `:reactions` route
  return the number type of the composition, which OptimaSolver's exact outer
  Jacobian differentiates through.
- **What is differentiated is no longer only the amounts.** The standard
  potentials of a species whose property functions capture dual numbers (a
  `NumericFunc` shifting a `ΔₐG⁰`, a parameter of a thermodynamic model), the
  parameters of an activity model, and the target of a constraint all carry
  their derivatives through the solve, checked against identities exact at
  every equilibrium: the derivative of a mass-action residual is `1/RT` with
  respect to a shift of one potential and zero with respect to anything else.
- **`equilibrate(state, solver)` lifts its answer the same way.** The back end
  chosen solves on the values; the answer is lifted at the point it returned,
  every level of a nested differentiation in turn, with whatever carries the
  duals: the state, the budget, the data or the model. Until 0.28.2 only the
  amounts of the state were seen, one level deep, through a sensitivity solved in
  `Float64`; a dual budget reached Ipopt and raised. Without OptimaSolver, that
  older sensitivity remains, and anything it cannot lift is refused by name.

### Changed: a kinetic run, a recipe and a surface differentiate

- **A kinetic run is differentiated with respect to its parameters.** The state
  and every buffer of the right-hand side take the number type of what the
  problem is given (a rate constant a closure captures, an initial amount, the
  temperature, the calorimeter, an explicit `heat_per_mol`); under partial
  equilibrium the partition of each step is the certified equilibrium lifted to
  those duals. The sensitivity of the partition that the heat balance uses, and
  its shift with temperature, come from the same tangent; until 0.28.2 they were
  a singular value decomposition in `Float64` with absent phases pinned by a
  threshold of its own. A differentiated run takes the plain run's route to the
  values of each partition (the interior point, escalated where its balance is
  poor) and only then lifts them, keeping, as the plain run does, the better
  balanced of the two answers: started cold, the certified solve on the duals
  did not always converge where the plain run did, and an accepted step of a
  differentiated hydration was left 3.4e9 mol out of balance.
- **The observables are exact derivatives of the solution.** `heat_flow` is the
  derivative of the interpolant, and the heat rate of `heat_release` is `−dH/dt`
  at each certified state, from the rates of the kinetic amounts and of the
  temperature with the partition's own sensitivity.
- **A recipe differentiates.** `Recipe`, its materials and constituents, the
  extents (a Parrott–Killoh extent tabulated in the number type of its rate
  constants, its temperature and its ceiling), `budget` and the readers of a
  `RecipeState` take the number type of what they are given: a water/binder
  ratio, a mass fraction, an oxide analysis or a rate constant carries its
  derivative into the equilibrium of the paste. Checked against closed forms: the
  water and the reacted mineral enter the budget by their columns, and scaling
  the rate constants of a Parrott–Killoh law scales its time.
- **A surface differentiates.** A site density reaches the budget of its family,
  and, when the sites follow their host, the conservation matrix and the
  equilibrium through it; the thickness of a Donnan layer reaches its potential
  and its contents, through the fixed point of `equilibrate_donnan`; and a
  `log K` given to `site_family` reaches the energy of its complex.

### Changed: identifiability, without differences

- `log_sensitivity` and `identifiability` differentiate exactly, by forward
  mode: one evaluation of the forward model on dual numbers per chunk of up to
  twelve parameters, where central differences cost two per parameter.
  `relstep` is accepted and ignored, with a deprecation warning. The 5 % step of
  the differences moved an exponent by 0.165 and gave a condition number of
  about eighty to an exact degeneracy, now 1e10; it crossed the kink of the
  Parrott–Killoh minimum and credited `k₂` with an eighth of `k₁`'s influence on
  a heat curve, where it has none; and a 20 K secant across a 15 K peak flattened
  the spectrum of a thermogram.
- **The rank is read off the standard errors when there is a noise level.**
  `identifiable_rank` counts the directions whose standard error in `log θ`, the
  noise level `σ` over the singular value, is below `tol` (default 1, a factor
  e); without a noise level, the gap rule stands. On exact spectra the gap had
  answered five directions for the six rate parameters of a heat curve, of which
  three are determined and the others known to within a factor of 89 at best, and
  three for the six parameters of a thermogram fitted exactly, all determined.
- **The noise level is given or estimated.** `identifiability(...; noise)` takes
  the standard deviation of the instrument; without it, `σ` is the residual
  standard deviation on the `n − p` degrees of freedom the parameters leave,
  never below `√eps` times the root-mean-square of the curve, so that an exact
  synthetic fit does not count a direction at the rounding. The residual RMSE
  over `n` that 0.28.2 scaled the standard errors by understated them by
  `√((n − p)/n)`: by 11 % for six parameters on thirty instants. `Identifiability`
  carries the level as `noise`.
- The covariance is formed from the singular value decomposition rather than by
  inverting `JᵀJ`, which squares the condition number: at an exact degeneracy
  that inverse had no sign left, and gave a trade-off a correlation of +1.
- The analyses stored in `scripts/hydration_calibration.jl` are regenerated with
  exact derivatives, and `main()` recomputes the stored ones on the same six
  parameters and thirty instants. The residual at the published parameters is
  25.94 J/g, where the stored analysis said 26.08: the stored value was older
  than 0.28.2, which gives 25.94 as well.

### Changed: requirements

- ChemistryLab requires **OptimaSolver 0.7.8** (`OptimaSolver = "0.7.8"`). Its
  `SolutionPhase(…; invert)` carries the inversion of the aqueous phase above;
  its `dual_newton_tangent` lifts every answer above; and its exact outer
  Jacobian is what the kinetic steps and the certified search run on. Earlier
  releases are not enough: 0.7.7 passes over every trial step without a
  composition, even from an iterate that has none, and so stops a solve at its
  first iterate where the start holds none, which the chloride binding of
  `test/chloride_binding_reference.jl` meets; 0.7.6 has no `invert`; and on
  cement pastes carrying a trace component (a trace of carbon nine orders of
  magnitude below the major elements) 0.7.5 could lose the trace and not bring
  it back, so that `equilibrate_certified` returned answers the certificate
  refused, the balance of that component wrong by its whole budget, where
  0.28.2 with OptimaSolver 0.7.3 certified.

### Fixed

- The explicit `heat_per_mol` of a kinetic reaction was converted to `Float64`,
  which stopped its derivative.
- Warnings and refusals that rounded a dual number for their message
  (`round(x; sigdigits)` on a `Dual` overflows the stack) print its value, and
  the hint that names the range of an activity model is no longer lost when the
  ionic strength is a dual number.
- The split seed of an automatic instance empties its receiver first: a
  rounding left from an earlier split had made one fail.

### Added

- **A recipe as a kinetic problem:** `KineticsProblem(recipe, system, rates, tspan)`.
  `rates` maps the name of a mineral constituent to its rate law. That
  constituent enters whole and unreacted, and dissolves into the primaries of
  the system by the reaction its column of the conservation matrix gives. Every
  other constituent is taken as `budget` takes it at the start of the run.
  The reacted part of an oxide constituent, such as the alkalis of a clinker,
  enters as the primaries that carry its elements: the protons a basic oxide
  consumes as hydroxide, the water an acidic one takes (SO₃) from the mixing
  water. The same recipe drives `hydrate`, where the extents are
  imposed, and `integrate`, where rate laws decide them. Until now a kinetic run
  was assembled by hand, species by species. A glass, known by its oxides only,
  has no formula to dissolve and is refused a rate.

## v0.28.2 — Kinetic steps seventeen times cheaper, and the volume fractions and windowed bound water of a recipe

### Changed: a step of a coupled kinetic run costs seventeen times less

The partial equilibrium of each step is solved from the previous one, and most
of its cost was spent where the inner fixed point of the dual Newton could not
converge. A solute on its way down moved by exactly 30 in log-amount per sweep,
the stall rule stopped the sweeps before it arrived, and the line search then
refused forty candidates per iteration for want of a converged inner solve that
the current point lacked as well. The solver of the steps now lets a solute fall
to its potential at once and asks a candidate only for what the current point
has. The options are OptimaSolver's (0.7.4), and nothing else uses them: static
equilibria, the certified replay of `speciated_states` and every certificate are
unchanged.

| three hours of a CEM I paste (82 steps) | before | now |
|:--|--:|--:|
| time | 176.1 s | 10.3 s |
| Newton iterations | 17 877 | 1 453 |

Where the rate laws read only degrees of reaction, as here, the trajectory is
the same to the last bit. Where the speciation feeds back, it moves within the
solver's tolerance: the semi-adiabatic run of its example page, whose
temperature follows the enthalpy of the speciation, takes 182 steps instead of
184 and its reference paste peaks at 40.3 °C at 1.02 d, as the page says, where
0.28.1 printed 1.04 d. The pages `coupled_hydration` and
`semiadiabatic_calorimetry` run in 269 s instead of 944 s. `DualEquilibriumSolver` takes the two options as keywords,
`inner_fall_bound` and `lenient_line_search`, with the defaults of before.

### Added

- `volume_fractions(rs::RecipeState)`: the share of the initial volume of a
  paste (reactants, water and residue, the reference of `porosity`) held by each
  species of the equilibrium and each unreacted constituent, closed by the void
  of chemical shrinkage, so the fractions sum to one and `void` is
  `porosity(rs).void`. A residue without a sourced density is refused by name,
  since its volume, and with it every fraction, is unknown. This is what a
  homogenization scheme reads off a recipe; until now it had to be assembled from
  `volume(rs)` by hand, and the residue was easy to forget.
- `bound_water(rs; window = (T₁, T₂), windows)`: the water the solids release
  between two temperatures, from a `DecompositionWindow` per solid holding water,
  unreacted minerals included. A source that reports bound water between 105 °C
  and 550 °C is compared over that range and not with the whole of the ignition
  water, which a thermogravimetric reading over a narrower range cannot weigh. A
  solid holding water without a window is refused by name. Without a window the
  function returns what it returned before.

## v0.28.1 — Equilibria two to four times cheaper, data read off figures that say so, and one notation throughout the documentation

### Added

- `literature_table_info(key, table)`: what a file of `data/literature` says
  about a table besides its values, its provenance kind, its location in the
  source and, for values read off a figure rather than printed as numbers, a
  `digitization` record: the method (a vector drawing or a raster image), the
  tool, and the error the reading adds, relative, absolute with its unit, or
  unquantified. The format accepts the field and checks it; the five tables read
  off figures (Deschner et al. 2012, Hirao et al. 2005 twice, Gruyaert et al.
  2010, Lavergne et al. 2018) carry it, and a test requires it of any table whose
  location says it was read off a figure. A validation that compares with such
  values can now add the reading to the uncertainty of the source.
  `LiteratureRecord` holds it in a new field, `table_info`; a record built with
  the constructor of 0.28.0, without it, is still accepted.
- A warning when a pure solid declared in a system repeats, up to a factor, an
  end-member of a declared non-ideal solid solution with the same Gibbs energy
  within 0.1 RT: ettringite declared pure beside `AFt_SO4_CO3`, whose SO4
  end-member is ettringite divided by three. It is the case already warned for
  two solid solutions, with one of them pure, and the same flat direction that
  can keep the certified search from concluding. Twelve such pairs exist in
  Cemdata18 among the four non-ideal AFm and AFt binaries.

### Changed: an equilibrium costs two to four times less

The certified equilibrium of a cement paste is the step a kinetic run repeats
thousands of times, so its cost was measured, on one machine, before and after,
in one process each, answers compared species by species:

| solve | 0.28.0 | now |
|:--|--:|--:|
| cement of the CEM IV page (107 species), cold | 0.065 s | 0.058 s |
| CEM I with the CNASH gel mixing ideally, cold | 0.179 s | 0.077 s |
| the same with the CNASH gel on its sites, cold | 0.433 s | 0.098 s |
| the same, ten warm restarts on neighboring budgets | 4.53 s | 1.35 s |

Every answer is certified and the same to 6e-11 in relative amount (bit for bit
where the route did not change). Where the time goes to continuations and
restarts rather than to the solves themselves the gain is small: 2 % over 32
blended-cement calculations, every printed result unchanged. Three causes, each
measured with a profiler:

- The coefficients `A` and `B` of the Debye–Hückel term depend on the
  temperature and the pressure only, through the water equation of state, and
  were recomputed at every evaluation of the activities: 55 % of the second
  solve. The last value computed from plain numbers is now kept, in an atomic
  field so that threads can share it; dual numbers, for derivatives in `T` or
  `P`, are computed each time.
- For a phase mixing on sites, the same problem under ideal mixing was solved
  first, as a starting point, even where the start of the linear program
  certifies at once: 47 % of the warm restarts. It is now solved only if the
  search reaches it, after that start and before the state as given.
- `pKw`, needed for the pH of every state a solver returns, was obtained by
  building the reaction H2O@ = H+ + OH- and combining its functions
  symbolically: 1.3 ms a call. It is now summed from the standard Gibbs energies
  of the three species, 10 µs, the same to 1.8e-15.

A certified answer whose route changed reports `route = :lp_start` where it
reported `:ideal_mixing`.

### Removed

- `scripts/blended_cement_kinetics.jl`. It represented a slag by the formula of
  anorthite and a metakaolin by a formula of its own, both with a placeholder
  Gibbs energy and assumed heats, and let the pozzolanic reaction consume
  portlandite without limit. None of it came from a source. The kinetics of
  blended cements return on published data.

### Fixed

- `data/NOTICE.md` did not list `cemdata18-cashplus.json`, and said the derived
  databases copy the Cemdata18 entries unchanged, while that one replaces
  `CaSiO3@`.
- Statements in the documentation that were not true:
  - the route table of the kinetics manual said the ODE route imposes the
    assemblage, which holds only without an equilibrium solver;
  - the docstring of `KineticsSolver` said the partition is re-speciated at each
    evaluation of the right-hand side, where it is once per accepted step;
  - the slag page announced a coupled slag calculation that does not exist;
  - the SIT docstring and the theory page gave an ionic strength of validity, 3
    to 4 mol/kg, that no source of this package states.
- Two comparisons with Reaktoro that nothing here had checked were removed: a
  comment on its kinetics, and a sentence of the kinetics manual.
- The notation of the documentation. The same quantity was written two ways
  (``C_p`` and ``C_P``; the conservation matrix as an italic ``A`` and a bold
  ``\mathbf{A}``), the hint over ``\varepsilon_r`` read it as a SIT coefficient
  of a reaction, and the multipliers ``y`` had opposite signs on two theory
  pages. A vector is now set in bold lower case, a matrix in bold capitals and
  their components in italic, the transpose as ``^\mathsf{T}``, and the
  multipliers keep the solver's sign throughout (``\mathbf{u} =
  -\mathbf{A}^\mathsf{T}\mathbf{y}``). The convention is stated at the top of
  the nomenclature, every symbol has one entry, and `test/docs_nomenclature.jl`
  refuses a duplicated name and the forms of a second convention.

## v0.28.0 — The CASH+ model of C-S-H, three pore-solution validations, and the linear-programming start without its slowdown

The C-S-H of a cement paste can now be described by the CASH+ model of Kulik,
Miron & Lothenbach (2022), with the sodium and potassium of Miron et al. (2022a,
b). The model is written in the compound energy formalism: the site mixing of a
sublattice model, plus the energy of the reciprocal reactions between its
end-members and regular interactions on each site. Until now ChemistryLab could
only mix ideally on the sites. Two checks against the papers:

- In the Ca-Si-H2O system at 25 °C, the model gives back the paper's invariant
  points to the digits it prints. Beside portlandite, the C-S-H has Ca/Si 1.640.
  Its bridging-tetrahedron sites are 77.2 % calcium, 20.3 % vacancy and 2.6 %
  silicate, and its interlayer sites 55.0 % calcium and 45.0 % vacancy. Beside
  amorphous silica, the C-S-H has Ca/Si 0.723.
- The 110 pseudocompounds of the discretized CASH+NK model that Miron et al.
  (2022a) publish are C-S-H compositions whose Gibbs energy the authors computed
  with their own implementation. Our model gives each of them to within
  0.1 kJ/mol, which is the effect of their formulas being printed to four
  decimals.

The model is also extended to Li, Rb, Cs, Mg, Sr, Ba and Ra (CASH+ext). Three
validations set the package against published pore solutions: a limestone
Portland cement from one day to 400 days, with CSHQ and with CASH+NK (Lothenbach
et al. 2008); 48 solutions of the first six hours (Schöler et al. 2017); and 55
solutions of fly-ash blends over 550 days (Deschner et al. 2012). The start the
linear program gives the certified search, which had made some cement
calculations 2 to 3.5 times slower since 0.26.0, now keeps its gain on cold
cements without that cost, and the activity models floor an amount at 1e-30 mol
instead of 1e-16, which a pore solution at pH 14 in a few grams of water needs.

### Breaking changes

- Below 1.0 a minor release is a breaking one for Julia's resolver: a package
  bounding `ChemistryLab = "0.27"` does not accept 0.28.0 and must widen its
  bound.
- A coefficient smaller than 1e-3 in a formula or a reaction is now kept, where
  it was read as zero (see Fixed). A formula or an equation that relied on it
  being dropped now carries that element or species.
- The certified search takes one candidate from the start the linear program
  gives, where it took three, and places each species of that start in its
  phase (see Fixed). A certified answer is the same to the tolerance of its
  certificate, but the start it comes from (`route`), `n_dual_solves` and the
  time can change.
- ChemistryLab requires **OptimaSolver 0.7.3** (`OptimaSolver = "0.7.3"`,
  which 0.7.1 and 0.7.2 do not meet). Its certificate holds a member of a present
  phase below the floor, whose potential the answer determines, to the inequality
  it had skipped, so an answer that certified with a solute left far below its
  equilibrium amount is refused (see Fixed).
- The activity models floor an amount at **1e-30 mol** (`_ACTIVITY_FLOOR`), where
  they floored it at `ϵ = 1e-16`. A species between the two now has the activity
  of its own amount, so trace amounts and a pH read at the old floor change. `ϵ`
  keeps its other roles (the bound of the interior-point back ends, the amount of
  an absent product in a cold state, the regularizations), and the dual solve and
  the implicit kinetic step now start with no amount below it, as the
  interior-point back ends did (see Fixed).

### Added

- `CompoundEnergyModel(lattice; interactions)` and
  `compound_energy_model("<key>:<model>", end_members)`. The first builds a
  model from a sublattice model and its site interactions. The second reads the
  same from `data/literature/<key>.json`.
  - The end-members must be every compound of the sites, each once; the
    constructor refuses any other set.
  - Their amounts are not unique, since the Gibbs energy depends on the site
    fractions alone. The solver is given the split in which the amounts are the
    product of the site fractions. To get it, the energy it minimizes adds `RT D`,
    where `D` is the divergence of the amounts from that product. `D` is never
    negative and vanishes, with its gradient, at the product, so the equilibrium
    and the chemical potentials are those of the model. Measured on twelve pastes
    of CASH+NK (Ca/Si 1 and 1.6; sodium, potassium or both), the equilibria
    certify with or without `D`, in under a second after compilation. What `D`
    adds is an answer whose member amounts do not depend on where the search
    started.
  - The activities depend on the standard Gibbs energies of the members. The
    solver passes those at the temperature of the solve, and the result is
    unchanged when each energy is shifted by that of its elements.
  - The convexity of such a model is decided in its site fractions, from the
    energies of its members at the temperature of the solve (`mixing_convexity(model,
    n; T, g)`, which the certificate calls). With two sites a bound proves it: each
    site's curvature against the coupling the reference surface puts between them.
    The CASH+ core is convex by that bound. CASH+NK is not decided by it, and no
    concave point is found on a lattice of site fractions, so its certificates are
    scoped `:kkt_point`. The verdict is kept, since it depends only on the model,
    the energies and the temperature.
- The derived database `cemdata18-cashplus.json`: Cemdata18 plus the twelve
  end-members of CASH+NK.
  - For the core end-members, G° and H° come from Table 8 of Kulik et al., and
    S°, Cp° and V° from Table 4. The H° of TSvh is the one Miron et al. (2022a)
    reprint.
  - It also carries the CaSiO3@ complex the model was fitted with (Table 9,
    accepted variant): its G° is −1514.14 kJ/mol, 3.42 kJ/mol above the value in
    Cemdata18.
  - It is the only base entry that is replaced, and the database records this.
- The `CASH+` (six end-members) and `CASH+NK` (twelve) phases in
  `data/solid_solutions.toml`, with `model = "compound_energy"`, and in
  `data/gel_models.toml` as one more model of the C-S-H gel. They are therefore
  refused beside CSHQ, CNASH_ss or the ECSH families.
- The `CASH+ext` phase: the model with the interlayer extended to Li, Rb, Cs,
  Mg, Sr, Ba and Ra (Miron et al. 2022a).
  - It has 33 end-members and the 55 interaction parameters of the interlayer
    site (Table A3).
  - Its energies give back the log K of the authors' Table A2 to 0.01, for every
    end-member.
  - The database adds the cations Cemdata18 lacks (Li+, Rb+, Cs+, Ba+2, Ra+2),
    with the properties the paper tabulates.
  - It also adds the Ca(OH)2@ complex the authors derived and kept when they
    fitted the alkali extension. The alkali systems of the CASH+ page and of the
    PC4 paste keep it, as the authors did. The core model leaves it out, as Kulik
    et al. fitted it.
  - The CASH+ page computes the gel with strontium and cesium beside sodium and
    potassium: at Ca/Si 1.2 it holds 78 % of the strontium and 4 % of the cesium.
- `data/literature/Miron2022a.json` and `Miron2022b.json`. They hold the six
  alkali end-members and the interaction parameters of the interlayer site
  (Tables A1 and A3), the discretized model, and the TCNh and TCKh fine-tuned
  for cement pore solutions (Table 5 of the second paper), which the database
  carries.
- A validation of CASH+NK on a hydrating cement: the Portland cement with 4 %
  limestone of Lothenbach, Le Saout, Gallucci & Scrivener (2008), from one day
  to 400 days, computed with CSHQ and with CASH+NK. It is on the CASH+ page, with
  its data in `data/literature/LothenbachLeSaout2008.json` and its recipe in
  `scripts/lothenbach_2008.jl`, and `test/validation_lothenbach2008.jl` checks it.
  - The paste certifies at every age with both models.
  - Miron et al. (2022b) reprint the cement and the kinetic constants. The test
    requires the two transcriptions to be the same numbers.
- A validation of the aqueous model on the 48 early pore solutions of Schöler et
  al. (2017): a CEM I 52.5 R alone and blended with slag, fly ash, limestone or
  quartz, analyzed during the first six hours. Each solution is speciated at its
  measured pH and its saturation indices set against the authors' (Table 7).
  - Portlandite and gypsum agree within 0.08, the Ca-rich C-S-H within 0.21.
  - Ettringite and monosulfate differ by up to 0.9. The differences obey
    `ΔE − ΔMs = 2 ΔGp` to 0.02, so what varies is a factor on aluminum, up to
    0.28 log units at 1.5–4.5 µmol/L of it. The page reports this as found.
  - Data in `data/literature/Scholer2017.json`, the calculation in
    `scripts/scholer_2017.jl`, the page `tutorials/validation_early_pore_solutions.md`,
    and the test `test/validation_scholer2017.jl`.
- A validation of the aqueous model on the pore solutions of Deschner et al.
  (2012), over 550 days: a CEM I 42.5 N alone and blended with 50 % of two
  siliceous fly ashes, of quartz, or of fly ash and limestone. Each of the 55
  solutions is speciated at its measured hydroxide, and its effective saturation
  indices set against the authors' (Table 3).
  - Ettringite and strätlingite agree within 0.03, monosulfate within 0.05,
    gypsum within 0.07 and portlandite within 0.10, the mean differences about
    0.01.
  - The paper prints the analyses only as plots. They are read from the vector
    drawing of its figures, marker by marker, through the major ticks of each
    panel, which adds less than 0.2 % to the values plotted.
  - Data in `data/literature/Deschner2012.json`, the calculation in
    `scripts/deschner_2012.jl`, the page `tutorials/validation_fly_ash_pore_solutions.md`,
    and the test `test/validation_deschner2012.jl`.
- The CASH+ page computes the gel at 50 and 90 °C. The pH falls by 1.8 units from
  25 to 90 °C, somewhat more than Kulik et al. state, and beside portlandite the
  silicon rises, as they state.
- `ParrottKillohExtent` takes `parameters`, constants of the law that replace
  those of Parrott and Killoh (1984), and `H`, the critical degree of hydration
  of its w/c factor, as Lothenbach et al. (2008) fit one per clinker phase. The
  rate constants are per day, as the papers print them, unless given with a
  unit. Left out, both give the law as before, bit for bit.
- `phase_list_system` takes `replace`, which declares another solid solution in
  place of one of the list's (another model of the same gel), and
  `exclude_aqueous`, which leaves out more aqueous species.
- `data/literature/Kulik2022.json`, with its notes. The end-members' G and H are
  those of Table 8, because the G° and H° columns of Table 10 are shifted by
  one row against its names. The printed H° of TSvh is off by 0.28 kJ/mol from
  its G° and S°, most probably through two transposed digits. The printed
  values are kept.

### Fixed

- The hint shown over an equation printed the HTML of a unit as text, as in
  `C m<sup>−2</sup> (kg/mol)<sup>½</sup>` for the Gouy–Chapman prefactor. It
  affected three entries of the nomenclature. The unit is now rendered like the
  name, and a test allows only `<sub>`, `<sup>` and `<b>` in these fields.
- A coefficient smaller than 1e-3 in a formula or a reaction was read as zero,
  so its element or species was dropped. In `Ca2.0993Si2.9298Na0.0004O11.0585H6.1988`
  the sodium disappeared, and in a reaction `0.0004Na+` went with it. Such a
  coefficient now stays, wherever a coefficient of a species is read, written or
  given: in formulas, in equations, in `Reaction` and in their printed forms.
  - Only a value at the level of round-off, 1e-12 and below, is still cleaned to
    zero.
  - `stoich_coef_round` itself is unchanged. It cleans the coefficients a
    computation produces, and a conservation matrix relies on it: a first version
    of this fix changed it, left entries of 1e-17 in the matrices, and made
    `reactions(cs.SM)` overflow the stack of the symbolic simplification.
  - A coefficient within 1e-3 of a simple fraction is still read as that
    fraction, as the database formulas need (`((CaO)1.25(SiO2)1(H2O)2.75)0.6667`
    has 5/6 Ca).
  - Measured on the 3391 species of the eight databases: no composition changes.
- The start the linear program gives the certified search (0.26.0) had made
  some cement calculations 2 to 3.5 times slower. Two causes:
  - Every species outside the vertex of the program started at its activity in
    moles, up to one mole, whatever its phase: pure phases the program found
    undersaturated were in it, and dozens of species near a mole. On four cement
    pastes the start was 11 to 66 mol off a budget of about 4 mol, and on two of
    them nothing certified from it. Each species now starts in its own phase. A
    solute starts at the molality the multipliers give it, per kilogram of the
    water at the vertex, and a member of a solid solution present at the vertex
    at its fraction of that phase, and neither below `ϵ`, the floor of the
    search. A pure phase, or a solid solution absent from the vertex, starts at
    zero. The start is then 1.3 to 1.8 mol off, and all four pastes certify from
    it at the first attempt.
  - Three candidates were drawn from that start: the answer of each back end
    from it, and the start itself. Measured on seven solves of five cements,
    only the default back end's answer ever paid. The interior point from the
    start never certified, and cost 1.5 to 10 s each time. The start itself
    certified only where a start from the state as given had already certified.
    The search now takes that one candidate and goes on to the state as given.
  - Thirty-two calculations of blended-cement pastes that took 3335 s with
    0.27.0 take 1277 s, with the activity floor below (see the next entry), to
    the same answers at the printed digits; one of them now certifies under the
    activity model it asks for first, where it fell back on another. Four of them
    had taken 230 s with 0.25.2 and 470 s with the start of 0.26.0. From the cast
    state of cement107 and of a CEM I with the CNASH gel, the search takes 0.05 s
    and 0.26 s, against 11.3 s and 10.1 s with `lp_start = false`, to the same
    composition (3e-10). A paste on which nothing certifies before the
    continuation pays one failed candidate more than with `lp_start = false`: 3 s
    of 45 s on the one measured.
- A certified answer could hold a solute far below its equilibrium amount, and
  the cause was the floor of the activities. Found on one of 32 cement
  calculations run to check the change of start above: H+ at 3e-100 mol in a
  paste whose potentials give it 1.2e-16, certified, and `pH(eq, model)`, which
  reads the amount, 0.09 high; every other quantity of the answer was right.
  - With a few grams of water per 100 g of binder, H+ at pH 13.5 to 14 is about
    1e-16 mol, the floor `ϵ` the activity models read an amount at. Below it the
    activity no longer follows the amount, and the dual solve, which recovers a
    solute from its own stationarity assuming it does, raised it by a fraction
    of a log unit per sweep and left it there.
  - The floor is now 1e-30 mol (`_ACTIVITY_FLOOR`). Leal, Kulik, Smith and Saar
    (2017) recall the recommendation of Leal, Kulik and Kosakowski (2016) that an
    unstable species hold less than one molecule (1.66e-24 mol) in a system of
    one mole. GEMS3K eliminates a solution species below 1e-30 mol, and PHREEQC a
    molality below 1e-30 mol/kg. Reaktoro bounds every amount by 1e-16 mol, a
    bound of its solver as `ϵ` is here.
  - A solute still found below the floor at a converged answer, where the
    potentials give it more, is given that amount (`floor·exp(r)`, the exact
    solution of the floored model).
  - The old floor had also been standing in for a start. From a state with
    species at exactly zero, the dual solve started them at `exp(−700)`, and only
    the flat activity below 1e-16 had kept an implicit kinetic step on calcite
    converging; at the new floor it did not (its extent came out at 2/3 of the
    rate times the step, uncertified). The dual solve and the kinetic step now
    start every amount at `ϵ` at least, the budget being that of the state, and
    the same step certifies with the element balance at 1e-13 where it was 2e-11.
    A step the suite recorded as a known fragility, C3A and gypsum with S-2 among
    the species, whose extents came out a factor 500 short, now certifies with the
    right ones.
  - The certificate missed it: it left a member of a present phase below its own
    floor (1e-25 mol) out of every test. OptimaSolver 0.7.3 holds such a member,
    where the answer determines its potential, to the one-sided form of the
    stationarity, reported as `stationarity_floored`.
- The small amounts the package works with are defined once, with the reason for
  each value, in `src/utils/numerical_floors.jl` (`_AMOUNT_FLOOR`,
  `_ACTIVITY_FLOOR`, `_CERTIFICATE_FLOOR`, `_LOG_UNDERFLOW`). They were literals
  in some forty signatures and formulas.
- Building species on several threads at once could crash Julia. Every
  thermodynamic factory memoizes the functions it compiles, the factories are
  shared, and the memo was a `Dict` without a lock, so independent calculations
  building their species on threads inserted into it concurrently. On four
  threads, 200 new parameter sets lost an entry in 2 runs out of 20, and a
  documentation build computing six coupled trajectories at once ended in a
  segmentation fault. The memo is now locked; after the first call for a set of
  parameters the lock guards a lookup only.
- Obtaining a database from several threads, or from two processes sharing a
  depot, is serialized. The checksum memo and the set of announced versions were
  unguarded, and two builds of the same derived database wrote the same `.part`
  file, one deleting it under the other. Each download or build now writes a
  temporary file of its own, moved into place when complete.

### Changed

- `_solid_solution_lna!` and `_ss_log_activities!` take the `ΔₐG⁰/RT` of the
  species as an optional last argument. Every existing model ignores it, so their
  results are unchanged bit for bit.
- The databases ChemistryLab derives from Cemdata18 (zeolites, chloride,
  CASH+) are built again on first use after the update, since the builder of
  one of them changed.

## v0.27.0 — Validated against measured pastes: phase lists, processes, a second instance on demand, and a nomenclature

Three published sets of pastes are computed from their papers' data and set
against their measurements and against GEMS3K run on the same budgets: a CEM I
42.5 N through its first year, four ternary cements with fly ash and limestone,
and four mortars of a CEM I 52.5 N with limestone and metakaolin, carbonated. The
two codes agree to 0.001 in pH before carbonation and 0.011 along it. A phase
declared `instances = :auto` is given its second composition only when the
certificate finds it wanting to split. The phases of a paste can be taken from a
phase list written after a paper, and the processes of the recipe layer run on
it. A nomenclature lists the symbols of the formulas, and hovering an equation
shows those it holds.

### Breaking changes

- Below 1.0 a minor release is a breaking one for Julia's resolver: a package
  bounding `ChemistryLab = "0.26"` does not accept 0.27.0 and must widen its
  bound.
- ChemistryLab requires **OptimaSolver 0.7.1** (`OptimaSolver = "0.7.1"`, which
  0.7.0 does not meet). Its test of the components a budget forces to zero is a
  fixed point; without it, a paste declaring phases of an element its budget
  lacks could lose the linear-programming start and take minutes.
- A solid solution with a Redlich-Kister or regular excess term is inverted by
  Newton's method instead of by successive substitution. A certified answer is
  the same to the tolerance of its certificate, but the route to it, the time it
  takes, and whether an answer certifies from a given start can change; the
  miscibility-gap page changed its conclusion with it.

### Added: a second instance when a phase wants it, `instances = :auto`

`SolidSolutionPhase(...; instances = :auto)`, or `instances = "auto"` in a
solid-solution file, declares one composition and allows a second. The phase is
solved with one, and when the certificate of `equilibrate_certified` finds it
wanting to split, on an answer that has otherwise reached the solution (stationary
and balanced to 1e-8, the split its worst violation), the system is rebuilt with a
second instance, the answer carried over, and the passes of `equilibrate_split`
look for the pair. A split read on an answer still far from the solution would
give a phase a second instance it does not need: on a slag cement paste, the AFt
of a first answer with an element balance of 0.17. The answer of the passes is
kept when its KKT error is smaller, and the certificate then names the phase in
`instances_added`. A caller whose composition stays outside the gap never pays for
the second copy, and one who did not know the phase would unmix no longer has to
declare two instances in advance. On a calcite–magnesite binary held inside its
gap by the budget, the pair found is the common tangent of the model to 1e-3, and
the answer of `instances = 2` to 1e-6. A phase whose composition stays outside
its gap keeps its single instance, as the AFm of the miscibility-gap page does.
`with_instances(cs, name => k)` rebuilds a system with `k` instances of a phase,
the primaries unchanged, and `with_instances(state, cs)` carries a state across.

### Added: validations against measured pastes, and against GEMS3K

The CEM I 42.5 N of Lothenbach and Winnefeld (2006) is computed from their
published data (composition, minor elements of the clinker phases, the rate law of
Parrott and Killoh with their water/cement factor) through its first year, and
its pore solution compared with the one they measured (their Table 3, transcribed
into `data/literature/LothenbachWinnefeld2006.json`, the values given only as
detection limits marked as such). The same budgets were run through GEMS3K, on
the Cemdata18 cement export of xGEMS with the same phases: the two codes agree to
0.001 in pH and 1.6 % on every element, so the differences with the paste belong
to the model. The largest is the alkalis: the `CSHQ` model of Cemdata18 holds
96 % of the sodium and 68 % of the potassium at 317 days, and the model's sodium
is a tenth of the measured one. The page is *Validation against a measured
paste*; `test/validation_lw2006.jl` holds the recipe to the budgets of the replay
and the certified equilibria to GEMS3K, with the fixture written by
`test/reference/xgems_replay.py`.

The ternary cements of De Weerdt et al. (2011) follow on *Validation against
measured blended pastes*: a CEM I, the same with 5 % limestone, a CEM II/B-V with
35 % siliceous fly ash and a CEM II/B-M (V-LL) with 30 % fly ash and 5 %
limestone, their clinker dissolving as their XRD measured it and their fly ash at
the rate their Fig. 7 prints, transcribed in `data/literature/DeWeerdt2011.json`.
GEMS3K on the twenty budgets agrees with ChemistryLab to 0.001 in pH and 1.6 % on
every element, the AFm sulfate and hydroxide declared in both as the Guggenheim
binary of Cemdata18. Against the pastes, the portlandite without fly ash is
within 1.5 wt.% from the seventh day; with fly ash the model consumes it far
faster (4.5 against 12.5 wt.% at 90 days in the CEM II/B-V), its C-S-H staying at
the Ca/Si of 1.58 that portlandite imposes, where the paper measures 1.4 and an
Al/Si of 0.13 that `CSHQ` cannot take. `test/validation_deweerdt2011.jl` checks
the transcription, the budgets and eight equilibria against GEMS3K.

The carbonation page is rewritten on the four mortars of Shi et al. (2016): a
white CEM I 52.5 N alone, with limestone, with metakaolin and with both, hydrated
91 days, then carbonated in steps up to 50 g of CO2 per 100 g of binder, their
materials and degrees of hydration transcribed in `data/literature/Shi2016.json`.
It replaces a page that carbonated an assumed composition with an uncertified
solver and static figures. The portlandite of the cement alone is within 1 % of
the thermogravimetric measurement, and the CO2 binding capacity within 2 % of the
authors' calculation in three pastes and 7 % below it in the metakaolin paste,
whose gel is richer in calcium. The pH holds above 13 while portlandite
carbonates, then falls to the plateau near 9.7 the authors computed, and the CO2
taken up when it passes 9.7 orders the pastes as theirs does. GEMS3K on the 44
carbonated budgets agrees to 0.011 in pH and 7 % on every dissolved element above
0.01 mmol/kg. `test/validation_shi2016.jl` checks the transcription of Table 5,
the budgets, and twelve steps of the four walks against GEMS3K.

### Added: phase lists

`data/phase_lists.toml` names the phases a paste may form, after the paper it is
written for, and `phase_list_system(name, substances; add, remove)` builds the
chemical system from it: its pure phases, its solid solutions as
`data/solid_solutions.toml` declares them, and the aqueous species their elements
allow. A calculation that departs from the list says so with `add` and `remove`.
The first list is the Portland paste of Lothenbach and Winnefeld (2006), which
the two validation pages now build their systems from instead of each writing
the same list; `phase_list(name)` and `phase_lists()` read the file.
`build_solid_solutions` takes `instances`, a `Dict` from a phase name to the
`instances` to declare it with in place of the file's, which is how the list
declares the AFm binary with `:auto`.

### Added

- A `ProcessResult` is indexed with `begin` and `end`, as `result[end]`, and has
  `keys`, so that `findfirst(predicate, result)` gives the step.
- `titrate` and `add_salt` take a formula that is not a species of the system,
  `add_salt(rs, "NaCl", amounts)`, entered through its primaries. A salt the
  database has no solid for could before be added only one ion at a time, which
  leaves the budget charged.
- `heat_release(rs1, rs2)` gives the heat a paste releases between two states of
  one recipe. The residue counts by what changed between the two, so an oxide the
  system cannot hold, set aside with the same mass in both, no longer makes the
  heat `NaN` because its enthalpy is unknown.
- Three material templates, the clinker, the siliceous fly ash and the limestone
  of De Weerdt et al. (2011), transcribed in `data/literature/DeWeerdt2011.json`.
  A template given by its phases can take `remainder = true`: what its analysis
  holds beyond the phases (the free lime, alkalis and magnesia of a clinker)
  becomes one oxide constituent, `"minor oxides"`. Hematite is among the Rietveld
  phases a template knows.
- Three material templates, the white Portland cement (CEM I 52.5 N), the
  limestone and the metakaolin of Shi et al. (2016), transcribed in
  `data/literature/Shi2016.json`. `remainder = "analysis"` keeps the oxides no
  phase holds at the amounts of the analysis, and scales the phases down when
  they leave less room than those oxides take. The phases of that cement, counted
  with their pure formulas, leave 1.3 % of it for the 3.9 % of magnesia, alkalis
  and sulfate its analysis gives, and `remainder = true` would have cut its
  potassium to a third. Free lime is among the phases a template knows.

### Changed: a phase with an excess term is inverted by Newton's method

The composition of a solid solution with a Redlich-Kister or regular excess term
is now recovered by Newton's method, as that of a sublattice phase already was,
instead of by successive substitution. On a CEM I paste declaring the published
AFm SO4/OH binary of Cemdata18 as one composition, the substitution held the
search for minutes and never certified; Newton certifies it in one route. With
the binary declared as the GEMS3K export of Cemdata18 declares it, the twenty
budgets of De Weerdt et al. (2011) give the same pore solution in both codes, to
0.001 in pH and 1.6 % on every element. Ideal mixing keeps the substitution, and
its results bit for bit.

The miscibility-gap page changes with it. Its paste's AFm lies outside the gap,
at a C4AH13 fraction of 0.28, and the published parameters with one composition
now certify; the page said that the certificate refused them because the phase
wanted to split, a verdict read off a search that had not converged.

### Fixed

- A paste declaring phases of an element its budget lacks could fall from the
  linear-programming start of the certified search to a search of minutes, and
  come back with a phase split it did not need. On a carbonated CEM I declaring
  Friedel's and Kuzel's salts without chlorine, one step took 121 s and gave the
  AFm a second instance, where the same paste without the two salts certified
  in 1.1 s with one. The cause was in the test that finds the components a
  budget forces to zero: the electron row, whose entries of one sign sit on
  perchlorate, was left free when chlorine was absent, and the dual solve
  stagnated on its multiplier. OptimaSolver 0.7.1 reads each row over the species
  still free, and ChemistryLab now requires it (`OptimaSolver = "0.7.1"`); the
  step takes 1.1 s, with one instance, as without the salts.
- A solve that builds a starting point (the ideal pre-solve, and the pre-solve
  with ideal mixing) could give a phase declared `instances = :auto` its second
  instance, and hand a state of the enlarged system to a search in the system it
  was meant for; the dual solve then indexed past the end of its conservation
  matrix. Met on a carbonated metakaolin blend. Those solves no longer split a
  phase: their answer stays in the system of the search.
- `leach` read every step in the system of the paste it started from, so once a
  step had given a phase declared `instances = :auto` its second instance, the
  next step failed on a size mismatch. Each step now reads the system of the
  answer before it.
- A glass or a remainder found by difference could not be built when two
  analyses of one material disagree: the limestone of De Weerdt et al. is 81 %
  CaCO3 by TGA, which holds more lime than its XRF analysis gives, and the oxides
  left then summed to more than their share, which `OxideConstituent` refuses.
  They are now taken as the whole of that share; a material whose analyses agree
  is built as before.

### Documentation

- `cemdata18_activity_model` keeps `b_γ` at its 25 °C value, and its docstring
  now says why beyond the paper giving no other: Helgeson et al. (1981) tabulate
  `b_γ` against temperature for six chlorides (their Table 26), not for KOH or
  NaOH.
- `manual/databases.md` has a table of the five models of the C-S-H gel the
  package ships, read from the files: members, mixing, source, and the elements
  the members hold besides calcium and silicon.
- The comparison of activity models has a fourth column, the B-dot model with
  the parameters Cemdata18 prescribes.
- A new page, *Mixing on sites: the CSH3T and CNASH gels*, solves CSH3T in
  water mixed by end-members and by sites (the silicon in solution differs by up
  to half), and reads the minimum chain length of Eq. (11) of Myers et al.
  (2014) from a C-(N-)A-S-H gel, the bridging-site vacancies of their Table 1
  transcribed into `data/literature/Myers2014.json`.
- The theory chapter and the miscibility-gap page no longer say that no code
  splits a phase by itself.
- A new page, *A CEM I paste replaced by fly ash, carbonated, salted and
  leached*, runs the four processes of the recipe layer on the CEM I of De Weerdt
  et al. (2011) at 90 days, its phases from a phase list: the pozzolanic
  consumption of portlandite as fly ash replaces the cement, the pH held by the
  alkalis while portlandite carbonates and then falling with the Ca/Si of the
  gel, chloride bound as Kuzel's salt and then Friedel's salt, the ionic
  strength leaving the range of the activity model, and a paste leached by
  renewals of its pore water.
- `manual/recipes.md` has a section on the bound water and the heat of a paste,
  and points to the phase lists.
- A **Nomenclature** page lists the symbols of the formulas, grouped by subject,
  with their units, the pages where a meaning holds when a symbol means
  different things in different chapters (A the conservation matrix or the
  Debye–Hückel parameter, φ the osmotic coefficient, a heat loss or the
  equilibrium map), and the value the package computes with for a physical
  constant, read from the library. **Hovering an equation** now shows the symbols
  it holds, with their meaning on that page: the symbols are found in the TeX
  source when the formula is typeset (`docs/src/.vitepress/nomenclature-match.mjs`),
  from `docs/nomenclature.toml`, the one file the page and the hints are built
  from. The home page and the theory chapter point to it.
- Section 8 of the surface-complexation theory defines each symbol where it is
  used (Ψ, σ, C, 𝒜, F, R, T). Formulas written as code on four pages are
  typeset, one of them a raw `\sqrt` shown as text; the heat-capacity polynomial
  of `:cp_ft_equation` has eleven terms, not ten, and its last is ln T.

## v0.26.0 — Cement modeling: databases from their publishers, the activity model of Cemdata18, sublattice mixing, a linear-programming start and recipes

The thermodynamic databases are no longer shipped: they are obtained from their
publishers on first use, checked against a checksum. The activity model that
Cemdata18 prescribes gets a name, and every certified answer states the ionic
strength it was reached at against the range of its model. The C-(N-)A-S-H gel
is mixed on its sites, as Myers et al. define it, and convexity is decided for
any number of end-members. The certified search starts from the linear program
over the pure phases, which refuses a budget no amounts can meet and makes cold
cement solves an order of magnitude faster. A layer of materials, extents,
recipes and processes poses a cement calculation in the terms it is described
in.

### Breaking changes

Below 1.0 the registry treats a minor bump as breaking whatever the API did, so
`[compat] = "0.25"` will not accept `0.26`, and a dependent must widen its bound.
The documentation environment of MeanFieldHomogenization.jl lists ChemistryLab
at `0.24, 0.25` and OptimaSolver at `0.5, 0.6`, and needs `0.26` and `0.7`
added; PoroMechanics.jl lists `0.15.2, 0.18, 0.22` and needs `0.23` to `0.26`.

Several behaviors change deliberately:

- **OptimaSolver 0.7 is required**, for its linear-programming start and its
  Newton inversion of a sublattice phase.
- **The databases are downloaded, not shipped.** `datapath` keeps working for
  every database name, but the first call needs the network, a directory named
  by `CHEMISTRYLAB_DATABASE_DIR`, or a copy installed with `install_database`.
  Code that opened the files under `data/` directly no longer finds them, and
  `cemdata18-merged.json` is gone (its substances are those of
  `cemdata18-thermofun.json`). The release obtained differs from the files
  shipped until now in ten Cemdata18 species and ten PSI/Nagra reactions, which
  change a result only where they enter it.
- **`CNASH_ss` is the published sublattice model**, so every CNASH answer
  changes.
- **A concave phase with more than two end-members is refused** at
  construction, as a concave binary was, and a certificate is scoped
  `:kkt_point` unless every mixing phase is proved convex.
- **A budget no amounts of the declared species can meet is refused at once**,
  with `budget_feasible = false` and `route = :infeasible`, and
  `STRICT_CONVERGENCE` raises; until now it went through the whole search.
- **The search starts from the linear program** (`lp_start = true`): the
  answers are the same, the route to them is not, and `lp_start = false`
  restores the former search.
- `parrot_killoh`, `parrot_killoh_avrami` and the literature key
  `ParrotKilloh1984` are deprecated in favor of the author's spelling, Parrott;
  the old names still work, with a warning.


### Thermodynamic databases are obtained from their publishers

`datapath("cemdata18-thermofun.json")`, and the same call for the PSI/Nagra,
aq17 and slop98 databases, resolves to the file its publisher distributes:
release v1.1.1 of ThermoHub, obtained on first use from GitHub (or the jsDelivr
mirror of the same commit), checked against the SHA-256 of the version
ChemistryLab is validated with, and kept in a cache of the Julia depot. Nothing
changes for code written with `datapath`. Four functions come with it:

- `database_info()` says where each database currently resolves;
- `fetch_databases()` obtains them all in advance, for a machine that will work
  offline, a documentation build or continuous integration;
- `install_database(path)` installs a copy downloaded by hand, after checking
  its checksum;
- `database_path(name)` is what `datapath` calls for a database name.

`ENV["CHEMISTRYLAB_DATABASE_DIR"]` names a directory of local copies, searched
first. When a file cannot be obtained, `DatabaseUnavailable` says which, what was
tried and the ways to provide it; the manual page *Database Interoperability*
explains every message.

The PHREEQC export of Cemdata18 (`CEMDATA18-31-03-2022-phaseVol.dat`) is
distributed by Empa through a page no program can use, so it is downloaded by
hand once and installed with `install_database`. Nothing ChemistryLab does by
default needs it: it is the input of `merge_json` and of the PHREEQC readers.

`cemdata18-zeolites.json` and `cemdata18-chloride.json` are built on first use
from the Cemdata18 file and ChemistryLab's own data (the published zeolite
tables, the fitted chloride end member in `data/chloride/cshq_cl.json`), and
rebuilt whenever either changes. Their added entries are identical to those of
0.25.2. `cemdata18-merged.json` is gone: its substances were those of
`cemdata18-thermofun.json`, byte for byte, and the pages that loaded it load the
latter, with the same results.

### The release of the databases

ThermoHub's release v1.1.1 carries values identical to the files ChemistryLab
read until now, except for three differences, which reach a calculation only
where the species concerned enter it:

- Cemdata18 has ten more species, with their reactions: `Fe(OH)3(am)`,
  `Fe(OH)3(mic)`, `FeCO3(pr)`, `S-2`, `CN-`, `MgSiO3@`, `AlHSiO3+2`,
  `FeHSiO3+2`, `Fe2(OH)2+4` and `Fe3(OH)4+5`;
- PSI/Nagra 12/07 writes ten reactions over different reactants;
- in slop98-organic, the enthalpy of `CH4@` differs by 1 J/mol.

Two consequences are measured in the test suite:

- The two amorphous iron hydroxides are now checked against Table 2 of
  Cemdata18 like the other rows, and close to the printed digit.
- `S-2`, in a sulfate paste without reductant, stays at the floor and leaves the
  equilibrium unchanged. Its presence nevertheless makes the implicit kinetic
  step fail to converge on the two C3A pathways of the test suite, with extents
  a factor 500 short. The outcome depends on the Gibbs energy of `S-2` in a way
  no chemistry explains, which marks it as numerical. The two tests leave `S-2`
  out, for that written reason, and a third records the failure as broken.

### Fixed

- **`merge_json` failed on every ThermoHub file.** It wrote its output by
  splicing text between the lines `"reactions": [` and `"elements": [`,
  assuming that order; ThermoHub lists `elements` first, and the splice ended
  in a `BoundsError`. The output is now written by the JSON writer, with the
  fields in the order the input gives them.
- **A ThermoFun record without its reference state made the whole database
  unreadable.** One substance of slop98-organic, `Eth@`, carries no `Tst`; it
  is now read at ThermoFun's reference state, 298.15 K and 1 bar.
- **The PHREEQC reader missed three options.** `parse_phases` looked for `-V⁰`,
  which no PHREEQC file writes, so the molar volumes (`-Vm`) were never read;
  it recognized `-log_K` but not `-log_k`, so a phase written in lowercase
  (`zeoliteP_Ca`) had no log K and was dropped from a merge; and it recognized
  `-analytical_expression` but not `-analytic`. Options are now recognized in
  every spelling the PHREEQC manual documents, whatever their case, and `-Vm`
  is kept as the molar volume of the phase, never as a volume of reaction.
- The chloride end member of CSHQ inherited, from the NaSiOH record it is built
  on, a `mass_per_mole` field that is NaSiOH's; ChemistryLab never read it, and
  it is no longer written.
- The documentation said the PHREEQC file brought the phase volumes to the
  merged database. The ThermoFun file already carries them for all 140
  crystalline phases; the note is removed.

### Corrected: the activity model of Cemdata18

The documentation, a docstring and two test comments said that CEMDATA18
carries no ion-size parameter, so that a GEM-Selektor run starts from `å = 0`.
The paper says otherwise: Appendix C, Eq. (C.1), prescribes the extended
Debye–Hückel equation with a common ion size of 3.67 Å and a B-dot of 0.123
for a KOH electrolyte (3.31 Å and 0.098 for NaOH), and the same B-dot on the
neutral species. The GEM-Selektor run reconstructed in the documentation was
configured with `å = 0`; that is a property of that run, now described as such.
The CEMDATA18 cement published with xGEMS runs Eq. (C.1), and GEMS3K on it is a
test oracle (`test/reference/xgems_cement.py`).

### Added: the activity model Cemdata18 prescribes, and its range

- `cemdata18_activity_model(:KOH)` or `(:NaOH)` returns that model, with its
  parameters read from the transcription of the paper.
- The certificate of `equilibrate_certified` reports `ionic_strength`,
  `activity_range` (the ionic strength up to which the model is stated valid,
  `nothing` when no range is stated) and `within_activity_range`.
- `equilibrate_certified(...; fallback_model, fallback_on)` solves again with a
  second activity model when the first does not certify, or, with
  `fallback_on = :out_of_range`, when its answer lies past its range.

### Added: the activity model of a PHREEQC database

`TruesdellJonesActivityModel` applies the WATEQ equation to every species a
PHREEQC database gives `-gamma` parameters, and PHREEQC's defaults to the others
(Davies for an ion, `0.1 I` for a neutral species). `phreeqc_gamma_parameters`
reads those parameters from a database file, master and secondary species alike.
A second `-gamma` for a species replaces the first, as PHREEQC reads it.

### Added: a budget no amounts can meet is refused, with its reason

`equilibrate_certified` first solves the linear program over the pure phases,
the equilibrium with every activity at one (OptimaSolver's `lp_start`). When no
non-negative amounts of the declared species meet the element budget, the
program proves it with a Farkas vector, a combination of the balances that every
species raises and the budget lowers, and the budget is refused at once: the
certificate has `budget_feasible = false`, `route = :infeasible` and
`unplaceable`, the reason in words ("the budget asks for −0.001 mol of Ca+2, and
every declared species holds it with a non-negative coefficient"), and
`STRICT_CONVERGENCE` raises. Until now such a budget went through the whole
search (continuation, restarts, repairs) and came back uncertified, with the
solver's non-convergence counter raised on the way.

### Changed: the search starts from the linear program

The vertex of that program, with every species it leaves out raised to the
amount its multipliers give it, is now the first start of the certified search.
Measured on two cement pastes against the same calls with `lp_start = false`:

| | cement107 | CEM I with CNASH_ss |
|:--|--:|--:|
| cold | 4.58 s → 0.21 s | 4.19 s → 0.09 s |
| warm, from the answer | 0.207 s → 0.208 s | 1.49 s → 0.09 s |
| a neighbor (1 % more water), from the answer | 0.156 s → 0.164 s | 2.40 s → 0.09 s |

every composition the same to 2e-8. `lp_start = false` restores the former
search. The certificate also says which start the answer came from (`route`:
`:lp_start`, `:state`, `:ideal_mixing`, `:ideal`, `:continuation`, `:restart`,
`:repair`) and how many dual solves the search ran (`n_dual_solves`).

### Added: ideal mixing on sublattices, and the CNASH gel as its authors define it

`SublatticeModel` is ideal mixing on several sites of a formula unit,
`ln aₖ = Σₛ mₛ ln y_{s,σₛ(k)}`, the model of Kulik (2011) for CSH3T and of Myers
et al. (2014) for the C-(N-)A-S-H gel. `sublattice_model("Myers2014:cnash",
members)` reads a published model from `data/literature`, and checks each
member's formula against the one the paper prints, rescaling the site
coefficients when the database stores another formula unit (the Cemdata18 records
of CSH3T are half Kulik's). `site_fractions` gives the site fractions. The
occupancy transcribed from Myers' Table 1 reproduces the activities of their
Appendix B, transcribed apart, and Kulik's Eq. (19) likewise; the tests hold the
two transcriptions against each other.

`CNASH_ss` in `data/solid_solutions.toml` is now that model. It was declared as
ideal mixing of its eight end-members, which is not the model of the paper that
defines it, so every CNASH answer of the package so far differs from the
published model's; results computed with it change. The solver inverts such a
phase by Newton's method (OptimaSolver 0.7, which this release requires): the
substitution it uses for other phases diverges on it. A system holding a
sublattice phase is first solved with that phase mixing ideally, and that answer
starts the sublattice solve (on a CEM I paste, measured before the start from
the linear program: 0.86 s, where the search alone took 20.5 s). Known limitation, measured: under the limiting law with `å = 0`, the
CEM I paste with the published gel does not certify, where it does under
`cemdata18_activity_model`.

### Added: the other C-S-H models of Cemdata18, and two Al/Fe binaries

- `CSH3T`, `ECSH1` and `ECSH2` are declared as ideal solid solutions, as
  Cemdata18 ships them (Lothenbach et al. 2019, Table 4), and listed in
  `data/gel_models.toml`, so that two models of the one gel stay refused
  together.
- `AFt_AlFe` and `AFm_AlFe`, the Guggenheim binaries of Cemdata18's Table 1
  (notes b and h) between ettringite and its iron analog and between
  monosulfate and its iron analog, with their published parameters and their
  miscibility gaps.

### Changed: convexity is decided for any number of end-members

`mixing_convexity(model, n)` says whether a mixing energy is convex over the
whole simplex (`:convex`, `:nonconvex` with a witness composition, or
`:undecided`), and why. A regular model is proved convex when the smallest
eigenvalue of `W/RT` on the tangent space of the simplex is at least −2, which
for two end-members is the familiar `W ≤ 2RT`, and shown not to be when a pair
exceeds `2RT`. Until now only binaries were checked: a concave ternary was
accepted without `instances = 2`, and its certified answers were scoped
`:global_minimum`, which is what they do not prove. Such a phase is now refused
at construction (or admitted with two instances), and a certificate is scoped
`:kkt_point` unless every mixing phase is proved convex.

### Added: a member that is a mixture of two others is reported

`SolidSolutionPhase` warns when a member is, in composition and to `0.1 RT` in
Gibbs energy, the mixture of two others: the phase then holds one substance
twice, and ideal mixing counts its configurations twice. The MgAl-OH-LDH ternary
of Cemdata18 is such a case (M6A is the average of M4A and M8A); it ships as
`MgAl_OH_LDH`, as published, with `acknowledge_degenerate = true`. Ordered
members whose Gibbs energy departs from the mixture, as T5C of CSH3T or the
siliceous hydrogarnets, are not reported.

### Added: recipes, extents and processes

A layer for what every cement calculation needs and each page wrote by hand:

- **materials** of mineral constituents (database phases) and oxide
  constituents (glasses, minor oxides), with `bogue` from the formulas of the
  library, `decompose` (non-negative least squares on an oxide analysis) and
  `reactive_part` (the glass by difference from a Rietveld analysis);
- **extents**: constant, tabulated, logistic in log time, the rate law of
  Parrott and Killoh integrated (optionally with their water/cement factor, as
  Lothenbach and Winnefeld (2006) apply it), and any of them capped by Powers'
  limit;
- **recipes** and their `budget`: the reacted parts enter the equilibrium, the
  unreacted ones are kept aside with their mass, volume and enthalpy, and an
  oxide whose element the system cannot hold is kept aside too, said so;
  a volume or a heat that needs a density or an enthalpy no source gives is
  reported missing, never estimated;
- `equilibrate_certified(recipe, system)` and a `RecipeState` read by
  `phase_masses`, `porosity`, `bound_water`, `pore_solution`, `volume` and
  `enthalpy`, the residue counted where it belongs;
- **processes** as sequences of certified equilibria: `hydrate`, `blend`,
  `titrate`, `carbonate`, `add_salt`, `leach`, tabulated by `process_table`;
- **templates** of published materials (`data/recipe_templates.toml`, references
  only; `material_template`), starting with the round robin of Durdziński et al.
  (2017): the Portland cement by its phases and by Bogue, the two slags, the
  siliceous fly ash with its glass by difference.

The layer poses exactly the problem the recipe of `scripts/gruyaert2010.jl`
built by hand (the element budget agrees to 1e-12), which the tests assert.

### Corrected: Parrott, not Parrot

The author of the 1984 hydration model is L. J. Parrott. The functions are now
`parrott_killoh` and `parrott_killoh_avrami`, the literature key
`ParrottKilloh1984`, and the table of Lavergne et al. (2018) `parrott_killoh_1984`.
The old names keep working, with a deprecation warning. The evidence is
indirect, since the 1984 paper has no DOI: the same author, with the same
initials and in the same years, signs as Parrott on Crossref, including two
papers on alite hydration (1981, 1983).

### Documentation: the application pages name their cement

The titles of the application pages give the EN 197-1 designation of the cement
they compute (CEM II/A-LL and CEM II/B-S, CEM III/A, CEM IV/A (V) and CEM IV/B
(V), CEM V/A (S-V), and the CEM I of each Portland page), where several said
only "a blastfurnace cement" or "the full Portland cement".

### Documentation: the cement pages run the activity model Cemdata18 prescribes

The pages on CEM I with its solid solutions, CEM II, CEM III/A, CEM IV, CEM V,
the miscibility gap and chloride binding now run `cemdata18_activity_model(:KOH)`
in place of the limiting law with `å = 0`, and each prints the molar K/Na ratio
of its alkalis (2.3 to 2.6) that makes the KOH parameters the right set. The CEM
IV and CEM V pages declare `CNASH_ss` as the published sublattice model. Every
number the prose quotes was taken again from the executed output:

- On the Portland and slag pastes the pH rises by 0.01 to 0.04 (CEM I 13.0994
  to 13.1156, CEM III/A 13.041 to 13.054), and the alkali members of the C-S-H
  hold a few percent more potassium and sodium.
- The CEM I cross-check against Reaktoro was run again with both codes on the
  Cemdata18 model: pH 13.1156 against 13.1426, total volume within 0.06 %. The
  two halves of `scripts/crosscheck` now read the database file and the activity
  parameters from the one file the Julia half writes.
- In the chloride page the paste binds more of the added chloride: 74 to 96 %
  instead of 65 to 94 % with the surface model, 85 to 97 % instead of 83 to 96 %
  with the chloride end member.
- On the CEM IV page the gel changes most. The published `CNASH_ss` takes a Ca/Si
  of 1.17 where `CSHQ` takes 1.60, so its paste keeps more portlandite. Every
  point of the portlandite sweep certifies, and so do both full-reaction limits
  of the CEM IV/B, which did not under `å = 0` (element balances of 0.26 and
  0.35 mol left). For `CSHQ` nothing else changed, so what failed there was the
  limiting law. The `CNASH_ss` limit lies at an ionic strength of 1.09 mol/kg,
  past the range Cemdata18 states, and its certificate says so. A cold solve now
  certifies the CEM IV/B at 28 days, and the continuation the page walked is
  gone.
- The zeolite block of the CEM IV page had left the alkalis of the clinker out
  of its budget. With the budget of the rest of the page no zeolite is stable in
  the full-reaction limit, where the page reported FAU-Y-K.
- On the CEM V page the fully reacted paste now certifies, at an ionic strength
  of 1.10 mol/kg, past the stated range.
- The miscibility-gap page stated three things its own output contradicted: the
  AFm composition of the single-composition case lies just outside the
  common-tangent pair, not inside; that case does not close its element balance
  to 2e-14; and the two instances of the last case do not share the amount
  lopsidedly. The text now says what the output shows, and the pH axis of the
  summary figure no longer hides two of its bars.

### License notices

Five source files carried no license header: `oxide_budget.jl`, `retention.jl`,
`implicit_step.jl`, `certified.jl` and `constraints.jl`. They now open with the
same `SPDX-License-Identifier: LGPL-2.1-or-later` line and copyright notice as
every other file of the package. `NOTICE` misspelled the name of the hydration
model's author (Parrott) and used a UK spelling; both are corrected.

## v0.25.2 — The options of the dual solver reach it, and the range of an activity model

Found while correcting a user's cement scripts, with OptimaSolver 0.6.2, which
fixes three defects of the dual Newton's search for the phases present.

### `equilibrate_certified` gave its dual solver no option

**`equilibrate_certified` built its `DualEquilibriumSolver` with the defaults,
whatever it was called with**, so a caller could not give the dual Newton more
iterations, a larger active-set budget or a tighter inner tolerance on the route
that certifies. It now takes `dual = (; maxit, max_active_updates, inner_tol,
inner_maxit, tol, si_tol)`, forwarded to that solver, to the ideal pre-solve and
to the route of a dual-number budget. `DualEquilibriumSolver` takes `inner_tol`
and `inner_maxit`, the tolerance and the sweep budget of the inner fixed point,
and forwards them to OptimaSolver's `DualNewtonOptions`, whose defaults it keeps.
The docstring says what the other keywords reach: `variable_space` sets the
formulation of the interior-point starts, Ipopt takes the common arguments of
Optimization.jl (`maxiters`, `reltol`, `maxtime`, `verbose`), and
`OptimaOptimizer` ignores the rest.

### Added — the range an activity model is stated for

`activity_model_range(model)` returns the ionic strength up to which the solving
manual states a model valid: 1 mol/kg for `HKFActivityModel`, 0.5 for
`DaviesActivityModel`, `nothing` where the manual gives no number. Neither formula
announces that it has left its range, and the Debye–Hückel limiting law of a
GEM-Selektor comparison (`å = 0`) is worse than inaccurate there: its `log γ`
keeps falling as the ionic strength rises, so a pore solution can run away and
the certified search not conclude. Comparing `ionic_strength(eq)` with this value
is how a caller knows which case it is in.

### Two messages that say what the composition or the species are doing

**A refusal of `equilibrate_certified` names an activity model used past its
range.** When the answer it gives up on has an ionic strength above
`activity_model_range(model)`, the message says so, with the value, and says
what to try: an ion size per ion if the model is the limiting law of a
GEM-Selektor comparison, a model stated for higher ionic strengths otherwise.
Measured on blended cement pastes taken to full reaction, that was the cause of
the refusal, and nothing in the message pointed at it.

**`ChemicalSystem` warns when two declared solid solutions, one of them
non-ideal, hold one substance twice**: an end-member whose composition is a
multiple of another's, with Gibbs energies that agree to that factor within
0.1 RT per formula unit. CEMDATA18's `ettringite03_ss`, the SO4 end of the AFt
binary, is ettringite divided by three to 5 J/mol, and GEM-Selektor's list
declares the binary beside the `ettringite` solid solution. With ideal mixing
that is harmless, and nothing is said. With the published Redlich–Kister model
on the binary, measured on four cement pastes, the certified search stopped
short of the solution on the flat direction between the two phases, and
certified once the substance was declared once. The warning names both phases
and both species, and says which to drop. `data/solid_solutions.toml` said it
in a comment; the code now says it to the caller.

### Data — the materials of the RILEM TC 238-SCM round robin

`data/literature/Durdzinski2017.json` carries Table 1 of the accepted
manuscript: the oxide analyses of the slags S1 and S2, the fly ashes SFA and
CFA, the Portland cement and the quartz filler, and the phase compositions of
the two fly ashes, with their amorphous fractions, and of the cement. The degrees of reaction the file
already held are those of the same materials, so the composition that reacts and
the fraction that reacts can now be taken from one source. The notes say where
the table was read and what it prints as it is: the S1 analysis sums to 100.4 and
gives no K2O.

### Documentation

- CEM IV example, section 7; CEM V example, section 6; database manual on the
  zeolite extension. The full-reaction limits of these pages do not certify with
  the limiting law the pages run, and the three pages gave three causes for it:
  a hard point for the solver, an unphysical question that no assemblage can
  answer, a gap in the phase list. Measured, the same limits certify with an ion
  size per ion, `HKFActivityModel()`, with nothing else changed. The pages now say
  that it is the limiting law that does not close there, and each runs the
  per-ion solve, which on the CEM IV page lands just past the range the manual
  states for that model too, and says so.
- CEM IV example, section 6: `CSHQ` also binds alkalis in the gel, through `KSiOH`
  and `NaSiOH`; the page said only `CNASH_ss` could.
- `equilibrate_certified` docstring: a section on the options of the solvers.
- Solving manual: `activity_model_range`, and what a refusal past the range
  means. Database manual: when declaring `Ettringite_ss` with `AFt_SO4_CO3` is
  harmless, and when it is not.
- `scripts/cem4_pozzolanic.jl` and `scripts/cem5_composite.jl` regenerated from
  their pages, which `test/scripts.jl` requires.

### Dependencies

`OptimaSolver = "0.6.2"`, in the root and the documentation projects: the
results the pages now show were measured with it.

### Continuous integration

The coverage job, the Julia 1.12 one, gets two hours. It runs without the
precompilation cache and took 49 and 55 minutes on the two runs of 0.25.0; on
0.25.1 it passed the one-hour limit and was canceled while the three other jobs
passed.

### Tests

- `test/test_dual_solver.jl`: the defaults of `DualEquilibriumSolver` are
  OptimaSolver's, the inner options are stored and reach OptimaSolver's
  `DualNewtonOptions`, an unknown option in `dual` is refused, and a
  certificate threshold no answer can meet refuses the answer the default
  certifies, whichever route produced the answer. On 0.25.1 the keyword was
  accepted and ignored: an unknown option raised nothing, and the unreachable
  threshold came back certified.
- `test/aqueous_properties.jl`: `activity_model_range` for each model, and the
  sentence the refusal adds past the range.
- `test/solid_solutions.jl`: one substance in two phases is silent with ideal
  mixing and said with the published model on the AFt binary.
- `test/literature.jl`: the two new tables close to their rounding, and the
  amorphous fractions and slag analyses they give.

## v0.25.1 — Two models of one C-S-H gel refused, and `equilibrate_split` under the strict flag

A system that declared `CSHQ` and `CNASH_ss` together was accepted without a
word, and counted its calcium silicate hydrate twice; `equilibrate_split` could
not be used with `STRICT_CONVERGENCE` set. Both surfaced while correcting a
user's cement scripts, where the first also kept the certified search from
concluding. The documentation said the pair was refused; it was not.

### `CSHQ` and `CNASH_ss` declared together were accepted

**A system declaring two models of the C-S-H that share no composition was
built, and the gel counted twice.** `ChemicalSystem` refuses two declared solid
solutions that share a composition, which catches `CSHQ` with an `ECSH` family,
their alkali end-members being one substance under two names. `CSHQ` and
`CNASH_ss` share none, so the pair passed, although `data/solid_solutions.toml`
warned against it and the CEM IV example stated that it was refused. On cement
pastes declared that way the certified search could fail altogether, and where
it concluded, the calcium, silicon and alkalis of the gel were shared between two
descriptions of one hydrate.

`data/gel_models.toml` now lists the end-member symbols of each model of one gel
— `CSHQ` (with `CSHQ-Cl`), `CNASH_ss`, `ECSH1`, `ECSH2` — and `ChemicalSystem`
refuses two declared solid solutions whose end-members belong to two models of
the same gel, naming both. Matching is by end-member symbol, whatever name a
phase is declared under. The two instances of one declaration stay exempt, one
model split over two declared phases is still accepted, and the composition test
is unchanged and runs first.

### `equilibrate_split` raised under `STRICT_CONVERGENCE` before seeding anything

**With the strict flag set, `equilibrate_split` raised at its first pass**, which
is by construction the pass that does not certify: the function exists to improve
on it. The passes now run with strictness suspended, as the starting routes of
`equilibrate_certified` already did, and the answer they end on is judged
strictly — an error if it does not certify, and the solvent check applied to it.
Without the flag nothing changes.

### Documentation

- CEM IV example: the "Never both at once" note said the pair was refused by
  name; it now says how each pair is refused, and since when.
- `data/solid_solutions.toml` and the manual page on choosing species: the notes
  on declaring two models of one gel name the new test.
- CEM I example with eight solid solutions: its last section said that the
  non-ideal AFm/AFt parameters could not be sourced and that a miscibility gap
  could not be represented. The parameters are in
  `data/literature/Lothenbach2019.json`, and two instances hold a gap; the page
  now points to the miscibility gap example.
- Solving manual: five aqueous activity models, not three (Pitzer and SIT), and
  `RegularSolutionModel` for a non-ideal ternary solid solution.
- `equilibrate` and `equilibrate_certified` docstrings: what `optimal = true`
  proves depends on `cert.scope` (0.25.0); a concave phase given two instances
  is, besides a waived convexity check, the case where the problem is not convex.
- README: `equilibrate` does not use Ipopt "under the hood"; with OptimaSolver
  loaded and an aqueous phase it takes the certified route, the default since
  0.14.
- README, database manual and the `build_solid_solutions` docstring: their
  examples passed the whole of `data/solid_solutions.toml` to `ChemicalSystem`,
  which declares `CSHQ` with `CNASH_ss` and is now refused; they take the phases
  by name, and say why the file is a catalog. The README's list of the shipped
  phases named five of the twelve.

### Tests

- `test/solid_solutions.jl`: `CSHQ` with `CNASH_ss` refused, whatever the phase
  names and the subset of end-members declared; `CNASH_ss` with `ECSH1` refused;
  one model over two phases accepted; the registry covers every C-S-H end-member
  of the database. Each refusal assertion fails on 0.25.0.
- `test/certified_equilibrium.jl`: under the strict flag, `equilibrate_split`
  returns the certified pair of the metastable binary instead of raising, and
  raises when no pass is allowed and the first answer does not certify.

## v0.25.0 — Sites that leave with their host, and certificates that say what they prove

A site family whose budget follows its host, the mechanism that lets a sorbent
dissolve with its sites, did not conserve charge; this release corrects it and
constrains every solver of the package with the corrected matrix. The
certificate of an equilibrium now states whether it proves a global minimum, a
point satisfying the optimality conditions, or a speciation consistent with its
own activities, since the extended activity models in general use are not the
gradient of one Gibbs energy. The derivatives of an equilibrium, which
PoroMechanics.jl takes per cell through dual element amounts, receive the
corrections that package carried as a local patch.

### Breaking changes

Below 1.0 the registry treats a minor bump as breaking whatever the API did, so
`[compat] = "0.24"` will not accept `0.25`, and a dependent must widen its bound.
The documentation environment of MeanFieldHomogenization.jl lists ChemistryLab
at `0.24` and needs `0.25` added; PoroMechanics.jl lists `0.15.2, 0.18, 0.22` and
needs `0.23`, `0.24` and `0.25`.

Several behaviors change deliberately:

- **A coupled site family gives other numbers.** The amounts of an equilibrium
  with a family under `SITES_FOLLOW_HOST` change, the charge now being
  conserved (see the erratum below). The free site of such a family sits at
  zero energy; a family whose free site is given the energy of the matter it
  carries is refused when the bias exceeds 0.05 log units, and so are a charged
  free site and a host whose formula lacks the atoms of its sites.
- **A site family that names no host is refused when it binds an element held
  by a solid solution of the system**, unless its support is declared
  `external = true`. Its sites could otherwise sit on the solid solution and
  count the same element twice, which the refusal of 0.24.0 caught only when a
  host was named.
- **A coupled family whose host is a kinetic species is refused**: a site budget
  that follows an amount moved by a rate law between two re-speciations is not
  accounted for.
- **A kinetic run re-speciates with the activity model of its equilibrium
  solver.** The certifying solver of the equilibrium partition was built with
  the model of the `KineticsProblem`, which defaults to the dilute one, while
  the interior point used the solver's; a run that gave the two different
  models, and was warned about it, re-speciated with the dilute model. An error
  in building that solver is now raised instead of being replaced by the
  interior point.
- **The derivatives of an equilibrium change where they were wrong**, and a
  sensitivity that fails its own stationarity or conservation check raises an
  error instead of being returned.
- **An entry with an entropy but no heat capacity is extrapolated in
  temperature.** It kept its tabulated Gibbs energy at every temperature, which
  contradicts the entropy it carries; it now takes a zero heat capacity. No
  entry of the databases shipped with the package is affected.

### Erratum: charge under host coupling, 0.22.0 to 0.24.0

The coupling subtracted `ν` from the site component in the host's column. With
a bare site component, which is charged, the total charge of the system then
drifted by `ν` times the host dissolved: on hydrous ferric oxide at Dzombak and
Morel's weak-site density titrated by hydrochloric acid, 1 × 10⁻⁴ mol for
3.2 × 10⁻³ mol of chloride. The host's formula in the database already contains
the surface groups its sites are made of, so the coupling now subtracts `ν`
times the column of the free site from the host's (Kulik 2002): every atom is
counted once, the neutral free site conserves the charge, to 10⁻¹³ in the
free-site and in the bare-component bases, and the two bases give one answer.

### Every solver on the same constraints

The interior-point back ends, the implicit-function sensitivities, the repaired
start of the certified route, the homotopy and the conservation matrix of a
kinetic partition now use one `_constraint_matrix`, which is `SM.A` itself when
nothing is coupled. Before, only the dual solver saw the coupling, and the other
routes held the site budget fixed; a test now forbids `SM.A` in the code of these
routes. A kinetic trajectory now keeps the sites at
`ν` times the dissolving host to 10⁻⁶ and the charge to 10⁻¹⁰.

`SurfaceSupport(...; external = true)` declares that a family's sites belong to
a solid outside the system, such as a gel frozen by `freeze_solid_solution`,
which the two-stage chloride route of 0.24.0 now does. An external support that
names a host is refused.

### What a certificate proves

`optimality_certificate` returns `scope` and `scope_reasons`. The optimality
conditions it audits prove a global minimum when the chemical potentials are the
gradient of one convex Gibbs energy. The first property is measured at the
audited composition from the symmetry of the Jacobian of the log activities,
which is 2 × 10⁻¹⁶ away from symmetric for the Debye–Hückel form with a common
ion size and no linear term, but 1.2 × 10⁻² for the ideal dilute model, 0.30 for
the B-dot setting of the GEMS runs of CEMDATA18 and 1.0 for the default B-dot and
Davies forms. The certificate reports `:global_minimum`; `:kkt_point` when the
conditions hold but their sufficiency is not established, for a solid solution
declared inside a miscibility gap or a constraint that shifts the activity of
water or makes the temperature an unknown; or `:self_consistent` when the
activities are not the gradient of one energy. `optimal` keeps its meaning.

### Derivatives of an equilibrium

The sensitivity system is equilibrated by rows and columns before a truncated
singular value decomposition: the curvature of a trace species reaches 10³⁰⁰
beside a conservation block of order one, and a pseudo-inverse of the unscaled
matrix discarded the conservation equations. Pure phases at their bound are
pinned, whereas a member of a mixing phase, an aqueous trace included, stays
free, since pinning it erased the perturbation of the component it carries. The
log activities are differentiated at a floor below the smallest amount present,
and a composition carrying dual numbers is certified at its values. The
regression cases of the PoroMechanics.jl patch are part of the tests, with the
call that package makes, checked against a centered difference to 10⁻⁵. The heat
of a partial-equilibrium run falls back on no sensitivity when one fails its
checks, and then carries the whole change of the partition in its jumps; on the
semi-adiabatic page, the peak of the run without feedback moves from 1.03 to 1.02
day, and the coupled run does not move.

### Thermodynamic data

A heat capacity given on several temperature intervals, as for quartz and
hematite in CEMDATA18, is taken on the interval that contains the reference
temperature by an explicit choice; it was the first interval listed, by the
order in which duplicate parameters overwrote each other, which gives the same
values for the databases shipped.

### Documentation

The first two pages of the Theory chapter defined the apparent and formation
Gibbs energies twice, one in temperature only and the other in temperature and
pressure. The argument is now made once: balanced reactions over the elements
and the charge, formation from the elements with the convention on the hydrogen
ion written as a charge term, the proof that the standard Gibbs energy of
reaction reduces to energies of formation, the same proof over primary species,
the standard potential integrated in temperature and pressure, and the relation
between the chemical potential and the apparent Gibbs energy the package stores.
The page on thermochemistry keeps what belongs to the code, including the models
by which that energy is evaluated, and no longer states that every certificate
proves a global minimum. The manual on stoichiometric matrices had their rows and
columns exchanged. The manual recommended the logarithmic variable space for
systems spanning many orders of magnitude; started from a recipe, whose absent
species sit at the regularization floor, an explicit solver in that space returns
the start (Ipopt) or stops without converging (OptimaSolver), and the manual now
has it refine a solved state. `equilibrate` without a solver, which certifies,
was never affected.

### Tests

The implementation of Pitzer's equations is checked to give activities whose
Jacobian is symmetric, the property on which the certificate rests. The Ipopt
extension, never loaded by the tests so far, is tested: its solve in
both variable spaces against the certified answer, a coupled family, the route
of dual numbers and the registration of its back end. `OptimizationIpopt` joins
the test dependencies, and the file runs last, since loading it adds a starting
point to the certified search.

### A correction to the notes of 0.21.0 and 0.22.0

Their breaking-change sections state that no package of the organization depends
on ChemistryLab. PoroMechanics.jl does, and so does the documentation of
MeanFieldHomogenization.jl.

## v0.24.0 — Chloride in the C-S-H of blended cements

The C-S-H of the CEMDATA18 pages is the CSHQ solid solution, which holds
calcium and alkalis but no chloride, so a salted paste put all of its bound
chloride in the AFm phases. The surface model of the chloride literature could
not be added to it: its sites bind calcium that CSHQ already counts. This
release gives the C-S-H its share in two ways, and says where each applies. The
first freezes the gel after a first equilibrium and puts the published surface
on it; it holds while portlandite buffers the gel. The second adds a chloride
end member to CSHQ, fitted on published sorption tests; it holds at any Ca/Si.
A new documentation page salts a CEM III/A along both routes and a CEM III/B,
without portlandite, along the second.

### Breaking changes

Below 1.0 the registry treats a minor bump as breaking whatever the API did, so
`[compat] = "0.23"` will not accept `0.24`, and a dependent must widen its bound.
The documentation environment of MeanFieldHomogenization.jl lists ChemistryLab
up to `0.22`, and needs `0.23` and `0.24` added.

Two behaviors change deliberately:

- **A site family hosted on an end member of a solid solution is refused when it
  binds an element that solid solution holds.** Silanol sites that bind calcium,
  hosted on an end member of CSHQ, counted the same calcium in the solid and on
  its surface, and nothing said so; `ChemicalSystem` now raises an error that
  names the family, the phase and the elements, and points to
  `freeze_solid_solution`. A family that binds nothing the phase holds, or that
  names no host, is unaffected.
- **A calorimeter under partial equilibrium takes its heat from the enthalpy of
  the states, and refuses a species without one.** The heat used to be that of
  the kinetic reactions alone, which under partial equilibrium only dissolve the
  anhydrous phases: the precipitation of the hydrates by the minimization was
  left out, `integrate` warned about it, and a semi-adiabatic cell overestimated
  its temperature. The heat of such a run changes accordingly, and a system in
  which a species carries no enthalpy of formation now raises an error naming
  it, since its term would drop out of the heat.

### Calorimetry under partial equilibrium

Enthalpy is a state function, so the heat an isothermal calorimeter records is
what the enthalpy of the whole composition loses, the hydrates precipitated by
the minimization included. Under partial equilibrium the calorimeters now
integrate `q̇ = −dH/dt` at the current temperature: the kinetic amounts as the
integrator moves them, and the equilibrium partition through its sensitivity to
the element amounts, from the optimality conditions of the last certified
partition. In a semi-adiabatic cell the partition is followed in temperature as
well, from the same conditions with the Gibbs–Helmholtz relation for the
temperature derivative of the potentials, and the heat it takes up as it shifts
joins the heat capacity of the cell. What a re-speciation changes beyond these
linear predictions, a phase appearing within a step, is added to the
calorimeter's state by the step callback, so the heat is conserved exactly: the
tests close the first law of a semi-adiabatic cell to 4 × 10⁻⁶ J out of 1.7 kJ,
and the isothermal heat agrees with `heat_release` over the certified states to
10⁻⁵. The partition is
re-speciated at the temperature of the cell, and `speciated_states` replays each
instant of a semi-adiabatic run at its own temperature rather than at the
initial one. The two warnings of `integrate` are removed.

### The temperature of a state kept through two fallbacks

`equilibrate_certified` falls back on a homotopy in the element amounts when a
cold start fails, and `equilibrate_split` seeds a second solve from the best of
several. Both built their first state at the default 298.15 K, whatever the
temperature of the state given, so the certificate certified the problem at the
wrong temperature. A paste at 20 °C that took either route came out at 25 °C:
the chloride test of this release failed on Windows only, where the 1 mol/L point
of Hirao's sorption test took the homotopy route and its pH came out lower by the
change of pK_w between the two temperatures. Both routes now keep the temperature
and pressure of the state.

### A solid solution, read and then frozen

`solid_solution_totals(state, name)` returns what a solid solution holds: the
amount of each end member, the moles of each element, and its mass, from the
molar masses the package computes. Ratios such as Ca/Si are quotients of the
element totals; the element type is the state's, so dual numbers give their
derivatives.

`freeze_solid_solution(state, name, target; release = (:Na, :K), buffer)` builds
the first state of a second system in which the solid solution no longer
reacts. Its elements are set aside, except the released ones, which return to
the solution as their cation with as much hydroxide: NaSiOH gives back 0.5 NaOH
and keeps its silica and water. Every element is conserved exactly across the two
stages. It refuses, by name, a species that holds matter but is missing from the
second system, a symbol whose composition differs between the two, an end member
that could form again, and an absent `buffer`: for a C-S-H that is portlandite,
without which a frozen composition models nothing.

The PHREEQC oracle of the C-S-H surface gains a two-stage case. PHREEQC
equilibrates Guo's inventory with CSHQ declared as an ideal `SOLID_SOLUTIONS` of
its own, then, without it, the surface on a C-S-H of the amount and composition
the first stage gave, swept in NaCl. ChemistryLab agrees at both stages, to
5.5 × 10⁻⁵ mol on an end member and 2.8 × 10⁻⁴ mol on a salt, so the first
stage is checked independently as well as the second.

### A chloride end member for CSHQ

`data/cemdata18-chloride.json` is CEMDATA18 unchanged, with `CSHQ-Cl` =
(CaCl₂)₀.₅ appended, and `CSHQ_Cl` in `data/solid_solutions.toml` is CSHQ with
it. Its Gibbs energy is the one fitted number, −4.9 ± 0.3 kJ/mol for its
formation from ½ Ca²⁺ + Cl⁻ at 20 °C, on the chloride bound by C-S-H in the
sorption tests of Hirao et al. (2005). Their Fig. 5 is a vector drawing; its
points were read from the coordinates into `data/literature/Hirao2005.json`, and
only the three up to 1 mol/L enter the fit, the limit of the B-dot activity
model. The model of the test reproduces the depletion measurement, including the
water the dried gel takes up as it rehydrates, which at 1 mol/L hides
0.12 mmol/g of the 0.41 the end member holds. `data/chloride/regenerate.jl`
builds the file and records the fit, its residuals and its uncertainty on the end
member.

The end member is effective. Plusquellec and Nonat (2016) found that chloride
does not adsorb specifically on C-S-H, and one parameter cannot follow the
measured points: the fit is within 0.05 mmol/g at 0.5 and 1 mol/L and four times
too high at 0.1 mol/L. A first candidate, written after NaSiOH with silica in its
formula, fitted as well but bound less chloride at a higher Ca/Si, against the
trend Zibara et al. (2008) measured; the generator fits both and records why the
first was rejected. The shipped one also makes a CaCl₂ solution bind twice as
much as a NaCl one at equal chloride, the ordering Tran et al. (2018) report,
which nothing was fitted or chosen on.

`build_solid_solutions` accepts a `database` key: an entry whose end members
exist in one database only is skipped without a warning when another is loaded.
With CEMDATA18 alone the shipped file therefore builds as before.

### Documentation

*Chloride binding in blended cements* salts the CEM III/A paste of the
blastfurnace cement page with up to 0.4 % chloride. The two routes agree on
Kuzel's salt, which holds most of the chloride, and differ on the C-S-H, which
holds half of the bound chloride at the lowest dose in the second route and a
fifth in the first. The page states two limits the first route inherits: the
surface takes 81 % of the portlandite's calcium before any chloride is added,
calcium the gel's Ca/Si already counted; and at the specific area of the model
the pore water is a film 0.57 nm thick, too thin for the ions of the diffuse
layer to be counted. The manual describes the new database.

A new section of the applications, *Outputs of a calculation*, holds what a
laboratory measures on a paste and can be read off computed states. *An
isothermal calorimeter, read off the states* computes the heat of the pastes of
Gruyaert et al. (2010) at the degrees of hydration their image analysis gives,
against their isothermal calorimetry, and the enthalpy their measurements imply
for the slag glass, which no database holds. *Bound water, and the thermogram it
integrates to*, moved out of the surface pages where it did not belong,
computes the bound water of the same pastes against their thermogravimetry; its
section on the decomposition windows is now labeled as the self-test it is.
*A semi-adiabatic calorimeter, inside the kinetics* puts the Portland cement of
the ionic hydration page in the calorimeter of Lavergne et al. (2018), with the
temperature among the unknowns of the kinetics, against the temperature they
measured: 56.4 °C at 0.82 day where they measured 52.1 °C at 0.75 day, with
nothing adjusted, and 40.3 °C when the same heat is integrated without its
feedback on the rates. It replaces the cell the ionic hydration page integrated
after the calculation, from a heat flow at 20 °C, and the precomputed heat table
no longer carries that temperature. The transcriptions are in
`data/literature/Gruyaert2010.json` and `data/literature/Lavergne2018.json`.

A new first page of the theory, *Energies, enthalpies and the chemical
potential*, rebuilds from the two laws the enthalpy a calorimeter measures, the
Gibbs energy the solver minimizes and the chemical potential between them; it
defines the absolute entropy and the entropy of formation, the reference of the
elements and that of the oxides, and the apparent quantities of the databases,
and shows on the slaking of lime why the Gibbs-Helmholtz relation holds for a
reaction but not for the apparent energy of one species.

Two snippets of the kinetics tutorial and manual counted the heat capacity of the
paste twice, once in the calorimeter's `Cp` and once in the sum the integrator
adds; `Cp` is now the vessel's alone. The license of `data/experimental/` names
all seven files it covers, where it named two.

## v0.23.0 — Published values out of the code, and a C-S-H surface checked against PHREEQC

0.22.2 gave published values a home in `data/literature/`. This release moves
the rest of them there: every value taken from an article, in the source, the
tests, the executed documentation and the scripts, is now read from a file that
names its source and where in that source it was found. Moving them meant
reading each one against its source again, which found wrong attributions, one
model error and several quiet defects. The release also adds a way to build a
surface from published reactions and to count the ions of its diffuse layer,
and uses both to check the C-S-H surface model of the chloride literature
against PHREEQC, alone and in a hydrated paste.

### Breaking changes

Below 1.0 the registry treats a minor bump as breaking whatever the API did, so
`[compat] = "0.22"` will not accept `0.23`, and a dependent must widen its bound.
The documentation environment of MeanFieldHomogenization.jl lists ChemistryLab
up to `0.22` and needs `0.23` added.

Three behaviors change deliberately:

- **`HKFActivityModel` computes different results.** It used the effective
  electrostatic radius of each ion (Table 3 of Helgeson et al., 1981) as the
  ion size of the Debye-Hückel denominator, 1.91 Å for Na⁺. The ion size is a
  property of an electrolyte, defined by their Eq. (125), which gives 3.72 Å for
  NaCl (their Table 2). The model now forms each ion's size from that equation
  for the electrolyte the ion makes with the NaCl background, as TOUGHREACT
  does. It reproduces Table 2 for the twenty-one salts of ions the package
  carries, and it follows Hamer and Wu's NaCl to 0.9 % at 0.1 mol/kg and 1.5 % at
  1 mol/kg, where it was 5 % and 19 % out. An explicit `sp[:a]`, the model's
  common `a` and `a_default` are ion sizes already and are unchanged.
- **A temperature or pressure given to a solve is refused.** `equilibrate`,
  `equilibrate_certified` and `EquilibriumSolver` passed the keywords they did
  not know to the optimizer, so `T = …` or `P = …` was dropped and the solve ran
  at the state's temperature without a word. At a fixed hydroxide amount that
  costs a pore solution the change in pKw, 0.17 between 20 and 25 °C. The four
  keywords `T`, `P`, `temperature` and `pressure` now raise an error that points
  to `set_temperature!`.
- **The shipped `data/solid_solutions.toml` follows Cemdata18 on the AFm and AFt
  binaries.** Its `AFm` entry mixed monosulfate and monocarbonate with
  Redlich-Kister parameters that have no published source, while Cemdata18
  treats the two as pure phases; it is removed, not made ideal. `AFm_SO4_OH`
  and `AFt_SO4_CO3`, which Cemdata18 publishes as non-ideal, were declared
  ideal because a phase that unmixes could not be expressed; they now carry its
  Guggenheim parameters with two instances each, and open the gaps the article
  prints. Code that built the whole file gets a different phase set.

### A surface from its published reactions

`site_family(name, reactions, aqueous; master, site, capacity, support, model)`
builds a `SiteFamily` from surface reactions written as PHREEQC writes them,
with their log K, or from the `SorptionReaction`s `read_sorption_model` returns.
Each complex receives the energies of the aqueous species its reaction consumes.
That one line of arithmetic was repeated by every hand-built surface, and when
it was left out, a cation's complex certified with nothing formed. Built from
`phreeqc.dat`, it reproduces the hydrous-ferric-oxide families the tests used to
build by hand.

### The ions of the diffuse layer

`DonnanLayer`, `diffuse_layer_contents` and `equilibrate_donnan` count the ions
that screen a charged surface, which a `DiffuseLayer` leaves in the solution. As
in PHREEQC's `SURFACE -Donnan`, read in its source, the surface keeps its
Gouy-Chapman potential and a layer of water of fixed thickness holds each
solute at the average Boltzmann enrichment whose charge balances the surface's.
`equilibrate_donnan` withdraws the layer's contents from the solution and
solves again until the two agree; the layer's water is added to the solution's,
as PHREEQC adds it, or taken from it, as a closed pore solution requires.
Against PHREEQC on the C-S-H surface in 18 solutions, the layer's chloride
agrees to 3 × 10⁻⁴ relative and its water exactly.

### The C-S-H surface against PHREEQC

The silanol surface of Elakneswaran et al. (2010), as Guo et al. (2018) use it,
is compared with PHREEQC on the same model: the same reactions, read by both
codes from one table of `Guo2018.json`, the Davies equation on both sides, and
Dzombak and Morel's diffuse layer. The fixture generator takes every energy,
molar mass and atomic mass it needs from the package itself.

- On eighteen closed NaOH–CaCl₂–NaCl systems every point certifies. Site
  fractions agree to 4.4 × 10⁻⁴, log a(H⁺) to 2.2 × 10⁻⁴ and the potential to
  0.04 mV.
- In Guo's paste (portlandite, monosulfate, ettringite, Friedel's and Kuzel's
  salts), swept in NaCl, the phases agree to 2.9 × 10⁻⁴ mol, the surface species
  to 4.5 × 10⁻⁴ mol and the bound chloride to 5 × 10⁻⁵ mol.
- A new page, *Chloride binding by C-S-H and Friedel's salt*, splits the bound
  chloride between the salts, the surface and the diffuse layer. The surface
  holds all of it before any salt forms and a third of it at 0.4 mol/kg; the
  layer, one Debye length thick, adds 2 to 8 %. The page states its scope first:
  a C-S-H of fixed composition, not the CSHQ solid solution, and a layer
  thickness that nothing published fixes.

Guo's deprotonation row is taken in its proton form. Guo, and Elakneswaran et
al. (2010), print it against OH⁻ with a constant that Elakneswaran et al. (2009)
compare with the proton-form values of Viallis-Terrisse (−12.3) and Pointeau et
al. (2006, −12.0).

### A prescribed pH on a diffuse layer certifies

`equilibrate_certified` under `FixedpH`, on a system with a `DiffuseLayer`,
threw a `DimensionMismatch` before measuring anything. The certificate received
the potential of each diffuse layer after the constraint's own unknown, and
`optimality_certificate` built the constraint's block alone. It now composes the
surface block as `solve` does, and audits the augmented problem. A `q` of any
other length is refused by name.

### Published values read from `data/literature/`

New records, each checked against its source (on the page image wherever the
text extraction loses signs or charges):
Atkins1992, BaroghelBouny1999, Blanc2012 (Thermoddem), Duan2016, Durdzinski2017,
DzombakMorel1990, Elakneswaran2009, Guo2018, HamerWu1972, Helgeson1981,
HongGlasser1999, IAPWS2014, Kettler1992, Kulik2002, Lothenbach2008 and
Lothenbach2010 (Cemdata07),
Lothenbach2019 (Cemdata18, Tables 2, 3, D.1 and D.2), MaLothenbach2020 and 2021
(the zeolites), PalandriKharaka2004,
ParkhurstAppelo2013, PlummerBusenberg1982, Pointeau2006 and Xu2012; Lavergne2018
and Powers1948 are extended. Every exported constant keeps its name and its
value, and a test compares it row by row with its file. The vendored
`cemdata18-zeolites.json` is rebuilt from the new records byte for byte.
Results of other codes (Reaktoro, GEM-Selektor, PHREEQC) are fixtures in
`test/reference/`, with versions and checksums.

- `literature_row(key, table, label)` returns a row by its label, and
  `literature_table` filters rows by column values, for tables in long format.
  The unit reader accepts an SI prefix the registry lacks, such as MPa, without
  modifying that registry.
- `water_surface_tension(T)` evaluates the IAPWS R1-76(2014) equation, tested
  against its Table 1.
- `build_solid_solutions` reads `guggenheim = "<key>:<pair>"`, published
  dimensionless parameters taken from `data/literature` rather than copied, and
  `instances`, for a solid solution that unmixes.

### Attributions corrected

- The Debye-Hückel A = 0.5114 and B = 0.3288 were credited to Table 1 of
  Helgeson et al. (1981), which gives 0.5091 and 0.3283. They are the 25 °C
  row of the LLNL table in the PHREEQC manual, which is now cited, as are
  PHREEQC's rule for uncharged species (the default `Kₙ` and `bₙ`) and
  Helgeson's 3.72 Å (the default `a_default`).
- The radii by charge are those of the TOUGHREACT V2 user's guide (Xu et al.,
  2012). The Computers & Geosciences article cited before has no such table.
- The clinker compositions 61.9/16.5/8.0/8.7 and 67.8/16.6/4.0/7.2/2.8 have no
  published source; they were credited to Lavergne et al. (2018) and are now
  labeled as assumed.
- The Redlich-Kister parameters of the AFm solid solution were credited to
  Cemdata18, which has no monosulfate–monocarbonate solid solution. They are
  kept, unchanged, as labeled placeholders.
- The degrees of reaction of Durdzinski et al. (2017) are in their Table 5,
  not Table 4.
- The quadratic heat-loss law with coefficients 0.3 and 0.003 was credited to
  Lavergne et al. (2018). The form is theirs, and the coefficients are
  illustrative; their calibrated device has 75 J/(h·K) and 0.26 J/(h·K²).
- Guo et al. (2018) used Cemdata07, not Cemdata18, and their chloride row is
  charge balanced. The validation chapter said otherwise on both counts.

### Fixed

- The self-desiccation tutorial took 0.0728 N/m for the surface tension of
  water at 25 °C, which is the 20 °C value. It now computes 0.0720 N/m, and the
  Kelvin radius at 80 % relative humidity becomes 4.7 nm.
- The warning of SciMLBase about parameters of mixed types survived 0.22.2, whose
  entry said it was gone: the vector of reactions had no concrete element type
  either. Both vectors are now held behind one wrapper, and a test checks the
  whole parameter tuple of a problem with two rate laws of different types.
- `ignition_loss`, `bound_water_per_phase` and the thermogram used typed molar
  masses for water and carbon dioxide when the system declares neither. They
  now weigh them from the formula, as every species is weighed. A script and a
  test gave a species written as CaAl₂Si₂O₈ a molar mass of 95 g/mol against the
  278.2 g/mol of its atoms, and the override is gone.
- The slag and metakaolin reactions of the blended-cement script created and
  destroyed elements; they are balanced from the formulas. Their heats per
  gram were credited to Gruyaert et al. (2010) and Lothenbach et al. (2011),
  neither of which gives them; they are stated as assumptions.

### Validation

New reference tests, each reading its published values from the files above:
every solubility product and HKF coefficient of Cemdata18 against its tables;
the HKF model away from the reference point against Duan et al. (2016), whose
pressure column is in bar although it is headed in pascals; Atkins et al.
(1992); a limestone blend; Guo's chloride binding; the three PHREEQC
comparisons above; alkali uptake by C-S-H against the 48 solutions of Hong and
Glasser (1999), where the pH agrees to 0.072 from 15 to 100 mM and the alkali
is over-bound, as expected of end members fitted to those data; and the
Cemdata07 generation, from Lothenbach (2010) and Lothenbach et al. (2008): which
of its solubility products and formation energies Cemdata18 kept, the molar
volumes of 28 solids (24 to half the printed digit), and the confirmation that
Guo et al. used it. The validation chapter pins every number it prints at the
precision it prints it, and opens with what has been checked and the traps
worth knowing.

## v0.22.2 — Published values in data files, and output that shows only the result

A value taken from an article is data: it has a source, a location in that
source, a unit and a status (measured, fitted, unverified), none of which
survives being typed into a source file, where a second copy then drifts from
the first unnoticed. This release gives such values a home, `data/literature/`,
one file per source, and moves there the first of them: the constants of the
rate laws. It also removes from the rendered documentation, and from any output
not written to a terminal, the noise that had come to outweigh the results.

### Published values are read from `data/literature/`

Each file, `data/literature/<key>.json`, is named after an entry of the
bibliography and records what was taken from it: named quantities with their
unit, kind and location, tables with a unit per column, how the transcription
was checked, and remarks. The new exported function `literature(key)` reads
and validates a file into a `LiteratureRecord`, whose quantities are `Traced`
values carrying their provenance; `literature_value` and `literature_table` are
the shortcuts a calculation needs, and `available_literature`,
`literature_path` and `LITERATURE_SCHEMA` complete the interface. A malformed
file is an error naming the file and the field, never a silently different
number, and a unit is unit arithmetic over the registry of `DynamicQuantities`,
evaluated without `eval`. Editing a file recompiles the package, so a stale
value cannot survive in the compiled image. A test checks every file against
`refs.bib`, key and DOI.

Coefficients that define an equation of state or a published model (the water
equation of state, the HKF constants, the Pitzer α and b) are the model rather
than data about it, and stay in the code with their source.

### The rate-law constants come from their sources

`PK84_PARAMS_*`, `WALLER_PARAMS_*`, `PK_BLAINE_REF`, the Powers ratios behind
`powers_alpha_max` and the parameters of the deprecated `PK_PARAMS_*` are now
read from `Lavergne2018.json`, `Waller1999.json`, `ParrottKilloh1984.json` and
`Powers1948.json`. The names are unchanged and every value is identical to the
literal it replaces. Moving them recorded what had not been recorded:

- Tables 3 and 4 and the fly-ash parameters of Lavergne et al. (2018) were
  checked against the article, page by page.
- The slag time of the Waller law, 100 days, has been attributed to Waller's
  thesis since the law was added. It does not appear in Lavergne et al., and
  the thesis has not been checked; the file records it as unstated, and the
  docstring of `WALLER_PARAMS_SLAG` says so.
- The parameters of the deprecated smoothed variant have no established
  source and are recorded as such.

### Output that shows only the result

- **Ipopt no longer warns on every solve.** OptimizationBase warned, at each
  solve, that Ipopt needs second derivatives and that it was building a
  `SecondOrder` backend itself; the certified search runs this back end many
  times per equilibrium, and the documentation printed the warning 209 times.
  The extension now passes that same pair explicitly, so the computation is
  unchanged.
- **Banners and progress bars only on a terminal.** `read_thermofun_database`,
  `build_species` and `build_reactions` drew boxed titles and progress bars on
  every call. Written to a file or captured into a document they are no longer
  drawn.
- A kinetics run no longer triggers the warning of SciMLBase about parameters
  of mixed types. The thermodynamic functions of the species are heterogeneous
  by nature; they are held behind a wrapper, and how they are called is
  unchanged.
- The warning of `pe` for a couple with a member at the solver floor was
  printed with runs of spaces inside its sentences.

### Documentation

- Pages that build their own `IpoptOptimizer` silence its iteration log, which
  ran to 1700 lines on the equilibrium tutorial. They also passed `abstol`,
  which Ipopt ignores, and the manual stated that it was forwarded; it now says
  that `reltol` becomes the Ipopt `tol` and that `abstol` has no counterpart.
- A warning that repeats what a printed table already reports, such as a
  refused certificate or a potential set by the solver floor, is no longer
  shown beside the table. The diagnostics of the kinetic runs, which no table
  repeats, are kept in folded sections under the results.
- The scripts shared by the ionic-hydration and calibration pages were
  included two or three times into the same module, redefining every
  documented name; their guards now test the module that includes them.

Nothing in the API or in any computed value changes, and the new functions
only add to the interface. The one change of behavior is the one above: the
database readers no longer print banners and progress bars when the output is
not a terminal.

## v0.22.1 — OptimaSolver again when Ipopt is loaded, and a documentation that reads in order

The documentation build had passed two hours, and looking for where the time
went found two defects of the solver that had nothing to do with the
documentation. In any session with Ipopt loaded, the OptimaSolver back end ran
on a degraded path; and a stage of the certified search had been unwired by
accident, a loss the first defect happened to hide. The rest of the release
answers the other complaint about the documentation: a reader new to chemical
thermodynamics found no way in.

### An `OptimaOptimizer` solve took the generic path whenever Ipopt was loaded

`OptimaSolverExt` solves an `EquilibriumSolver` built on an `OptimaOptimizer`
with the exact conservation matrix and the exact gradient, and its comments
measure what happens without them. `OptimizationIpoptExt` defines the same
`solve` for every back end, and the OptimaSolver method was meant to outrank it.
It did not: its signature, `EquilibriumSolver{F, <:OptimaOptimizer, V} where
{F, V}`, left the struct's parameters free of the bounds the struct declares,
and Julia does not then rank it above the bare `EquilibriumSolver`. In any
session that loaded `Optimization` and `OptimizationIpopt` besides
`OptimaSolver` -- the equilibrium tutorial does, and so does the documentation
build -- every `OptimaOptimizer` solve went through the generic
`OptimizationProblem`, with the conservation matrix rebuilt by finite
differences and the gradient of `dot(n, μ(n))` taken by automatic
differentiation. On the 109-species CEM IV paste its start lay 2e-3 mol away
from the right one, and a certified solve that costs 0.7 s cost 34 to 78 s.

The signature now bounds its parameters, `EquilibriumSolver{<:Function,
<:OptimaOptimizer}`. With both extensions loaded `which` names
`OptimaSolverExt`, and since the test suite does not load Ipopt, a test asserts
the ranking against the generic signature directly.

### The certified search pre-solves under ideal activities again

The stage that solves the problem under ideal activities and restarts the
non-ideal search from that answer -- the fix of the ulp sensitivity of a CEM I
with eight solid solutions, recorded in the docstring of `_ideal_start` -- had
been removed in v0.18.0 by the commit that removed the linear-programming
start, although that commit argued for removing the LP alone. The degraded
OptimaSolver start above happened to give the dual solve a start it could use,
and so masked the loss. With that path corrected, a cold CEM IV paste without
ash no longer certified from either back end, and the points of the page that
depend on it were refused in turn. The stage is restored where it stood, runs
only when nothing has certified yet, and from its answer the paste certifies at
3.6e-15. Nothing in the test suite had noticed the removal, since the only test
reaching the stage used a budget no stage can rescue; that paste, solved from
its cast state, is now a test of its own. Replayed with the packages the documentation loads, every certified
and every refused point of the cement pages is again the one their text
describes.

One page described a refusal that was the defect's and not the chemistry's: on
the CEM I with three solid solutions declared and the other end-members left as
pure phases, the search now certifies the answer of the complete declaration,
and the section says so, with the correct reason -- a pure end-member cannot
outcompete the solution it belongs to, and a partial phase list certifies
nothing about the solutions it leaves out.

### The certified search no longer solves the same start twice

`equilibrate_certified` offers its cached starting points again after the ideal
pre-solve, after the continuation, after each restart and in each repair round,
and every offer ran the dual Newton and the certificate once more from a start
already solved. `solve_certified` now takes a `memo`, an `IdDict` keyed on the
start object, and `equilibrate_certified` keeps one for the duration of a call.
The dual solve and the certificate are deterministic once the solver, the
budget, `ϵ` and the constraint are fixed, so the record is the result: the test
suite checks that a repeated call returns the very objects of the first, and an
in-process comparison on the CEM IV paste finds the same composition bit for bit
with and without the memo. The saving is modest, about a tenth of the time of a
refused point.

### `ΔₐG⁰` is an apparent energy of formation, and now says so

The builder documentation and the Manual described `ΔₐG⁰` as "the Gibbs energy
of formation". It is the apparent one: the elements are held at the reference
temperature (Benson-Helgeson convention), so that `ΔₐG⁰(Tr)` equals the `ΔfG⁰`
supplied, and away from `Tr` the two differ by terms of the elements that cancel
from any balanced reaction and from nothing else. The code was right; the
glosses were not.

### The documentation

- **Reading paths.** The home page names three entry points by what the reader
  already knows, and *Getting started* is rewritten around one calculation,
  with the solubility constant obtained from the database before any solver
  runs.
- **Prerequisites and exits.** Every page outside the API opens with what it
  assumes and ends with where to go next.
- **Tutorials that are tutorials.** `tutorials/equilibrium.md` and
  `tutorials/kinetics.md` keep the calculation; the options -- back ends,
  constraints, derivatives, activity models, solid solutions, rate functions,
  rate constants, kinetic reactions, calorimeters -- move verbatim to two new
  Manual pages, *Solving an equilibrium* and *Writing a kinetic model*. Anchors
  keep their names, and the URL cited by an earlier section of this file,
  `tutorials/equilibrium/#sec-aqueous-properties`, still resolves.
- **Theory.** A new page, *Standard states*, gives the convention of each class
  of species as the code implements it, the passage between the molality and
  mole-fraction scales (a constant `RT ln(m°M_w)`, −9.96 kJ/mol at 25 °C), the
  pressure the heat-capacity polynomials do not carry, and the case where a
  reference energy stops being a convention. *Thermochemistry* gains the
  apparent energies, checked against the function the package builds, Euler's
  theorem and the Gibbs-Duhem relation, and the structure of the HKF model;
  *Proving that an answer is the answer* gains stable, metastable, partial and
  local equilibrium, and the comparison between mass-action and minimization
  formulations. Anderson and Crerar (1993) is cited by section, and Tanger and
  Helgeson (1988) is added to the bibliography, its DOI resolved on Crossref.

### The documentation build

Measured with the packages `docs/make.jl` loads, in one process and without
threads, the executed blocks of the whole site take 45 minutes on the author's
machine. The calibration page is the heaviest at 14; CEM IV, which paid a full
cascade for each refusal the ideal stage now avoids, went from 468 s to 335 s
under the same conditions. The continuous-integration runners are slower, but
integrate the coupled trajectories on four threads.

A pull request no longer waits for the whole site. The draft pre-flight of
`docs/make.jl` resolves every cross reference of the whole tree in minutes, and
five partial builds, run side by side, execute every page between them; the
groups are in `docs/shards.jl`, balanced on measured cost, and a page no group
names falls into the last. The full build, the only one that deploys, runs on
`main` and on demand. `scripts/docs_timing.jl` produces the per-block cost map
the groups were balanced on, loading the packages `docs/make.jl` loads and in
the same order, since which extensions are active decides the route a solve
takes.

## v0.22.0 — a sorbent that appears and disappears

The site budget can now follow the phase that carries it. Until this release a
capacity was a number posed once; a C-S-H that precipitates as a paste hydrates
carries its sites with it, and a budget that cannot move is a statement about a
quantity that no longer exists.

Three defects found by an external audit are fixed first, because the feature
leans on all three.

### A named host is the host that was named

`SurfaceSupport.host` was resolved by symbol, correctly, and the symbol was then
thrown away in favor of the formula. The kinetics index registers both, so a
formula shared by two polymorphs is written twice and keeps whichever came last.
Calcite and aragonite are both `CaCO3`; holding 1 and 100 mol, a law asking for
`Cal` got Arg's amount and a reactive area a hundred times too large.

Nothing changes for an ordinary species, which takes its formula as its symbol.
It differs exactly where the two differ, and a species with no symbol at all is
now required to have an unambiguous formula rather than resolving to whichever
sibling was declared last.

### A family member is matched by identity, not by its label

`_resolve_site_families` looked its members up by symbol against the system's
own vector and checked nothing about what it found. The file carried a comment
asserting that a shared member was unreachable. It is reachable: the species a
label lands on need not be the species the family validated, so two families
with **different** site symbols and the same member labels resolved to the same
indices and one family's conservation row was simply absent from the matrix.

The same gap let a family hold the `AS_SURFACE` copies it qualifies while the
system kept the caller's unqualified originals, leaving `idx_surface` empty with
the site mixing still running.

### A declared capacity is compared with the state that uses it

`site_moles` had no caller anywhere in `src`. The site row took its right-hand
side from the initial amounts like every other row, so declaring `1e-6 mol` of
sites and initializing `1e-3` gave an equilibrium holding `1e-3`, with
`optimality_certificate` reporting `optimal = true`. The certificate was not
wrong — it checks the budget it is given. Changing the capacity alone changed no
number, which also makes calibrating one impossible.

`site_budget_residual`, `check_site_budget` and `host_consistent_state` close
that. The tolerance is relative and its default comes from a measurement: a
state initialized with its occupied sites at the `1e-12` solver floor is
correctly initialized and still `2e-9` off, so `1e-9` would fire on the ordinary
case.

### The site budget follows its host

`SITES_FOLLOW_HOST` on the support, a capacity measured per unit mass or per
unit specific area, and the budget becomes `ν` moles of sites per mole of host —
equation (30) of Kulik (2002) read as a coefficient. Which capacities qualify is
**measured** rather than listed: `sites_per_host` evaluates the capacity at two
scaled host amounts and requires the budget to scale with them, with `n₀` held
fixed, since scaling it too makes a `ShrinkingCoreArea` ratio equal one
everywhere and hides the nonlinearity the probe exists to find.

The coupling is one entry of the constraint matrix, and which entry is not a
free choice. Two plausible routes are wrong, both measurably:

Subtracting from the row of the **free site** subtracts that primary's whole
composition, because the rows are indexed by species and `XsOH` carries an
oxygen and a hydrogen. At Dzombak and Morel's weak-site density that invents
seven percent of the oxygen of `Fe(OH)₃`. Repairing it needs a preimage of the
pure site pseudo-element, and there is none — least-squares residuals of `0.378`
on an amphoteric oxide and `0.500` on a cation exchanger, structural rather than
a quirk of one basis.

Appending a row **instead** of replacing one leaves two equations on one
quantity, which together say the host may not dissolve at all. Measured: the
solve returned `MaxIters`, the host moved regardless, and the site total stayed
at its initial value.

The **bare** site as the component removes the obstruction instead of working
around it. It is a component and not a substance, so it need not be among the
species; the residual is then zero, the preimage is the unit vector, and element
conservation is exact by construction. It carries the charge the free site
carries with its site symbol — `XsOH` is `Xs⁺ + OH⁻`. A neutral one leaves the
charge row among the primaries, the two appear in one ratio everywhere, and only
their sum is identifiable: measured, the multipliers ran to `±2.3e5` while their
sum stayed at `−60`, the dual Newton stalled, and the host came out thirteen
percent wrong.

Measured, on portlandite carrying sites at three host amounts: converges as well
as the uncoupled solve, holds the constraint to `10⁻⁷`, conserves its elements
to `10⁻¹⁵`, and reports every present phase at `log SI = 0` to `10⁻¹⁴`.

### A species at the activity floor was worth `ln 2` more than the floor

Molalities are built from `max.(n, ϵ)` and were then passed to `log(mᵢ + ϵ)`.
Two regularizations stacked, so any species sitting at the floor came back at
`log(2ϵ)` instead of `log(ϵ)` — in the HKF, Davies and Pitzer closures, while
the dilute model, which takes the log bare, was right.

`ln 2` on a species nobody looks at would be harmless. A **primary** can sit at
the floor, and `saturation_indices` reads each element potential off its primary
species, so the offset reached every phase carrying that element. Measured on
amorphous ferric hydroxide, where `Fe³⁺` at pH 7 is a `10⁻¹⁶` species: the
solid, present and at equilibrium, reported `log SI = 0.298` — which is
`ln 2 / ln 10` — under HKF and Davies and `0` under the dilute model.
`optimality_certificate` reported `optimal` and was right; it reads the solver's
own multipliers. The index was the thing that lied, on every non-dilute model
the package ships.

The second `ϵ` is gone. No number moves anywhere else: the molality is already
strictly positive, so the term only ever did something at the floor.

### What a coupled family costs, and the one thing it does not yet settle

`SITES_FOLLOW_HOST` is a constraint, and as a constraint it is exact. It is not
yet a complete thermodynamic model, and the gap has a size.

With a fixed budget the free site's `ΔₐG⁰` cancels out of every surface
reaction, so zero is free: shifting a whole family by 20 kJ/mol moves nothing by
more than `3e-11`. Coupled it does not cancel, because the host carries `−ν` of
the site component — and `ΔₐG⁰ = 0` on a free site is not a gauge but a claim,
namely that a surface hydroxyl forms from the elements for nothing. `XsOH`
carries a real oxygen and a real hydrogen. At Dzombak and Morel's weak-site
density, `ν = 0.2`, that claim is worth 8.3 log units on the host's own
solubility, and measured, it dissolves an amorphous ferric hydroxide outright
where the same system with a fixed budget holds its solid.

The reference is not a convention to choose: it is the energy of the matter the
free site carries, `μ°(H₂O) − μ°(H⁺) = −237.2 kJ/mol` for an oxide, read off the
same matrix the constraint is built from — so `host_coupling_bias` computes it
for any free site, an exchanger's included. Set it, and the coupling costs
nothing measurable: the host keeps `9.999993e-4 mol` against `9.999693e-4` with
a fixed budget, the site total is `ν` times the host amount to seven digits, the
solve certifies, and across the whole band the guard allows the answer moves by
`7e-8`. Leave it at zero and the family is refused at construction, with the
value to use in the message.

Kulik (2002) reaches the same place from the other side, and the theory page now
says so: he keeps the free site out of the balance entirely, as a *surface
monolayer solvent* of fixed activity with `μ_n = 0`, and carries the capacity in
a surface activity term. That formulation needs no reference energy at all. This
one does, and now states it.

### A charge component is refused on a measurement, not on its presence

The first form of this guard refused a coupled family whenever `Zz` survived
among the primaries. That is a symptom and not the defect. Charge stays an
independent component whenever it is independent of the element rows — one
element in two oxidation states is enough — and such a system is perfectly well
posed: measured on hydrous ferric oxide in a mixed-valence iron chloride
solution, the coupled matrix has full identifiable rank and its charge row is
nonzero on a ferrous complex, so it is not the site row at all. Worse, the
charge the refusal then suggested was itself refused on the next call, so the
two suggestions pointed at each other.

The decision is now the identifiable rank of the matrix the solve will be
constrained with, read off its singular values, and the message reports what it
measured. It refuses everything the old one correctly refused — the tightest of
those has a spectral gap of 13 — and admits the redox systems it wrongly did,
the tightest of which sits at 2.2.

### A saturation index that agrees with the stationarity the solver reached

Two rules that were right for the charge row and wrong for a coupled site row.
A primary absent from the species was given zero potential — which shifted every
surface index by more than forty log units, `y_site` being `+95.7`. And the
index was formed with `SM.A` rather than the matrix the solve was constrained
with, leaving the host's own `−ν` out of its own index: portlandite came back at
`log SI = 0.0015` while present and at equilibrium, where a present phase is
zero by definition.

The missing potential is not a convention: the free site is a species and is
always present, so its stationarity determines it exactly.

### The site-density scale a constant was fitted at

An intrinsic adsorption constant is fitted at some total site density and its
value depends on that choice, so two constants fitted at different densities are
not comparable — including two from the same paper.
`convert_logk_site_density` is Kulik's equation (21), and a doctest reproduces
his published `−0.73` and `−2.33` on Dzombak and Morel's own two densities.
`REFERENCE_SITE_DENSITY_NM2` is the `12.05 nm⁻²` he writes; the SI form is
derived from it through Avogadro rather than written twice.

### A worked example of the thing itself

`examples/evolving_sorbent.md` titrates hydrous ferric oxide carrying Dzombak
and Morel's own weak-site density until it is gone, with manganese on it. It is
written for a reader who has never done surface complexation: what a site is,
why its budget is a conservation row and not an element, and why posting the
budget as a number fails on a solid that dissolves — you would need the answer
in order to set up the question.

Two things happen in one sweep and telling them apart is the point. Up to about
`0.2 mol` of acid per mole of oxide nothing dissolves and the acid simply takes
the manganese off the surface, which a fixed budget describes perfectly well.
Past that the oxide goes and its sites go with it, the ratio holding at `ν` to
six decimals. Every number on the page comes out of a block the build executes,
including the three checks it closes on: the relation to `1.5e-7`, element
conservation to `5.1e-12` and every present phase at `log SI = 0` to `1.0e-11`.

### PHREEQC as the oracle

PHREEQC has coupled a `SURFACE` to an `EQUILIBRIUM_PHASES` mineral since v2.
`test/reference/phreeqc_evolving_surface.py` titrates a sorbent to exhaustion
and the site totals are the declared coefficient times the phase amount to
`2.0e-10` across five partially dissolved states.

Two traps were caught by assertions rather than by reading, both leaving a
plausible table: PHREEQC punches the initial solution too, so the output carries
one row more than there are reaction steps; and `-molalities` are per kilogram
of water, which is not 1 kg once acid has been added. The generator also records
the PHREEQC **engine** version, which the older surface generator does not.

### The eliminating kinetic route keeps what the parent declared

`coupling = :species` rebuilt the free side from its species and the parent's
primary names, dropping the solid solutions and the site families. A surface
therefore could not be combined with that route at all, before this release and
independently of it.

### Breaking changes

Below 1.0 the registry treats a minor bump as breaking whatever the API did, so
`[compat] = "0.21"` will not accept `0.22` and any dependent must widen its
bound. No package in this organization depends on ChemistryLab, so there is
nothing else to change.

Beyond that: eighteen exported names are new and none was removed, but four
behaviors change deliberately. A rate law's host is now found by symbol rather
than by formula, which is a different species exactly when two share a formula
and one carries a distinct symbol. A site family whose members do not match the
system's species by formula, site symbol and aggregate state is now refused
where it was accepted. And a state whose site amounts contradict its declared
capacity is refused where it was solved. A site family that follows its host and
whose free site is left without a reference energy is refused where it was
solved, and a saturation index computed on a primary that has collapsed to the
activity floor now reports the equilibrium instead of the floor.

## v0.21.0 — charged surfaces, published models, and numbers that say where they came from

Two threads, and they meet.

**A surface acquires charge**, so this release adds the electrical work of
putting one more charge on a charged object — in both forms the literature uses
— along with cation exchange, and it reads two published models from the files
their authors released rather than rebuilding them here.

**And doing that ran into the other thread.** A published model is a body of
fitted constants, and what makes one usable by somebody else is not the numbers
but what travels with them: which measurement each came from, how well it is
known, and — once a parameter is identified rather than looked up — whether the
data determined it at all. So the release also carries `Traced`,
`identifiability`, and the two observables that make an identification possible.

### Cation exchange, in the convention it was fitted in

`VanselowMixing` and `GainesThomasMixing` differ in what fraction an exchanger's
activity is: a fraction of *particles* or a fraction of *charge equivalents*.
For a homovalent exchange they coincide; for Ca²⁺ against Na⁺ they do not, and
the two are related by the exchanger's own composition — which is what the
calculation is solving for, so no constant converts between them. This package
therefore **declares** the convention and converts nothing implicitly. A
constant fitted under one and used under the other is a different model.

Checked against Reaktoro, which carries both: agreement 4.1e-9 on each.

### A charged surface, two models, and the unknown only one of them needs

`ConstantCapacitance` and `DiffuseLayer` both write the electrical work as
`z_k ψ̃` and differ only in the closure that gives the potential.

For a constant capacitance, `σ = CΨ` inverts to something linear in the
composition, the charging work is a genuine quadratic potential with a positive
semi-definite Hessian, and **the certificate covers it unchanged** — a proof,
not an expectation.

For a diffuse layer, Gouy-Chapman's relation is transcendental in `Ψ` but
monotone in it, so it inverts too: `ψ̃ = 2 asinh(σ/κ√I)`. Writing it that way is
exact and the solver cannot always follow it — the inner loop that recovers a
mixing phase reads `lnγ` at the previous iterate, a fixed point that contracts
only while `electrostatic_stiffness` stays below
`ELECTROSTATIC_STIFFNESS_LIMIT`. So the potential is carried as an **unknown of
the outer Newton** instead, with the closure as its equation. Over eighteen
points spanning three ionic strengths and six pH values, the eliminated route
certifies three and the unknown certifies all eighteen, agreeing with PHREEQC
to 2.1e-4 on a site fraction.

The threshold is measured on both sides rather than chosen: 3.38 certified,
6.18 not, nothing between. It also explains a number v0.20.0 recorded as a
conditioning artifact — a capacitance below about 3 F/m² on ferrihydrite fails
for exactly this reason, not that one.

### A surface carries one potential, not one per family

Ferrihydrite is two families — strong sites and weak sites — on one oxide, so a
proton bound to a weak site charges the same surface a proton bound to a strong
site does, and both feel the same `Ψ`. The potential therefore belongs to the
**support**, and the charge that raises it is summed over every family standing
on it. `support_group` is that grouping; a support whose families disagree about
having a potential, or about the area of the surface they share, is refused.

This is what makes the published calibration reachable. Dzombak & Morel's zinc
edge, run with the layer on and its ionic strength matched to PHREEQC's:

| pH | no layer | with layer | PHREEQC, with layer |
|---:|---------:|-----------:|--------------------:|
| 6.0 | 0.3732 | **0.1825** | 0.1816 |
| 7.0 | 0.9044 | **0.7094** | 0.7069 |

The layer moves the edge by about half a pH unit, which is why constants fitted
with it cannot be used without it, and the agreement is 0.47 % — slightly better
than the 0.94 % of the same comparison without the layer.

### What the certificate means for each, stated rather than implied

Measuring two criteria that one word had been hiding separated them:

| | is it a gradient? | is that gradient extensive? |
|:--|:--|:--|
| ideal site mixing | yes | yes |
| constant capacitance | **yes** | no |
| diffuse layer | **no** | no |

Only the first column is the certificate's. A constant capacitance fails the
second because its charging work is homogeneous of degree two in the amounts
while a Gibbs energy is homogeneous of degree one — the area is a parameter,
like a volume, and extensivity in the amounts alone is not what a surface of
given area has. A diffuse layer fails the first, and that one is a real cost:
its potential depends on a bulk ionic strength that does not depend in return
on the surface, so the activity Jacobian is asymmetric and is the Hessian of
nothing. That is the Dzombak-Morel approximation itself, which PHREEQC's
default `SURFACE` block makes too; what comes back is a self-consistent
speciation, not a certified minimum, and `is_gradient_consistent` says which.
`site_gradient_asymmetry` measures it rather than asserting it.

Stacking two electrostatic models is refused. Adding two potentials is how a
Stern model is *drawn* and not how it works: the capacitances belong to
different charge planes, and summing them puts every species on both.

### SIT, the model the published compilations are written in

`SITActivityModel` implements the Specific ion Interaction Theory: Debye-Hückel
with a deviation term summed over ion **pairs** rather than one coefficient
times the ionic strength.

That difference is the point. Davies says the correction depends on how much
salt there is; SIT says it depends on which salt, and
`ε(Na⁺,Cl⁻) = +0.03` against `ε(Na⁺,SO₄²⁻) = −0.12` are of opposite sign. It
also means `γ` turns around instead of falling forever — in NaCl at 25 °C,
`log₁₀γ(Na⁺)` reaches −0.174 near `I = 1` and is back to −0.156 at `I = 3`, which
a one-parameter deviation term cannot reproduce while matching the dilute end.

It matters here beyond being one more option: the NEA reviews and ANDRA's
ThermoChimie database are **calibrated in SIT**, so a `log K` taken from them and
used under Davies is not the constant that was fitted — the same category of
error as using a surface constant fitted with a diffuse layer in a model without
one.

Checked against PHREEQC on the formula alone: the fixture carries the molality
*and* the activity coefficient for every species, so the comparison involves no
speciation, no database of `log K` and no convergence path. Agreement
1.4 × 10⁻³ in `log₁₀γ` over four decades of ionic strength, the residual being
the Debye-Hückel slope rather than the coefficients, which are the same numbers
on both sides.

**No compilation ships with the package.** `build_sit_parameters` reads the
`SIT` block of a PHREEQC-format database the caller already has; the one
PHREEQC distributes is the ANDRA/RWM ThermoChimie-TDB, which is not USGS-
authored and is not ours to redistribute. Reading rather than transcribing is
also what keeps five hundred coefficients from becoming five hundred chances to
mistype one.

`Traced` carries **where each coefficient came from**, and
`missing_epsilon_pairs` reports which of a system's pairs are resting on the
literature's convention that an unlisted `ε` is zero. A convention silently
applied is indistinguishable from a coefficient somebody determined.

And it differentiates **with respect to the coefficients**, not only the
composition: `∂ln a(Na⁺)/∂ε(Na⁺,Cl⁻) = ln10 · m(Cl⁻)` exactly. That is what
makes identifying an interaction coefficient from data a well-posed question
rather than a finite-difference exercise, and getting there needed the element
type to follow the parameters — the promotion trap that lets every piece pass
its own test while the chain breaks.

**A coefficient for a neutral solute now reaches the calculation.** The matrix
of `ε` excluded every pair involving a neutral, so the closure's neutral branch
could never receive one — a compilation that carried a neutral-ion coefficient
had it dropped on the way in, and the result was ideal with nothing said. SIT
still leaves a neutral ideal *by default*, which is the absence of a
coefficient and not a rule imposed here; what changed is that "unless the
compilation says otherwise" is now something that can happen. Found by a
coverage report: the branch was unreachable, which is why no test covered it.

A side effect worth recording: SIT gives the diffuse-layer comparison a third
aqueous model, and it **rules the aqueous model out** as the cause of that
comparison's residual. Ideal 2.12 × 10⁻⁴, Davies 2.07 × 10⁻⁴, SIT
2.09 × 10⁻⁴ — two per cent apart. Together with the dielectric constant and the
surface species set, that is three explanations measured and discarded; the
residual stays small and unattributed, which is what the measurements support.

### Reading a published sorption model

`read_sorption_model` reads the `SURFACE_MASTER_SPECIES`, `SURFACE_SPECIES`,
`EXCHANGE_MASTER_SPECIES` and `EXCHANGE_SPECIES` blocks of a PHREEQC-format
database into site families, exchangers and reactions — **keeping what travels
with each constant**. A published sorption compilation writes it into every
line:

```
Am+3 + Ilt_sOH = Ilt_sOAm+2 + H+   # … error: 0.26 ref: Marinich_ea:2024:rep:
```

and a reader that took the `-log_K` and dropped the rest would be discarding the
part that says how much to believe it. Each `log K` therefore arrives as a
`Traced`, its source the `ref:` tag and its uncertainty the `error:` one.

On ClaySor 2023 that is 187 reactions across six edge-site families and four
exchangers, all published — and **60 of the 187 state no uncertainty at all**,
which `provenance_report` says in one line and which no amount of reading the
numbers would reveal.

`Traced` gains an `uncertainty` field for this, which is what a second real
customer is for. `nothing` there means the source says nothing, which is not the
same as saying the value is exact.

**Nothing ships.** ClaySor 2023 is CC-BY-4.0 and freely available from its
Zenodo deposit; this reads the copy a user has. And reading a model is not being
able to solve it: one is written against a particular *aqueous* database — ClaySor
names PSI/Nagra TDB 2020 in its own first lines — and its constants are that
database's, in exactly the way a surface constant fitted with a diffuse layer is
not the same constant as one fitted without.

One parsing trap, because it produces a plausible wrong answer rather than an
error: `+` is a charge as well as a separator, so splitting
`Ca+2 + 2 MntxNa = Mntx2Ca + 2 Na+` on the character yields three terms, none of
them the calcium ion. The separator is a plus with whitespace on both sides.

### A published clay model, run end to end

ClaySor 2023 — the 2SPNE SC/CE model of Bradbury and Baeyens, two-site
protolysis plus cation exchange on illite and montmorillonite — is the first
model in this package that somebody else published and that is read from their
own file rather than rebuilt here.

It exercises something no previous case did: **four site budgets on one solid**,
three edge families counting particles and an exchanger counting charge
equivalents, competing for the same solution. Against PHREEQC on the
Na-montmorillonite subset:

| what | worst relative gap |
|:--|--:|
| edge sites (protolysis) | **2.5 × 10⁻⁶** |
| exchanger (Na/Ca) | **2.0 × 10⁻²** |

Four orders of magnitude apart, and structurally so. The edge sites see only the
proton, whose activity is prescribed on both sides, so no aqueous model enters
and what is compared is the surface model alone. The exchanger sees the sodium
and calcium activities — and that is a measurement, not an excuse: the identical
system under ideal activities is off by 23.8 % and under Davies by 1.95 %, a
factor of twelve from changing nothing but the solution.

**It is not a reproduction of ClaySor, and the fixture says so as a field.**
ClaySor's constants were fitted against PSI/Nagra TDB 2020, which this
repository does not have; both sides run the sorption model over `phreeqc.dat`'s
aqueous chemistry instead. Using published constants over a different aqueous
database is a different model, in the same way a surface constant fitted with a
diffuse layer is a different constant from one fitted without — and the whole
point of saying which is that both look like the same number.

### A number that says where it came from

`Traced` attaches a [`ProvenanceKind`](@ref) and a source to a value, so what a
data file knew about a number survives into a table, a figure or a fitted
result. This is not a new idea in this repository — `data/pitzer-reardon1990.toml`
already carries `origin = "estimated:<analog>"` on every coefficient Reardon had
to borrow, and explains why in its own header — it is that idea promoted from
data into a type.

Six kinds, ordered from the weakest claim to the strongest, and two of them are
the point:

  - **`PROV_UNSTATED` is the weakest, not the middle.** A number that forgot to
    say where it came from must never strengthen a result.
  - **`PROV_FITTED` is not evidence.** A fitted value may be excellent, and
    whether it is depends on the data, the model and whether the parameter was
    identifiable at all — three questions a predicate cannot answer. So
    `is_evidence` is true for `PROV_MEASURED` and `PROV_PUBLISHED` only.

`Traced` is deliberately **not** a `Real`. A value that flowed silently into
arithmetic would arrive at the far end with its history gone; unwrapping it is
an act a reader can see. `weakest` is what a derived quantity can honestly claim
about itself, and `provenance_report` is what to print beside a result — a table
where nine coefficients are published and one is a placeholder is a different
object from one where all ten are, and the difference shows in none of the
numbers.

`SITActivityModel` is the first model built on it: a borrowed `ε` keeps its own
standing instead of inheriting its compilation's. `docs/src/manual/where_the_numbers_come_from.md`
gains the section that says when to reach for it.

### Which parameters a measurement can actually determine

`identifiability` answers, for any forward model and any parameter vector, the
question that has to be settled before a fitted number is quoted. It is the
reasoning `scripts/hydration_calibration.jl` worked out for calorimetry, taken
out of that one script and made to work on anything: the singular values of
`∂y/∂log θ`, the parameter correlation matrix, and the linearized standard
errors.

Three deliberate choices in it.

**Against the logarithm**, so a rate constant and a dimensionless exponent are
comparable — a matrix mixing `∂y/∂k` with `∂y/∂n` has a spectrum that says more
about the units than about the data.

**A rank read off a gap, not a threshold**, because a threshold has units and a
gap does not. The spectrum that calibration measured, `[420, 100, 60, 6.3, 1.4,
0.20]`, has its largest ratio between the third and the fourth — which is why it
fits three parameters and not six. That case also **set the default**: the ratio
there is 9.5, so a round threshold of 10 would have answered six on the very
spectrum the rule exists for, and the default is 5.

**And it closes onto `Traced`.** `as_traced` returns the fitted parameters
carrying `PROV_FITTED`, the dataset, and the standard error as their
uncertainty — and marks `PROV_PLACEHOLDER` any parameter lying mostly in the
directions the data do not constrain, because a number the data did not
constrain is one the optimization had to put somewhere, not one it determined.

Which parameters those are is a **subspace** question — `null_participation` —
and not the parameter's position relative to the rank. The rank counts
directions; the order is the caller's packing. For two parameters entering only
as their product the participation is `[0.5, 0, 0.5]`, so neither is determined
alone, whereas reading the rank by index would have cleared the first and
flagged the second: the wrong answer by a plausible route.

Tested on a model with a collinearity **put in on purpose** — two parameters
entering only as their product — so the method has a known right answer rather
than a plausible one. It finds rank 2 of 3, names the trade-off direction as
equal and opposite in exactly those two, and reports a correlation of 1.

And then it was pointed at a model in this package rather than at a synthetic
one, which changed two things.

**The cut goes at the largest qualifying gap, not the first.** Measuring a rate
law against its own shrinking-core exponent gives `[1.6e-5, 2.7e-6, 5.1e-11]`: a
ratio of six, then one of fifty thousand. Cutting at the first answered *one*
determined direction — the amplitude alone — when the amplitude and one exponent
combination are both determined and only their split is not. A factor of six is
ordinary conditioning; a factor of fifty thousand is a structure.

**And a coarse difference step can hide an exact degeneracy.** The default 5 %
moves that law's exponent by 0.165, far enough that the second-order differencing
error differs between two parameters that are exactly collinear: the condition
number comes out 80 instead of 3 × 10⁵, and the correlation −0.97 instead of
−1.000. Refining to 1 % restores both. `log_sensitivity` now says so, with the
habit it asks for — refine the step and see whether the answer moves — and the
manual shows the sweep rather than describing it.

`where_the_numbers_come_from.md` gains the section on what to do when the number
does not exist yet, both the synthetic collinearity and that real one, ending on
the caution that a fit is not a mechanism: adjusting a site density can absorb a
denticity and still fit.

### What a solid assemblage loses on heating

`ignition_loss` and `bound_water` compute, from the formulas the database
already carries, the total a thermogram integrates to. Thermogravimetry is the
second observable the calibration example asks for by name, and the reason is
identifiability: calorimetry constrains three combinations of six kinetic
parameters, and a measurement that sees the *phases* breaks correlations heat
cannot.

It counts **hydrogen**, not formula water, and the difference is not pedantry —
portlandite is `Ca(OH)₂`, has no `H₂O` written in it, and loses one water per
formula unit. A rule searching for `H₂O` would report zero for the second most
abundant hydrate in a paste. The aqueous phase is excluded, which is the
distinction the water-budget page exists for: pore solution is water and is not
bound water.

**What it is not is a thermogram, and that is stated rather than approximated.**
Turning the total into a curve needs the temperature window each phase releases
in, and those are literature values rather than consequences of a formula. This
package does not carry them; `bound_water_per_phase` is the quantity they would
attach to, which is why it is exposed per phase. Supplying them completes a TGA
observation operator, and inventing them would fabricate exactly the part of the
measurement that does the identifying.

### A thermogram, and the windows it takes to have one

`ignition_loss` gives the total. `thermogram` gives the curve, which is what
identifies phases — C-S-H, AFt and AFm all release below 200 °C and are told
apart by the shape of the release, not by its size.

A window is **not** a consequence of a formula, so `DecompositionWindow` carries
both its parameters as [`Traced`](@ref) values and makes a curve say where they
came from: a publication (`PROV_PUBLISHED`), a measurement, or a placeholder
that keeps printing as one. A bare number is `PROV_UNSTATED`, the weakest claim
there is.

**And they are recoverable from a curve**, which is what closes the chain.
`thermogram` is written as a smooth function of its windows so that
`window_parameters` hands them to an optimizer and `identifiability` says
afterwards which of them the curve determined. On three separated peaks, a
Gauss-Newton started 40 K and 40 % away recovers them exactly and the rank is
6 of 6.

**Where it stops is the useful part.** Two phases releasing 5 K apart — the
ordinary case in a paste — give a rank of **1 of 4** and a condition number of
332: the curve sees one peak with a position and a width, not two with four
parameters between them. A fit would still return four numbers, and
`as_traced` marks as `PROV_PLACEHOLDER` every parameter lying mostly in the
unconstrained directions — for exactly that reason.

Two things are stated rather than smoothed over. A phase with no window
contributes to the starting mass and never leaves, so `phases_without_windows`
reports it rather than letting a curve integrate quietly to the wrong total. And
a logistic has infinite tails, so the curve does not start at zero; no
renormalization happens, because rescaling a curve to start at zero would put
the discrepancy somewhere a reader cannot see.

### An area that follows the grains, and the measurement that says what it is worth

The rate laws here scale with the binder's fineness, and that factor has always
been a **number**: computed once from the fineness given at construction and
carried unchanged through the integration. A dissolving grain does not keep its
area, so `parrott_killoh_avrami` and `waller` now also accept a
`ShrinkingCoreArea`, and then the factor follows the amount left.

What is new is not only the mechanism but what comes with it. Multiplying the
shell-formation branch `k₃(1-ξ)^n₃` by `(1-ξ)^p` gives `k₃(1-ξ)^(n₃+p)`: on that
branch the new exponent and the old one are **the same parameter written
twice**, and fitting both fits a sum. The Jander branch carries no such
exponent, so there `p` is a genuine new shape; and the Waller sigmoid moves its
two exponents together while `p` moves one, so there `p` is distinguishable from
`n`. Which of these a given dataset sees is a measurement, and `identifiability`
— added in this same release — is what makes it.

The published constants were fitted with the factor frozen, so an evolving area
leaves the calibration it came with. That is said in the docstrings, in the
theory chapter and here, rather than left for a reader to discover from
parameters that no longer mean what their source said.

Nothing changes for anyone who does not ask for it: at `n = n₀` the evolving
factor is exactly one, and a frozen factor now takes a branch that returns the
same object rather than multiplying by one. Both are asserted bit-for-bit.

### Reproducibility, which turned out to need the most work

Three habits were removed from this repository, all of them silent and all of
them found by review rather than by a test.

Physical constants were written as literals in sixteen files while
`src/utils/constants.jl` exists to take them from `DynamicQuantities`. Standard
Gibbs energies were retyped from databases the package ships, and had already
drifted — water read `-237181` against slop98's `-237183.0`, portlandite
`-897010` against CEMDATA18's `-897013.0`. And `55.5 mol` of water, which is
0.99983 kg rather than 1 kg, was the basis of every molality.

`docs/src/manual/where_the_numbers_come_from.md` is new, written for a reader
who has not met these traps: how to ask the package for a constant, a molar
mass, a kilogram of water or a standard energy, which databases exist and why
not to mix them, and the three kinds of number it *is* right to type.

The cross-code fixtures move from Julia blocks pasted into test files to
`test/reference/*.json` with their provenance as fields. PHREEQC's `phreeqc.dat`
is vendored — unmodified, from tag v3.7.3, with the USGS User Rights Notice
beside it as that notice requires — so the md5 the fixtures record now pins a
file in this repository instead of one inside somebody's conda environment.

### What this release does not do

**The support is still fixed.** A site budget is a number set once, and the
destination — sorption on a C-S-H that precipitates as the paste hydrates —
needs it to follow the phase. The data model was written for that from the
start: `site_moles` takes the host's amount rather than returning a constant.
What is missing is not the data model but two things that only show up when the
support is allowed to move.

The first is where the sites come from. A capacity written as a site density
times an area, or times a dry mass, is **linear** in the host's amount, so the
site conservation row stays linear and folds into the budget exactly — no
bilinear term, contrary to what one expects. But a free site is a chemical
species: `XsOH` carries an oxygen and a hydrogen. Growing the support therefore
creates surface hydroxyls, and those have to be debited from the water rather
than appearing from nothing, or the element balance is violated by exactly the
number of sites added. Any outer iteration that updates the budget between
solves has to carry that bookkeeping, and an implementation that only updates
the site row would balance the sites and silently unbalance the oxygen.

The second is whether sorption feeds back on the support's own stability. Giving
the host a coefficient in the site row shifts its saturation index by the site
potential, which is thermodynamically the statement that a sorbing surface is
more stable than a bare one. That is a real effect and a real modeling decision,
and it is not one to make implicitly by writing a matrix entry.

Both are named here rather than left as an absence, because the architecture
makes them look like small steps and they are not.

**No multidentate species, no lateral interactions.** Unchanged from v0.20.0:
the site balance holds for any denticity, ideal mixing of occupied and free
sites is exact only for one, and anything else is refused at construction.
Frumkin and quasi-chemical mixing are the route, and neither is here.

**The evolving fineness is a mechanism, not a calibration.** It is delivered
with the measurement that says when its exponent is identifiable and when it is
a rewriting of a parameter the law already had — and with no fitted value,
because fitting one against constants calibrated at a frozen area would be
making up a number.

### Breaking changes

**Below 1.0 the registry treats a minor bump as breaking whatever the API did,
so `[compat] = "0.20"` will not accept `0.21` and any dependent must widen its
bound.** No package in this organization depends on ChemistryLab, so there is
nothing else to change.

Beyond that the release is **additive**: sixty-one exported names are new, none
was removed, and no existing signature or number moves. A system that declares
no surface is bit-for-bit unchanged, and so is a rate law given a plain fineness
— the evolving area is reached only by asking for it, and equals the frozen
factor exactly at the initial amount.

## v0.20.0 — chemistry that happens on a surface

An equilibrium can now have part of its chemistry on an interface. That is the
one thing the formulation could not express at all: alkali held by a C-S-H,
chloride bound without forming Friedel's salt, a radionuclide retained by an
oxide — all of them described until now, when they were described, by a bulk
proxy.

Nothing about a system without a surface changes. No exported name, signature or
number moves for it.

### Breaking changes

`ChemicalSystem` gains a sixth type parameter for its surface declarations, so
code spelling its parameters out — `ChemicalSystem{T, R, C, S, SS}` — must add
one. Partial parameterizations such as `ChemicalSystem{C, S}` are unaffected.

`AggregateState` gains `AS_SURFACE` and `Class` gains `SC_SURFCOMPLEX`, both
**appended**, so no existing member changes its integer value. Code that
enumerates either exhaustively now sees one more.

`molar_mass(::KineticReaction)` raises when no reactant carries a molar mass,
where it used to return **0.1 kg/mol**. That default was silent and
multiplicative — it scales the reactive area and therefore the whole rate, and
could be wrong by an order of magnitude. The same default is gone from the
internal helper the rate factories use.

And below 1.0 the registry treats a minor bump as breaking whatever the API did,
so a downstream bound of `"0.19"` will not accept this release.

### Added — sites, and the machinery that was already there

A **site** is a conserved quantity that is not a chemical element. The package
already had one of those: electric charge, carried as a pseudo-element in a
species' formula and turned into a conservation row by the ordinary matrix
assembly. Site families are carried the same way, each with its own
pseudo-element from the twenty-four reserved `SITE_SYMBOLS`.

The consequence is that **not one line of `stoich_matrices.jl` changed**. On an
oxide with one family and its protolysis, the assembly produces

```
  H2O@   [1, 0,  1, 0, 0, 0,  0]
  H+     [0, 1, -1, 0, 0, 1, -1]
  Ca+2   [0, 0,  0, 1, 0, 0,  0]
  XsOH   [0, 0,  0, 0, 1, 1,  1]      <- the site balance
```

and the rest of the matrix is the protolysis itself.

`SiteFamily` groups the species sharing one budget, free site included, with an
`AbstractSiteCapacity` in one of three forms — a site density per square meter,
a capacity per kilogram of dry support, or a prescribed total. Three because
published data comes in three shapes and converting between them needs a number
nobody measured: turning a clay's exchange capacity into an area density would
mean inventing a BET area to divide by.

### Added — Langmuir, as a consequence

The activity of a species on a site is its **site fraction**, and that is the
whole model. Competitive Langmuir is not implemented; it falls out:

```math
n_j / N = β_j / (1 + Σ_k β_k),   β_j = K_j a_j
```

with the shared denominator that is the signature of a finite capacity. No
isotherm is written anywhere in the package.

Verified against that derived form to a relative `1e-6` across four pH values,
and against **PHREEQC** computing Dzombak & Morel's model from its own database
to **1.7e-8** on the protolysis of hydrous ferric oxide.

Kulik's surface activity coefficient `1/(1 − θ)` is deliberately **absent**: it
exists so that an *eliminated* free site can reproduce what an explicit one
already does, and applying it here would count saturation twice. The identity is
asserted instead, to `1e-8`.

### Added — one area abstraction for the kinetics and the surface

Two notions of surface lived in the kinetics and never met: `AbstractSurfaceModel`,
a total area with no caller at all, and `blaine_factor`, a dimensionless ratio
frozen at construction and used by every cement run.

They are one family now, and the **measurement is the type**. The package warned
three times in prose that a BET area is not a Blaine fineness — silica fume is
about 20 000 m²/kg by the first and 2 000 by the second, a factor of ten on the
hydration rate — and nothing enforced it. `area_ratio` is defined only between
two areas of the same kind, so that substitution now raises instead of returning
a plausible number. `blaine_factor` is that ratio, and its numbers do not move.

`total_area(model, n, n₀, M)` carries the initial amount from the start, because
an area that cannot see where it began cannot follow a microstructure that
evolves. `ShrinkingCoreArea` is the consequence — area as `(n/n₀)^p`, with
`p = 2/3` for spheres and `p = 1` recovering the mass-proportional law exactly.
The power law is regularized so its derivative stays finite at exhaustion, which
an integrator does reach.

`SurfaceSupport` names the host solid once, which is what let the 0.1 kg/mol
default go.

### Added — coupled to the kinetics

The equilibrium is re-solved on a rebuilt system at every accepted step, and
that rebuild takes its declarations by keyword: anything not passed is dropped
silently, as solid solutions were until 0.8.2. Site families now travel with it,
on the same all-or-nothing rule, and a family split between the partitions is
refused by name rather than dropped.

A surface species is not aqueous, so the automatic partition would have made it
kinetic — taking its family's conservation row with it. It no longer does.

Checked on a solid releasing calcium at a constant rate onto an inert sorbent:
the site budget holds at every reported instant, the released calcium is exactly
`k·t` and all of it is accounted for, the binding reaction's law of mass action
holds on the **integrated** state, and tightening the integrator by two decades
moves the answer by under `1e-5`.

### Added — two families, and what a cross-code difference measures

Dzombak & Morel's ferrihydrite has strong sites and weak ones, forty times more
numerous and a thousand times less avid for a metal. That contrast is what makes
a sorption edge bend, and it needs two families on one support.

Against PHREEQC on the zinc edge, `-no_edl` on both sides, the worst departure
is **0.94 %** — a hundred million times worse than the protolysis comparison.
The difference is **not** in the surface model, and the suite measures that
rather than claiming it: the identical system under ideal aqueous activities is
49 % off, under Davies 0.94 %. A factor of fifty, from changing nothing but the
aqueous phase. Zinc is divalent and its activity coefficient leaves 1 quickly.

That is also why the protolysis comparison carries no metal: with only protons,
and the proton activity prescribed on both sides, no activity coefficient enters
and the surface chemistry is compared alone.

### Documentation

A theory chapter written for someone who has never done surface chemistry, with
the derivation of Langmuir in full and an explicit account of what the model
does **not** cover — no surface potential, no evolving support, no multidentate
species, no exchange conventions.

Two worked cases, both executed at build time: one family titrated by pH against
the closed form, and the two-family zinc edge against PHREEQC. Their figures
carry numbers the build computed.

An oracle bench under `test/reference/`, with PHREEQC and GEMS3K alongside the
Reaktoro generators that were already there. GEMS3K **loads and cannot run**: it
needs a system export that the ClaySor deposit, which ships a GEM-Selektor
database, does not contain. That is recorded rather than worked around — a
tolerance floor invented from a comparison nobody ran would be worse than none.

### What this release does not do

No surface potential, so a parameter set fitted **with** an electrostatic term —
including Dzombak & Morel's own published calibration — is being used outside
its model if used here. The comparisons above run `-no_edl` on both sides for
that reason, and are cross-code checks rather than reproductions.

No evolving support: the site budget is fixed. No multidentate species: the site
balance holds for any denticity, but ideal mixing of occupied and free sites is
exact only for one, so anything else is refused at construction rather than
returned as a number that looks like an isotherm.

## v0.19.0 — two polymorphs are two species, and a dimensioned argument is an error again

Three defects a test suite could not see because the assertion guarding each one
tested something adjacent to it; a global that leaked into every package loaded
beside this one; and a measurement of which activity models actually come from an
excess Gibbs energy.

### Breaking changes

- **`isequal` on species now compares the symbol as well**, so two polymorphs are
  two species. Calcite and aragonite are both `CaCO3`, both `AS_CRYSTAL`, both
  `SC_COMPONENT`; in CEMDATA18 their standard Gibbs energies differ by 821 J/mol,
  which at 298 K is 0.33 in `ln K` — the whole difference in solubility between
  them. They compared equal. Code relying on two differently-named species of one
  formula being `==` now sees them as distinct, which is the point.
- **`AggregateState` gains `AS_LIQUID`**, appended, so no member is renumbered.
  Metallic mercury in `slop98-inorganic-thermofun.json` now comes back
  `AS_LIQUID` instead of `AS_UNDEF`, and `instances(AggregateState)` has one more
  entry.
- **`log`, `exp` and some thirty other `Base` functions no longer accept a
  `DynamicQuantities.Quantity`.** Loading ChemistryLab used to add those methods
  globally. Code that passed a dimensioned value to one of them was getting an
  answer; it now gets a `DimensionError`, which is what it should always have
  got. Form the dimensionless ratio and take the logarithm of that.
- **Below 1.0 the registry treats a minor bump as breaking whatever the API did**,
  so a downstream bound pinned to `"0.18"` will not resolve `0.19` and must be
  widened.

### Fixed — a `Dict` could not tell two species apart, and said so inconsistently

`hash` mixed in the symbol while `isequal` ignored it, so `isequal ⟹ hash` was
broken in both directions at once: `calcite == aragonite` returned `true` while
`Dict(calcite => 1)` raised `KeyError` on `aragonite`.

That disagreement had already cost a day. CI was green on Julia 1.13 and red on
1.12 from the same commit, because the `Reaction` constructor strips
pseudo-species with `delete!`, which finds a key by `hash` and confirms with
`isequal`; whether it reached `ELECTRON` through `Species("e")` depended on where
the table put them, hence on the Julia version. The root cause is gone.

Identity is now formula, aggregate state, class, and whatever the symbol adds to
those. A symbol that merely spells the species' own formula adds nothing and is
dropped, so `Species("H2O")` and `Species("H₂O")` remain one species — and so do
`ELECTRON` and `Species("e")`, and `Species("Ca+2")` and `Species("Ca⁺²")`. The
stored symbol is untouched: it is also the lookup key, and `cs["H2O"]` must keep
working.

Two designs were tried first. Canonicalizing the stored symbol broke that lookup
key outright. Canonicalizing inside the comparison, to `unicode(formula(s))`,
looked right and is not: `unicode` is **not constant** on the equality classes of
`Formula`, which are composition and charge — `Formula("e") == Formula("e-")`
while their spellings are `"e"` and `"e⁻"`. Both were found by running the suite.

### Fixed — a dimensioned argument to a dimensionless function

`thermo_factories.jl` extended `Base` math functions to `Quantity` by stripping
the unit first. Type piracy by construction, so the methods were global. And it
did not only remove a check: `ustrip` returns the value in SI **base** units, so
`log(1u"mol/L")` gave `6.908` and `log(1u"mmol/L")` gave `0.0` — the same
concentration scale a thousand apart, two different answers, nothing in the
source to show it.

`_adim` restores the check by dispatch, so `Float64`, `ForwardDiff.Dual` and
`Symbolics.Num` take the identity method untouched and only a genuine `Quantity`
pays for the test. A dimensionless quantity passes, which is the case the
extensions existed for. The Aqua exemption that declared the piracy is gone with
the methods, which is what stops them coming back.

Measured: Aqua 10/10 with no `piracies` entry, and a warm cement equilibrium at
2.4 ms against 2.5 ms before, same session and same script.

### Fixed — a printing keyword that was assigned and never read

`pprint` built its row labels from `col_label`; `row_label` was assigned on the
line above and never used, so asking for one accessor on the rows and another on
the columns labeled both by the column's. Two `eval`s go with it — one reading a
numeric literal the term regex had already constrained, one turning a printing
keyword into a function — neither a vulnerability, both replaced by code that
cannot be anything else.

### Measured — which activity models come from one excess Gibbs energy

The consolidation roadmap asked for an analysis of reported discrepancies between
chemical potentials and the gradient of an energy built from them. A sharper test
than the one planned answers it: if `μᵢ = ∂G/∂nᵢ` for any `G`, second derivatives
commute, so `∂ln aᵢ/∂nⱼ = ∂ln aⱼ/∂nᵢ`. That is read off the Jacobian of
`ln_activities` by automatic differentiation — no `G` is built, no solver runs,
nothing is differenced against a second solve.

Worst relative asymmetry over pairs, ion/ion and solvent/ion split:

| model | ion/ion | solvent/ion |
|:--|--:|--:|
| Pitzer | 1.3e-15 | 1.3e-15 |
| B-dot, `å` per ion, `Ḃ ≠ 0` (the default) | 1.8e-01 | 2.9e-02 |
| B-dot, common `å` **and** `Ḃ = 0` | 0 | 1.8e-14 |
| Debye-Hückel limiting law | 0 | 3.4e-13 |
| `å = 0`, `Ḃ = 0.0976` (the GEMS setting) | 1.2e-01 | 3.6e-01 |

The limiting law is exact, which is the classical result, and Pitzer is exact,
which is what it claims for itself. It takes **both** corrections to break it:
dropping `Ḃ` while keeping per-ion radii makes the ion/ion asymmetry worse. The
reason is one line — `∂ln γᵢ/∂nⱼ = ln10 · f′(I; zᵢ, åᵢ) · zⱼ²/(2 kg)`, so symmetry
demands `f′/zᵢ²` not depend on `i`; a per-ion `å` puts `i` in the denominator, and
`Ḃ·I` is added with the same coefficient to every ion, so it contributes `Ḃ` and
not `Ḃ·zᵢ²`.

Neither is a defect of this implementation: both are the published extended form,
which every geochemical code uses. The theory page states the consequence where a
reader meets it — when the activities are not the gradient of a single `G`,
`optimal = true` is a stationary point of the residuals rather than a minimum of a
potential — and `test/activities.jl` pins the zeros.

### Changed — a solid-solution parameter says which convention it is in

The literature writes Redlich-Kister coefficients two ways: CEMDATA18 gives the
AFm sulfate/hydroxide binary as `A₀ = 0.188, A₁ = 2.49` in RT units, while this
package takes J/mol and divides by `RT` internally. The factor is `RT ≈ 2478`, and
a bare `Float64` cannot say which was meant.

`RedlichKisterModel` and `RegularSolutionModel` now accept a `Quantity`:
`u"J/mol"` passes, `u"kJ/mol"` converts, `u"J"` raises. A bare number still means
J/mol, so no existing call changes. `excess_ln_gamma_expression` returns the
symbolic `ln γₖ`, and a test asserts it agrees with the compiled path — which
makes the formula in the docstring true by construction rather than by
proofreading. The compiled path itself is untouched.

`log10_gamma_expression` does the same for the activity kernels, and buys more
than inspection: the derivation behind the table above — that Maxwell demands
`f′(I; zᵢ, åᵢ)/zᵢ²` not depend on the ion — is now **checked** rather than only
its consequence measured. Differentiated symbolically, the ratio agrees between
charges to `1e-14` with a common `å` and no `Ḃ`, and disagrees as soon as either
correction is present. `_log10γ_ion`, the compiled half, gains the docstring it
never had.

### Fixed — pure water, where the assertion that passed was not the one that mattered

An audit reported the Schur-complement route failing the autoprotolysis test on
macOS ARM64 while CI stayed green on Linux. It is not a platform difference.
Measured on Linux: both ions come out about twelve times too large, so the ion
product is off by a factor of 140 — while their **ratio** stays within half a
percent of one, and the ratio was all this route asserted. The element balance
closes exactly either way, so the route returns a point that conserves matter and
is not an equilibrium, and it reports `MaxIters` while doing so. The missing half
is a `@test_broken`, which turns red the day it is fixed.

### Added — Windows in the CI matrix, and what had to change for it

`windows-latest` joins the matrix on the current release. Three things had to
change first, and one was a defect: `display_data_path` decided a file was outside
the package by asking `relpath` and reading `..` off its answer, which on two
different drives it does not produce — so a path plainly outside could have come
back rewritten, into the banner every documentation page showing a database load
captures verbatim. It is a prefix test now. Line endings are normalized to LF by a
`.gitattributes`, and the partial documentation build falls back to a copy where a
symlink cannot be created.

### Changed — the shipped databases, measured rather than assumed

PR #59 left two questions open by name. Both are settled. Of seven distinct unit
strings in `data/*.json`, two do not parse — `cal/(mol*bar)` with 1768
occurrences and `kbar` with 109 — and **neither reaches the parser**:
`eos_hkf_coeffs` is converted by an explicit SUPCRT-to-SI table, written because
the JSON's own unit metadata is wrong for `a3` and `a4`, and `m_expansivity` is
not read by this package at all. Every unit on a field that *is* read parses,
fractional exponents included, and a test now walks the files and says so.

The same guard is applied to the classifications. `extract_classification`
matches a label by name and falls back when there is no match, and a fallback is
a valid value, so a label the enum does not carry is data the import throws away
without saying so -- `AS_LIQUID` was one. A test now requires every
`aggregate_state` and `class_` on the `substances` of every shipped database to
resolve to a member, and both enums document the correspondence: the four
aggregate states ThermoFun uses here, and the four substance classes. `ELEMENT`
and `CHARGE` appear only in the `elements` section, so `Class` is right not to
carry them.

Adding a member also turns every hand-written list of "all of them" into a
filter. `idx_speciation` carried one -- four aggregate states named when there
were four -- and reads `instances(AggregateState)` now. Its class list stays
explicit, because `SC_SSENDMEMBER` is excluded on purpose.

### Fixed — a docstring that documented nothing, and the guard that let it through

`raw"""` is a string **macro**, not a string literal, so Julia does not attach
it: a definition below one has no documentation at all, silently. Two docstrings
in this release were written that way, to escape the backslashes of a LaTeX
block. `@doc raw"""` attaches; they use it.

The cost was not the mistake but the delay in seeing it. `check_docrefs.py`
stripped the `raw` prefix and credited the docstring, so the local guard said
nothing, and Documenter reports an unresolvable `@ref` at its `CrossReferences`
stage — **after** every example on the site has run. Ninety minutes of CI to
learn it.

Both ends are closed. The script refuses a bare `raw"""` outright and names it.
And `CHEMLAB_DOCS_PREFLIGHT_ONLY=1` runs the static checks and a draft build and
stops: 61 s when the tree is clean, 84 s to reject that same `@ref`, against the
ninety minutes it took before.

### Changed — the documentation can be checked in minutes

Not a change to the package, and the reason it is here: verifying that an example
still runs used to cost over two hours. Doctests have their own CI job and finish
in 28 s, and `CHEMLAB_DOCS_ONLY` builds only the pages a change touched — one
manual page in 1 min 49 s. A partial build states what it does not check and
refuses to deploy; the full build is untouched and remains the gate.

### Compatibility

No compat bound moves. `OptimaSolver` stays at `"0.5.5, 0.6"`.

## v0.18.1 — a database file was a program, and one solve could silence the next

Two defects that had nothing to do with chemistry and everything to do with
whether an answer can be trusted: **a database file was executed rather than
read**, and **a flag one solve set could be left behind for every solve after
it**. Around them, the work needed to check any of this without paying two hours
for it.

### Fixed — imported metadata is interpreted as data, not as code

Reading a ThermoFun database executed what the file said. Not as a corner case of
malformed input: the ordinary path through `read_thermofun_database` evaluated
two of its fields as Julia source.

**Classification labels** — `aggregate_state` and `class_` — went through
`eval(Meta.parse(...))`. They are now matched by name against the corresponding
enum, which is what they always were.

**Unit strings** went straight to `uparse`. DynamicQuantities maps the
*arguments* of a call expression into its unit registry but leaves the call
**head** alone, and the head is then evaluated inside a module which — like every
Julia module — sees `Base`. So a unit field reading `write("…", "…")` was not
rejected before it ran, and did not even raise: the call executed, its return
value was a number, and a `Quantity` came back as though the field had been an
ordinary unit.

The whole expression tree is now validated before `uparse` is reached. Numbers,
registered unit symbols, the arithmetic operators, roots and named constants are
accepted; every other call, and every macro, assignment, index and block, is
refused and takes the default-unit fallback that unparsable fields already took.
`safe_uparse`, an unguarded wrapper on the same function that nothing but its own
test called, is deleted rather than hardened.

**Nothing about the bundled databases changes.** They read identically — 876
assertions across the import, zeolite, thermodynamics and HKF suites — and the
scientific questions this deliberately does not touch stay open: a unit the
installed DynamicQuantities does not know, such as `cal`, still falls back
silently rather than warning, and the label `AS_LIQUID` (one species in
`slop98-inorganic-thermofun.json`, metallic mercury) is still not in the
`AggregateState` enum and still resolves to `AS_UNDEF`. Both belong to a separate
piece of work on the validity of imported values, not to a fix whose whole point
is to be narrow.

### Fixed — a solve's policy no longer escapes into the next one

Two module-level flags were saved, written and restored in a `finally`. That is
correct while no two of them overlap, and two solves on two tasks do overlap.
Measured, with A entering the scope first and leaving first while B was still
inside its own: B lost the flag inside its own region, and — the half that
matters — B's restore put back the value it had saved, which was A's `true`, so
the flag was left **set** after both had finished. `_EXPLORING_STARTS` exists to
silence the non-convergence warnings of a start search, so from that point on
every solve in the session swallowed its own. A real failure stopped announcing
itself.

`_EXPLORING_STARTS` is internal and is now a `ScopedValue`: dynamic extent,
inherited by child tasks, invisible to siblings, and nothing to restore.

**`STRICT_CONVERGENCE` keeps its `Ref`**, because it is public and documented and
a `ScopedValue` cannot be assigned. What the package needed was never to change
that setting but to *suspend* it while computing a starting point, and it now
does exactly that, scoped to one task. A caller's `ChemistryLab.STRICT_CONVERGENCE[] = true`
behaves as before and is never written by the package.

`NONCONVERGED` becomes a `Threads.Atomic{Int}` — `[]` reads and writes unchanged,
which is its whole documented interface, but concurrent increments can no longer
be lost. Its docstring claimed `integrate` reset it and reported the total;
nothing in `src/` resets or reads that counter, and the docstring now says so.

### Fixed — Ipopt no longer buries the output it was called from

`_default_ipopt_solver` set no print level, so Ipopt ran at its default 5: a
banner, a full iteration table and a summary **per solve**. This back end is one
of several in a multi-start cascade that can run it hundreds of times in a single
call. On the quickstart page, one equilibrium produced 110 lines of output, about
80 of them Ipopt's, around the three the example exists to show; it now produces
28, with the same pH. A caller who wants the trace builds their own optimizer and
passes it, which makes seeing it a choice rather than an accident.

### Fixed — a printing keyword that was assigned and never read

`pprint` built its row labels from `col_label`. `row_label` was assigned on the
line above and then never used, so asking for one accessor on the rows and
another on the columns labeled both by the column's — silently, a wrong label
printing as happily as a right one.

Two `eval`s go with it, neither of which was a vulnerability and neither of which
is claimed to have been. `eval(col_label)` becomes a choice among four accessors:
a printing keyword needs that, not the power to evaluate an expression in this
module. `eval(Meta.parse(coeff_str))` read a coefficient the term regex had
already constrained to a rational, a decimal or an integer — so it called `eval`
at run time, inside a function, to read a numeric literal, and the type check
guarding its result could not fire for the same reason. The replacement returns
the same values *and* the same types on every form the regex admits.

### Development — the documentation can be checked in minutes

Not a change to the package, and the reason it is here: verifying that an example
still runs used to cost over two hours, so it was not done often enough to be a
check.

**Doctests have their own CI job.** `makedocs` ran them by default, so a
docstring example printing one digit differently failed a two-hour build. The
suite now runs in 27.6 s, and the preamble it shares with `docs/make.jl` lives in
`docs/doctest_setup.jl` — it had been written out twice and had drifted, the
commented-out job naming a package this one does not depend on while missing one
several doctests use.

**`CHEMLAB_DOCS_ONLY` builds only the pages you touched**, one manual page in
1 min 49 s. Filtering `pages` is not enough and was tried: Documenter walks the
*source* directory and builds every `.md` it finds, so a filtered page tree was
still inside `examples/` fifty minutes later. `docs/partial.jl` stands up a farm
of symlinks to `docs/src` without the pruned pages instead. Such a build states
what it does not check — links and docstring coverage — and refuses to deploy.
The full build is untouched and remains the gate.

### Compatibility

No exported name, signature or numerical result changes, and no compat bound
moves. A downstream package bounded on `"0.18"` resolves this release without any
edit.

## v0.18.0 — every cement in EN 197-1, and a documentation that builds

Four things a blended cement needs and this package did not have: a **redox
variable**, the **C-S-H that carries aluminum and alkalis**, a way to enter a
**glass** that has no phases to name, and an honest answer to *how much of it has
reacted*. With those, CEM I through CEM V all certify.

The last of the four is the one that changed the most results, and it is not a
refinement. A Gibbs minimization reacts whatever budget it is handed, so what it
is handed is the modeling decision: hand it the whole binder and it answers what
the paste becomes after every grain has dissolved, which no paste does. Two
related omissions had the same shape — a Bogue clinker carries no alkalis, so a
budget built from it returned a portlandite floor and called it a pore solution.

Around them, the documentation build: the coupled trajectories now run
concurrently, and the three stages that used to kill a three-hour build at its
very last step now check themselves at its first.

And underneath all of it, three defects in the **conservation laws themselves** —
a species projected onto components that cannot carry it, a rank decided by a
singular-value threshold, and a null space one vector short. Those are below the
cement chemistry rather than beside it: they decide what the certificate's
element balance is a balance *of*.

### Breaking changes

- **`ChemicalSystem` refuses two declared solid solutions that share a
  composition**, naming the pair. CEMDATA18 carries three descriptions of one
  C-S-H gel — `CSHQ`, `CNASH_ss` and the `ECSH` family — and declaring two of
  them counts the same hydrate twice. The overlap is exact, not approximate:
  `KSiOH`, `ECSH1-KSH` and `ECSH2-KSH` all carry `((KOH)2.5SiO2H2O)0.2`. A
  script that declared two of these families got an answer before and now
  raises; the answer it got was wrong.
- **`data/solid_solutions.toml` gains `CNASH_ss`**, so a script that loads the
  whole file and declares everything in it now also declares that phase — and,
  by the rule above, can no longer also declare `CSHQ`. Loading the file has
  always been an explicit act; which phases to declare remains the caller's.
- **`[compat] OptimaSolver` moves to `"0.5.5"`.** 0.5.3 carries
  `phase_split_measure`, without which a mixing phase that is present is
  certified on the stationarity of its members alone — and stationarity cannot
  see that the Gibbs minimum for a non-ideal phase is two coexisting
  compositions. The bound is **0.5.4 and not 0.5.3** because 0.5.3 shipped that
  test with a defect this package's own documentation exposed: the measure
  probed compositions the element balance forbids, reading a sentinel
  multiplier as a chemical potential. On the CEM III/A page it turned a
  converged equilibrium — element balance 3.8e-14, pH 12.489, assemblage
  unchanged — into `optimal = false` with a violation of +54.06, the same +54.06
  on every unrelated system carrying the same declaration. Worse than a wrong
  verdict: `equilibrate_certified` ranks its routes on that flag, so a false
  negative sent it through the whole cascade to return a *worse composition*.
  The bound is 0.5.5 rather than 0.5.4 for a second reason: **the documentation
  cannot be built without it.** Under 0.5.4 a warm cement equilibrium costs
  583 ms, and the site performs thousands of them; the build exceeded two hours
  and was killed by its own timeout. 0.5.5 takes that equilibrium to 17 ms.
  `OptimaSolver` 0.5.5 was registered ahead of this release for that reason.
  The bound now reads `"0.5.5, 0.6"`: 0.6.0 adds `split_trials` to its
  certificate and `split_starts` to `SolutionPhase`, which is what
  [`equilibrate_split`](@ref) needs to seed a second instance at all. Both are
  read through `hasproperty`/`hasfield`, so the package works under either — but
  the split loop has nothing to act on under 0.5.
- **A species that no declared component can carry is now refused**, where it
  used to be projected onto the components by the least-squares decomposition.
  A system that was built before and made moles of that species out of nothing
  now raises at construction, naming the species and the components. The answer
  it gave before was not a conservative approximation; it conserved the wrong
  thing. The smallest case is one line long, and this package's own test suite
  asserted it: `ChemicalSystem([H2O, H+, OH-], [H2O])` — water as the only
  component. Those three species span a **two**-dimensional space, so one
  component cannot express the other two, and the matrix came out reading
  `H+ = 0.4 H2O` and `OH- = 0.6 H2O`, which balances arithmetically and lets a
  solver make `H+` out of water with no `OH-` and no charge to pay for it. Two
  components are needed and the refusal says so. See the fix below.
- **The registry treats a minor bump below 1.0 as breaking whatever the API
  did**, so `[compat] ChemistryLab = "0.17"` will not accept `0.18` and
  downstream bounds must be widened. `MeanFieldHomogenization.jl` depends on this
  package only in `docs/Project.toml`.

### Fixed — three holes in the conservation laws themselves

The stoichiometric matrix is what every other answer rests on: it states which
quantities are conserved, and the certificate measures the element balance
against *it*. Three defects in how it was built are fixed here. None of them
announced itself, and two of them made the solver report a perfectly balanced
answer to the wrong problem.

**A species no component could carry was projected onto them instead of
refused.** The decomposition runs through `pinv`, and a projection never fails:
asked to write a species over components that cannot express it, it returns the
least-squares answer and says nothing. Measured — magnetite declared in a system
with no iron component came back as `4 H₂O@ − 8 H⁺`, a column that satisfies
every row of the matrix, and the equilibrium then made **2.4 mol of magnetite out
of a budget holding no iron at all**, with the certificate confirming the element
balance to 1e-11. It was right to: that balance was the one the matrix stated.
Membership in the span is now checked and a species outside it **refused by
name**, with the components listed and the reason given. The check is a rank
comparison in exact rational arithmetic, not a tolerance on the least-squares
residual, and that distinction is not stylistic: several CEMDATA18 formulas carry
decimal stoichiometry — jennite is `(SiO2)1(CaO)1.666667(H2O)2.1` — and the
parser keeps `5//3` on the calcium row while the oxygen row sums the decimal, so
a perfectly expressible species shows a numerical residual of 1e-6. A threshold
placed above that is a threshold, with everything that follows from one; a rank
comparison has none.

**A rank decided by a singular-value threshold.** Two booleans were read off
`rank(A; rtol = 1e-6)`: whether the charge row survives as a conservation law
independent of the elements, and which species are independent enough to be
components. A LAPACK SVD deciding a combinatorial question about integer element
counts is a threshold deciding something exact, and it decides differently on
different CPUs and LAPACK builds. The package's CI showed precisely that — the
same commit green on one Julia and red on another, with the sulfate/sulfide
half-reaction "balancing with no electron" because the charge row had been
dropped and the electron's column was identically zero. The rank is now computed
by exact rational row reduction and is the same everywhere.

**The null space was one vector short.** The exact elimination behind it was
fraction-free (Bareiss), on the strength of the theorem that its division is
exact — but the division is `÷` on `BigInt`, which **truncates** rather than
throwing, and it is not exact once pivots are skipped. On a 5×6 integer matrix of
rank 4 (smallest singular value exactly zero), two divisions left a remainder, an
eliminated column kept a stray entry, the pivot count came out 5, and the null
space came back with one vector instead of two. A null space one vector short is
a **missing conservation law**, not a slow answer. Both it and the rank now come
from one exact rational row reduction; validated on 1800 matrices — random
shapes, constructed rank deficiencies and rational entries — where the exact rank
matches the numerical one, the null space has the complementary dimension, and it
annihilates the matrix exactly.

### Fixed — `equilibrate_split` had never been executed

Two defects, both invisible to reading and both fatal on the first call, in a
function this release introduces. `optimality_certificate` rebuilds its return
value field by field and **dropped `split_trials`** — the incipient composition
Michelsen's analysis computes — so the loop read `nothing` and returned on its
first pass, silently doing nothing. Past that, the pass was accepted on
`cert.worst_violation`, a field this package's certificate does not have, so
reaching the line at all raised a `FieldError`.

The certificate now forwards `worst_violation_split`, `split_phases` and
`split_trials`, guarded by `hasproperty` so a back end that does not run the test
costs nothing. A pass is kept on the **KKT error** — the worst of stationarity,
element balance, supersaturation and the constraint residual — and not on one
residual of it, which is the same ranking mistake `_kkt_error` exists to prevent.
The seed also moves material **from the fuller instance into the emptier one**:
moving a share of an instance holding 1.5e-4 mol while its twin holds 5.1e-2 is a
perturbation of three parts in a thousand, and a seed that cannot move the answer
is indistinguishable from no seed at all. There is now a test that runs the loop.

### Added — a miscibility gap, detected and then represented

Detection and representation are separate problems, and 0.18.0 closes both.

**Detected.** With `OptimaSolver` 0.5.3, the certificate tests a **present**
mixing phase for wanting to split. Measured on the AFm sulfate/hydroxide binary
with the published Redlich-Kister parameters, whose spinodal is
x ∈ [0.631, 0.914]:

| model | certificate | worst violation |
|:--|:--|--:|
| ideal mixing | `optimal = true` | +1.0e-10 |
| published Redlich-Kister | **`optimal = false`** | +8.7e-03 |

The second used to certify. It was a KKT point and not a minimum, and nothing in
the output said which: stationarity is blind to the one failure that matters for
a non-ideal phase — that the minimum is two coexisting compositions rather than
the one reported.

**Represented.** `SolidSolutionPhase(...; instances = 2)` asks `ChemicalSystem`
for a second copy of each end-member, under a derived symbol (`monosulphate12#2`)
sharing the same thermodynamic record, so a composition vector carrying one
amount per species can hold material in either lobe of the gap or in both. Two
coexisting compositions of one substance need the substance to appear twice.

The duplication is done by copying the species rather than by letting two phases
share them, which keeps `ss_groups` **disjoint** — the activity assembly, the
mole-fraction fill and the certificate are unchanged. The duplicated column is a
copy of one already present, so the row rank of the conservation matrix is
unchanged and the copies start empty, leaving a budget computed as `A n`
unchanged too.

`instances > 1` is **refused for a convex model**, and not as a formality: two
instances of a convex phase are degenerate, every split of the amount between
them having the same energy, so the minimum becomes a flat manifold. Inside a
spinodal the common-tangent pair is unique. So the second instance is admitted
exactly where it is needed.

This is how GEM-Selektor represents the same thing — CEMDATA18 ships the AFm and
AFt binaries under two names each, so its users declare the binary twice. The
criterion is the same object in both codes: GEMS' phase stability index
Λ_k = log₁₀ Ω_k is, term for term, what `phase_split_measure` computes, derived
independently from the same KKT conditions [Kulik et al. 2013]. Neither code
splits a phase by itself; the difference is only that the duplication is asked
for here by a keyword rather than carried in the database.

**Located.** `common_tangent(model; T)` returns the pair itself — the two
compositions at which one straight line is tangent twice to the molar Gibbs
energy of mixing, equivalently where both end-members have equal chemical
potentials in the two phases. It is two equations in two unknowns, solved by
Newton with `ForwardDiff` derivatives, and it involves **no part of the chemical
system but the mixing model**, which is what makes it cost microseconds. This is
the construction PHREEQC uses for binary solid solutions, after
[Glynn & Reardon 1990]. Validated against an analytic oracle: for a symmetric
model the pair comes out at (0.070720, 0.929280), residual 2.3e-13, symmetry
exact to the last bit.

`miscibility_split(model, x̄; T)` then applies the lever rule: inside the gap the
two compositions are **fixed** and only their proportions move with the overall
composition, which is the property that makes the gap flat in a phase diagram.

**Found by minimizing, when the budget pins the composition — and not otherwise.**
This is the part that had to be measured rather than argued, and the first
measurement pointed the wrong way.

Where the element balance *fixes* the overall composition inside the gap, two
declared instances separate onto the common-tangent pair by themselves, and the
certificate proves it. Measured on a calcite/magnesite binary — two different
substances, so 0.025 mol of each pins x̄ = 1/2 whatever the energetics say —
with a Redlich-Kister gap, against `common_tangent` computed from the mixing
model alone and told to nothing in the solve:

| `a₀` | binodal | the two instances | certificate |
|--:|:--|:--|:--|
| 8 kJ/mol | (0.052846, 0.947154) | 0.0528 and 0.9472 | `optimal = true`, 1.5e-10 |
| 14 kJ/mol | (0.003662, 0.996338) | 0.0037 and 0.9963 | `optimal = true`, 1.0e-10 |
| 20 kJ/mol | (0.000315, 0.999685) | 0.0003 and 0.9997 | `optimal = true`, 1.0e-10 |

with the amounts in the proportions the lever rule asks for — 0.5010 against
0.4990 at x̄ = 1/2. That is the whole construction recovered by the
minimization, from a start that had everything in one instance.

The AFm binary of a real CEM I is **not** that case, and the difference is
chemical rather than numerical: nothing pins the phase's composition there. The
sulfate has somewhere else to go — ettringite — and the hydroxide is abundant,
so the composition is free, the two instances start at the same x and stay at it.
The symmetric state satisfies every first-order condition jointly, so it *is* a
stationary point and no descent direction leads away from it.

Which is why `common_tangent` exists as well: for a phase whose composition is
free, the pair is **computed** from the model rather than discovered by the
solver, and that is not a workaround but the standard construction — neither
GEM-Selektor nor Reaktoro asks a global minimization to find a binodal either.
`examples/miscibility_gap.md` shows the whole chain with the numbers: the
certificate refusing a single composition, the pair computed, and the lever rule
applied.

Also new: **`with_symbol`**, the same species under a different label, which is
what builds those copies.

### Added — how far a binder actually reacts, and under which cure

A Gibbs minimization reacts whatever budget it is handed, without comment. Every
cement page of this documentation used to hand it the **whole** binder, which
asks what the paste becomes after every grain has dissolved — a question about
geological time, not about a specimen at 28 days. On a CEM I at w/c 0.5 the
difference is small enough to ignore. On a CEM V at 48 % replacement it is not:
the full-reaction budget puts in alkalis and aluminum that no real paste
releases, and the answer comes back with a pH of 14.4, an element balance stuck
at 3·10⁻¹ and supersaturated layered double hydroxides it has no room to
precipitate. The failure was in the question.

- **`powers_alpha_max(w_c; curing = :sealed | :saturated)`** — the water ceiling
  now carries its curing convention. Powers gives both coefficients: **0.42**
  sealed and **0.36** under water, the 0.06 g/g between them being the chemical
  shrinkage, which a sealed specimen pays out of its own water and an immersed
  one draws from the bath. The one-argument call is unchanged and still sealed,
  so nothing that used it moves.

- **The ceiling applies to every constituent, not only to the clinker.** The
  coefficient is built from the water the hydrates bind and the water the gel
  holds at arrest; it is a property of the **pore space**, and water in a gel
  pore two molecules wide is no more able to reach a slag particle than an
  unhydrated alite core. Which is to say: nothing in the argument mentions
  clinker, and the pages no longer pretend it does.

- **For a glass the kinetic ceiling is the lower one, and it is now cited.** The
  RILEM TC 238-SCM round robin [Durdzinski2017] measured the degree of reaction
  of two slags and a siliceous fly ash in seven laboratories at exactly the
  geometry these pages use — Portland cement blended at 40 % and 30 %, w/b
  0.40 — and reports, at 28 days by SEM image analysis, 38–49 % for the slags and
  about 20 % for the ash. Against a water ceiling of 0.95 at that w/b, the glass
  is limited by its own dissolution and not by the water. The reacted fraction is
  the **smaller** of the two ceilings.

The CEM II, CEM III, CEM IV and CEM V pages now state a reacted fraction per
constituent and sweep it, and `theory/cement_water_budget.md` gains two sections
saying why the mechanism transposes and what a cure is. Two consequences are
visible in the results and are the point of the change: a CEM IV/A at 28 days
keeps about half of its portlandite where the full-reaction answer had exhausted
it, and **a CEM V/A certifies** where it previously had no admissible assemblage
at all.

A third is sharper still, and it is the aluminum question of the CEM IV page made
decisive. At 45 % replacement and a 28-day reacted fraction, `CNASH_ss` certifies
with an element balance of 1·10⁻¹² and **`CSHQ` does not**, stopping at 8·10⁻².
At 23 % the two models merely disagreed about where the aluminum sat; at 45 % the
model with no aluminum end-member has nowhere to put it.

- **`SaturatedCuring(V_ref)`** — and the specimen is genuinely opened, not only
  its ceiling moved. The mirror image of `CapillaryWater`: where the sealed
  constraint lets the saturation fall and lowers the water activity by the Kelvin
  term, this one holds the total volume at the fresh paste's and draws water in
  to make up what chemical shrinkage empties. One unknown, the water imbibed; one
  column in the conservation rows, so the system is open to water and closed to
  everything else; one equation, `Σ V̄ᵢ nᵢ = V_ref`, which is **linear** in the
  composition — so unlike `CapillaryWater` the problem stays a convex
  minimization on an affine set and `optimal` keeps its full meaning.

  The amount drawn in is not a numerical device: it is the **chemical
  shrinkage**, which is what a chemical-shrinkage test measures by watching a
  specimen drink — and it is computed from standard molar volumes rather than
  supplied. That makes it a check the package was never fitted to pass, and
  `examples/cement_wc_ratio.md` runs it. On a w/c = 0.40 paste:

  | | per gram of reacted cement |
  |:--|--:|
  | shrinkage from the standard molar volumes | 0.0606 cm³ |
  | water the cured specimen drew in | 0.0604 g |
  | Powers, as the gap between his two ratios 0.42 − 0.36 | 0.0600 g |

  The two internal routes agree with each other to 0.3 % — one is a difference of
  molar volumes, the other is the titrant the constraint had to admit — and both
  land within 1 % of a coefficient measured on pastes in 1948, with nothing
  fitted anywhere. The same page shows the `void` column at 0.0803 sealed and
  zero under water.

  The obvious objection is answered on the page rather than argued: Powers' 0.42
  does enter, through `α`, so is the agreement his number coming back out? The
  page moves `α` from 0.60 to 1.00 and the shrinkage per gram of reacted cement
  moves by **0.3 %** — flat, because the normalization divides by the same
  reacted mass. And it lists every hypothesis the number rests on: a closed
  14-species list, ideal molar volumes with no excess term, a dilute pore
  solution, the fresh volume as reference, and the fact that Powers' 0.06 is
  itself a difference of two separately measured averages. The obvious alternative, `FixedActivity("H2O@", 1.0)`, is wrong
  and the docstring says why — a cement pore solution sits near a_w = 0.98 from
  its salts alone, so prescribing 1 would imbibe without bound. A bath does not
  fix the activity inside the specimen; it fixes the availability, and that is a
  volume statement.

  What is still closed is a **coupled** run: the kinetic integration remains a
  closed system, so a cure enters it through `α_max` and not through its water
  balance.

### Fixed — a Bogue clinker has no alkalis, and the pH said so

Every blended-binder page entered its clinker through a Bogue composition, which
returns four phases — C₃S, C₂S, C₃A, C₄AF — and **no sodium or potassium**: they
are minor oxides, outside the four-phase decomposition. The pages nonetheless
declared `KSiOH` and `NaSiOH` as C-S-H end-members. Phases present, elements
absent: the package's own trap 4, taken from the wrong side.

Nothing failed. Every solve certified, with element balances at the level of
rounding. What came out was a pH of **12.510 on three CEM II pastes and a
CEM III alike**, to three decimals — because with no alkalis in the budget,
portlandite is the only thing setting the pH, and a portlandite buffer is by
construction insensitive to everything else. A real cement pore solution sits
above 13 for exactly the complementary reason: the alkalis dissolve almost
entirely and stay in solution while the calcium is held at that floor.

The diagnostic worth keeping is the **invariance, not the value**: a quantity
that does not move when the inputs move is either buffered or not computed from
them. The corroboration was the CEM V, whose fly ash carries 2.5 % K₂O — the only
page with alkalis in its budget, and the only one that did not return 12.510.

The clinker's alkalis are now part of the budget on all four blended pages, at a
usual industrial level (Na₂O equivalent 0.73 %), labeled `ASSUMED` like every
other composition not in the deposit, and carried at the clinker's own reacted
fraction since they leave the grain with it. The four pastes that all returned
12.510 now return 13.280, 13.261, 13.182 and 13.041, **ordered by how much
clinker was replaced** — the dilution comes out of the element balance with
nothing prescribing it.

Because that one input now carries the pH almost by itself, the CEM III page
**sweeps it** across the industrial range. A factor of three on the alkali
content moves the pH by 0.40 unit and the portlandite by under 5 %, all three
points certified — which is the separation stated as a measurement: the calcium
is held by portlandite and cannot follow, so the alkalis *are* the pH. The page
then says what to do about it: a pH quoted from these pages is worth what the
assumed alkali content is worth, so substitute a real analysis, and a durability
argument that turns on pore-solution pH cannot be settled by a calculation whose
alkali input was assumed.
`manual/choosing_species.md` records the trap and its signature — the invariance,
not the value.

### Added — oxidation state

Nothing in the package could hold sulfur at two valences, and no test anywhere
solved a multi-valence system. A slag-blended cement is exactly that problem: the
slag brings S(-II), the pore solution carries S(+VI).

Half of it was already there and never exercised — `StoichMatrix` keeps the charge
row as an independent component when an element appears at several valences, so
the oxidation state is conserved separately from the elements. What was missing
is the intensive variable conjugate to it.

- **`ELECTRON`** — the electron at the conventional standard state, zero for
  every thermodynamic function, as `H+` is. A convention, not a measurement, and
  the docstring says so: every potential computed from it inherits it.
- **`half_reaction(state, oxidized, reduced)`** — the couple balanced over `H+`,
  water and the electron. No coefficient is transcribed; they come from the
  element and charge balance.
- **`pe`** and **`Eh`** — the electron activity inferred from that
  half-reaction's `log K`, and the same number through Nernst.
- **`FixedpE`** and **`FixedEh`** — equilibrium at a prescribed potential, for a
  system genuinely open to a redox buffer. The titrant mechanism of `FixedpH`
  could not express it (there is no electron species to prescribe an activity
  for), so `_titrant_blocks` now takes a **linear combination** of
  log-activities; `FixedActivity` and `FixedpH` are its one-term case.

Validated against published half-reaction constants, computed here from
CEMDATA18's own Gibbs energies — so the agreement also checks that the two
datasets share a reference state:

| half-reaction | computed | published |
|:--|--:|--:|
| `SO4-2 + 9 H+ + 8 e- = HS- + 4 H2O` | 33.69 | 33.66 |
| `Fe+3 + e- = Fe+2` | 13.02 | 13.03 |

And with a prescribed potential, every point certified, the sulfur partition
moves three orders of magnitude over six pe units — `pe` read back through an
accessor that knows nothing of the constraint agreeing with the prescribed value
to 1e-3.

The documentation is explicit about the limit. Different couples need not agree,
and on one solution carrying both, iron reports pe = +13.0 while sulfur reports
−3.7. A paste has a single redox state only if its couples are at mutual
equilibrium, which on the time scale of hydration they are not.

### Added — the C-S-H a blended cement actually forms

`CNASH_ss` is declared, with its eight end-members. They were in the shipped
database all along; the phase was simply never declared. `CSHQ` has no aluminum
end-member at all, so with `CSHQ` alone the Al released by a slag or a calcined
clay has nowhere to go but the AFm/AFt phases and the aluminum balance comes out
wrong.

Its known degeneracy is documented rather than hidden: the eight end-members span
a space of rank 5, because Myers' model carries site constraints an ideal
eight-component mixture does not. The feasible set stays bounded, so a solve is
well posed, but the individual amounts are not determined by the element balance
alone — read the total and the ratios, not the eight numbers.

### Added — an oxide analysis as an element budget

`oxide_budget` is the entry route for a material with no phases. A clinker phase
has a formula; ground granulated slag, a fly ash and a natural pozzolana are
glasses, reported by their oxide analysis and by nothing else. Bogue does not
help — it inverts a decomposition over phases that exist.

`primary_decomposition` does the algebra and **refuses** above a 1e-8 residual: an
oxide outside the span of the primaries has no decomposition, and a least-squares
approximation of one would put elements into the budget that the oxide does not
carry. The analysis is **not renormalized** — a datasheet summing to 0.96 is
missing its loss on ignition, and scaling it to 1 invents material.

The docstring is equally explicit that a budget says what a glass *contains* and
nothing about what it does: a slag and a quartz sand of the same analysis give
the same `b`, and the degree of reaction is a kinetic quantity supplied from
outside.

### Added — zeolites, in a database of their own

`data/cemdata18-zeolites.json` is CEMDATA18 with 28 zeolites appended. It is
**generated** by `data/zeolites/regenerate.jl` and shipped, so it can be
reproduced and audited rather than trusted.

A pozzolanic or an alkali-activated binder at high alkalinity precipitates
zeolites. Without them in the species list the alkalis have nowhere to go but the
pore solution and the calculated pH comes out too high — an error in the phase
list that looks like an error in the solver. CEMDATA18 carries five zeolites;
clinoptilolite, heulandite, mordenite, phillipsite, analcime, stilbite and the
gismondine/faujasite/LTA series, in both their Na and their K forms, are not
among them.

The data are transcribed number by number from two open-access papers by the
laboratory that produced CEMDATA18 itself — Ma & Lothenbach, *Cement and
Concrete Research* **135** (2020) 106111 and **148** (2021) 106537, both DOIs
resolved against Crossref. Nothing is estimated, interpolated or adjusted.

The generator refuses on three grounds rather than warning: a symbol that would
overwrite a CEMDATA18 entry, a dissolution that does not balance in elements and
charge when re-derived from the formula string, and a `log Ksp` that does not
close to within 0.05 log units when recomputed from `ΔfG⁰` through CEMDATA18's
own aqueous Gibbs energies. That last one is what makes the merge defensible at
all: two thermodynamic datasets may only be merged if they share a reference
state, and the usual failure is silent — an offset of a few kJ/mol on `Na+` moves
every dissolution equilibrium by an order of magnitude with no solver
complaining. All 28 phases agree to within 0.026. The same three checks are
asserted in the test suite, because a generator can only refuse at the moment it
runs.

Three other candidate datasets were examined and rejected on measurement rather
than on preference; `data/zeolites/README.md` records which and why.

### Fixed — the documentation build, and then un-fixed the fix

It had reached 3 h 20 and was being canceled by its own timeout. The cause was
not diffuse: four pages ran coupled hydration trajectories, and one coupled
forward solve cost 364 s while everything else on the site together cost about
ten minutes.

The first answer was to compute them **once** and store the result, which worked
and brought the build to 23 minutes. It also brought a staleness guard comparing
the commit and the resolved solver version, a documented refresh procedure, and a
list of ways it could silently go wrong.

The second answer removed all of that. `OptimaSolver` 0.5.5 takes a warm cement
equilibrium from 583 ms to 17 ms, so the trajectories are affordable at build
time again and **every number in the manual is computed by the build that shows
it**. A stored result is a claim about code that may since have changed; nothing
here makes that claim any more.

`scripts/precomputed.jl` holds those runs and memoizes them per process — which
matters, because Documenter runs the whole site in one process and a page asks
for a trajectory's phase history and its calorimetry as two tables. Two calls,
one integration. The shared process is a hazard everywhere else in this file and
here it is what makes the arrangement work.

Two further changes came out of the diagnosis, which took four wrong turns before
it took the right one:

- the build **reports which block is slow**. Under `JULIA_DEBUG=Documenter` it
  named each block and timed none of them, so a three-hour build identified
  nothing. A logger now reports the duration of each block that exceeds a
  threshold, with a running total.
- `equilibrate_certified` **solves a starting point only when the search asks for
  it**. It offered every registered back end's answer as a start and returned at
  the first that certified, so the later ones were computed and thrown away. The
  honest measurement is in the docstring: this buys nothing where the first start
  does not certify, which includes the cement case.

A third change closes the arithmetic. Eight 28-day trajectories stand behind the
site, they share nothing at all — separate systems, separate solver buffers,
separate ODE states — and Documenter expands pages one at a time, so they ran
strictly in series. `warm_precomputed` now computes each page's own trajectories
**together**, on whatever threads the session was started with; the workflow asks
for four. Measured on the two runs behind the calibration target, on two threads:
1093.4 s one after the other against **593.2 s together, a factor 1.84** — and
both vectors came back identical bit for bit, maximum deviation exactly 0. The
only globals the solve path mutates are a warning gate and a warning counter,
neither of which enters the numerics. With one thread the call is an ordinary
`map`, so a default session behaves as before.

Two smaller economies, both of them fixes rather than trims:

- the alite-sensitivity table no longer integrates its own reference curve. The
  δ = 0 row **is** the calibrated fit, so it is read from the fit the page
  already plots — which also removes a real inconsistency, that run having
  silently omitted the fitted induction period and so measured the sensitivity of
  a different model from the one it reported.
- three tables whose *point* is a configuration that does not certify now pass
  `autostart = false`. The multi-start cascade — every back end, then the ideal
  pre-solve, then the homotopy continuation — cannot help where no admissible
  single-composition minimum exists; it spends minutes arriving at the answer the
  first route already gave. The certificate reported is the same one.

### Fixed — a one-character bibliography entry killed a three-hour build, and the build now checks itself first

`ExpandBibliography` is one of the **last** stages of `makedocs`: it runs after
every executed block of the site. So an entry DocumenterCitations cannot parse
does not fail the build in seconds — it fails it once every hour of computation
has already been spent, and nothing is written.

That is what happened. A title carrying `CNASH\_ss`, the LaTeX escape for an
underscore, threw `ArgumentError: Invalid command: \_ss` from the TeX parser and
took a three-hour run down at the very end.

The entry is fixed — the underscore is written bare inside braces — but the entry
was never the real problem. **`docs/make.jl` now formats every bibliography entry
before `makedocs` is called**, with exactly the call the late stage makes, so this
class of failure takes milliseconds and names the offending key. It immediately
earned itself: a second entry, `SmilauerKrejci2009`, threw
`ArgumentError: Premature end of tex string` because the `:authoryear` style
abbreviates a given name by slicing the raw string and cut inside `{\v s}`. That
one would have killed the next build. Its names are now written in UTF-8, which
survives any slice — the names themselves unchanged.

`CrossReferences` is a late stage for the same reason, so the same treatment is
applied to it: every `@ref` written as an explicit anchor is now resolved against
the `(@id ...)` anchors and the header slugs **before** `makedocs`, naming the
file and the target. The bare form is checked too — `[Some Heading](@ref)`
resolves against the *heading text*, slugified and case-sensitively, so a heading
renamed or merely recapitalized silently breaks every link to it, and three were
broken that way.

**And then a draft pass, which is the one that closes the class.** The static
checks above catch what can be caught by reading the markdown. `missing_docs`
cannot be: it needs the module loaded and the `@autodocs` filters applied, and it
is decided in `CheckDocument`, which runs *after* `ExpandTemplates`. A build died
there at minute 70 over five undocumented constants, having executed every
example on the site to reach the check and then terminating **before rendering**,
so the seventy minutes bought nothing. A **draft** build runs the same pipeline
with the `@example` blocks skipped and reaches the same checks in **24.2 s**,
measured, so `docs/make.jl` now runs one first — into a temporary directory, with
its own `CitationBibliography` since the plugin carries state across a build. It
is Documenter's own check, run early, rather than a second implementation of it.

Two things had to be turned off in that pass, and both for the same reason —
draft mode breaks them by construction, and each was found by running it rather
than by reasoning about it. `cross_references`: the figures on this site are
written by the `@example` blocks themselves, so with the blocks skipped every
`![](...)` on fourteen pages is an invalid local link. `size_threshold`: an
HTML-renderer limit that the real build, rendered by DocumenterVitepress, does
not have. `missing_docs` — the one the pass exists for — stays strict.

### Documentation

- **`theory/redox.md`** — why charge is a conservation law independent of the
  elements and when the rank test keeps it, why the electron activity must be
  inferred rather than read, and what a slag cement does and does not get from
  the calculation.
- **`manual/cement_notation.md`** — the oxide alphabet as a table, checked at
  build time against `CEMENT_TO_MENDELEEV` (formulas *and* molar masses), the bar
  convention, and the trap that `Species("C3S")` and `CemSpecies("C3S")` both
  print `C₃S` and differ threefold in molar mass.
- **`examples/example_stoich_matrix.md`** — a page titled "Stoichiometric
  Matrix" that displayed no stoichiometric matrix now shows three, and restores
  the parameterized decomposition of Chen & Brouwers: a C-S-H written
  `C_a S A_b H_g`, inverted symbolically over the anhydrous oxides, then
  collapsed to numbers and differentiated. One error in the old script is
  deliberately not carried over — it wrote ettringite without its alumina, which
  parses, weighs 1153 g/mol and is not a cement phase.
- `manual/chemical_system_state.md` — the volume and porosity example reported
  `0.0 m³` and `NaN`, because species built from formulas carry no molar volume.
  It now uses database species, and the rescaling of a state is documented.
- **`manual/binder_families.md`** — the map of EN 197-1: the families and their
  composition ranges as a table, what each constituent brings to the element
  budget and therefore which of this package's models the calculation needs, and
  the measured heats of the seven shipped records in the order of their
  replacement level, 376 J/g down to 234 J/g.
- **`examples/miscibility_gap.md`** — the same CEM I run three ways on the AFm
  sulfate/hydroxide binary: ideal mixing (certified, to a question that was
  changed), the published parameters with one composition (refused at
  construction, and uncertifiable when the refusal is waived), and the published
  parameters with two instances — which the page shows returning the **same**
  composition twice, then computes the pair the minimization did not find with
  `common_tangent` and applies the lever rule to it. With the mixing-energy curve
  and its spinodal drawn, and an explicit statement of what `optimal` does and
  does not prove once the problem is no longer convex.
- **One executed page per blended family**, each an element budget in and a
  certified assemblage out, beside the measured calorimetry of a real specimen
  where one exists:
  - `examples/cem2_blended.md` — a CEM II/A-LL against a CEM II/B-S, with the
    limestone removed from the first to isolate the carbonate effect from the
    dilution. The carbonate takes the AFm site, the sulfate stays in ettringite,
    and the mechanism comes out of the element budget with nothing fitted.
  - `examples/cem3_slag.md` — the glass entry route, hydrotalcite, and the
    sulfur ladder that makes charge a component of its own.
  - `examples/cem4_pozzolanic.md` — the same paste solved with `CSHQ` and with
    `CNASH_ss`, which is where the aluminum question becomes visible, and
    portlandite along a replacement sweep run **twice**: at the reacted fraction
    a specimen has at 28 days, and in the limit. The two answer different
    questions and are routinely confused — in the limit a CEM IV/A exhausts its
    portlandite inside the EN 197-1 range, at 28 days the same binder keeps most
    of it. The sweep advances by continuation, each point starting from its
    neighbor's certified answer, and it carries the **certificate into the
    figure**: a point the certificate refused is drawn hollow rather than
    dropped, because a smooth curve through a hole in the evidence is worse than
    the hole. **This is the one page with no measured specimen behind it** — the
    deposit carries no CEM IV record — and it says so at its head.
  - `examples/cem5_composite.md` — slag and fly ash at once, on one additive
    budget, with every element traced to the constituent that brought it, and a
    sweep over the round robin's own 7-, 28- and 90-day columns so that the
    sensitivity to the reacted fraction is measured on published ground rather
    than on an abstract parameter. Like the CEM IV sweep it advances by
    continuation in the reacted fraction — a younger paste has released less of
    everything, so it is a smaller perturbation of pure water and a good start
    for an older one. On a convex problem that cannot change what is found, only
    whether it is found, which on a 135-species cement is the whole difficulty.

  Every composition not in the deposit is labeled `ASSUMED` at the point of use,
  at the midpoint of the EN 197-1 range for its designation. The deposit reports
  fineness, water/binder ratio and calorimetry, and reports neither the clinker
  phase composition nor the replacement level of any blend. Each page also states
  a **reacted fraction** per constituent — the water ceiling for the clinker, the
  measured 28-day value for a glass — because handing the whole binder to a Gibbs
  minimization asks a question about geological time.

- **`theory/kinetics.md`** — every rate law the package ships, the physics each
  one encodes, its parameters with their units, and a table saying which numbers
  are published and which were fitted here. The saturation-ratio form and its
  catalysts; the canonical Parrott–Killoh with its three competing mechanisms,
  and why the attribution of the older variant was withdrawn rather than
  repaired; the Waller sigmoid for supplementary materials; the three
  multiplicative corrections, including the genuine discontinuity at 80 % RH; and
  `α_max`, which is the one parameter that is not about speed. It also sets the
  package's own Waller kinetics against the RILEM round robin's direct
  measurements and reports that they disagree, in opposite directions, rather
  than quoting whichever is convenient.
- **`manual/choosing_species.md` gains a seventh trap** — a solid solution whose
  declared range cannot reach where the answer is. CEMDATA18's siliceous
  hydrogarnet binary spans x(Al) ∈ [0, 0.5] only, which no CEM I ever needs and
  every aluminum-rich blend does; the minimization stops at the edge of the range
  with aluminum left over and says nothing. Including why the tempting fix —
  declaring all three members as an ideal ternary — is worse, not better.

### Fixed — the API section had vanished from the navigation bar

A navbar curation added in August folded `API` and `References` into a dropdown
labeled `Reference`. On v0.17.0 that read, correctly, as the docstring reference
having been removed from the manual. `API` is a top-level entry again, and the
mechanism carries a comment saying why it is not what gets folded.


## v0.17.0 — the water that is there, the water that counts, and every solid solution declared

Three halves, which is one too many for the metaphor and an honest count of
the release. An **ion-interaction activity model**, the second of the two
ingredients the roadmap named for predicting where a sealed paste stops. The
**coupled kinetic run** that 0.16.0 described and never performed, which turns
out to have been impossible for a reason worth recording. And a cement whose
**phase list is no longer chosen by hand**: declaring every solid solution the
database defines and letting the minimization decide was not usable before this
release, and now is.

Around them, the documentation was reorganized into four chapters, and every
activity model in the package states the physics behind its formulas and the
provenance of every default it carries.

### Breaking changes

- **`ΔₐG⁰overT` is renamed `ΔₐG⁰overRT`.** It was a typo: the quantity is
  `ΔₐG⁰/(RT)`, dimensionless, and it has to be, since it is added directly to
  `ln aᵢ` to form `μᵢ/RT`. The parameter tuple is part of the documented
  interface — an `activity_model` closure reads `p.ΔₐG⁰overRT` — so a custom
  activity model written against the old name now fails with a
  `NamedTuple has no field` error rather than silently. 45 occurrences renamed.
- **The registry treats a minor bump below 1.0 as breaking whatever the API
  did**, so a downstream bound pinned to the previous minor will not accept
  `0.17` and must be widened. `MeanFieldHomogenization.jl` depends on this
  package only in `docs/Project.toml`, which already reads
  `"0.14, 0.15, 0.16, 0.17"` and needs no change.
- **`PoreHumidity` no longer differentiates its retention law at full
  saturation.** Its value there is unchanged; its derivative is now zero instead
  of infinite. Anything that read the derivative at `S = 1` — nothing could,
  since the integration it exists for did not run — changes.
- **`SolidSolutionPhase` refuses a model whose mixing energy is concave.** A
  phase with a spinodal unmixes: the Gibbs minimum there is two coexisting
  compositions, and this formulation has one amount per species to describe it
  with. Code that built such a phase and got an answer now raises at
  construction, with the interval named. `check_convexity = false` restores the
  old behavior, and the optimality certificate — whose sufficiency rests on
  convexity — is then explicitly given up. Nothing in
  `data/solid_solutions.toml` is affected: every shipped model is convex.
- **The composition a solve returns may differ** where it previously could not
  certify. An answer whose solvent had been taken by the solids is no longer
  ranked against ones that conserve mass, and a route that failed before may now
  certify, so a script that recorded an uncertified result will see a different
  one.

### Added — `PitzerActivityModel`

An ion-interaction model, and a different kind of object from everything else in
`activities.jl`. Those are corrected Debye-Hückel laws: one screening term, one
size correction, one empirical term for the rest. This one expands the excess
Gibbs energy as a **virial series in the molalities** — a coefficient per ion
pair, one per triplet, no per-species radius — which Anderson & Crerar derive as
a cluster expansion with osmotic pressure in place of pressure.

The consequence is the reason to have it: `γ` and the osmotic coefficient are
partial derivatives of **one** function, so the Gibbs-Duhem relation between
solutes and solvent is an identity of the algebra. Measured, with the derivative
taken analytically: the residual is **exactly zero** along three composition
directions at 0.1, 1 and 3 mol/kg, where the B-dot model's is not small. A
finite difference cannot see this — its own truncation error at 0.1 mol/kg is
larger than the quantity — which is why the test uses AD.

Also validated: the dilute limit reproduces `log₁₀ γ± → −A|z₊z₋|√I`, and the
model agrees with the independent B-dot implementation at a millimolal, as two
models sharing a limit must.

**Nothing is parameterized by default, and completeness is checked against the
species list rather than against the set.** A table is complete or not
*relative to a system*, so `PitzerParameters` takes every table as a keyword
without a default — `UndefKeywordError` before a number is computed — and the
check that matters happens when the model meets a `ChemicalSystem`: every
cation-anion pair present must have a `β⁰`, and the error names those that do
not. A missing `θ`, `ψ` or `λ` is zero, which is the convention of the
literature the tables come from; a missing `β⁰` cannot be, because falling back
on ideal behavior for one pair of a Pitzer calculation is not an answer.

Two limitations are in the docstring rather than left to be discovered: the
higher-order electrostatic terms are not implemented, which is exact for a
symmetrical pair and an omission in a Na/Ca mixture; and no temperature
dependence of the interaction parameters themselves.

### Added — a cited cement parameter set, with its provenance per entry

`data/pitzer-reardon1990.toml`, read through `build_pitzer_parameters`, from
Reardon (1990), *Cement and Concrete Research* **20**, 175–192, whose tables are
those of Harvie, Møller & Weare (1984) with three exceptions he names. Nothing
about it is automatic: the file has to be named.

Two properties of the set change what a user should conclude from it, so they
are recorded in the file and readable from code:

- the silicate, aluminate and ferrate parameters are **estimates, not
  measurements**. Reardon says which analog each borrows — HSO₄⁻ for Fe(OH)₄⁻,
  Al(OH)₄⁻ and H₃SiO₄⁻, SO₄²⁻ for H₂SiO₄²⁻, H₂CO₃⁰ for H₄SiO₄⁰ — and those are
  exactly the ions a cement assemblage needs. Every affected entry carries
  `origin = "estimated:<analog>"`, read back by `pitzer_origin`;
- the set assumes a **fully dissociated** speciation. The association of Ca with
  SO₄, of Na with OH, is already inside these β coefficients, so a species list
  that also carries `Ca(SO4)@` or `CaOH+` — as CEMDATA18 does — counts each
  association twice, and CEMDATA18 also names the silica species differently.
  The completeness check refuses such a system rather than returning a number,
  and that refusal is the correct outcome in any code, not a shortcoming here.

Transcription was verified rather than trusted: Table 2's 216 values were
checked against the PDF text layer and all matched, the only surplus tokens on
the PDF side being OCR fragments of the Greek symbols. Tables 5 and 6 could not
be checked that way — that text layer turns one −0.0677 into −0.0077 and drops
two others — so they were read from the pages rendered at six times
magnification, and the data file's header says so.

### Fixed — `PoreHumidity` could not be integrated

0.16.0 built it, wired the dispatch, and said in the self-desiccation page that
handing it to `parrott_killoh_avrami` and integrating was how to obtain the
arrest in time. Nothing ever did, and doing it found why: **the integration did
not advance past `t = 0`.**

Every retention law of the van Genuchten family has an unbounded `dp_c/dS` at
full saturation, and a sealed paste starts saturated. An implicit solver
differentiates the rate law at its first step, so it received an infinite
Jacobian entry — the gradient of the humidity with respect to the composition
came back with an infinite norm — and returned immediately with every degree of
hydration at zero. No unit test on the object could see this: the object is well
behaved, and only an integration reaches the singular point.

The value at saturation is not in doubt, so it is now returned as a constant and
the derivative there is zero. The code says plainly that this is a
**regularization and not an identity**, that it applies within `1e-10` of
saturation — the initial condition alone, in practice — and that the solver
leaves that point on its first successful step, after which the real and stiff
derivative applies.

What the coupled run then gives, with `α_max = 1.0` and no Powers bound: the
arrest is a **result**. The uncoupled rate reaches the same α(C3S) = 0.9164 at
every w/c, because it knows nothing about how much water there is; coupled, α
runs from 0.556 at w/c 0.25 to 0.837 at 0.50. The drier mixes stop at an
internal humidity of 0.800 and a pore saturation of **0.786** — which is the
`S*` the static water budget reads off the same measured retention curve at
RH 0.80. The trajectory arrives at the number the budget assumes, by a different
route, with nothing arranged to make it so.

It also shows that the proportionality `α_max ∝ w/c`, exact in the static
construction, does **not** survive integration: `k = (w/c)/α_max` is near 0.45
at the dry end and rises with w/c. The static page is careful to say that the
proportionality is structural and no evidence; this is what it looks like when
the structure is removed.

### Added — `water_activity`, and a `kelvin_shift` on `log_activities`

0.16.0's own release notes recorded that `log_activities` did not know about the
capillary shift, so a solve posed at `a_w = 0.90` read back as 0.999995. Both
accessors now take a `kelvin_shift` keyword, `0.0` by default so that nothing
changes for a caller who does not ask, and `water_activity(state, model;
kelvin_shift)` returns the composed value. The chemical and capillary lowerings
are independent, their chemical potentials add, and so their activities
multiply — which is what the keyword implements and what the docstring derives.

### Documentation — four chapters, and the physics behind every formula

The page tree was four flat lists in the order the pages were written. It is now
grouped, and each chapter answers one question: **Theory** why this is the right
calculation, **Manual** how an object is written, **Tutorials** how to drive a
calculation, **Applications** what a real case looks like and what the choices
cost in numbers. Nine pages that were object syntax filed under Tutorials moved
to the Manual. **Cementitious media is a subsection of four chapters rather than
a chapter of its own**, because the material is the subject of the package and
belongs wherever its question is being asked. No page was removed; three paths
that released notes link by URL were deliberately left where they are.

A **Theory** chapter, with the rule written on its own index that a theory page
may show code where the correspondence with the implementation is the point, and
does not build systems, solve, or print tables:

- **Thermochemistry** — the definitions and identities in the code's own
  notation, where `μ°(T,P)` comes from, equilibrium as a constrained
  minimization with the component potentials as its multipliers, and the
  saturation index as a difference of those potentials. That identity is now
  asserted in the test suite, which recomputes `LogSI` by hand from the formula
  the page states.
- **Proving that an answer is the answer** — the convexity proof, the KKT
  conditions and the dual solver, moved out of the tutorial they were living in.
- **Activity models** — where the `√I` comes from, in three steps with the
  assumptions each makes, since those assumptions are what later fails; why `A`
  and `B` are properties of water rather than fitting constants; the Debye
  length, which at a cement pore solution's ionic strength is 0.55 nm — the
  width of the gel pores where a paste keeps its last water, so the
  continuum-dielectric and mean-field assumptions are strained exactly where the
  arrest happens; and `Ḃ` in its real status as a deviation function fitted to
  one salt.
- **Solid solutions** — mixing entropy, the excess Gibbs energy, the Margules
  and Redlich-Kister forms as implemented, and what the sign of `W` means. Two
  identities are checked rather than claimed: that Redlich-Kister reduces to a
  regular solution, and that `W = 2RT` is the critical point.
- **The water budget of a hydrating paste** — why a Gibbs minimization predicts
  that clinker survives below `w/c ≈ 0.30` while Powers reports 0.42, and why
  the gap is not a thermodynamic quantity: the capillary term is two orders of
  magnitude too weak, and what stops a real paste is a broken liquid path
  across four orders of magnitude of length, which a 0D model has no
  representation of.

Every default of every activity model is now classified as derived, tabulated,
or **a convention with no source recorded in this package** — `Ḃ = 0.041`,
`Kₙ = 0.1` and `å_default = 3.72` are in the third class, and that is said out
loud rather than given a plausible citation.

Long solver output is folded into collapsed blocks rather than unrolled: seven
pages ended on an assignment whose value Documenter then displayed in full, up
to 48 species and a whole conservation matrix between two paragraphs.

### Validated — against measurement, not only against itself

Every earlier claim about the Pitzer model was internal consistency, and
internal consistency cannot tell a correct parameter set from a self-consistent
wrong one. Hamer & Wu (1972), Table 16 — a critical compilation of the osmotic
*and* mean activity coefficients of NaCl at 25 °C — settles it, and each half of
the model is checked separately:

**Pitzer follows the measurement to better than half a percent from 0.001 to
6 mol/kg**, through the minimum near 1 mol/kg and the climb back to 0.99 at six
molal, neither of which a Debye-Hückel form can produce. The osmotic
coefficient agrees to the same order. On the same points the B-dot model is 5 %
out at a tenth molal, 19 % at one and 44 % at six — its stated range is real,
and nothing in its output announces the exit.

That is also a check on the transcription: the Na/Cl coefficients were read off
a scanned table, and nothing mistyped reproduces a measured curve over four
decades.

### Added — a worked cement, from the clinker up

`examples/cem1_from_clinker.md`: four anhydrous phases, a w/c, and everything
after that computed — the degree of hydration of each phase, the hydrates that
appear, the porosity, the chemical shrinkage, the internal humidity. Three
clinkers, of which one is the measured CEM I of Baroghel-Bouny et al. and two
are constructed to trade alite for belite at constant silicate, which isolates
one variable rather than comparing three cements nobody has made.

Writing it corrected three statements that reading could not have caught, and
one of them is a trap worth knowing: the aluminate reaction
`C3A + 3 Gp + 26 H2O → ettringite`, driven by a Parrott-Killoh rate, **violates
mass conservation**. That rate follows its own clinker phase and does not watch
its co-reactants, so the extent keeps advancing after the gypsum runs out —
demanding 0.28164 mol against 0.25499 present, with the gypsum floored at zero
rather than going negative, so sulfate is created. Nothing in the package
objects: `extent_residual` measures integrator drift, and the feasibility
machinery guards an equilibrium sub-solve that is not running. A
fixed-stoichiometry kinetic reaction is only safe when its co-reactants cannot
run out, and the page says so with those numbers.

Six figures there, and four more in the Applications pages: the three activity
models against molality with the limiting law, their water activities, the
Gibbs-Duhem residual on log-log axes, and the mixing free energy of a regular
solution across the critical point, which makes the miscibility gap visible
rather than tabulated.

### Fixed — a page could change the plot font for every page after it

`examples/cem1_solid_solutions.md` called `default(fontfamily = "Computer
Modern")`. Documenter runs every `@example` block in one process, so that
applied to every page built afterwards — and Computer Modern has no Unicode
sub- or superscripts, which the axis labels here use freely (`ΔₐG⁰ [J.mol⁻¹]`,
`Ca²⁺`, `CO₃²⁻`). GR then prints `GKS: glyph missing from current font` for each
one and looks for a fallback: instant on a developer machine, where fontconfig
is warm, slower in a CI container, and in either case hundreds of log lines that
bury what the build is actually doing.

The labels keep their typography. What changed is that the font is now chosen
**once**, in `docs/make.jl`, as one that carries those glyphs, and `make.jl`
refuses to build if any page sets `fontfamily` itself — naming the file and
line, in two seconds rather than after an hour. A page may still set
`framestyle`, `grid` and the rest.

The workflow also sets `JULIA_DEBUG=Documenter`, which makes Documenter name
each page as it expands it. Without it the stage that executes all 299
`@example` blocks prints nothing at all, so a build working normally is
indistinguishable in the log from one that has died.

### Coverage

`src/equilibrium/pitzer.jl` and `src/databases/pitzer_toml.jl` are covered
completely. The gaps were not scattered lines but three untested behaviors — a
neutral solute reaching the λ terms, a gas phase, and solid-solution
end-members, the last of which is the silent failure where an aqueous model that
forgets the solid-solution branch leaves them as pure phases.

### Fixed — the search could settle outside the model's domain

`_keep_better` ranked candidates on the KKT error alone. A composition in which
the solids have taken the water — solvent mole fraction 0.033, against the
`SOLVENT_FRACTION_FLOOR` of 0.5 — is not an answer this model can describe, and
letting it win on a smaller residual hid every candidate that conserved mass
behind it. Admissibility is now compared first, in both directions, as the
optimality flag already was. The package diagnosed that state already; what it
could not do was stop one from being chosen.

### Fixed — a lost solid solution could not be recovered

`_repair_start` skipped solid-solution end-members, on the correct observation
that the saturation index of a member at the bound reports a small mole fraction
rather than a phase that should form. The phase has its own criterion, and it is
now used: for ideal mixing a solid solution exists only at `xᵢ = 10^SIᵢ` with
`SIᵢ` the index of the **pure** end-member, so it is saturated exactly when
`Ω = Σᵢ 10^SIᵢ = 1` and should form above it. The reported index gives `SIᵢ` back
once `ln aᵢ / ln 10` is added, which also removes the `0/0` a phase sitting
entirely at the floor would produce.

### Added — the ideal model as a stepping stone

When no route certifies, the same minimization is solved first under
`DiluteSolutionModel` — no activity coefficients, well conditioned, and it
certifies — and its answer becomes the start for the non-ideal solve, which then
begins with the correct active set instead of discovering it. Only when nothing
else certified, so the ordinary case pays nothing.

### Added — `spinodal_interval`, and a refusal that names the interval

A solid solution exists as one homogeneous phase only where its mixing energy is
convex. Where `d²g/dx² < 0` the Gibbs minimum is **two coexisting compositions**,
and a formulation with one amount per species cannot hold them.
`SolidSolutionPhase` now refuses such a model at construction and says where the
interval is; `check_convexity = false` proceeds anyway, with the certificate's
sufficiency — which rests on convexity — explicitly given up.

The classical symmetric threshold is recovered as a check: a regular solution
unmixes above `W = 2RT`. The parameters this matters for are real ones: the AFm
and AFt Redlich-Kister sets used for a CEM II are concave over `x ∈ [0.63, 0.91]`
and `[0.56, 0.83]`, which is why GEM-Selektor declares each of those binaries
twice — ten phases where the distinct chemistry is eight. Run with them, the
solve stopped at an element balance of 2.6e-3 with nothing reported missing;
with ideal mixing the same eight phases certify to 4.8e-13 and reproduce the
reference to 0.005 units of pH.

Nothing in `data/solid_solutions.toml` is refused: every shipped model is convex.

### Compatibility

Tested on Julia 1.12 and 1.13. `[compat] julia = "1.12"` is unchanged and
already admitted 1.13.


## v0.16.0 — the water a paste cannot use

A sealed cement paste stops hydrating before it runs out of cement, and until now
nothing in the package represented why. `powers_alpha_max` supplied the empirical
cap and said in its docstring that the bound is "a statement about transport and
access, not about thermodynamics". This release supplies the missing state
variable — the internal relative humidity of a partially saturated pore space —
and the tutorial that validates it recovers Powers' coefficient from independent
data.

### Breaking changes

Nothing in the API breaks and no default behavior changes. Two consequences are
breaking in practice:

- **The registry treats a minor bump below 1.0 as breaking whatever the API did**,
  so `[compat] ChemistryLab = "0.15"` will not accept `0.16`. Downstream packages
  must widen. `MeanFieldHomogenization.jl` depends on this package only in
  `docs/Project.toml` and needs `"0.14, 0.15, 0.16"`.
- **`solve_certified` now writes the *winning* candidate's parameters** into the
  caller's `parameters` Ref. It wrote the last candidate's while returning the
  best-by-error one, so a prescribed-pH scan reported a titrant amount belonging
  to a different composition. Code that relied on the old value was reading a
  mismatch.

### Added — water retention, and the Kelvin relation

`WaterRetention` and its subtypes describe how tightly a pore space holds the
water still in it. That is a constitutive input, measured on a particular material
at a particular age, so **nothing ships with a default value**: every named law
takes its parameters as keywords without defaults, and Julia raises
`UndefKeywordError` before any number is computed. The enforcement is the
language's, not a check that can be forgotten.

- `TabulatedRetention(; S, a_w)` takes a measured isotherm and interpolates
  linearly in `ln a_w`, the variable the coupling reads. Its inner constructor
  refuses a non-monotone table rather than interpolating nonsense, and leaving the
  measured range clamps and warns once instead of extrapolating a logarithm into
  confident absurdity.
- `VanGenuchten(; a, m)`, with the published parameters of Baroghel-Bouny et al.
  (1999) quoted in its docstring together with the convention they use — they
  write the same expression with `b = 1/m`, and getting that inversion wrong is
  silent.
- `FunctionRetention(f)`, the escape hatch: a bare function of the saturation,
  returning the activity directly. It is the one law nothing validates, since
  `f` is opaque, and its docstring says so — `CapillaryWater` still tests it at
  `S = 1`. `CapillaryWater` and `PoreHumidity` wrap a bare function in it
  themselves, so no caller has to name it.
- `kelvin_activity`, `kelvin_radius`, `capillary_pressure`, `water_activity`.
  RH 80 % is a meniscus of radius 4.8 nm, the gel-pore scale, which is why the
  water Powers assigns to gel pores and the water a sealed paste cannot use are
  the same water.

### Added — `CapillaryWater`, a constraint with a certificate

One unknown, the Kelvin shift of the solvent's chemical potential; one equation,
the retention law at the current saturation `S = V_liquid/(V_ref − V_solid)`. No
titrant column, because a sealed specimen exchanges water with nothing — the whole
difference from `FixedActivity`. It refuses eagerly when a species present carries
no standard molar volume, since such a species contributes zero to the pore volume
in silence and would corrupt the saturation with nothing to show for it.

Measured on a CEM I paste: certified at every w/c from 0.30 to 0.50, constraint
residual between 1e-16 and 1e-13, internal humidity falling from 0.87 to 0.32 as
the mix gets drier.

Two caveats are in the docstring rather than left to be discovered. The
certificate proves a **KKT point of the constrained problem** — stationarity, mass
balance, no absent phase supersaturated, and the capillary closure satisfied — but
not global optimality, because a composition-dependent shift of `ln a_w` is not
derived from a convex `G` for an arbitrary retention law. And **`log_activities`
does not know about the shift**: it evaluates the activity model, so with the
shift at `ln 0.90` it returns `a_w = 0.999995` while the solve was posed at 0.90.

### Added — `PoreHumidity`, and where the arrest actually comes from

The Kelvin term is **two orders of magnitude too weak to arrest hydration**.
Measured, with a certificate on every answer: imposing a water activity from
saturation down to 0.80 leaves the equilibrium assemblage of a CEM I paste
unchanged to six digits. Below that the constrained problem stops certifying, so
nothing is claimed there — and the arithmetic says nothing is expected there
either. At `a_w = 0.80` the shift is
553 J per mole of water, worth about 1.8 kJ per mole of alite against a hydration
Gibbs energy of order −100 kJ/mol; nulling it would need `a_w ≈ 5e-6`, a Kelvin
radius smaller than a water molecule. A real paste stops at 75–80 % RH because
transport and nucleation stop.

So the humidity belongs in the rate law, and `humidity_factor` has implemented
that cut all along. What was missing is that its argument was exogenous:
`_humidity_at` took a constant or a function of **time**. The rate closure already
receives the full composition, so `PoreHumidity(retention, system; reference)`
computes the saturation and returns the pore humidity, and a third method reads
it. The other two ignore the new argument, so nothing a caller wrote before
changes.

### Fixed — the certificate now audits the problem that was solved

`optimality_certificate` rebuilt a `FixedTP` problem whatever constraint had been
applied. That misses two things at once: the conservation rows lose their `Aq q`
term, so a prescribed activity or pH is measured against a budget short by exactly
the titrant amount; and `hq` is skipped, so a constraint that shifts a chemical
potential is measured against the unshifted one and can never certify however
right it is. It takes `constraint` and `q` keywords, both defaulting to the old
behavior, and returns `param_residual`.

`_kkt_error` counts that residual, so the multi-start search can no longer prefer
a candidate that minimizes the Gibbs energy while violating the very equation that
makes it a constrained answer.

### Compatibility

Tested on Julia 1.13, which the CI matrix now runs explicitly alongside the
1.12 floor. `[compat] julia = "1.12"` already admitted it and did not change.

### Documentation — Powers' 0.42, taken apart

A new tutorial, [Self-desiccation: where Powers' 0.42 comes from](https://micropochemomechanics.github.io/ChemistryLab.jl/stable/tutorials/self_desiccation/),
with `scripts/self_desiccation_powers.jl` as its runnable companion.

It derives `α_max = (w/c)/k` with `k = b + s S*/(1−S*)`, then fills each term from
its own source: `b = 0.3095` g/g and `s = 0.0639` cm³/g from certified equilibria
at imposed degree of hydration (both constant in α to four digits, which is the
check on the linear budget the derivation assumes), and `S* = 0.7862` at RH 0.80
from a measured desorption isotherm.

**Inverted, Powers' 0.42 corresponds to an arrest at 77.5 % relative humidity** —
assembled from a chemical shrinkage computed out of the thermodynamic database, an
isotherm fitted for a drying study two decades earlier, and Powers' own water
split. That is the window sealed pastes are independently reported to stop in.

The page is explicit that the *proportionality* `α_max ∝ w/c` is structural and
therefore no evidence, that only the coefficient is predicted, and that the
residual 11 % is attributable: the model's formula water exceeds Powers'
non-evaporable water by 0.0795 g/g, which is the interlayer water CEMDATA18 writes
into the C-S-H formula and D-drying removes. Not the same quantity, so not an
error — and the `CSHQ` solid solution binds more, 0.3684 g/g, so it moves away
from Powers rather than toward him.

It also **executes the two controls that could have refuted it**, rather than
asserting them: the imposed-activity sweep above, which is why the arrest is read
from the rate law and not from the Gibbs energy; and the same construction run on
a second measured curve from the same table, which returns a ratio
`α_max/(w/c)` constant to every digit printed in both cases and differing between
them by 18 %. The constancy is an identity of the construction and carries no
information; the value, `1/k`, carries all of it.


## v0.15.2 — a solution that is no longer a solution

A certificate proves that a composition minimizes the Gibbs energy of the problem
**as posed**. It says nothing about whether the problem was posed inside the
model's domain, and there is one way to leave that domain while producing an
ordinary-looking `ChemicalState`: let the solids take all the water.

### Added — `solvent_fraction`, and a guard on the answer

`solvent_fraction(state)` is the mole fraction of the aqueous solvent within the
aqueous phase, `n_w / Σ_aqueous n`. Every quantity built on that phase —
molality, ionic strength, activity, pH — is defined *per kilogram of solvent*,
and `DualEquilibriumSolver` parameterizes its interior variables by the solvent's
chemical potential. All of it presumes the solvent **is** the phase. A real
electrolyte keeps `x_w` above about 0.9: seawater is 0.99, a saturated NaCl brine
0.90.

Measured on a sealed cement paste taken below its stoichiometric water demand —
w/c = 0.28 on the mix of the w/c example — the Gibbs minimum consumes the free
water down to **6e-9 mol**, the solver's floor. The solvent then holds a **fifth**
of its own aqueous phase and the ionic strength is reported as **409 mol/kg** by
a Debye-Huckel model valid to about one. Nothing in the certificate objects,
because nothing is wrong with the minimization: forming more hydrate always
lowers the energy, and nothing in the model penalizes a solution concentrated
past any physical meaning. The water activity of a real paste collapses as the
pores empty and stops the reaction; an activity model extrapolated that far goes
on returning finite numbers.

`equilibrate_certified` now checks `solvent_fraction` on the answer it returns
and says so, naming what happened, when it falls below `SOLVENT_FRACTION_FLOOR`.
Warned rather than raised under the default flag, because the composition of the
*solids* still carries the mass-balance information a caller may legitimately
want; under `STRICT_CONVERGENCE[]` it raises, like any other answer that is not
one.

`SOLVENT_FRACTION_FLOOR = 0.5` is not a modeling choice but a floor no real
solution comes near — about 28 mol of solute per kilogram of water, past
saturation for anything. A value below it means the solve drove the water into
the solids, not that the solution is concentrated.

### Documentation — a usable answer below the stoichiometric demand

A high-performance concrete is mixed at w/c between 0.25 and 0.35, which is an
ordinary regime and must be computable, with a certificate and with all three of
the things one wants from it: how much clinker stays unhydrated, which hydrates
form, and what the pore solution contains. The w/c example now works that case.

The construction is to stop the reaction where the physics stops it rather than
ask the minimizer to discover an arrest point it has no term for: react a fraction
α of the clinker with **all** the water, and leave the rest unhydrated. The
equilibrium is then computed on a system that still has a solution in it, and α
is what `powers_alpha_max` supplies. Imposing the reacted fraction is the standard
construction of cement thermodynamic modeling — it is how Lothenbach & Winnefeld
(2006) compute a hydrating paste.

Measured: at w/c = 0.25 with α = 0.595 the solve certifies, the solvent holds
0.9993 of its phase, the ionic strength is 0.0351 mol/kg and the pH 12.39, with
40.5 % of the clinker unhydrated and a total porosity of 0.2056 — against the
409 mol/kg and vanished solvent of the unconstrained solve at the same w/c.

Predicting the arrest point instead of imposing it is a well-posed thermodynamic
question the package cannot answer yet, and the page says which two ingredients
are missing. An **activity model valid at very high concentration**, Pitzer-class,
because what physically stops hydration is the collapse of the water activity as
the last of the pore solution is consumed, and an extended Debye-Huckel model
extrapolated to 409 mol/kg goes on returning finite numbers instead of collapsing.
And a **coupling between pore structure and water activity**, the Kelvin term,
because water in a fine pore is held at a reduced activity whatever its
composition — which is what self-desiccation is, and it is poromechanics rather
than solution chemistry.

### Documentation — the water-limited regime, and Powers' 0.42

`docs/src/examples/cement_wc_ratio.md` claimed that "no clinker survives at any
w/c", explained by the minimum "always forming a less hydrous assemblage rather
than leave alite standing". The first half is true only over the 0.30-0.60 range
it scans, and the second is false: the least hydrous assemblage available still
binds water, and when there is not enough, alite stays.

Measured on the page's own species list, clinker left as a fraction of the
cement: 54.6 % at w/c = 0.15, 34.3 % at 0.20, 14.1 % at 0.25, 2.0 % at 0.28, and
zero from 0.30 up. The crossover is a property of the hydrate assemblage, not a
constant, so the page measures it in an executed block instead of quoting one —
and reads it for one thing only, the *existence* of the regime, since past it the
aqueous phase has vanished and the amounts are not equilibrium values.

That does not put the page at odds with Powers, and the reason is worth stating
because the two numbers measure different things. Powers' 0.42 g of water per
gram of cement is **not** a stoichiometric demand: it is about 0.23 g of
non-evaporable water, written into the hydrate formulae and a mass balance any
Gibbs minimization must respect, plus about 0.19 g of **gel water** held in the
C-S-H gel pores, physically present and chemically unavailable. A sealed paste
stops by self-desiccation with water still in the specimen. So between roughly
0.23 and 0.42 the equilibrium and Powers disagree and **both are right**, and
below the stoichiometric demand they agree. `powers_alpha_max`'s docstring said
"hydrating one gram of cement binds about 0.42 g of water", which is the
misreading that produced the wrong claim; it now separates the two halves and
notes that curing water moves the bound to about 0.36.

The page also solves with `equilibrate_certified` rather than through Ipopt, and
prints whether every point certified. Ipopt is kept as a documented alternative
in an admonition, not executed: not a worse optimizer, but a bare interior point
that returns an iterate and no statement about it, and an extra binary dependency
for a calculation the package can already prove. The visible difference is small
and telling — absent phases come back at exactly zero instead of sitting at the
solver's 1e-8 lower bound. Every number on the page was re-measured on the
certified route and is unchanged: pH 12.3924 flat across the scan, porosity
12.2 % to 41.0 %, and at w/c = 0.50 the one-argument porosity 0.2846 against a
two-argument total of 0.3376 with the volume shrinking 7.41 %.


## v0.15.1 — a starting point that conserves mass, and a strict flag that is read

The automatic initial approximation shipped in 0.15.0 could hand the solver a
starting point that is not on the constraint surface, and then build every later
step on it. Reported from a Windows run of the same cement paste that certifies
here: `optimal = false`, an element balance off by **6.7 mol**, every hydrate at
zero — and a table of amounts that reads like a result.

### Added — the certificate names the missing phase, and the route puts it in

`saturation_indices(state, model)` returns `LogSI` for every species: zero for a
phase at equilibrium with the solution, negative for an undersaturated one,
positive for one that **should have precipitated**. It is what GEM-Selektor
prints as `LogSI`, and what `optimality_certificate` had been compressing into a
single worst violation without saying which phase it was.

No fitting is involved: the row labels of the conservation matrix are the primary
species, so a component's element potential is that primary's chemical potential
and `LogSI_s = [Σ_c A_cs μ_c/RT − μ_s/RT] / ln 10`. The check comes with it —
every phase actually present at an equilibrium must come out at zero, and on a
CEM I paste the twelve present solids land within 1.2e-12.

`equilibrate_certified` now **acts** on that. When its answer is uncertified
because a phase sits at the lower bound while the solution is supersaturated with
respect to it, the route puts that phase in and solves again, at most
`_MAX_RESTARTS` times.

The failure this exists for is a **phase swap**, which an active-set loop that
admits one phase at a time cannot perform. Reported from one machine while
another certified the same source: all 0.02515 mol of magnesium sat in brucite
with `hydrotalcite` absent and supersaturated by 5.58 log units, and admitting
the hydrotalcite requires dissolving the brucite entirely and taking aluminum
back from the hydrogarnet in the same step. The solve was otherwise impeccable —
stationarity 1.5e-16, element balance 1.8e-14 — which is exactly what an answer
converged onto the wrong active set looks like.

It is reproduced exactly, and now fixed: solving the paste with `hydrotalcite`
out of the phase list gives a certified 74.1514 cm3 with all the magnesium in
brucite, which is the reported state to seven digits; putting hydrotalcite back
and starting from there,`equilibrate_certified` certifies at 74.1899 cm3 with the
magnesium where GEM-Selektor puts it.

Each missing phase is given what the recipe could make of it,
`min_c b_c / A_cs` over the components it consumes, scaled by a tenth — a
chemical bound, not a guess at the answer. A carbonate in a system holding
1e-9 mol of carbon is offered 1e-9 mol and no more. What matters is only that
the phase starts well away from the boundary, since being *at* the boundary is
what the active-set loop cannot recover from.

### Fixed — a rung is accepted only if it conserves mass

Two things can go wrong at a rung of the continuation, and only one of them
raised. A back end can throw, which was handled; or it can *return* a
composition that violates the element balance, which was accepted and carried
forward as the start of every rung after it. The walk then walks away from the
problem it was posed.

`homotopy_initial_state` now tests each rung before accepting it, and a refused
rung is retaken by halving the distance back to the last `λ` that worked, up to
`max_bisections` times per target — the fixed `steps` ladder becomes a
suggestion, and a rung that cannot be taken in one jump is taken in two.

The test is `|rᵢ| ≤ balance_atol + balance_rtol · scaleᵢ` on the element-balance
residual, row by row, with `scaleᵢ = max(|bᵢ|, Σⱼ |Aᵢⱼ| nⱼ)`, and it **has** to
be mixed rather than relative. Two conservation rows of this problem carry a
legitimately negligible budget: electroneutrality is exactly zero, and a cement
recipe is routinely given a carbon trace of 1e-9 mol. Measured at λ = 0.01 on
the paste, a residual of 1.7e-10 mol on the charge row scores 168 against its
own budget and 1.1e-10 mol on the carbon row scores 10.6 — both physically
nothing. A purely relative criterion rejected 35 of 66 rungs on that walk and
never reached its first three targets; the mixed one accepts all ten rungs and
reaches every target. A row holding 1e-11 mol cannot be balanced better than the
solver's absolute floor, and asking it to be is a category error.

Neither tolerance certifies anything, and both are exposed as keywords rather
than buried as constants. `balance_atol` sits above the accuracy the interior
point itself reaches — about 3e-6 mol on this class of problem — because a rung
is a guess and not an answer. The certificate judges the result, afterwards, and
it is unchanged.

Measured on the CEM I paste at w/c = 0.5, the answer is the same and its balance
is better: certified at 74.1899 cm3 against GEM-Selektor's 74.2136 and
pH 13.0994 against 13.0957, with an element balance of 4.4e-15 where it was
6.3e-14.

The coupled kinetic step is untouched by construction: `implicit_step` passes
`autostart = false`, so it never enters the continuation.

### Fixed — a successful solve printed warnings about its own candidates

A call ending `optimal = true` still printed "equilibrium solve returned
`MaxIters`" and "the dual equilibrium solve did not certify optimality", which
reads as a failed solve and is not one. Those come from **candidates**:
`equilibrate_certified` runs every back end from several compositions precisely
because none of them works on every problem, and keeps whichever answer the
certificate proves, so a candidate falling short is what the search is for.

An internal scope now marks the stretches where starting points are computed or
tried — the back-end solves, the continuation's rungs, the multi-start search
itself, and the coupled kinetic step's warm start — and the two diagnostics stay
quiet inside it. Nothing is hidden: the verdict on the *answer* is still
pronounced once, on its certificate, and `verbose = true` still reports every
rung and every rejected start.

The same distinction fixes something worse than noise. The back-end start solves
are wrapped in a `try`, so under `STRICT_CONVERGENCE[] = true` a `MaxIters`
candidate **raised**, was swallowed, and the search silently lost it — leaving a
caller who asked for strict results with a worse search than one who did not.
Measured on a CEM I paste where the interior point ends on `MaxIters`: with the
flag set the route returned an element balance of 27.6 mol, and with it clear the
very same call returned 1.8e-14. A start is not a result, and the flag now
applies only to results.

### Fixed — the conservation matrix reached the solver with an abstract type

`ChemicalSystem` stores its stoichiometry as `Matrix{Real}` whenever integer and
rational coefficients coexist, which a cement's does — `C3AFS0.84H4.32` and its
kind. That is right for the chemistry and wrong for the solver: an abstract
element type boxes every entry and turns `mul!(res, A, x)` into the generic
fallback with a dynamic dispatch per element, on a product evaluated at every
objective and constraint call. SciMLBase had been saying so on every run of such
a system, warning that "arrays or dicts to store parameters of different types
can hurt performance" — a warning that looked like noise about the library's
internals and was in fact a correct report of a type instability in the hot loop.

`EquilibriumProblem` now narrows the conservation matrix, and the default
`b = A * u0` which inherits the same abstract element type, to a concrete
floating-point array. `DualEquilibriumSolver` had always converted; the
interior-point path had not. A caller who passes a concrete array keeps exactly
what they passed — an exact `Matrix{Rational{Int}}`, a `Matrix{Int}`, or a
`Matrix{<:Dual}` for someone differentiating through it — and the exact
stoichiometry is untouched in `system.SM.A`, where it belongs.

### Fixed — the multi-start search could discard the answer it just computed

`equilibrate_certified` ranks the answers of its multi-start search, and the
comparison was on the **stationarity alone**. A composition can be stationary to
1e-3 while violating mass conservation by moles, and that is not a near-answer:
it is not an answer to the problem posed. Reported from that Windows run, the
first route came back stationary to 2.4e-3 with an element balance off by 6.7 mol
and a phase supersaturated by 45 — and it beat every candidate the continuation
produced, because those were stationary to only 1e-2 while conserving mass. The
search computed a usable answer and threw it away, which is why the certificate
came back **bit-identical** across two successive library fixes.

The ranking is now `_kkt_error` — the worst of stationarity, element balance and
supersaturation, with a negative supersaturation counted as zero since every
absent phase undersaturated is what optimality requires. That is the ranking
`solve_certified` already used internally over its own starts, so the two agree
where they previously disagreed, and there is one definition of the quantity
instead of two.

The certificate of a failed solve now also reports what the automatic initial
approximation did — not reached, declined, produced no usable start, ran without
improving, improved without certifying, or certified. On a solve that fails on
one machine and not another, that is the first thing anyone needs to know, and
it should not require a second run with `verbose = true`.

### Changed — `STRICT_CONVERGENCE[] = true` now raises on an uncertified answer

This is what should have made that Windows run loud instead of plausible. The
flag was read only on the interior-point return code, so a caller who set it —
asking that a non-converged solve never pass as a result — still got an
uncertified answer out of `equilibrate_certified` with nothing but a `@warn`.
An uncertified answer from that route is precisely what the flag exists to
refuse: it can violate mass conservation by moles and still be an ordinary
`ChemicalState`.

**This changes behavior for code that sets the flag**, which is why it is called
out rather than filed under fixes. The default is untouched: with
`STRICT_CONVERGENCE[]` at its default `false`, an uncertified answer is still
returned with a warning, and `optimality_certificate` is still the way to audit
it. Only the opt-in path is affected, and only in the direction the flag asks
for.

One consequence had to be handled inside the package. The coupled kinetic step
computes its warm start with `equilibrate_certified`; that answer is a *starting
point*, not a result, so the flag is cleared around it and restored in a
`finally` — the same distinction the continuation already makes for its own
rungs. Left set, the strict flag would have turned that guess into a raise, the
surrounding `catch` would have swallowed it, and a caller asking for strict
results would have silently got a worse start than a caller who did not.


## v0.15.0 — the alkali end-members of the C-S-H, and a readable aqueous state

The shipped `data/solid_solutions.toml` described a C-S-H that could not hold
alkalis. `CSHQ` was declared with four end-members, and CEMDATA18's `KSiOH` and
`NaSiOH` — the alkali-uptake end-members of that same phase — were reachable
from no file in the package, although both are fully parameterized in both
shipped databases (Cemdata18, Lothenbach et al. 2019; the CSHQ model itself is
Kulik 2011). Uptake of alkalis by the C-S-H is not a refinement: it is what sets
the pore-solution pH of a real paste. `docs/src/examples/cement_wc_ratio.md` said so already, in its
list of limitations — "the alkalis themselves, which in a real paste raise the
pore-solution pH to 13 or above".

Measured against a GEM-Selektor reference on a CEM I at w/c = 0.5 (100 g of
oxides, 50 g of water, the CEMDATA18 phase list, and the extended
Debye-Huckel model that run used): GEMS puts **74 % of the total potassium and
90 % of the total sodium into the C-S-H**. With the six-end-member phase
declared and `HKFActivityModel(Ḃ = 0.097637, Kₙ = 0.0)` on an ion size of zero,
this package now returns pH 13.0994 against GEMS' 13.0957, a total volume of
74.1888 cm3 against 74.2136, `KSiOH` to -0.46 % and `NaSiOH` to +0.06 %, and
activity coefficients within 0.25 % (monovalent) and 1.16 % (divalent) of the
ones GEMS printed. Reaktoro 2.13 on the same problem agrees to 0.08 % on the
volume.

### Added

- **`CSHQ` now has six end-members**: `CSHQ-TobD`, `CSHQ-TobH`, `CSHQ-JenH`,
  `CSHQ-JenD`, `KSiOH`, `NaSiOH`.
- **`C3(AF)S0.84H`**, the CEMDATA18 Fe-siliceous hydrogarnet
  (`C3AFS0.84H4.32` + `C3FS0.84H4.32`), ideal mixing. It is where the iron of a
  ferrite phase ends up, and it took 0.0507 mol in the reference run — second
  only to the C-S-H and the portlandite among the aluminate and ferrite
  hydrates. `build_solid_solutions` therefore returns six phases where it
  returned five.
- The `"C-S-H"` reporting group of `volume_fractions` covers the two alkali
  end-members, which would otherwise have been counted under `"other"`.

### Added — the aqueous state is readable

The activity closures computed the molalities, the ionic strength and the
activity coefficients on their way to the log-activities, and kept all three to
themselves. Anyone comparing a state against GEM-Selektor, PHREEQC or Reaktoro
needs them species by species, and had to reach into
`activity_model(cs, model)` and `ChemistryLab._build_params(state)` to get them.
They are public now, together with the two things that make them meaningful:

- **`molalities(state)`** — `mᵢ = nᵢ / (n_w Mw)` for every solute, mol/kg.
- **`ionic_strength(state)`** — `I = ½ Σ mⱼ zⱼ²`, mol/kg. Model-independent: it
  is a property of the composition, and it is the first thing to compare
  against another code, because an ionic strength that disagrees means the two
  are not describing the same solution whatever their volumes agree on.
- **`activity_coefficients(state, model)`** — γᵢ of every aqueous species,
  evaluated from **the model's own formula** rather than as a ratio `a/m`. The
  ratio agrees for an abundant solute, and the tests check that it does, but a
  species parked at the solver's 1e-16 mol lower bound has its log-activity
  dominated by the closures' `+ ϵ` regularization, and the ratio then returns
  values of order 1e300 for a charge class whose only members are trace. The
  formula depends on the ionic strength and the charge alone, so it is exact at
  any amount. Measured on a CEM I pore solution: the charge classes |z| = 3 and
  4 now come back at -2.7 % and -4.7 % of the coefficients GEM-Selektor
  reports, where the ratio gave 1e300.
- **`log_activities(state, model)`** and **`activities(state, model)`** — the
  vector the Gibbs energy is built from, for every species.
- **`concentration_scale(model)`** — `:molality` or `:molarity`. An activity is
  a number on a scale and nothing in the number says which:
  `DiluteSolutionModel` is on the molarity scale but takes ρ = 1 kg/L, so its
  activities coincide numerically with molalities, while the other two are on
  the molality scale. A custom model should declare its own.
- **`pH(state, model)`** and **`pOH(state, model)`** — `−log₁₀ a(H⁺)` and
  `−log₁₀ a(OH⁻)`, the **activity** convention that GEM-Selektor, PHREEQC and
  Reaktoro report. The existing one-argument `pH(state)` is a different
  quantity: `−log₁₀ c(H⁺)` in mol/L over the computed liquid volume,
  reconstructed through `pKw` when the solution is basic. On a Portland cement
  pore solution at I ≈ 0.2 mol/kg, with γ(H⁺) ≈ 0.61, the two differ by
  **0.21 units** (13.31 against 13.10) — enough to be mistaken for a modeling
  error. Both are now documented side by side.

### Added — the starting point is found, not asked for

`equilibrate_certified` computes an initial approximation when its ordinary
starting points fail to certify, so a realistic cement is solvable without the
caller knowing anything about the answer.

This was the practical obstacle. From the cold state of a CEM I paste of 135
species — all the mass in the reactants, every product at the `ϵ` floor — **no
back end reached the optimum**: the answer came back `optimal = false` with a
worst supersaturation of order 1e1 and a total volume 13 % wrong, and
`equilibrate` exited on `MaxIters` in both `Val(:linear)` and `Val(:log)`. The
problem is convex, so this was never a local minimum; it is the conditioning of
an interior-point method started against the boundary. With 120 of 135 species
at 1e-16 and nine at ~1 mol the barrier gradients span sixteen orders of
magnitude, and a solid-solution end-member has `ln a = ln x → −∞` as its mole
fraction goes to zero, so the objective's gradient is unbounded on exactly the
face where a mixing phase vanishes.

`homotopy_initial_state(state)` walks the solute amount up from a dilute system:
everything but the aqueous solvent is scaled by `λ`, and `λ` goes to 1 with each
step started from the answer to the previous one. At `λ = 1` the composition is
the one given, so the element balance is unchanged; it is a starting point, not
a certified equilibrium. Nothing in it is chemical — "everything but the
solvent" needs no guess at which hydrates will form, and `λ` is not a physical
parameter.

Measured on that paste, `equilibrate_certified(state)` from the cold state now
certifies and returns a total volume of 74.1888 cm3 against GEM-Selektor's
74.2136 (-0.033 %) and pH 13.0994 against 13.0957, identical to what a
chemically informed seed gives. With `autostart = false` the same call fails, at
83.7 cm3.

The route also **restarts from its own answer**. The continuation ends on a
composition that is nearly the equilibrium but not certifiably so, and returning
it as the least-bad answer left the last step to the caller: measured on that
same paste under the per-species Debye-Huckel model, stationarity 9.9e-7 and no
certificate, where one more solve started from that answer gives 1.5e-16 with
the worst absent phase 1.4e-5 below saturation. It is the same observation that
motivates the continuation, applied once more — a start near the answer is what
this problem needs, and the best one available is the answer already in hand.
Bounded, and it stops as soon as a round buys nothing.

Doing that exposed a flaw in how rounds were compared: ranking on the KKT error
whenever the new answer was uncertified let an uncertified point with a smaller
stationarity displace a certified one, trading a proof for a residual. It could
not fire while the incumbent was always uncertified; the restart loop reaches
that rule from a certified state, so the optimality flag is now compared first,
in both directions.

Two design points:

- **It costs nothing in the ordinary case**, because it only runs when nothing
  else certified, and its answer is kept only if it is actually better —
  certified beats uncertified, and among uncertified the smaller KKT error wins.
- **`autostart = false` declines it**, and the coupled kinetic step passes that.
  There the caller already supplies a starting point — the previous instant of
  the integration — and a handful of extra solves inside an implicit ODE step
  would be paid at every step. This is the same principle as a coupled Reaktoro
  run, where the solver is carried across instants.

Differentiability is unaffected: under `ForwardDiff` the certified route strips
to the primal state, solves in `Float64` and attaches the sensitivity through
the implicit function theorem, so the continuation never sees a `Dual` and the
derivative does not depend on how the starting point was found.

The walk is done under `DiluteSolutionModel` whatever the target model is, and
deliberately: walking under the extended Debye-Hückel model with a common ion
size of zero runs away to an ionic strength of 18 mol/kg, its coefficients
falling with `I` raising solubility raising `I`. The ideal model has no such
feedback, and its endpoint is a good start for the non-ideal one.

For the record, the approach this replaced: GEM-Selektor computes its initial
approximation by linear programming — `AutoInitialApproximation` in GEMS3K's
`ipm_simplex.cpp`, an "LPP-based automatic initial approximation of the primal
vector x" using a "modified simplex method with two-side constraints" (Kulik et
al., *Comput. Geosci.* 2013). Reproducing that needs a genuine LP solver. Posing
the same linear program and handing it to the barrier method here does **not**
work: measured, it does not move off the cold state, leaving the linear
objective at -723.8 where -755.9 was feasible. A simplex-based initial
approximation remains the principled option and would need an LP dependency.

### Added — the rest of the CEMDATA18 solid solutions, and a model of any arity

`data/solid_solutions.toml` went from six phases to **eleven**. The five added
complete the set of multi-end-member phases a GEM-Selektor CEMDATA18 run of a
Portland cement is given: `Straetlingite_ss` (straetlingite / straetlingite7),
`AFm_SO4_OH` (C4AH13 / monosulphate12), `AFt_SO4_CO3` (tricarboalu03 /
ettringite03_ss), `Hydrotalcite_AlFe` (Mg3AlC0.5OH / Mg3FeC0.5OH) and `MSH`
(M075SH / M15SH). Every end-member is in both shipped databases.

Ideal mixing is an **assumption** there, and the file says so. CEMDATA18
documents non-ideal parameters for some of these; they are not reproduced
because they could not be sourced with confidence, and a mixing parameter
written from memory is worse than an ideal model honestly labeled. The file also
records what was measured on a CEM I paste: all five are undersaturated
(LogSI −0.032 to −9.981, and a binary ideal solid solution gains at most
log10 2 = 0.301 over its best pure end-member, so mixing cannot bring them in),
and a solid solution whose every end-member sits at the solver's lower bound is
numerically awkward — Reaktoro fails to converge when any single one of the five
is declared on that paste.

**`RegularSolutionModel`** fills a gap in arity, not in chemistry.
`RedlichKisterModel` is the general binary form and is restricted to two
end-members; `IdealSolidSolutionModel` takes any number but no interaction at
all. A C-S-H with six end-members, or the CNASH and ECSH families of CEMDATA18,
had no non-ideal option. The new model is the symmetric multi-component
Margules form, `G^ex = Σ_{i<j} W_ij x_i x_j`, with

    ln γ_k = (1/RT) [ Σ_{j≠k} W_kj x_j − Σ_{i<j} W_ij x_i x_j ]

and `W` in J/mol, symmetric. The tests check the three things that matter: it
reduces **exactly** to `RedlichKisterModel(a0 = W₁₂)` for two end-members, it
satisfies `Σ x_k ln γ_k = G^ex/RT` in a ternary — the Gibbs-Duhem consistency a
hand-written `ln γ` usually gets wrong — and a six-end-member phase accepts it
where Redlich-Kister raises. `build_solid_solutions` reads it from
`model = "regular"` with either `w` (a binary) or `W` (a full matrix).

### Added — two ionic strengths, and a Setschenow coefficient per species

`ionic_strength(state; kind = :effective | :stoichiometric)`. The effective one
(the default, and the previous behavior) sums over the speciated free ions, so a
neutral pair such as `Ca(SO4)@` contributes nothing; it is the one every activity
model here is a function of. The stoichiometric one sums as if every complex were
fully dissociated over the system's primaries, which is the analytical ionic
strength of the recipe and what some salting-out and diffusivity correlations are
fitted against. The gap between them measures how much salt is associated: 0.7 %
on a Portland cement pore solution, far more on a sulfate brine. The
stoichiometric sum reuses the decomposition the mass balance already uses
(`SM.A`), so it needs no separate table of dissociation reactions.

Neither is a difference of *formula* — all three codes compute `½ Σ mⱼ zⱼ²`, and
the spread between GEM-Selektor (0.2097), this package (0.2121) and Reaktoro
(0.2178) on that pore solution comes from their converged compositions, not from
their definitions.

The Setschenow (salting-out) coefficient of a neutral aqueous species is now read
from `sp[:Kₙ]` when present, falling back on the model's global `Kₙ`. CO₂(aq),
the noble gases and the neutral silicates are not equally salted out, and one
coefficient for all of them was the B-dot literature's simplification rather than
a fact. No table of coefficients is shipped: a table is data, and data belongs in
a database or in the caller's hands.

### Verified — the water activity, against Gibbs-Duhem

A benchmark was added for the property that separates this package from its
neighbors. For a 1:1 electrolyte the osmotic coefficient of the extended
Debye-Hückel model has a closed form, and integrating Gibbs-Duhem over the
model's own activity coefficients must reproduce it. `HKFActivityModel` does, and
the test pins NaCl(aq) at 25 °C to a_w = 0.996657, 0.983603 and 0.966898 at 0.1,
0.5 and 1.0 mol/kg — the Gibbs-Duhem-consistent values, which coincide with the
tabulated water activities of NaCl.

Reaktoro's `ActivityModelDebyeHuckel` does not. Measured on the same three
molalities, its reported water activity is 0.20 %, 1.70 % and 3.92 % below what
its **own** activity coefficients imply, so the inconsistency is internal and
needs no external data to demonstrate. On a CEM I pore solution it reports
a_w = 0.8824 where GEM-Selektor gives 0.992588 and this package 0.993723, and the
consequence is quantitative: `ettringite` and `ettringite30` differ by two water
molecules, so their ratio goes as `1/a_w²`, and `(0.9926/0.8824)² = 1.265`
against a measured ratio of ratios of 1.266. That single discrepancy accounts for
the whole difference in the AFt split. Anything whose stoichiometry differs by
water — the C-S-H hydration states, the AFm and AFt series, the hydrogarnets — is
sensitive to it, which is most of a cement.

### Added — a common ion size, settable on the model

`HKFActivityModel(; å)` imposes **one** effective radius on every charged
aqueous species, short-circuiting the per-species tables. That is what
GEM-Selektor, PHREEQC's `-gamma` and most published cement models actually use,
and it was previously reachable only by mutating `sp[:å]` on every species
before building the `ChemicalSystem`. `å_default` looks like the knob for it and
is not: it is the last resort of the lookup chain and is never reached for an
ion either table covers, so setting it changes essentially nothing. Both facts
are now stated in the docstring and asserted in the tests.

`å = 0` collapses the Debye-Hückel denominator to 1, giving the limiting law
plus the B-dot term. With `HKFActivityModel(å = 0.0, Ḃ = 0.097637, Kₙ = 0.0)`
the package reproduces the activity coefficients of a GEM-Selektor CEMDATA18
run to 0.25 % on the monovalent ions and 1.2 % on the divalent ones — a run
whose model can be recovered from its own output, since CEMDATA18 carries no
ion-size parameter at all.

### Fixed

- `STRICT_CONVERGENCE[] = true` disabled the new initial approximation instead
  of hardening it. The continuation's early rungs are *expected* to fall short —
  they are starting points on the way to `λ = 1`, not results — so under the
  strict flag each one raised, was caught, and every later rung started cold
  again, leaving the walk useless in exactly the mode a careful caller turns on.
  The flag is now saved, cleared for the duration of the walk and restored in a
  `finally`, so a non-converged *result* still raises while an intermediate rung
  does not. Found by running a caller's script with the flag on, and pinned by a
  regression test that checks both the outcome and that the flag comes back.
- `docs/src/tutorials/equilibrium.md` and `README.md` built an AFm solid
  solution from `dict["Ms"]` and `dict["Mc"]`. Neither symbol exists in any
  shipped database, so those examples could never have run. They now use
  `monosulphate12` and `monocarbonate`, which is what the TOML has always used.
- The TOML-format example in `docs/src/tutorials/databases.md` showed the same
  two symbols, and a four-member `CSHQ` that no longer matches the file.
- `docs/src/examples/simplified_clinker_dissolution.md` said "the only built-in
  model is `DiluteSolutionModel`". There are three.

### Documentation

- A new tutorial section,
  [Reading the aqueous properties back](https://micropochemomechanics.github.io/ChemistryLab.jl/stable/tutorials/equilibrium/#sec-aqueous-properties),
  covering the accessors above, the two conventions of pH and why γ is not a
  ratio; and an `Aqueous properties` section in the equilibrium API reference.
- The bibliography gained Helgeson (1969), Helgeson, Kirkham & Flowers (1981),
  Davies (1962), Kulik (2011), Parkhurst & Appelo (2013), Robie & Hemingway
  (1995) and Xu et al. (2011), which the activity-model docstrings had been
  citing in prose without entries. Every DOI was resolved against the Crossref
  REST API; Davies (1962) is a Butterworths monograph that Crossref does not
  index, and its author, title and publisher are confirmed by the
  contemporary review in *Science* (doi:10.1126/science.143.3601.37), while the
  authorship of USGS Bulletin 2131, absent from its Crossref record, is
  confirmed by the USGS publications catalog.
- `data/solid_solutions.toml` now says at its head what it is **not**. `CSHQ`
  and `C3(AF)S0.84H` match the GEM-Selektor phases of those names, but `AFm`,
  `Hydrogarnet` and `Hydrotalcite` are deliberate alternatives to the CEMDATA18
  phase model: GEMS treats `monocarbonate`, `C3AH6`, `C3FH6` and `hydrotalcite`
  as *pure* phases, its AFm solid solution is `C4AH13` + `monosulphate12`, and
  its hydrotalcite solid solution is `Mg3AlC0.5OH` + `Mg3FeC0.5OH` at
  Mg:Al = 3. Reproducing a GEMS result means declaring the phases in the script,
  not taking this file wholesale — and the file no longer lets a reader assume
  otherwise.


## v0.14.2 — holding AMD at a version that still has `SS_Int`

A release with no change to ChemistryLab itself. It exists because an upstream
patch release broke every environment that resolves this package's test target
or its documentation.

AMD 0.5.4, published on 2026-09-04, removed `SS_Int`. `SparseColumnPivotedQR`'s
AMD extension calls it and bounds AMD only at `"0.5.1"`, so any fresh resolve
picks 0.5.4, and `LinearSolve` — hence `OrdinaryDiffEq`, hence this package's
test and docs environments — fails to precompile. Nothing in ChemistryLab is at
fault and nothing in it can avoid the clash.

`AMD = "0.5.1 - 0.5.3"` is therefore declared in the test target and in
`docs/Project.toml`. AMD is **not** a dependency of the package: the bound
constrains nothing downstream, only the environments built here. Remove it once
AMD restores the binding or SparseColumnPivotedQR tightens its own bound.


## v0.14.1 — the ambiguities a test suite cannot see

`Aqua.test_all` now runs as part of the suite. It checks eight properties no
domain test looks at — method ambiguities, unbound type parameters, undefined
exports, dependency hygiene, type piracy, tasks left running at load — and it
found three real defects, all fixed here. No exported name, signature or result
changes; the calls below simply used to fail.

### Fixed

- **Fourteen method ambiguities from one signature.** `convert(::Type{T},
  f::Formula{T}) where {T}` left `T` unconstrained, so dispatch weighed it
  against every `convert(::Type{T}, x)` in the ecosystem — `Missing`,
  `Nothing`, `Ref`, `FunctionWrapper`, `VecElement`, polynomial types and more.
  `Formula{T <: Number}` already implies the bound; saying it in the signature
  costs nothing and removes all fourteen.
- **`convert(::Type{<:ForwardDiff.Dual}, ::Formula)`** was ambiguous with
  ForwardDiff's own catch-all and therefore unusable. It now resolves, spelled
  with the same type parameters ForwardDiff uses so that it is genuinely more
  specific.
- **`promote_rule(Species, Species)`** was the missing diagonal of the
  `Species`/`AbstractSpecies` pair, and errored.
- **`Species(...)` and `CemSpecies(...)`** declared a `where {T}` bound only by
  the varargs, so `T` was unbound when no pair was passed. `T` was never used in
  the body; the bound is gone.
- **`gather_species`** on a dictionary whose keys *and* values are species was
  ambiguous in dispatch and in meaning. It now says which two readings are
  possible instead of raising an unexplained `MethodError`.

### Declared, not fixed

`thermo_factories.jl` extends some thirty `Base` math functions to
`DynamicQuantities.Quantity`, stripping the unit first. That is type piracy —
function and type both belong elsewhere — and it is deliberate. It is now
declared to Aqua as such rather than left unexamined. It remains global: any
code loading ChemistryLab gets these methods whether it asked or not.

## v0.14.0 — an equilibrium that comes with a proof, and a kinetic step that is one problem

### Breaking changes

- `equilibrate(state)` now solves through **every** loaded back end and returns
  the answer `optimality_certificate` proves optimal. It is slower and, on some
  systems, gives a different — correct — answer. `equilibrate(state; certify =
  false)` restores the single-back-end behavior.
- `OptimaSolver = "0.5"` is required. The parameter block the new constraints and
  the kinetic step ride on arrived there; against 0.5's predecessor the extension
  fails with a `MethodError` on `DualNewtonProblem`.
- The `equilibrate` keyword `model` is now named explicitly in the one-argument
  form, and `certify` and `constraint` are new keywords. Code passing positional
  arguments is unaffected.

### An interior-point answer whose charge balance was wrong in the second digit

On 1 mmol of calcite dissolving in 1 kg of water, the interior-point path returned
`‖An − b‖∞ = 3×10⁻⁶`, which looks like round-off. It is not: the rows of a
conservation matrix do not share a scale, and against its own budget that residual
is **3 % on the charge row** and 0.3 % on the carbonate, while the water row at
55.5 mol sits at `1.8×10⁻⁸`. The pH came out 9.5210 against a correct 9.5274.

The measure is fixed in `OptimaSolver` 0.5.0. The answer is fixed here, and the
iteration is not the place: traced over twenty-three iterations, `α` is pinned at
its ceiling of 0.15 and `‖dn‖` decays at `1 − α` with the residual frozen — the
fraction-to-boundary rule lets through 15 % of a correction the next iteration
re-poses. Thirty thousand iterations change nothing, nor does a tolerance of
`10⁻⁶`.

What does is that the problem is convex, so its KKT conditions are sufficient and
`optimality_certificate` **decides** rather than ranks. `equilibrate` therefore
solves by every route and keeps a proved answer:

| case | interior point | dual Newton | certified |
|:--|--:|--:|--:|
| calcite in water, 10–40 °C | 3.0e-2 | 3.0e-12 | proved |
| calcite + 1 mmol CO₂ | 1.0e-3 | 6.3e-13 | proved |
| calcite + 50 mmol CO₂ | 1.3e-16 | 8.5e-14 | proved |
| pure water | 1.3e-16 | 1.3e-16 | proved |
| CEM I, w/c 0.45 and 0.60 | 1.5e-14 | 2.0e-14 | proved |
| CEM I, w/c 0.30 | 7.2e-16 | 3.3e-10, **not** proved | proved |

Ten of ten, against seven through the interior point alone and nine through the
dual Newton alone.

One bug found in the process, and it mattered: `solve_certified` recomputed the
component totals from **each** starting point. A start that violates the balance
therefore posed a different problem, and the dual solve certified the answer to
the shifted one — two "certified" compositions 0.2 % apart on dissolved calcium,
which a convex problem with one minimum cannot have. `b` is now fixed once, from
the state as given.

### Equilibrium under something other than fixed (T, P)

`Adiabatic`, `FixedEnthalpy`, `FixedVolume`, `SealedVolume`, `FixedpH` and
`FixedActivity`, passed as `equilibrate(state; constraint = ...)`.

The prescribed quantity is an unknown of the solver's own square system, not an
outer loop: an adiabatic solve costs one equation, not a solve per trial
temperature. Two vehicles, and they are not interchangeable — a prescribed
**property** adds a parameter and its residual, a prescribed **chemical
potential** adds a *column* to the conservation matrix, an unknown amount of a
titrant the system may draw on, and that amount comes back as part of the answer.

Validated against a published number: adiabatic `H⁺ + OH⁻ → H₂O` gives
**−55.85 kJ/mol** at every amount tested, against the accepted −55.8, with
nothing fitted to it. `ΔT = 6.62 K` for 0.5 mol in a kilogram of water, against
27.9 kJ over 4.18 kJ/K. Enthalpy conserved to `10⁻¹²` relative, and solving at a
*fixed* temperature equal to the one found returns the same composition to
`10⁻¹²`.

A volume constraint on a condensed system is **refused**, with the lever it
measured named in the error. The molar volumes of water and of the minerals in
these databases are exactly pressure-independent, so the relative lever
`(∂V/∂P)·P/V` is about `10⁻⁶` — some 9 600 bar to change the volume by one
percent. That is the physics: the volume of an incompressible condensed system is
fixed by its composition. Declare a gas phase, or use `porosity(state, reference)`
for a sealed specimen at fixed pressure.

### A kinetic step as Leal poses it: one problem, fully implicit

`KineticStepSolver` and `kinetic_step`. The reaction extents are unknowns of the
same Gibbs minimization as the amounts and the element potentials, and the rate is
evaluated at the **end-of-step** composition. No frozen speciation in a right-hand
side, and no lag between the kinetic and the equilibrium species. The reactivity
constraint being linear, it joins the conservation block, so the algebraic cost of
kinetics is the number of **reactions**, not of species.

Measured against an analytic answer: at a constant 1 µmol/s, `Δξ = Δt·r` to
`10⁻¹⁰` relative and the calcite left is `n₀ − kΔt` to the last bit, with the
element balance at `10⁻¹⁴`.

Measured on the feedback path nothing exercised before: with `r = k(1 − Ω)` at
`k = 10⁻⁵ mol/s`, an explicit step of `10⁶ s` would dissolve 10 mol — a thousand
times the calcite present. The implicit step lands at `Ω = 0.999989`, approaching
saturation from below and never crossing it, at any step size.

A **solid solution** takes part: C₃S dissolving into a four-member CSHQ solution
comes out certified at a stationarity of `2×10⁻¹³`, with all four end-members
present and ten times the step giving ten times the C-S-H. Two settings are
decided by the certificate rather than guessed, because neither answer works on
both kinds of problem — `warm_start` for the admission of a mixing phase from a
cold start, `pin_minerals` for whether the kinetic minerals are held in the active
set. Both are documented with the measurements that decide them.

The certificate of a kinetic step is taken on the **augmented** problem. A mineral
held back by a rate law is supersaturated by construction, so tested against the
unconstrained equilibrium it reports a residual of 7 `RT` and a correct answer is
called wrong.

`coupling = :species`, which would impose the assemblage instead of proposing it,
**raises**. It was implemented and does not converge: pinning a pure phase by a
linear row leaves its stationarity row in the system as well, made trivially
satisfiable by that row's own multiplier, and the Jacobian degenerates — the
extents came out 4.3 times short. For a model whose assemblage IS the
stoichiometry, the ODE route already does exactly that, and the tutorial now opens
with a table saying which route answers which question.

### Choosing the step length, and a case the stiff ODE route gets wrong

`kinetic_step_adaptive` chooses `Δt` from a Richardson estimate on the extents —
one step against two half-steps, which for a first-order method is the error of
the coarse one. Reaktoro has no equivalent: the step is the caller's and nothing
reports what taking it cost.

Measured on calcite under `r = k(1 − Ω)`, `k = 10⁻⁴ mol/s`, over `10⁵ s`:

| route | steps | result |
|:--|--:|:--|
| `Rodas5P`, and OrdinaryDiffEq's default polyalgorithm | 6 | **extent −457 mol**, `retcode = Success` |
| `Tsit5`, explicit | 85 626 | correct, 519 s |
| one `kinetic_step` of `10⁵ s` | 1 | wrong, and its certificate says so |
| `kinetic_step_adaptive` | **7** | correct to eight digits |

The cause is not a missing Jacobian term, which is the natural guess. The ODE
route's residual reads the speciation **frozen** at the last accepted step, so
`∂(du)/∂bₑ = 0` is exact for the system being integrated. What goes wrong is that
the frozen speciation makes the right-hand side inconsistent with the state within
a step: an implicit method steps past the point where `Ω` crosses one, the rate
changes sign, and the run enters a branch it never leaves. Three measurements
settle it — bounding `dtmax` to `10³ s`, 103 steps against 6, returns the
identical wrong value to seven digits; removing the re-speciation returns a sane
answer; and routing the re-speciation through the certified route moves −457 to
−383. So the fix is not to freeze, which is what the implicit step does.

The `OrdinaryDiffEq` extension now warns when a trajectory ends on amounts no
chemistry can produce, naming the species and the budget it exceeded. That is
what the ODE route can offer; a rate law that reads the solution belongs on the
implicit one.

### Re-speciation escalates to the certified route when it fails the balance

`respeciate!` went through the interior point alone, and where that is wrong it is
wrong by percent. It now escalates: when the plain solve's absolute element-balance
residual exceeds the retry tolerance, the partition is re-solved through every
back end and the answer the KKT certificate proves optimal is kept. The cost of a
normal run is unchanged, because the escalation fires only on the solves that need
it.

### `ode_solver = :auto`

Hands the choice of integrator to `OrdinaryDiffEq`'s default polyalgorithm, which
detects stiffness at run time. It is offered rather than made the default: on
these problems it lands on the same stiff method as `Rodas5P` and returns the same
answers, so switching would change nothing while removing a reproducible choice.

A step is accepted on the estimate **and** on the certificate. Richardson's
difference is blind to an error the two resolutions share: on calcite over
`10⁵ s` the coarse step and both half-steps each dissolve the entire mineral,
their extents agree to `5×10⁻¹¹`, and the estimator grades that step excellent.
The certificate is not fooled, since such a composition violates
`Δξ − Δt·M·r(n) = 0` by the whole extent. This also closed a hole in
`kkt_certificate`, which checked stationarity, the linear rows and phase
stability but not the **nonlinear** parameter residual, so a step violating its
own rate equation could be proved optimal (`OptimaSolver` 0.5.0).

One further design point, because the obvious choice fails. The tolerance
is relative to the amount each reaction acts on, **not** to the extent. Scaling by
`Δξ` makes the tolerance vanish with the step while the equilibrium solve's noise
does not, so the measured error grows as the step shrinks: the controller halved
to the floor without advancing, and the answer got worse as the tolerance was
tightened — 2.5e-5 at `reltol = 1e-3` against 9.2e-5 at `1e-5`.

### `coupling = :species`: prescribed products inside the strong coupling

The third combination, which neither Reaktoro nor the ODE route offers. One
constraint per kinetic species, `nᵢ = nᵢ(0) + Σⱼ νᵢⱼ Δξⱼ`, so a solid-to-solid
scheme in the form of [Lavergne2018] imposes its products while the aqueous phase
still minimizes `G` under element conservation. Measured on two C₃A pathways, the
extents come out at `Δt·r` exactly and every species follows the stoichiometry to
eight or nine digits, ettringite and monosulfoaluminate included, where
`:reactions` leaves the AFm at zero because the minimization is free to rearrange
the solids.

It was refused in an earlier draft of this release, and finding out why turned up
a genuine defect in `OptimaSolver`: a pinning row for a product that starts absent
has one positive entry and a zero budget, exactly the shape `degenerate_components`
reads as "this component is absent from the system". It pinned that row's
multiplier and declared the species dead — stationarity residual 458, extents 4.3
times short. `conservation_rows` now restricts that criterion to rows that
conserve something.

The pinned species are **eliminated**, not constrained, and that is what makes it
work. Holding them by a linear row — the obvious implementation — stalls: the
species' stationarity row stays in the system, satisfied by a multiplier that must
reach the mineral's own chemical potential, of order 10²–10³ in `RT` units, and
the Newton sits at a fixed point its line search cannot leave. Measured, an
element balance of `6.1×10⁻⁷` mol whatever `maxit`, `tol` or the number of
active-set updates. Their amounts being an explicit affine function of the
extents, they are removed from the system instead, their element content is
subtracted from the budget, and a plain equilibrium is solved over what remains
with a Newton on the `nr` extents around it: the balance goes to `2×10⁻¹⁴` and the
step is certified.

### A single step far beyond the relaxation time can find the other root

For a rate law that vanishes at equilibrium the implicit step has two solutions,
and the second is the composition with the mineral wholly dissolved: it satisfies
the element balance and the reactivity row while violating `Δξ − Δt·M·r(n) = 0`
by the entire extent. The certificate refuses it, and `kinetic_step` reports the
step uncertified rather than passing it off.

Which root the Newton finds is a property of the build, not of the chemistry.
Measured on calcite under `r = k(1 − Ω)` with `k = 10⁻⁵ mol/s`: steps of `10⁴ s`
and `10⁶ s` converge to the right root on Julia 1.12 and to the other one on
1.13.0-rc4 — which is how this was found, a suite green on one and five failures
on the other. `kinetic_step_adaptive` is build-independent, because it refuses an
uncertified step and halves until one certifies; it reaches the equilibrium
values to eight digits in seven steps on either.

The test now asserts the contract rather than the outcome of a Newton on a
particular build: a certified step never crosses saturation, an uncertified one is
flagged and its rate-equation residual is of order one, steps well inside the
relaxation time certify on any build, and the adaptive march reaches the right
answer. Asserting the physics unconditionally on a step of `10⁶ s` was asserting
which root a given LLVM finds.

### `pin_minerals = :auto` ranks the two formulations by margin

It used to return the first one whose certificate passed. That makes the choice a
**threshold** decision, and a threshold decision between two answers that straddle
the tolerance is settled by the last bits: the same commit then behaves differently
on a different machine. Found exactly that way — the suite green locally and five
failures in CI on the identical commit, in the test asserting that a saturating
rate law never crosses equilibrium.

Both formulations are now computed and ranked by a single margin —
stationarity, element balance, worst supersaturation among absent phases, and the
residual of the rate equation. Ranked, the two are not close: on the case that
exposed it the free branch has run to equilibrium and violates its own rate
equation by the whole extent while the pinned one satisfies it to `10⁻¹³`, fifteen
orders of magnitude apart. The adaptive march resolves `:auto` the same way.

### ForwardDiff through the certified route

A composition carrying dual numbers takes the implicit-function route: the
certified solve runs in real arithmetic and the derivative is attached at the
answer. This had to be built rather than inherited — the certified path converts
the component totals to `Float64` to fix them once, which silently dropped every
partial. `∂pH/∂n(CO₂)` now agrees with a central difference to `4.3e-9` relative,
and the primal value is the certified one.

Pushing duals through the search would be wrong in any case: an active set has a
discrete component, so the map `b ↦ n*(b)` is smooth only piecewise and the
derivative belongs at the solution with the active set frozen.

### Correctness of published values

- **The "maleic acid" titration is malonate.** The SLOP98 entries are
  `MALONIC-ACID,AQ` / `H-MALONATE,AQ` / `MALONATE,AQ`, the three-carbon diacid.
  Their `ΔₐG⁰` give pKa **2.851 / 5.696** against a tabulated malonic 2.83 / 5.69,
  and nowhere near maleic's 1.92 / 6.27. The database was right and the name was
  wrong. The script and the page are renamed, and both pKa are now derived rather
  than hard-coded — the old figures were plotted as dashed lines that crossed the
  curve at no half-equivalence, so the defect was visible in the published figure.
- **`parrott_killoh` is deprecated and no longer attributed to Parrott & Killoh.**
  Its nucleation term carries no Avrami logarithm, `K₃` sits where the canonical
  form has `k₂`, `N₁ = 3.3` is the canonical `n₃`, and `k₃ = 1.1` has no
  counterpart: two different models, not two parameterizations. With `PK_PARAMS_*`
  the diffusion branch takes over at α ≈ 0.003 (C₂S), 0.013 (C₃S) and 0.057
  (C₃A), and those three then land on **α(7 d) = 0.2386 whatever their `K₁`**,
  while C₄AF is limited by its own nucleation branch at 0.193. A CEM I at
  w/c = 0.40 is reported near 0.61. All demos and doc pages now use
  `parrott_killoh_avrami` with `PK84_PARAMS_*`: the CEM I paste moves from
  ᾱ(7 d) = 0.234 to **0.628**, ΔT from 2.0 to **14.2 °C**, and the heat released
  from 115 to **308 kJ/kg**. The slag and metakaolin laws move to `waller`.
- **The calorimeter's `Cp` was double-counted.** The denominator is
  `Cp + Σᵢ nᵢ Cp°ᵢ(T)` with the second term recomputed at every step, so adding
  the sample to `Cp` counts it twice and understates ΔT by a factor of about 1.75.
  Three scripts and one doc page corrected; the docstring, which contradicted
  itself, now carries the warning.
- **A documented temperature sweep ran on an empty state.** `ChemicalState(cs)`
  holds no matter, so every element balance was zero and all 21 points returned
  the same trivial answer — flat lines, with prose describing a curve that was not
  there. Repaired, it also shows the note was backwards: `Kₛₚ` is retrograde
  (`10^-8.411` at 10 °C to `10^-8.517` at 30 °C, matching the database to 0.003
  log units, and `−8.48` at 25 °C is the accepted value), and yet dissolved
  calcium **rises**, because the pH falls 0.4 and shifts carbonate to bicarbonate.
  The familiar statement holds for a system buffered at fixed pCO₂, not a sealed
  one.
- **The w/c page was out of navigation, and not by accident.** Its heavy blocks
  were never executed and the scan contradicts its analysis on three points:
  ettringite never forms (the sulfate all goes to `monosulphate12` at
  equilibrium), no clinker survives at any w/c including 0.30, and the porosity has
  neither minimum nor inflection — so no Powers threshold, which is kinetic and
  not thermodynamic. Rewritten with executed blocks, the sealed-curing porosity
  convention (`porosity(eq, fresh)`, 5.3 points and a 7.4 % volume shrinkage apart
  from the one-argument form at w/c = 0.50), its assumptions written out, and put
  back in the navigation.
- Stale documentation removed: the `@test_broken` on `CaOH⁺` at ×2.47 has been
  gone since the convergence test moved to the true KKT error at `μ = 0`, and the
  `pKw = 13.979` that came with it was computed from the pre-fix `OH⁻` — it is
  13.9994 against Reaktoro's 14.0001.

### New tests

`test/published_values.jl` pins the malonate pKa and the retrograde calcite sweep
against the database's own `Kₛₚ`; `test/certified_equilibrium.jl` the certified
route, including that `b` must be fixed once; `test/equilibrium_constraints.jl` the
adiabatic enthalpy of neutralization, the refused volume constraint and the
implicit titrant; `test/kinetics/test_implicit_step.jl` the analytic kinetic step,
the impossibility of crossing saturation, several reactions on one mineral, and a
solid solution inside the step.


## v0.13.0 — a data file is named, not located

### Breaking changes

No name was removed or renamed and no signature changed. Below 1.0 the resolver
treats a minor bump as breaking regardless, so a downstream
`[compat] ChemistryLab = "0.12"` must be widened to `"0.13"`. In this repository
group that is `MeanFieldHomogenization.jl`'s `docs/Project.toml`, which pins
`ChemistryLab = "0.12"` and will otherwise refuse the new version.

One behavior did change, and it can only turn a former failure into a success:
the database readers now fall back to the bundled `data/` directory when a path
does not resolve against the working directory. A call that already worked keeps
resolving to exactly the same file — see below.

### `datapath`, because a script was only runnable from one directory

`build_species("data/slop98-inorganic-thermofun.json")` resolves its argument
against the working directory. That is the repository root when the script is
launched from a shell sitting there, and almost never otherwise: running the
same file from an editor whose REPL started elsewhere failed with

    SystemError: opening file "data/slop98-inorganic-thermofun.json"

Half the scripts in `scripts/` were written that way and half were not; the test
suite had independently converged on `joinpath(pkgdir(ChemistryLab), "data", …)`
in eighteen places, and the documentation pages carried
`"../../../data/…"`, a depth that silently depended on where the page sat in
`docs/src/`.

The exported [`datapath`](@ref) names a bundled file without reference to any
working directory:

```julia
substances = build_species(datapath("cemdata18-thermofun.json"))
ss_phases  = build_solid_solutions(datapath("solid_solutions.toml"), dict)
readdir(datapath())                     # what ships with the package
```

Every script, documentation block and test now uses it, so the code a reader
copies out of a page is code that runs.

### The readers resolve a bundled name, so older code keeps working

`read_thermofun_database`, `build_solid_solutions`, `extract_primary_species` and
`merge_json` (for its two inputs, never its output) resolve their path argument
in this order: as given, then under the bundled `data/` directory, then relative
to the package root. The working directory comes **first**, which is what makes
this safe: a call that resolves today resolves to the very same file tomorrow, a
local database still takes precedence over a bundled one of the same name, and
only a name that *is* one of the shipped files can reach the fallback — so a
mistyped path to a file of your own still fails, now with an error listing what
is available.

The banner these functions print shows the path relative to the package
(`data/cemdata18-thermofun.json`) rather than the absolute one. Documenter
captures that banner into the built page, and an absolute path would have baked
a build machine's directories into the documentation.

### Also

- `scripts/README.md`, which did not exist, records how to run a script and the
  one rule that keeps the collection working: a script consumed by the
  documentation or the test suite (`ionic_hydration.jl`,
  `hydration_calibration.jl`) never calls `Pkg.activate`, because the active
  project is global process state and switching it mid-build breaks every
  `@example` block that follows, on every page.
- The `Usage:` headers of the scripts said `julia --project` while the scripts
  activated `scripts/`; they now agree, and mention that running the file
  directly from an editor works.
- `using Revise` is gone from the three `tutorial_*.jl` scripts — it belongs in
  a `startup.jl`, and it resolved only through the machine's global environment.
  `SymPy`, which `tutorial_symbolic_reactions.jl` imports, is now declared in
  `scripts/Project.toml`, so that script runs in the environment it activates.
- The `merge_json` snippet in `docs/src/man/advanced.md` named a `.dat` file that
  does not exist (`cemdata18.dat`, where the shipped one is
  `CEMDATA18-31-03-2022-phaseVol.dat`) and passed a bundled path as the output
  argument. Both corrected.

### A red test suite, and the version pin that hid it

`test/kinetics/test_calibration.jl` errored with
`UndefVarError: NelderMead not defined in Main`, taking the whole suite down
(496 passed, 1 errored). The cause is upstream: **OptimizationOptimJL 0.4.19
stopped re-exporting Optim's algorithm names**, and
`scripts/hydration_calibration.jl` reached `NelderMead` through that re-export.

What made it look like a local mystery is that the script kept working: a
`scripts/Manifest.toml` is gitignored, and the one on the author's machine still
pinned 0.4.18, where the re-export was there. The test environment resolves fresh
from the registry, gets 0.4.19, and fails — so the suite was red on continuous
integration while every local run of the script was fine.

The script now does `import OptimizationOptimJL: NelderMead`. The binding exists
in that module either way, exported or not, so the explicit import works on both
versions and no longer depends on which one an environment happens to resolve.

## v0.12.0 — calibrated against a measured heat curve, and an isothermal heat that is a heat

### Breaking changes

**The ODE state layout gained one slot for `IsothermalCalorimeter`.** Code that
reads `sol.u` directly, or carries its own offset arithmetic into the state
vector, has to be revisited. `cumulative_heat(sol, ::IsothermalCalorimeter)` and
`heat_flow(sol, ::IsothermalCalorimeter)` also change what they return — from a
reaction extent in moles to a heat in joules. See below.

No name was removed or renamed and no signature changed. Below 1.0 the resolver
treats a minor bump as breaking regardless, so a downstream
`[compat] ChemistryLab = "0.11"` must be widened to `"0.12"`.

### `cumulative_heat` was returning a number of moles

`IsothermalCalorimeter`'s docstring promised that "cumulative heat
`Q(t) = ∫₀ᵗ q̇(τ) dτ` [J]" was "tracked as an extra ODE state". It was not.
`build_u0` appended a trailing slot only for `SemiAdiabaticCalorimeter`, and the
right-hand side wrote `du[end]` under the same condition, so with an isothermal
device `u[end]` was the **last reaction extent ξ_M, in moles**. `cumulative_heat`
returned it as a heat and `heat_flow` finite-differenced it. The function's other
branch, `Q(t) = H(0) − H(t)` from `p.saved_H`, was unreachable: `saved_H` and
`saved_t` are read in exactly two places in the package and written in none.

The isothermal calorimeter now contributes a real state integrating `dQ/dt = q̇`,
so `n_extra_states`, `extend_u0` and `extend_ode!` — until now reachable only from
the test suite — describe what the integrator actually does. The dead branch is
gone. Every index keyed on the state's length was audited; the semi-adiabatic path
already went through `n_extra_states(cal)` and is unaffected, and `p.has_T`
remains false for an isothermal device, so nothing collides on `u[end]`.

A new test compares both routes to the heat on the stoichiometric formulation,
where both are valid, and requires them to agree — the assertion that would have
caught this. The state-layout docstrings were also wrong in a second way, omitting
the reaction-extent block entirely; both now describe all four segments.

One caveat is stated where it belongs rather than left to be discovered. `q̇` here
is `heat_rate`, the heat of the **kinetic** reactions. That is the heat of
hydration when those reactions produce the hydrates. Under partial equilibrium
they only dissolve the anhydrous phases into ions, the hydrates are precipitated
by the Gibbs minimization, and this sum cannot see their heat — worth hundreds of
joules per gram on an ordinary Portland cement. `integrate` now warns on that
combination, as it already did for the semi-adiabatic device, and points at
`heat_release`.

### Calibrating hydration kinetics against measured calorimetry

`scripts/ionic_hydration.jl` has always run a complete CEM I forward with the
*published* Parrott & Killoh parameters, and its `IONIC_CALIBRATION` dictionary
has always been an explicit, deliberately unused hook. Nothing in the repository
closed that loop: there was no experimental dataset, no loss function and no
inverse problem anywhere in `src/`, `scripts/` or `test/`.

`scripts/hydration_calibration.jl` and
[its documentation page](https://micropochemomechanics.github.io/ChemistryLab.jl/stable/examples/hydration_calibration/)
close it, on real measurements, fitting **kinetic parameters only** — the CEMDATA18
thermodynamics are measured quantities, and the Blaine fineness, w/b ratio and
temperature are reported by the experiment.

**The data.** Two CEM I records now ship under `data/experimental/`, subsets of the
CC-BY-4.0 Zenodo deposit of Šmilauer & Reiterman
([10.5281/zenodo.15212785](https://doi.org/10.5281/zenodo.15212785)): a CEM I
52.5 R at w/b 0.50 and Blaine 397 m²/kg over 262 h, fitted, and a CEM I 52.5 R at
w/b 0.45 and Blaine 415 m²/kg over 617 h, held out. They are the only source found
that is at once openly licensed for redistribution, real CEM I measured by
isothermal calorimetry, and documented well enough to *fix* the non-kinetic inputs
— every file states the fineness, the w/b ratio, the temperature and the
normalizing binder mass.

`regenerate.jl` in that directory rebuilds both files **byte for byte** from the
deposit, so nothing about them is unverifiable, and `--check` is an exact
comparison. The subset is the rows nearest to 500 log-spaced times, by nearest
index rather than by interpolation: every number in the vendored files is one the
calorimeter reported. A fourth column carries the depositors' own fitted curve,
read from their tabulated data rather than reimplemented from a functional form
recalled from memory, so a fit can be judged against somebody else's on the same
record. `data/experimental/LICENSE` states the terms — CC-BY-4.0, separate from
the package's LGPL — and the companion article, published CC BY-NC-ND 3.0, is
cited with nothing reproduced from it.

**The published parameters are already close.** With no adjustment at all the
coupled model gives Q(262 h) = 381.9 J/g against 376.0 J/g measured, +1.6 %. The
depositors' four-parameter affinity model, fitted to this very record, gives
388.9 J/g. So parameters fitted in 1984 to other cements, carried through
CEMDATA18 thermodynamics and a Gibbs minimization, land closer at the endpoint
than a purpose-fitted empirical curve. What a calibration can earn here is the
*timing*.

**Three parameters, not six, and the number is measured.** Six candidates were
written down, one per mechanism with a distinct signature on the curve. The
singular values of `∂Q/∂log θ` on the coupled model come out
`[422, 101, 60, 6.3, 1.4, 0.20]` — a factor of nine between the third and the
fourth — so the measurement determines three *combinations*. The correlation
matrix says which: `k₁_C3S` with `n₁_C3S` at −0.985, and `k₃_C3S` with `n₃_C3S`
at −0.977, the latter pair being essentially the whole leading singular direction.
One per pair is all that can be fitted, and the rate constant is kept in each; the
third is `k₁_C3A`. Belite's `k₃_C2S` carries a weight of 0.001 in the three
directions the data see, which is a quantitative way of saying that eleven days of
heat cannot see belite.

Activation energies are not fitted at all, and the reason is not caution: all these
records are at 20 °C, an activation energy is a temperature sensitivity, and a fit
that reported one would be reporting a number the data cannot contain.

**Two failures worth having in the record.** A cheap stand-in for the coupled
model — same rate laws, hydrate assemblage written down by hand — was tried and
does not work, in two distinct ways. First, `C₃A + 3 Gp + 26 H₂O → ettringite` is
stoichiometrically impossible in this mix: 11 % C₃A demands about 1.1 mol of gypsum
per kilogram of binder and 4.6 % gypsum supplies 0.27, so a rate law reading only
C₃A drives the gypsum **negative** and the released heat to 1609 J/g against a
measured 376. Second, with the aluminate sent to hydrogarnet instead, the level is
right to some 9 % but a fit on the *normalized* curve saturates whatever parameter
bounds it is given, whichever subset is fitted — the imposed assemblage produces a
curve shape the Parrott–Killoh family cannot make. Both are now documented sections
of the page and assertions in the test suite, and together they are the
quantitative argument for the model `ionic_hydration.jl` already implements.

**A docstring claim that does not survive measurement.** `parrott_killoh_avrami`
states, after Parrott & Killoh, that "C₃S has no diffusion-controlled stage". The
sensitivity of the released heat to `k₂` for alite was expected to be exactly zero
and is not: the Jander term `k₂(1-ξ)^{2/3}/(1-(1-ξ)^{1/3})` falls as ξ grows, so
past a high degree of hydration it does become the minimum of the three branches,
and `k₂` moves the heat by up to about an eighth of what `k₁` does — all of it
after the first day. The published statement describes the 1984 fit's intent; this
implementation of it has a diffusion-limited tail. Pinned by a test so the claim
cannot drift back.

**Parameter-space automatic differentiation does not work, and the manual said it
did.** `docs/src/man/kinetics.md` claimed "the entire chain is
ForwardDiff-compatible" and showed a derivative through `integrate`. The ODE
right-hand side is AD-clean in the state and in time — deliberately, because
`Rodas5P` needs the time gradient — but not in the parameters: `build_u0` returns
a `Vector{Float64}`; `build_kinetics_params` casts the temperature, the initial
amounts, the stoichiometry and the calorimeter constants to `Float64`; the
surface-area constructors cast to `Float64` despite being declared `{T<:Real}`;
`system_enthalpy` accumulates into a `0.0`; and `respeciate!` is `Float64`-only by
construction. There was no test through `integrate`. The claim is corrected, the
blocking lines are named, and the calibration uses `AutoFiniteDiff()` with
`NelderMead` — still the SciML interface, `Optimization.jl`, and honest about why.

**Reusable, not general.** The helpers live in the script, not in `src/`, on
purpose: the interface should be settled by use before it is exported. The seam for
somebody else's data is `read_calorimetry`, which reads a documented CSV and takes
the experimental conditions from its comment header, and `CALIB_SPEC`, a list of
(phase, rate-law field, prior, box) rows.

### The dormant period is now on by default in `scripts/ionic_hydration.jl`

`run_ionic_hydration` and `ionic_reactions` take `induction` and
`induction_phases`, and the default is **on**: `τ = 5 h`, `m = 2.5`, applied to
the two silicates. A CEM I has a dormant period; Parrott–Killoh does not.

The numbers are round on purpose. The calibration returns 5.6 h and 3.56 on one
record and then shows `τ` correlated with `k₁_C3S` at 0.994, so the data barely
determine either separately. Carrying 20 228 s into an uncalibrated example would
lend a precision the measurement does not support. `induction = nothing` recovers
the previous behavior exactly.

What it changes is the early heat and little else. Measured over 28 days on the
example's own formulation, `Q(6 h)` falls from 74.1 to 40.2 J/g — the point of the
exercise — while the 28-day pore solution (pH 12.754) and porosity (0.3635 against
0.3634) do not move, and the integrator takes 216 accepted steps instead of 200.

Re-measuring both configurations caught an unrelated drift.
`IONIC_DEFAULT_SYSTEM`'s docstring claimed 202 steps and pH 12.58, where the model
*without* any dormant period now gives 200 and 12.754 — that figure had gone stale
for some earlier reason, and the docstring now says so rather than letting the new
default take the blame.

### Downstream: MeanFieldHomogenization.jl

Two things follow for MFH, and the second is not optional.

1. **`docs/Project.toml` pins `ChemistryLab = "0.11"` and must be widened to
   `"0.12"`.** Below 1.0 the resolver treats a minor bump as breaking whatever the
   API did, so the MFH documentation will not resolve against this release until
   that bound moves. MFH's own `Project.toml` does not depend on ChemistryLab at
   all — the coupling lives only in its docs environment.
2. **MFH keeps its own copy of the ionic setup**, in
   `scripts/common/ionic_hydration.jl`, with its own `IONIC_CALIBRATION` and its
   own `parrott_killoh_avrami` call. The dormant period has to be added there, and
   in `scripts/common/stoichiometric_hydration.jl` for the stoichiometric route — on the
   clinker silicates only, since the silica fume there already goes through
   `waller`, whose sigmoid carries its own onset delay.

   It will move a *reported* number, not just an input.
   `scripts/common/paste_micromechanics.jl` returns `E = 0` until the hydrate foam
   percolates — the setting transition as a genuine zero of the self-consistent
   fixed point — and `45_ionic_hydration_micromechanics.jl` reports a setting time
   in hours from the first non-zero modulus. Without a dormant period the model
   accumulates hydrates from `t = 0`, so percolation, and the setting time with it,
   arrives too early.

### `heat_release`'s `states` keyword never worked

`heat_release` accepts `states` so a caller who already holds the certified replay
does not pay for it twice, and its docstring says so: *"it is the expensive part,
and a caller that needs the compositions anyway should not pay for it twice"*. The
line that implemented it was

```julia
states = something(states, speciated_states(sol, kp; times = times))
```

and `something` is an ordinary function, so Julia evaluated **both** arguments and
the replay ran whether or not the states were supplied. Measured on the ordinary
Portland cement of `scripts/ionic_hydration.jl` over 28 days, 60 instants:

| call | before | after |
|--- |--- |--- |
| `heat_release(...; states = sts)` | 121 s | **0.017 s** |
| `heat_release(...)`, replay included | 118 s | 118 s |

`docs/src/examples/ionic_hydration.md` computes its replay once and hands it to
both `ionic_phase_history` and `heat_release` for exactly this reason; it has been
doing the replay twice. A timing assertion in `test/kinetics/test_postprocessing.jl`
now guards it, deliberately — a redundant computation whose result is discarded
cannot be detected any other way.

Note what this does **not** speed up: a calibration loop calling `heat_release`
without `states` already paid one replay and still does. The replay is the price of
correctness, and the same docstring records why — reading the running composition
instead gave 1174 J/g where the certified answer was a few hundred.

Found by benchmarking, not by reading. The hypothesis under test was that the
enthalpy sum dominated; it costs 0.029 s. An optimization written for it was
reverted, because fifty lines and a fallback branch to save 26 ms is not an
improvement, and the measurement that refuted it is what exposed the real defect.

### Proper names out of the identifiers

Author names belong in the bibliography, not in constants. In the two example
scripts and the page that includes one of them:

| was | is |
|--- |--- |
| `LAVERGNE_MIX_C100`, `_VESSEL_CP`, `_LOSS_A`, `_LOSS_B` | `CALORIMETRY_MIX_C100`, `_VESSEL_CP`, `_LOSS_A`, `_LOSS_B` |
| `lavergne_semiadiabatic` | `semiadiabatic_cell` |
| `PK_L2018_C3S` and siblings | `PK_SMOOTHED_C3S` and siblings |

Nothing in `src/` is affected and no exported name changes. `Lavergne2018` and
every citation in prose stay: that is where the attribution belongs.

The names of *models* stay too, and the distinction is deliberate.
`parrott_killoh`, `parrott_killoh_avrami` and `waller` are how the literature
designates those rate laws, in the same way it says Arrhenius, Avrami, Jander,
Powers, Langavant or Blaine — they are technical terms, not credits, and they are
part of the public API. What has been removed is the name attached to a *source of
parameter values*, which is a citation wearing an identifier's clothes.

### Also

  - `ionic_reactions` and `run_ionic_hydration` take a `pk_params` keyword that
    overrides the published Parrott & Killoh sets, so the calibration varies the
    rate laws of the existing model instead of restating it. The default sets are
    now a single `IONIC_PK84` constant rather than a dictionary literal repeated
    per call.
  - `scripts/opc_semiadiabatic_calorimetry.jl` cited Lavergne et al. (2018) with
    the DOI `10.1016/j.cemconres.2017.11.007`, which Crossref resolves to a
    different paper (Machner et al., *CCR* **105**, 1–17). Corrected to
    `10.1016/j.cemconres.2017.10.018`. The same defect was found and fixed in
    `docs/src/refs.bib` in v0.4.0; the script was missed then.
  - Nine bibliography entries added, every DOI resolved against Crossref and the
    Šmilauer entries also against the publisher's own landing page. The dataset and
    the article have different author lists — two creators and three respectively —
    and each is credited for its own artifact.

## v0.11.0 — a proved answer, whichever route found it

### Breaking changes

One new exported name, `solve_certified`. Nothing was removed or renamed and no
existing signature changed, but below 1.0 the resolver treats a minor bump as
breaking regardless, so a downstream `[compat] ChemistryLab = "0.10"` must be
widened to `"0.11"`.

### `solve_certified`: several routes, one proof

`solve_certified(des, starts; b)` solves from each starting composition in turn and
returns the first answer `optimality_certificate` **proves** optimal, or — if none
is proved — the one with the smallest KKT error together with its certificate, so
the caller always sees what it is getting.

Offering several routes is rigorous here rather than opportunistic, and the reason
is the certificate. For a convex problem it does not rank answers, it decides them:
the proof is the same proof whichever start produced the point, and nothing about it
depends on having predicted the winner. What would be a fudge is choosing a route by
taste and reporting its output unproved.

It is also necessary, because no back end dominates. Measured on an LC³ equilibrium
solved COLD at four degrees of reaction, with the dual Newton started from each
interior-point back end in turn:

| degree of reaction | from `OptimaOptimizer` | from `IpoptOptimizer` |
|---|---|---|
| 0.05 | certified, 2.7e-12 | certified, 9.1e-13 |
| 0.25 | **not certified**, 7.2 | certified, 4.2e-12 |
| 0.50 | **not certified** | certified, 4.6e-12 |
| 1.00 | certified, 1.8e-11 | **not certified**, 3.5e-3 |

Each solves a case the other misses. With both offered, every degree of reaction is
certified from a cold start — where before, an intermediate one could only be
reached by continuation from a nearby solution. The starts are supplied by the
caller, so this adds no dependency: whichever back ends are loaded are the ones
available.

### Two assertions corrected because they were wrong

Neither was in the way; both stated something untrue.

The negative control on water autoprotolysis asserted that forcing the Schur
complement back on gives `[H⁺]/[OH⁻] > 2`, pinning a direction of error. Its point
is that the null-space step is what carries the correct ratio, so the assertion is
now that the answer is NOT one — it used to come out at 3.78 and now at 8.7e-5, both
far from unity, and pinning the sign made the test fail for a change that did not
touch what it is about.

The `solve_certified` testset was written asserting `9.5 < pH < 10.3` for a system
that carries 0.01 mol of CO₂. That is the pH of calcite in PURE water, from a
different test; with the CO₂ the answer is 6.43, which is Reaktoro's own value for
this system (`H⁺ = 3.69e-7` in the reference table). The certificate was right and
the assertion was not.

### `speciated_states` walks up to the first instant it is asked for

Every replayed speciation warm-starts from the previous one — and the FIRST one
requested has no previous one, so it started from the cast composition, which
carries no active set at all. The chain is now walked up to it through a few
earlier times whose compositions are discarded and whose only purpose is to hand
the guess an active set.

The consequence was not cosmetic. On the reference OPC replayed at forty instants,
the first of them came back with **56 interior species where the answer has 25** —
every candidate hydrate present, four of them at 1e-5 to 1e-6 mol, the signature of
an interior-point iterate that never reached a vertex — and neither the certifying
Newton nor its continuation recovers from that, because both inherit the start.
That instant is now certified, with stationarity 9.1e-13 and element balance
6.1e-15 mol: the whole replay is proved optimal where before one instant in forty
fell back silently to an uncertified composition. The 28-day heat of that paste
moves from 420.2 to 420.3 J/g, which is the size of the error the fallback was
hiding.

Two other routes to the same instant were implemented and measured, and neither is
in the code: retrying it from each of the thirty-nine certified neighbors, nearest
first, left it unproved with every one of them; and the continuation between
certified instants was rewritten to step forward adaptively rather than bisect —
correct in itself, since bisection lowers the upper end on failure and so abandons
the target for good, but it does not reach this instant either. What was missing was
the start, and only the run-up supplies it.

### The full Portland cement, through its pore solution

New example page and script, `scripts/ionic_hydration.jl`: a complete CEM I —
alite, belite, aluminate, ferrite, gypsum, limestone filler — dissolving into ions,
with the whole hydrate assemblage decided by Gibbs minimization at every accepted
step. The aluminate cascade that a stoichiometric model has to encode by hand comes
out as a result: run with 3.5 % limestone the ettringite survives to 28 days, run
without it the ettringite is depleted to zero, and no line of code distinguishes
the two cases.

The page carries the calorimetry with it, isothermal and semi-adiabatic
(NF EN 196-9), with every number of the cell written out: 420 J/g against 405 J/g
at 28 days, a rise of 19 K against an adiabatic 75 K. The heat is taken from the
enthalpy of the certified compositions, not from the kinetic reactions — those only
dissolve the clinker, so the sum `Σ rᵢ(−Δ_r H⁰ᵢ)` cannot see the precipitation heat
at all, and driving the cell from it gave a rise of 207 K.

## v0.10.0 — the heat of hydration, and solid solutions the solver can hold

### Breaking changes

Five new exported names — `enthalpy`, `heat_capacity`, `missing_enthalpy`,
`system_enthalpy` and `heat_release`. Nothing was removed or renamed and no
existing signature changed, but below 1.0 the resolver treats a minor bump as
breaking regardless, so a downstream `[compat] ChemistryLab = "0.9"` must be
widened to `"0.10"`.

`[compat] OptimaSolver` moves to `"0.4"`, which is required: the solid-solution
work below depends on `SolutionPhase(...; mole_fraction = true)`, introduced
there.

### Calorimetry that a coupled model can actually use

`enthalpy(state)` and `heat_capacity(state)` sum `nᵢ ΔₐH⁰ᵢ(T,P)` and
`nᵢ Cp⁰ᵢ(T,P)` over the composition, and `heat_release(sol, kp; times)` returns
the cumulative heat and the heat rate along a trajectory. This is Eqs. (17)–(21)
of Lavergne et al. (2018): enthalpy is a state function, so its drop between two
states at the same temperature is the heat given off, with reactants, ions and
hydrates each counted once and no reaction stoichiometry to write down.

That last point is what makes it necessary. `heat_rate` sums `rᵢ(−ΔᵣH⁰ᵢ)` over
the KINETIC reactions, which is right when those reactions produce the hydrates —
and wrong under partial equilibrium, where they only dissolve the anhydrous
phases into ions and the hydrates are precipitated by the Gibbs minimization.
Driving a semi-adiabatic cell from it put an ordinary Portland cement at a
temperature rise of 207 K. That combination now warns instead of returning the
number in silence.

`heat_release` reads the **certified** speciations of `speciated_states`, not the
composition the integrator carries. The in-run minimization is warm-started and
uncertified, and a single hydrate is worth hundreds of kilojoules: read that way
the curve came out at 12.7, 145, 1174, 936 and 631 J/g at 1 h, 6 h, 12 h, 1 d and
2 d — heat that rises and then falls. Certified, the same cement gives 185, 288,
347 and 420 J/g at 1, 3, 7 and 28 days, and the curve is monotone.

`missing_enthalpy(state)` lists the species that carry no `ΔₐH⁰` and are
therefore absent from the balance, because a heat curve missing one hydrate is
not visibly wrong.

### Every species now matches Reaktoro, `CaOH+` included

`equilibrium_reference.jl` recorded `CaOH+` as `@test_broken`: out by a factor 2.5
at 4e-9 mol, and left alone because Ipopt landed on the same point, which made it
look like a property of the problem rather than of one solver. It was neither. The
interior-point stage was reporting points it had reached without ever satisfying
the KKT conditions at `μ = 0` — see `OptimaSolver` v0.4.0 — and the trace species
are exactly where that shows. The assertion is now a plain `@test`.

### `s[:Cp⁰]` returned zero on a species that had one

The thermodynamic functions are built on demand. `getproperty` knew that;
`getindex` did not, and returned its not-found value `0` — an `Int64` — for any
of `:Cp⁰`, `:ΔₐH⁰`, `:S⁰`, `:ΔₐG⁰`, `:V⁰` on a species whose functions had not
been forced yet. Callers found out only when they tried to evaluate it. Returning
0 is right for a missing ATOM, which is what that fallback is for; it was never
right for a property the species can produce.

### A solid solution now has a reference worth the name

`DualEquilibriumSolver` took the first end-member of every solid solution as the
phase reference. The reference is the member whose stationarity the outer system
carries rather than inverting, so it has to be the one the phase is mostly made
of — for the aqueous phase that is the solvent, and not by convention. Picking a
minor end-member sends the outer unknown `ln x_ref` towards −∞ and the Jacobian
with it. It is now chosen by magnitude from the caller's own composition, and
solid solutions are declared to `OptimaSolver` as mole-fraction phases, which is
what lets them be solved at all.

## v0.9.0 — equilibria that are proved, not hoped for

### Breaking changes

Two new exported names, `DualEquilibriumSolver` and `optimality_certificate`.
Nothing was removed or renamed and no existing signature changed, but below 1.0
the resolver treats a minor bump as breaking regardless, so a downstream
`[compat] ChemistryLab = "0.8"` must be widened.

### The interior-point solver could not say whether its answer was the answer

`EquilibriumSolver` minimizes `G` by an interior-point method and, on a cement,
it stops on `MaxIters` — at any tolerance. That was known. What was not known is
how much it costs, because nothing in the package could check.

The Gibbs problem is **convex**: an ideal mixing entropy, which is convex, plus
terms that are LINEAR in the amounts of the pure phases, whose potential does not
depend on how much of them there is, over the polyhedron `{A n = b, n ≥ 0}`. The
minimizer is therefore unique and every stationary point is global — so a solver
returning different answers from different starting points is not finding local
minima, it is stopping short of stationarity, and the KKT conditions are
*sufficient*: they can be checked, and checking them is a proof.

**`optimality_certificate`** does that check, on any composition and whatever
produced it. It reports the stationarity of the interior species, the component
balance, and the worst supersaturation among absent phases. A species below the
floor is at its bound, where the condition is the inequality `μ + Aᵀy ≥ 0` and
not the equality — imposing the equality on an amount held at `1e-16` whose
mass-action value is `e⁻³⁰⁰` misstates its log-activity by 263 RT units.

### `DualEquilibriumSolver`

Newton on the KKT system in element-potential space — the Brinkley–Karpov
formulation of the geochemical Gibbs-minimization codes. Writing `u = −Aᵀy`, an
aqueous species obeys `aᵢ = exp(uᵢ − gᵢ)` and a pure phase is present exactly
when `uᵢ = gᵢ` and absent when undersaturated, which is the classical
phase-stability criterion.

**The algorithm is not in this package.** It is
`OptimaSolver.dual_newton_solve`, where it belongs: the active set, the
degeneracy test on the rows of `A`, the two-level Newton and the certificate are
statements about a convex program, not about chemistry. What this package
supplies is the four things that *are* chemistry — the conservation matrix and
the reference potentials `Δ_aG⁰/RT`; the activity model as the callback `h` with
`∇f = g + h(n)`; which species are strictly positive at any solution and which
may vanish; and which species is the solvent, whose activity is a mole fraction
and therefore bounded above, so that its stationarity cannot be inverted. This
release requires OptimaSolver 0.3.

Two levels. The inner one inverts the **solutes'** mass-action laws at fixed
potentials and fixed solvent amount; the outer is a Newton on `1 + m + |P|`
unknowns — the solvent, the `m` element potentials, the amounts of the active
phases — some fourteen numbers for a cement, against forty-seven species in the
interior-point route. Parameterizing the solutes by `ln n` makes their positivity
automatic, so the fraction-to-boundary rule that capped the interior-point step
at *every* iteration has nothing left to act on.

Six things had to be right, and each was found the hard way:

  - **the solvent cannot be inverted through its own mass-action law.** Its
    activity is a mole fraction, so `ln a_w ≤ 0` always, and an arbitrary `y` can
    demand `ln a_w > 0`, for which no finite composition exists. It belongs to the
    outer system, where the balance determines it;
  - **two phases declared saturated over-determine `y`.** Their stationarity rows
    are then jointly infeasible. With portlandite wrongly admitted at 1e-6 mol
    from the starting guess, the residual sat at 29.5 and fell by 1e-3 an
    iteration; dropping it *during* the Newton rather than after, the same solve
    reaches 5e-12 in eleven steps;
  - **an amount must not be clamped to `ϵ` on the way out**, for the reason given
    above;
  - **the active set must be updated during the Newton, not after it**, and one
    phase admitted at a time. Two phases both declared saturated over-determine
    `y` and their stationarity rows are jointly infeasible; admitting a batch
    feeds a cycle in which a phase is admitted, driven negative, dropped, and
    readmitted. On a cement without limestone that produced a solve converged to
    2e-12 — of the wrong subproblem, with an absent phase supersaturated by 10.9.
    Visited active sets are recorded, which makes the loop finite in practice as
    well as in theory;
  - **a component whose total has vanished must be removed from the unknowns.**
    Its element potential is determined by nothing, the Jacobian is singular in
    that direction, and this is the ordinary state early in a hydration run, where
    the iron, sulfur and aluminum totals are at machine zero. The test is not
    `bₖ ≈ 0` but `bₖ ≈ 0` **with the non-zero entries of row `k` sharing a sign** —
    only then does `Σ Aₖᵢnᵢ = 0` with `n ≥ 0` force each term to vanish. The `H⁺`
    row carries `+1` for `H⁺` and `−1` for `OH⁻`, so its zero total is the
    ordinary state of pure water; treating it as degenerate kills the entire
    acid–base system and returns pH 7.000 with the calcite undissolved;
  - **the replay must be seeded and, if need be, continued.** `speciated_states`
    tries a ladder of early instants until one certifies, then walks forward; when
    an instant still resists it bisects the interval in `bₑ` from the last
    certified one, which is a homotopy in the component totals and terminates.
    Without the ladder the first requested instant of a cement could not be
    proved — neither the interior-point answer nor the cast composition puts the
    Newton close enough.

### Measured

On the package's own Reaktoro reference — calcite, CO₂ and water, ideal
activities — the certified answer matches **every species to 1 %**, including
`CaOH+`, which `equilibrium_reference.jl` records as `@test_broken` because the
interior-point answer is 147 % high. On calcite in pure water the certified pH is
**9.90**; the interior-point answer is 6.96, which is not an imprecision but a
wrong answer, and nothing in that solver's output reveals it.

`speciated_states` now certifies each instant it replays, and **names any it
cannot**: falling back silently would leave the caller unable to tell a certified
trajectory from an uncertified one, and those are not the same object. On a full
ordinary Portland cement over 28 days, with and without limestone, **all forty
replayed instants of both mixes are certified**, with element balances between
1e-11 and 1e-13 mol — where the six-hour instant, the AFt peak at which the
assemblage switches, stood at 6.9e-2 mol before this work.

### What is not certifiable yet

**Solid solutions.** `DualEquilibriumSolver` refuses a system that declares one,
and does so explicitly rather than returning a number. A pure phase has unit
activity, hence `gᵢ = uᵢ` and a bound-constrained variable with an active set. An
end-member of a solid solution has `μᵢ = gᵢ + ln aᵢ(n)`, whose activity goes to
`−∞` as its mole fraction goes to zero: it is never exactly absent while the
phase exists, so the active set belongs at the level of the PHASE and not of the
species. That is a different algorithm. Treating an end-member as a pure phase
would not fail loudly — it would return a composition that looks reasonable and
is not the minimum — which is why the refusal is explicit and
[`speciated_states`](@ref) reports such instants as uncertified.

This matters for what comes next: slag and calcined-clay binders form C-A-S-H,
hydrotalcite and AFm solid solutions, and certifying those needs the phase-level
active set.

**The in-run speciation.** `respeciate!` still uses the interior-point solver, so
the compositions *inside* the integration are not certified; only the replay is.
Measured, that makes no difference here, because `bₑ` is integrated from the
rates alone and a Parrott–Killoh or Waller law reads only its own degree of
reaction. A rate law reading log-activities — a saturation ratio, or a
pH-dependent dissolution law — would feed the speciation back into the
trajectory, and would need this closed first.

### Still true

The interior-point solver is unchanged and still does not report convergence on a
cement. It is what gets the certifying Newton into the neighborhood — on 74 of
80 measured cement equilibria the certified answer is reached from its guess —
and that division of labor is deliberate.

## v0.8.2 — a replay that conserves matter, solid solutions that reach the solve

Bug fixes. No API changed and nothing was removed, so a downstream
`[compat] ChemistryLab = "0.8"` accepts this release unchanged.

The `docs/src/examples/coupled_hydration.md` page is rewritten on the real API:
it described four calls that never existed and had no executed block, so nothing
had ever checked it. Every block now runs when the documentation is built.

The `OptimaSolver` bound is relaxed from `"0.2.7"` to `"0.2"`. The patch-level
lower bound guarded against 0.2.6, where `OptimaOptimizer` discarded the
caller's `u0`; a fresh resolution always takes the newest patch, so the bound
only ever mattered when reviving an old manifest.

### The relative residual saturated on near-empty rows

`_row_residual` scaled each row by the larger of its own element total and the
matter flowing through it. For a row whose total is a millionth of the largest —
an element barely present, or the charge row early on — that divisor is
essentially zero, and a rounding error came out as an alarm. A full ordinary
Portland cement balanced to **1.4e-10 mol** at 28 days was reported at `3.2e-2`,
and the near-empty rows of the first steps saturated the measure at `1.0`, which
read as a 100 % violation and was nothing of the sort. Every row is now floored
at a millionth of the system scale, and the same run reports `3.9e-6`.

### A replayed speciation could stop on a broken balance

Two guards were added to [`speciated_states`](@ref), both judged on conservation
of matter rather than on a retcode:

- **A retry from a guess carrying no active set.** The warm start is the previous
  speciation, and it brings the previous *active set* with it. Where the
  assemblage switches — an ordinary Portland cement around six hours, when the
  ettringite peaks and the sulfate starts moving into AFm — that set is the wrong
  one and an interior-point method started inside it does not cross over. When
  the balance exceeds `1e-6 mol` the instant is solved again from the cast
  composition, and whichever result closes better is kept.

- **A feasibility restoration that actually converges.** Alternating projection
  converges linearly, at a rate set by the angle between the box and the affine
  set, and on a cement that angle can be small: at that same six-hour instant,
  where the iron row carries 0.013 mol across thirteen species, 200 sweeps left
  a residual of 8.4e-1 and 20 000 were needed to reach 6.7e-9. A replay runs a
  handful of times and now asks for the budget it needs.

  Inside the ODE right-hand side the budget stays where it was: raising it there
  made the run *worse* — the worst in-run balance went from 1.1 mol to 8.5 at
  2000 sweeps and to 41 at 100 000 — which is the same unpredictability the
  back-end shows elsewhere and is not understood.

Measured on a full OPC replay, with OptimaSolver 0.3: the six-hour instant —
the AFt peak, where the assemblage switches — goes from **6.9e-2 mol to
7.1e-15**, and every other instant lands between 1.8e-15 and 3.6e-9. Over the
201 accepted steps of the coupled run the worst balance is 6.0e-6 mol and the
median 2.7e-9. On the 40-instant grid the chapters use, the porosity and the
elastic modulus are both monotone across every interval and the pore solution
holds at pH 12.52–12.59.

### Solid solutions were dropped from the equilibrium partition

`_equilibrium_subsystem` rebuilds a `ChemicalSystem` for the partition the
equilibrium is actually solved on, and it did not carry `solid_solutions` over.
The loss was silent and total: a coupled run could declare CSHQ, AFm or
Hydrogarnet and the solve would treat their end-members as separate pure phases,
the mixing entropy never entering the Gibbs energy. Measured on alite and belite
with the four CSHQ end-members, the run was bit-identical with and without the
declaration. With the fix the two differ — the element balance goes from 1.5e-12
to 3.6e-15 and the pore solution from pH 12.48 to 12.37.

A solution survives into the partition only if **all** its end-members are in it;
one split across the kinetic and equilibrium sides is not a phase the equilibrium
can mix, and is dropped rather than passed truncated.

This makes the code path work. Whether a given solid solution is *stable* in a
given system is a separate, thermodynamic question: with these Cemdata18
end-members the small silicate system above still precipitates no C-S-H.

### The warning could not tell the trajectory from a probe

The right-hand side is evaluated far more often than the solution advances —
Jacobian finite differences and rejected steps included — and a poor speciation
on a *perturbed* `bₑ` never enters the result. Reporting the worst over all
evaluations alarmed about something that does not affect the answer.

`integrate` now reports the worst on the **accepted steps** first, and the
all-evaluations figure second. It also states when the distinction matters:
`bₑ` is integrated from the rates alone, so a rate law reading only its own
degree of reaction (Parrott–Killoh, Waller) gives a trajectory independent of the
speciation, while a law reading log-activities feeds it back in.

### The warning now gives moles

`integrate` reported only the relative figure, which cannot be judged without
knowing what it was divided by. It now leads with the absolute worst in moles —
`1e-10 mol` is machine precision whatever the system, `1e-2 mol` against a
0.3 mol sulfate budget is not — and says where such cases come from: the first
steps, where the paste has barely reacted. On the OPC above the honest reading
is 1.05 mol at the worst early step and machine precision from three days on.

## v0.8.0 — a public replay of the equilibrium partition

### Breaking changes

- **New exported name: `speciated_states`.** Below 1.0 Julia's resolver treats a
  minor bump as breaking whatever the API did, so a downstream
  `[compat] ChemistryLab = "0.7"` will **not** accept 0.8 and must be widened.
  Nothing was removed or renamed, and no existing behavior changed, so widening
  the bound is the only adjustment required.

### The composition at a given time was not recoverable

`state_at` returns the purely kinetic reconstruction `n(0) + νᵀξ` and says so:
the redistribution performed by the equilibrium solve is not recoverable from
the stoichiometry. For a cement that reconstruction is meaningless — every
dissolved element in solution, not one hydrate. The composition left in the
solver's buffers is no better: it is rewritten at every right-hand-side
evaluation, Jacobian differences and rejected steps included, so it is not the
accepted composition at any time. Reading it is what made a full OPC look as
though it held 0.244 mol of ettringite against 0.267 mol of sulfate.

**`speciated_states(sol, kp; times)`** does the recovery properly: the kinetic
species from the ODE state, and the equilibrium partition re-solved from the
element totals the run carried, walking the instants in order.

Two things it must do, both found the hard way:

- **cap and restore the guess.** The warm start is the equilibrium of the
  *previous* `bₑ`, so once an element has been spent it demands more of it than
  now exists and the interior-point solve begins outside its own feasible set.
- **use a back-end instance the integration has not touched.** Replaying through
  `p.eq_solver.solver` — the object that drove thousands of solves during the
  run — returned pH 14.2 with 0.31 mol of ettringite and no AFm, where a clean
  instance of the same type and settings gives pH 12.58 with the sulfate
  entirely in AFm, on the same trajectory, the same `bₑ` and the same guess.
  **The back-end carries state across solves**; this is a defect upstream, and
  until it is fixed a replay must not inherit it.

The first requested instant is the loose one, having only the cast composition
to start from; every instant after it lands at machine precision. Ask for a
first point within the first hours.

### The back-end defect is fixed upstream

The state the replay had to avoid was `OptimaOptimizer._cache`: with
`warm_start = true` the algorithm object started every solve from its previous
solution, discarding an explicit `u0` — including the guess this package builds
in `respeciate!`, so the caps and projections above were partly defeated during
the run itself. **OptimaSolver 0.2.7** consults the cache only when the caller
supplies no interior point, and `[compat]` now requires it. `speciated_states`
still builds a clean back-end instance, which costs nothing and keeps the replay
correct against an older back-end.

## v0.7.1 — the element balance of a re-speciation, measured and enforced

Bug fixes only. No API was removed or renamed, and no documented behavior
changed, so a downstream `[compat] ChemistryLab = "0.7"` accepts this release
unchanged.

### The reported element-balance residual could not see a violated element

`_respeciate_solve!` scaled the whole residual `|Aₑnₑ − bₑ|` by
`maximum(abs, bₑ)` — in a cement paste, by the water budget. Water is 34 mol
against 0.27 mol of sulfate, so a violation of 0.465 mol of sulfate, i.e. 174 %
of the sulfate present, was reported as `1.4e-2` and read as a converged solve.

The residual is now taken row by row and scaled by each element's own budget,
or by the matter flowing through that row when the budget is zero — as it is for
the charge row, where dividing by `bᵢ` alone turned a rounding error into a
reported residual of 2·10³.

This is a diagnostic change, but not a cosmetic one: every convergence study
run against the old number was measuring the wrong quantity.

### A non-converged speciation was handed on as the next warm start

`eq_warm` was set after any solve that did not throw, non-converged ones
included, so a single bad point seeded every step after it and the error was
locked in for the rest of the run. A speciation is now handed on only when its
per-element residual is within `EQ_RESIDUAL_TOL`; otherwise the step falls back
to the cold stoichiometric reconstruction.

### Re-speciation could start outside its own feasible set

The Gibbs minimization is posed with `Aₑn = bₑ` as a hard equality, and the
warm start is the equilibrium of the *previous* `bₑ`. Once an element has been
spent — the sulfate of an ordinary Portland cement, after the gypsum is gone —
that guess demands more of it than now exists, and the interior-point solve
starts infeasible.

Two guards now run before the solve, both of which leave a feasible guess
untouched:

- `_budget_clip!` caps each species at what the totals can supply,
  `nⱼ ≤ minᵢ bᵢ/Aᵢⱼ`, over the rows it consumes. Rows with a negative total —
  H⁺ in a hydrating cement, met by the hydroxides — bound nothing.
- `_restore_feasibility!` then lands the guess in `{Aₑn = bₑ, n ≥ 0}` by
  alternating projection, ending on the positivity clamp so no small negative
  amount survives for the barrier to reject.

On an isolated ordinary-Portland-cement equilibrium — alite, belite, aluminate
and gypsum dissolved, with ettringite, monosulphate, katoite, portlandite and
C-S-H free to precipitate — this takes the achieved element balance from 174 %
to `2e-15`, machine precision, and returns the textbook assemblage: AFm 0.2323,
ettringite 0.0116, C-S-H 1.8556, portlandite 2.3021, at pH 12.58. Both budgets
close on the last digit — sulfate `0.2323 + 3x0.0116 = 0.2672` and aluminum
`2x(0.2323 + 0.0116) = 0.4879` — against the 0.2672 and 0.4879 available.

### The cold start was supersaturated in every phase at once

The cold-start guess added the stoichiometric reconstruction `Ne^T xi`, which
places every dissolved element in solution with no hydrates at all. For an
aqueous-only system that is harmless; for a cement it is close to the worst
possible start, supersaturated in every phase simultaneously, with an H+ entry
so negative that it is clamped to the floor and loses the acidity the hydroxides
must balance. The guess is now the composition the specimen was cast with, which
`_restore_feasibility!` then carries onto the current `be`.

### What this unblocks

A full ordinary Portland cement now runs end to end: alite, belite, aluminate,
ferrite and gypsum, over 28 days, 202 accepted steps, `retcode = Success`. The
pore solution holds at pH 12.58, and the aluminate sequence comes out of the
thermodynamics with no sequencing rule written anywhere -- ettringite forms
early, peaks at 6 hours and converts to monosulphate once the sulfate is spent,
the AFm settling at 0.26719 mol, exactly the sulfate budget at one SO4 per
formula.

Note that `p.n_full` is a scratch buffer rewritten at every right-hand-side
evaluation, Jacobian differences and rejected steps included; it is not the
accepted composition at `t_end`. Recover a speciation from a solution by taking
`be` from the ODE state and re-equilibrating, as `state_at` already documents.

## v0.7.0 — the real porosity of a setting binder, and a warm-started coupling

### Breaking

- Being a minor bump below 1.0, a downstream `[compat] ChemistryLab = "0.6"`
  will not accept 0.7 and must be widened. No API was removed and no documented
  behavior changed; the additions below are new methods.

### The porosity of a cement was not available

`porosity(state)` returns `(V_liquid + V_gas) / V_total`, and both ends of that
ratio are wrong for a hydrating binder: the denominator is the *current* volume,
which shrinks as the reactions proceed, while a sealed specimen keeps the volume
it was cast with; and the numerator has no gas term, so the empty porosity left
by the chemical shrinkage is structurally invisible. The errors compound — on a
w/c = 0.5 paste at 28 days the method returns 0.327 where the porosity referred
to the specimen is 0.375, the volume having shrunk 7.2 %.

The method is unchanged, and correct for a fixed-volume aqueous system, but now
carries that warning. Two new methods give the right calculation:

- **`porosity(state, reference)`** → `(; liquid, void, total)`, referred to the
  fresh material and counting the Le Chatelier contraction as empty porosity.
- **`saturation(state, reference)`** — a sealed paste desaturates as it hydrates
  though no water ever leaves it.

`porosity(state, ref).void` is exactly the `"void"` entry of `volume_fractions`
called with the same reference, so a transport or micromechanical model fed by
either sees the same thing.

### Fixed — the equilibrium sub-solve was restarted cold at every step

`respeciate!` rebuilt its starting guess from the reaction extents each time,
placing every dissolved element in solution with zero hydrates — a wildly
supersaturated composition. Harmless for an aqueous-only system, nearly the worst
possible start for a cement. It now warm-starts from the speciation left by the
previous accepted step, which is an equilibrium for a nearby `bₑ`.

## v0.6.0 — the coupled path made trustworthy

The v0.5.0 release opened the door to coupled dissolution/precipitation modeling.
Building a cement model on it exposed five defects, four of which were silent.

### Breaking

- **`EquilibriumSolver` gained a `model` field**, so its signature is now
  `EquilibriumSolver{F, S, V, M}`. Code constructing the raw struct positionally
  must be updated; the documented `(cs, model, solver)` constructor is unchanged.
  A new accessor `activity_model(::EquilibriumSolver)` returns it.
- Being a minor bump below 1.0, a downstream `[compat] ChemistryLab = "0.5"`
  will not accept 0.6 and must be widened.

### Fixed — the activity model was silently discarded in a coupled run

A coupled run must rebuild the equilibrium solver for the equilibrium
*sub-system*, because the compiled potential does not carry over. Having no
record of the model it was built from, it rebuilt with `kp.activity_model` —
`DiluteSolutionModel()` by default. Asking for `HKFActivityModel` on a cement
pore solution at I ≈ 0.5 mol/kg therefore gave an infinitely dilute solve, with
no warning. The solver now remembers its model, the rebuild uses it, and a
disagreement with the problem's own model is reported.

### Fixed — the equilibrium sub-solve started on its own bound

`respeciate!` floored its starting guess at `p.ϵ` (1e-30), which
`EquilibriumProblem` raises to exactly the `1e-16` lower bound. An interior-point
method started on its bound stalls: the package's own Reaktoro reference case
reported six non-converged solves out of eight steps for this reason alone. The
guess is now floored strictly inside the box.

**Do not "fix" this class of stall by loosening the optimizer tolerance.**
Measured against Reaktoro 2.13 on that same case, the worst species error grows
from 4.3 % at the default `tol = 1e-10` to 13 % at `1e-9`, 38 % at `1e-8` and
252 % at `1e-7`. A green retcode bought that way costs a factor of sixty in
accuracy.

### Fixed — a bad speciation was invisible

A non-success retcode was a `@warn` at `maxlog = 1` whose result was used anyway,
and it never reached the failure count reported by `integrate`, which only saw
solves that *threw*. Over thousands of steps that is one warning for any number
of bad speciations.

- `NONCONVERGED` counts them, and `integrate` reports the total.
- More usefully, `integrate` also reports the worst `|Aₑn − bₑ|∞` over the run.
  That is the criterion with physical meaning here: a solve can stop short of the
  optimizer's tolerance and still satisfy the element balance to machine
  precision, which is exactly what the remaining stalls do (9e-12 on the
  reference case).

### Fixed — `integrate(kp; reltol = …)` threw a `MethodError`

The documented shortcut forwarded its keywords to a concrete method that accepted
none. Precedence is now explicit: call site beats the solver's own settings,
which beat the defaults.

### Fixed — two solid solutions were dead entries

`data/solid_solutions.toml` declared AFm on `Ms`/`Mc` and Hydrotalcite on
`Ht_OH`/`Ht_CO3`. None of those four symbols exists in either shipped database,
so `build_solid_solutions` skipped both with a warning. AFm being the only
Redlich-Kister entry, the non-ideal mixing path had no live case at all. They are
now `monosulphate12`/`monocarbonate` and, for hydrotalcite, the OH/CO₃ couple at
matching Mg:Al = 2 (`hydrotalcite`/`Mg2AlC0.5OH`). All five phases build.

### Documentation

- `saturation_ratio` stated `ln Ω = Σνᵢlnaᵢ + ln K`, contradicting both the next
  line and the code. The sign is corrected.
- `examples/cement_carbonation.md` claimed `cemdata18-thermofun.json` does not
  contain calcite. It does; the merged database is needed for the phase volumes.
- `man/kinetics.md` used `cs.dict_reactions` as though it were a catalog of the
  database, when it holds only the declared kinetic reactions.

## v0.5.0 — the bridge from a kinetics run to a microstructure

### Breaking

- **The ODE state vector gained the extents of reaction.** Its layout is now
  `[bₑ, nₖ, ξ, (calorimeter)]`, where `ξ` integrates `dξ/dt = r`, one entry per
  kinetic reaction. Code indexing `sol.u` positionally for anything other than
  `bₑ` or `nₖ` must be updated; the calorimeter's slot is addressed from the end
  of the vector and is unaffected, as are `temperature_profile`,
  `cumulative_heat` and `heat_flow`.
- Being a minor bump below 1.0, a downstream `[compat] ChemistryLab = "0.4"`
  will not accept 0.5 and must be widened regardless.

### A rate law may now depend on any species

This is what the state change buys, and it is the point of the release.

Previously only the kinetic species evolved inside the ODE residual: every other
amount was pinned to its initial value for the whole run, because the buffer
holding them was refreshed only by `respeciate!`, which runs only when an
equilibrium solver is present. A rate closure gating on a consumed reactant
therefore never saw it move, and the failure was silent — the kinetic mass
balance stayed exact and the solver reported success while the gated species went
negative. In the cement model that motivated this release, 0.44 mol of ettringite
formed out of 0.27 mol of gypsum.

The residual now reconstructs every species from `n = n(0) + νᵀ ξ` before
evaluating the rates, so reaction sequencing — sulfate available for ettringite,
portlandite available for a pozzolanic reaction, product inhibition — is written
directly. With an equilibrium solver the equilibrium partition is still owned by
the equilibrium solve and refreshed once per accepted step.

As a side effect, [`reaction_extents`](@ref) and [`state_at`](@ref) read the
extents from the solution instead of re-integrating the rates by quadrature: both
are now exact to the solver's tolerance, and the `nsub` keyword is gone.

### From moles to volume fractions

The pieces to turn a hydration run into the input of a mean-field homogenization
scheme were all present — `V⁰` from CEMDATA18, `volume(state, species)`,
`porosity` — but nothing joined them, and the sealed-volume balance was missing
entirely.

- **`volume_fractions(state)`** and **`volume_fractions(state, groups)`**, the
  latter aggregating species into the phase families a homogenization scheme
  consumes (`"C-S-H"`, `"AFt"`, `"anhydrous"`, …). A species listed in two groups
  is an error; species in no group are collected rather than dropped.
- **`chemical_shrinkage(state, state₀)`** and the `"void"` phase. Passing
  `reference` to `volume_fractions` selects the sealed-curing convention: the
  fractions are referred to the initial volume, held fixed, and the deficit left
  by the reactions becomes an explicit gas-filled phase. Without it the fractions
  sum to less than one and the microstructure is wrong.
- **`missing_molar_volumes(state)`** — a species with no `V⁰` is silently absent
  from every volume computation in the package; this makes that visible.

Aqueous solutes have negative standard partial molar volumes, so an individual
fraction can be negative. Those contributions are kept, which is what makes
`volume_fractions` and `volume` agree exactly.

### Post-processing a kinetics solution

`n_full` is a buffer mutated in place — after a run it holds the last accepted
step and nothing else — so there was no way to recover the composition at an
arbitrary time.

- **`reaction_extents(sol, kp)`** — the extent of each reaction, read from the
  ODE state and therefore exact at any instant.
- **`extent_residual`** — the drift between `nₖ` and `νₖᵀ ξ`, which are redundant
  by construction, so the integrator's own consistency is measurable.
- **`state_at(sol, kp, t)`** — the full `ChemicalState` at any instant, every
  species rebuilt from `n = n₀ + νᵀξ` so the result conserves exactly what the
  reactions conserve.
- **`degrees_of_hydration`** and **`mean_degree_of_hydration`**, replacing the
  `phase_alpha` closure copy-pasted into three shipped scripts.

### Parrott & Killoh, canonical formulation

The shipped `parrott_killoh` implements a smoothed variant whose parameters are
not those of the 1984 paper — its diffusion branch uses the same `K₃` and `N₃`
for all four clinker phases. It is unchanged, and now documented as one of two
variants.

- **`parrott_killoh_avrami`** with **`PK84_PARAMS_C3S/C2S/C3A/C4AF`** — the
  canonical Avrami / Jander / power-law form, `α̇ = min(α̇₁, α̇₂, α̇₃)`. With these
  parameters C₂S has no nucleation–growth stage and C₃S no diffusion-controlled
  stage, which is a sharp check on a transcription.
- **`waller`** with **`WALLER_PARAMS_FLY_ASH/SILICA_FUME/SLAG`** — supplementary
  cementitious materials do not follow Parrott & Killoh; `blended_cement_kinetics.jl`
  had to invent PK parameters for slag and metakaolin for want of this.
- **`blaine_factor`**, **`humidity_factor`**, **`powers_alpha_max`** — the three
  rate corrections, previously either absent or retyped inline in every script.

The Avrami branch vanishes at `α = 0`, so `α ≡ 0` solves the ODE and hydration
never starts. Parrott & Killoh's own discrete scheme escapes this by integrating
over the first time step; a continuous solver cannot, so the argument is floored
at `PK_AVRAMI_SEED`.

### Fixed

- **A time-dependent rate law broke every Rosenbrock solver.** The ODE right-hand
  side typed its rate vector from `eltype(u)` alone, but Rosenbrock methods
  (`Rodas5P`, `Rodas4`, `Rosenbrock23`, …) need a *time* gradient, which they take
  by calling the residual with a dual `t` and a plain `u`. Any rate depending on
  `t` then failed with "First call to automatic differentiation for time gradient
  failed". `parrott_killoh` ignores `t`, so nothing exposed it until `waller`. The
  type is now promoted with `typeof(t)`.
- **Non-kinetic amounts were frozen inside the residual** — see the section above,
  which is the substance of this release.

### Bibliography

Every entry of `docs/src/refs.bib` was checked against Crossref (six DOIs) or,
for the four entries that have no DOI, against the publisher or an authoritative
catalog record. Three defects were corrected:

- The DOI of `Lavergne2018` pointed at `10.1016/j.cemconres.2017.11.007`, which
  Crossref resolves to a **different** paper (Machner et al., *CCR* **105**, 1–17).
  Corrected to `10.1016/j.cemconres.2017.10.018`.
- `Lothenbach2015` held the Cemdata18 paper, published in **2019**, while
  `data/solid_solutions.toml` uses the tag `Lothenbach2015` for the genuinely
  different Lothenbach & Nonat (2015) paper and `Lothenbach2019` for Cemdata18.
  The entry is now keyed `Lothenbach2019`, and the real
  Lothenbach & Nonat (2015), *CCR* **78**, 57–70, is added.
- The author of `ParrottKilloh1984` is **Parrott**, with two t's. The citation key
  is unchanged, being referenced throughout the sources and documentation.

## v0.4.0 — Equilibrium actually coupled to kinetics, and dual numbers everywhere

### Breaking

- **`ChemicalState` gained a fourth type parameter**, `R <: Real`, carrying the
  number type of the dimensionless diagnostics: the signature is now
  `ChemicalState{C, S, Q, R}`. Code that spelled the type out in full has to be
  updated; code that used `ChemicalState` without parameters is unaffected.
- `pH`, `pOH`, `porosity` and `saturation` no longer return `Float64`
  unconditionally. They return whatever number type the composition carries, so
  that `∂pH/∂composition` and `∂porosity/∂(w/c)` are obtainable.

### Fixed — the equilibrium was never actually solved during a kinetics run

Four independent defects, each sufficient on its own to disable re-speciation.
They combined to make the whole coupling dead code, and the errors were hidden.

- **`equilibrate(state, ::EquilibriumSolver)` throws.** That is the call
  `build_kinetics_ode` made. `equilibrate` expects a SciML *algorithm* and wraps
  it into an `EquilibriumSolver`; handing it one that is already built wraps it
  twice, and the solve has no applicable method. The entry point for a prebuilt
  solver is `solve(esolver, state)`, which is what is called now.
- **The failure was swallowed by a bare `catch`**, so a run in which
  re-speciation never once succeeded looked exactly like a healthy one.
  Failures are now counted, the first is reported with its exception, and the
  total is reported at the end of `integrate`.
- **`KineticsSolver(; equilibrium_solver = …)` was never read.** `integrate`
  built its parameters from the `KineticsProblem` alone, so the documented way
  of passing a solver did nothing. It is honored now; a solver set on the
  problem still wins, and the conflict is reported.
- **The integrated element amounts `bₑ` were discarded.** The ODE integrates
  `dbₑ/dt = Aₑ νₑᵀ r`, but the re-speciation rebuilt its input from the
  persistent buffer instead, so the constraint of the equilibrium sub-problem
  (Leal et al. 2017, Eq. 54) was not the quantity being integrated. The previous
  speciation is now projected onto `bₑ` through `pinv(Aₑ)` before the solve.

### Fixed — three defects that made the `kinetic_species` API produce nonsense

Found while implementing the partitioned coupling, all on the 3-argument
`KineticsProblem(cs, state, tspan)` path:

- **The controlling mineral was mis-identified.** `_find_mineral_idx` takes the
  first crystalline reactant and, failing that, the first reactant of any kind.
  On a reaction generated from the nullspace — where the kinetic mineral may sit
  on either side — that is the solvent, so the ODE integrated **water** as the
  dissolving phase. A species the system declares kinetic now wins over any
  guess.
- **The stoichiometry had an arbitrary orientation**, so the kinetic mineral
  could carry `ν = +1` and the ODE grew the clinker: α(C₃S) came out at −4.2
  over seven days. Reactions are now normalized to `ν = −1` on their controlling
  mineral, which also delivers the `1/|νₖ|` scaling the `ChemicalSystem`
  docstring already promised.
- **The conservation basis was inconsistent.** `Aₑ` was cut from the canonical
  element matrix `CSM.A` while the solve is posed on `SM.A`, the matrix with
  respect to the primaries; and the equilibrium sub-system re-derived its own
  primaries instead of inheriting the parent's. Either mismatch makes every
  equilibrium solve fail — 89 failures out of 89 steps in one intermediate
  state. `Aₑ` now comes from the sub-system itself, which shares the parent's
  primaries.

With the three fixed, seven days of C₃S/C₂S hydration at `w/c = 0.4` gives
α(C₃S) = 0.277, portlandite and C-S-H in comparable amounts, an alkaline pore
solution at millimolar calcium, and `‖Aₑnₑ − bₑ‖∞ = 8.7e-8`.

### Changed — the coupling is now the partitioned formulation of Leal et al.

`respeciate!` solves `nₑ = φ(bₑ)` over the **equilibrium partition only**, on a
sub-system built for the purpose, with the integrated element amounts `bₑ`
handed to the solver as the constraint rather than derived from a composition.
Running the minimization over the whole system, as the previous code did, would
equilibrate the kinetic minerals instantaneously.

**Validated against Reaktoro.** Calcite dissolving at a constant rate makes the
kinetic half analytic, so every difference of rate-law convention drops out and
only `nₑ(t) = φ(bₑ(t))` is under test. The two codes agree to 0.3 % at `t = 0`
and 4.3 % at 3600 s, the largest deviation always on `CO₂(aq)` at `2.3e-8` mol.
Assertions in `test/coupling_reference.jl`, oracle in
`test/reference/reaktoro_coupling.py`.

### Changed — operator splitting instead of a solve inside the right-hand side

The equilibrium sub-problem is no longer solved inside the ODE right-hand side.
It is solved once per accepted step, by the new `respeciate!`, wired as a
`DiscreteCallback`. The right-hand side is evaluated many times per step and
differentiated to build the Jacobian; solving an optimization problem in there
made the cost unpredictable and the Jacobian inconsistent with the model being
integrated. With the solve outside, residual and Jacobian see the same
speciation.

### Fixed — the conservation matrix was rebuilt by finite differences

`OptimaSolverExt` did not pass `A` and `b` through the problem parameters, so
`OptimaSolver` fell back to reconstructing the constraint Jacobian by finite
differences with a `1e-7` step. The element balance was capped at ~2e-6 mol and
did **not** respond to the requested tolerance. `A` is known exactly; it is
handed over now. Measured on calcite/water: feasibility `2.17e-06 → 4.02e-14`.

This also settles the choice of default back-end. Once the matrix is exact,
OptimaSolver reaches machine-level feasibility *and* is 3 to 26 times faster
than Ipopt on the cement equilibria, so it remains the high-priority default.

### New — differentiating *through* the equilibrium solve

`ForwardDiff` now crosses an `equilibrate`. No solver is asked to iterate on
dual numbers — Ipopt is a C library and never could. The equilibrium is solved
once at the primal values and the sensitivities come from the optimality
conditions of [Leal2017](@cite), the problem Reaktoro solves: one saddle-point
factorization serves every partial derivative, and the answer is exact rather
than a finite difference.

The complementarity block is what partitions the species, and skipping it does
not degrade the answer gently. On calcite + CO₂ in water with a gas phase
declared, the unreduced system puts the whole response into the **absent** gas
species (`n = 5e-11` mol) — a sensitivity that satisfies the element balance to
`4e-16` and means nothing. No back-end returns the stability multipliers `z`, so
the active set is recovered internally: a species negligible on the scale of the
system that nonetheless takes a leading share of the response is pinned and the
system re-solved.

Verified against the package's own finite differences (`9e-5`, the truncation
error) and against **Reaktoro 2.13 reading the same Cemdata18 file**, over the
same eleven species and under the same activity model: the two agree to
`4.4e-4` relative or better on every species, below Reaktoro's own `7.2e-4`
spread across finite-difference step sizes. The equilibrium amounts agree to
the same order. The absent gas species gets exactly zero from the active-set
treatment, against `2e-9` by finite differences.

A cross-code comparison has three knobs — database, species list, activity
model — and each is worth tens of percent. Leaving Reaktoro on HKF against this
package's default `DiluteSolutionModel()` moves `∂Ca²⁺/∂(CO₂)` from `+0.1520`
to `+0.2179`; dropping the aqueous calcium complexes makes `∂Ca²⁺/∂(CO₂)` and
`∂calcite/∂(CO₂)` mirror each other, and restoring them breaks that mirror in
*both* codes alike. The element balance closes to `2e-16` throughout.

### Fixed — trace species now agree with Reaktoro

On calcite + CO₂, every species except `CaOH⁺` agrees with Reaktoro to 5 % or
better, and `pKw` comes out at 13.979 instead of 13.40. Before, `CaOH⁺` was
20.8× high, `OH⁻` 3.19× and `H⁺` 1.24×.

The cause was in the back-end's convergence test, which excluded variables
judged "at their bound" using a threshold scaled by the *largest* amount in the
problem. With a solvent at 55 mol that threshold was `5.5e-7`, so every trace
ion below it was declared to sit on a bound of `1e-16` and its stationarity was
never enforced. `OptimaSolver` 0.2.5 judges it against the variable's own bound.

`CaOH⁺` remains 2.5× high at `4e-9` mol. Ipopt lands on the same value (×2.46),
so it is attributable to neither back-end, and at that amount it is chemically
inconsequential. The assertion is `@test_broken`.

### Changed — the analytic gradient is handed to the back-end

`∂G/∂nᵢ = μᵢ` exactly, the remaining term vanishing by Gibbs–Duhem. The
extension now supplies it instead of leaving the back-end to differentiate
`dot(n, μ(n))`, where that cancellation is only exact in theory: evaluated, the
solvent term alone is `55 × 0.018 ≈ 1` against a `μ` of order 10. This does not
resolve the defect above, but a Gibbs energy is stationary exactly where its
gradient is, and steering on a gradient wrong by ten percent is not defensible
whatever else is true.

### Added — documentation

- **Coupling kinetics and equilibrium** — the partitioned formulation derived,
  including why the state is `(bₑ, nₖ)` rather than `(nₑ, nₖ)`, and what to
  check when a run finishes.
- **The hydrating paste, end to end** — a worked example computing the hydrate
  assemblage instead of imposing it, with measured output.
- **Validation against Reaktoro** — what agrees, what does not, and the three
  knobs a cross-code comparison has to match.

### Fixed — the solver's return code was ignored

Neither extension looked at the optimizer's `retcode`; a non-converged iterate
was written into the state as though it were the equilibrium. It is now checked
and warned about. The warning, rather than an error, is deliberate: the flag is
unreliable in both directions on these problems — pure water comes back flagged
`MaxIters` while giving `[H⁺]/[OH⁻] = 1.000003`. Set
`ChemistryLab.STRICT_CONVERGENCE[] = true` to raise instead.

### Fixed — dual numbers could not enter a `ChemicalState`

The constructor took its element type from the *temperature*
(`Q = typeof(T_q)`), which silently cast the composition to `Float64`: nothing
downstream of a `ChemicalState` was differentiable. The type is now the
promotion of `T`, `P` **and the amounts**. The dimensionless diagnostics follow
(see Breaking above); only their *display* strips the dual part.


## v0.3.1 — Maintenance

- `[compat]` upper bounds raised: `OrdinaryDiffEq` to `"6, 7"`, `TimerOutputs`
  to `"0.5, 1"`.
- CI badge restored; Runic badge.
- Installation instructions updated for registration in Julia's General
  registry; documented the optional optimization backend required to solve
  equilibria (`Optimization`+`OptimizationIpopt`, or `OptimaSolver`).
- Confirmed each GitHub Release keeps archiving automatically to Zenodo's
  existing concept DOI `10.5281/zenodo.17756074` via the native
  GitHub↔Zenodo integration (no workflow or token needed).
- MPCM-Registry deprecated in favor of the General registry.

## v0.3.0 — Chemical kinetics module

### New features

**`src/kinetics/` — mineral dissolution/precipitation kinetics**

- `KineticReaction` — couples a reaction to a rate law; supports
  `transition_state`, `first_order_rate`, and the empirical
  Parrott–Killoh (1984) model for cement clinker hydration
- `KineticsProblem` / `KineticsSolver` — ODE problem formulation
  following the SciML `(u, p, t)` convention; integrated via
  `KineticsOrdinaryDiffEqExt` (weakdep, activated by `using OrdinaryDiffEq`)
- `SemiAdiabaticCalorimeter` / `IsothermalCalorimeter` — coupled
  heat-balance ODE; `temperature_profile`, `heat_flow`,
  `cumulative_heat` post-processing helpers
- `StateView{T,I}` — O(1) named access to species data vectors in the
  ODE hot path (no per-step allocation)
- `RateMechanism`, `RateModelCatalyst`, `BETSurfaceArea`,
  `FixedSurfaceArea` — building blocks for custom rate closures
- `parrott_killoh(params, mineral_name)` factory with built-in
  Schindler & Follliard (2005) Arrhenius correction and default
  parameters for C₃S, C₂S, C₃A, C₄AF
- ForwardDiff-compatible throughout; `KineticFunc` and `transition_state`
  closures accept `Dual` numbers

### Infrastructure
- Relicensed to LGPL-2.1-or-later
- Registered in MPCM-Registry (OptimaSolver resolved via registry)
- GitHub Actions workflows: CI, Documentation, Register, CompatHelper, Format, TagBot
- Multi-version documentation deployment (`docs/deploy_docs.jl`)

## v0.2.3 — Activity models & solid solutions

- Extended Debye–Hückel and Davies activity models
- Redlich–Kister solid solution phases
- `SolidSolutionPhase`, `build_solid_solutions` from TOML
- `with_class` for end-member requalification

## v0.2.0 — Equilibrium solver integration

- `EquilibriumProblem` / `EquilibriumSolver` / `equilibrate` API
- `OptimizationIpoptExt` and `OptimaSolverExt` extensions
- Implicit-differentiation sensitivity via OptimaSolver
- HKF aqueous solute thermodynamic model (`NumericFunc`)
- ThermoFun JSON and PHREEQC `.dat` database parsers
