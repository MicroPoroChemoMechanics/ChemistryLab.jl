# Oracle for test/equilibrium_reference.jl.
#
# Reaktoro 2.13 reading the *same* Cemdata18 ThermoFun file this package ships,
# over the *same* species, under the *same* (ideal) activity model — the three
# knobs that have to match for a cross-code comparison to mean anything.
#
#   conda create -n reaktoro-env -c conda-forge reaktoro thermofun
#   conda run -n reaktoro-env python test/reference/reaktoro_calcite_co2.py
#
# Writes `test/reference/reaktoro_calcite_co2.json`, which the test reads.

import hashlib
import json
import os
import sys

import reaktoro as rkt
import numpy as np

DB = "data/cemdata18-thermofun.json"

SPECIES = "H2O@ H+ OH- CO2@ HCO3- CO3-2 Ca+2 CaOH+ Ca(CO3)@ Ca(HCO3)+".split()
N_H2O, N_CAL, N_CO2 = 55.5, 0.05, 0.01

db = rkt.ThermoFunDatabase.fromFile(DB)
solution = rkt.AqueousPhase(" ".join(SPECIES))
solution.set(rkt.ActivityModelIdealAqueous())   # matches DiluteSolutionModel()
system = rkt.ChemicalSystem(db, solution, rkt.MineralPhase("Cal"))
names = [s.name() for s in system.species()]


def speciate(n_co2):
    state = rkt.ChemicalState(system)
    state.temperature(25, "celsius")
    state.pressure(1, "bar")
    state.set("H2O@", N_H2O, "mol")
    state.set("Cal", N_CAL, "mol")
    state.set("CO2@", n_co2, "mol")
    result = rkt.equilibrate(state)
    assert result.succeeded(), "Reaktoro failed to equilibrate"
    return np.array([state.speciesAmount(n) for n in names])


amounts = speciate(N_CO2)

# Central differences, and the spread across step sizes that bounds how well
# this oracle can be trusted in the first place.
derivs = {}
for h in (1e-3, 1e-4, 1e-5):
    derivs[h] = (speciate(N_CO2 + h) - speciate(N_CO2 - h)) / (2 * h)
spread = np.max(np.abs(derivs[1e-3] - derivs[1e-5]))
d = derivs[1e-5]

with open(DB, "rb") as handle:
    md5 = hashlib.md5(handle.read()).hexdigest()
payload = {
    "generator": "test/reference/reaktoro_calcite_co2.py",
    "python": sys.version.split()[0],
    "reaktoro": rkt.__version__,
    "database": DB,
    "database_md5": md5,
    "aqueous_model": "ideal, to match DiluteSolutionModel()",
    "n_H2O": N_H2O,
    "n_Cal": N_CAL,
    "n_CO2": N_CO2,
    # the finite-difference spread across h in {1e-3, 1e-4, 1e-5}: nothing below
    # it is meaningful in the sensitivities
    "fd_spread": float(spread),
    "species": {
        name: {"n": float(amounts[k]), "dn_dCO2": float(d[k])}
        for k, name in enumerate(names)
    },
}
here = os.path.dirname(os.path.abspath(__file__))
path = os.path.join(here, "reaktoro_calcite_co2.json")
with open(path, "w", encoding="utf-8") as handle:
    json.dump(payload, handle, indent=2, sort_keys=False)
    handle.write("\n")
print(f"wrote {os.path.basename(path)}")
