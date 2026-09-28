# GEMS3K on the paste of Lothenbach & Winnefeld (2006), as an oracle for the
# certified equilibria ChemistryLab computes on it.
#
#   julia --project=docs test/reference/xgems_lw2006.jl        # the budgets first
#   conda run -n mpcm-oracles python test/reference/xgems_lw2006.py
#
# The Julia half writes, for each age of their Table 3, the element budget the
# recipe puts into the equilibrium and the species of ChemistryLab's system. This
# half runs the CEMDATA18 cement export of xGEMS (the one `xgems_cement.py`
# fetches, at the same pinned commit, into a cache outside this repository) on
# those budgets, at 20 degrees C, with the SAME phases: every solid and gas the
# system does not declare is removed by raising its standard Gibbs energy by
# 1e6 J/mol, and every aqueous species is kept, so that the traces of elements
# the budget lacks have somewhere to be. (An upper bound of zero on those species
# makes GEMS3K stop without converging.) The aqueous model is the export's own,
# the extended Debye-Hueckel of Cemdata18. What is written, `xgems_lw2006.json`,
# is numbers computed by running the code, with the identity of the code, of the
# export and of the budgets recorded as fields.

import hashlib
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from xgems_cement import FILES, PREFIX, REPO, fetch, resolve_commit  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.environ.get(
    "CHEMISTRYLAB_ORACLE_CACHE", os.path.join(os.path.expanduser("~"), ".cache", "chemistrylab-oracles")
)
T_K = 293.15
ELEMENTS_REPORTED = ("K", "Na", "Ca", "S", "Si", "Al")
REMOVAL_J_PER_MOL = 1.0e6


def main():
    import xgems

    path = os.path.join(CACHE, "lw2006", "budgets.json")
    if not os.path.isfile(path):
        sys.exit(f"missing {path}: run test/reference/xgems_lw2006.jl first")
    raw = open(path, "rb").read()
    inputs = json.loads(raw)
    ours = set(inputs["species"])

    commit = resolve_commit()
    d, digests = fetch(commit)
    cwd = os.getcwd()
    os.chdir(d)
    try:
        e = xgems.ChemicalEngine("CemHyds-dat.lst")
    finally:
        os.chdir(cwd)
    elements = [e.elementName(i) for i in range(e.numElements())]
    species = [e.speciesName(i) for i in range(e.numSpecies())]
    missing = sorted(ours - set(species))
    if missing:
        sys.exit(f"the export lacks species the system declares: {missing}")
    aq = e.indexPhase(e.aqueousPhaseName())
    e.setPT(T_K, 1.0e5)
    removed = []
    for i, s in enumerate(species):
        if s not in ours and e.indexPhaseWithSpecies(i) != aq:
            e.setStandardMolarGibbsEnergy(i, e.standardMolarGibbsEnergy(i) + REMOVAL_J_PER_MOL)
            removed.append(s)
    iw = species.index("H2O@")
    M_w = e.speciesMolarMasses()[iw]

    rows = []
    for row in inputs["rows"]:
        b = np.array([max(row["elements"].get(el, 0.0), 1.0e-12) for el in elements])
        # A cold start first; GEMS3K sometimes stops on one (code 4), and a warm
        # start from the previous age then converges. Which one served is recorded.
        for start in ("cold", "warm"):
            e.setB(b)
            if start == "cold":
                e.setColdStart()
                code = e.reequilibrate(False)
            else:
                e.setWarmStart()
                code = e.reequilibrate(True)
            if e.converged():
                break
        n = np.array(e.speciesAmounts())
        kg = n[iw] * M_w
        in_solution = np.array(e.elementAmountsInPhase(aq))
        mmol = {el: 1000 * float(in_solution[elements.index(el)]) / kg for el in ELEMENTS_REPORTED}
        mmol["OH-"] = 1000 * float(n[species.index("OH-")]) / kg
        rows.append(
            {
                "time_h": row["time_h"],
                "elements": row["elements"],
                "converged": bool(e.converged()),
                "code": int(code),
                "start": start,
                "pH": float(e.pH()),
                "mmol_per_kg_water": mmol,
            }
        )
        print(f"{row['time_h']:8.2f} h  converged={e.converged()} ({start} start)  pH={e.pH():.4f}")

    out = {
        "generator": "test/reference/xgems_lw2006.py, after test/reference/xgems_lw2006.jl",
        "code": "xgems " + getattr(xgems, "__version__", "unknown"),
        "export": {"repository": REPO, "commit": commit, "path": PREFIX, "files": FILES, "sha256": digests,
                   "license": "LGPL-3.0 (xGEMS); the export files are not redistributed here"},
        "database": "CEMDATA18, as exported in the CemGEMS system of xGEMS",
        "budgets_sha256": hashlib.sha256(raw).hexdigest(),
        "temperature_K": T_K,
        "species": sorted(ours),
        "removed_by_raising_G0": {"J_per_mol": REMOVAL_J_PER_MOL, "count": len(removed)},
        "rows": rows,
    }
    target = os.path.join(HERE, "xgems_lw2006.json")
    with open(target, "w") as fh:
        json.dump(out, fh, indent=1)
    print("wrote", target, len(rows), "rows")
    return 0


if __name__ == "__main__":
    sys.exit(main())
