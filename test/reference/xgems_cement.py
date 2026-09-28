# GEMS3K on a CEMDATA18 cement, as an oracle for the activity model.
#
#   conda run -n mpcm-oracles python test/reference/xgems_cement.py
#
# The xGEMS repository (LGPL-3.0) publishes GEMS3K exports of several systems,
# among them a CEMDATA18 cement paste ("CemHyds"). This script downloads that
# export at a pinned commit into a cache OUTSIDE this repository -- the export
# files are never committed here -- runs it with xgems, and writes
# `xgems_cement.json` beside itself: numbers computed by running the code, with
# the identity of the code, of the export and of the database recorded as
# fields.
#
# What it records:
#   * the parameters of the aqueous model the export runs, read from its IPM
#     file (the common ion size and the B-dot of the extended Debye-Hueckel);
#   * for the export's own paste, diluted with water in steps, the ionic
#     strength and the activity coefficient of each charge class at 20 and
#     25 degrees C: what `cemdata18_activity_model` must reproduce;
#   * the standard Gibbs energies at 25 degrees C of about twenty species, to
#     compare with the Cemdata18 file ChemistryLab reads before trusting that
#     the two describe the same database.

import hashlib
import json
import os
import re
import sys
import urllib.request

COMMIT = "7ee3d4e101b37ee1ab293dcec2403fd8515893ca"  # the export this fixture was computed from
REPO = "gemshub/xGEMS"
PREFIX = "demos/resources/CemGEMS-keyvalue"
REFERENCE_IONS = ("H+", "OH-", "K+", "Na+", "Ca+2", "SO4-2", "CO3-2", "AlO2-", "HSiO3-")
FILES = ["CemHyds-dat.lst", "CemHyds-dch.dat", "CemHyds-ipm.dat", "CemHyds-dbr-0-0000.dat"]
HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.environ.get(
    "XGEMS_EXPORT_CACHE", os.path.join(os.path.expanduser("~"), ".cache", "chemistrylab-oracles", "xgems")
)


def resolve_commit():
    """The pinned commit, unless XGEMS_COMMIT names another one."""
    return os.environ.get("XGEMS_COMMIT", "").strip() or COMMIT


def fetch(commit):
    d = os.path.join(CACHE, commit)
    os.makedirs(d, exist_ok=True)
    digests = {}
    for f in FILES:
        path = os.path.join(d, f)
        if not os.path.isfile(path):
            url = f"https://raw.githubusercontent.com/{REPO}/{commit}/{PREFIX}/{f}"
            urllib.request.urlretrieve(url, path)
        with open(path, "rb") as fh:
            digests[f] = hashlib.sha256(fh.read()).hexdigest()
    return d, digests


def aqueous_parameters(ipm_path):
    """b_gamma and the common ion size of the aqueous phase, from <PMc>."""
    text = open(ipm_path).read()
    m = re.search(r"^<PMc>\s*([^<#]+)", text, re.M)
    values = [float(v) for v in m.group(1).split()]
    return {"b_gamma": values[0], "ion_size_angstrom": values[1]}


def main():
    import numpy as np
    import xgems

    commit = resolve_commit()
    d, digests = fetch(commit)
    cwd = os.getcwd()
    os.chdir(d)
    try:
        e = xgems.ChemicalEngine("CemHyds-dat.lst")
    finally:
        os.chdir(cwd)
    names = [e.elementName(i) for i in range(e.numElements())]
    b0 = np.array(e.elementAmounts())
    iH, iO = names.index("H"), names.index("O")
    species = [e.speciesName(i) for i in range(e.numSpecies())]
    charges = [e.speciesCharge(i) for i in range(e.numSpecies())]

    rows = []
    for T in (293.15, 298.15):
        for water_factor in (1.0, 1.5, 2.0, 3.0, 5.0, 10.0):
            b = b0.copy()
            # Add water: the paste's own water times (factor - 1), counted on its
            # oxygen, which the export's recipe holds mostly as H2O.
            extra = (water_factor - 1.0) * (b0[iO] / 2.0)
            b[iH] += 2 * extra
            b[iO] += extra
            code = e.equilibrate(T, 1.0e5, b)
            lng = np.array(e.lnActivityCoefficients())
            gamma = {s: float(np.exp(lng[species.index(s)])) for s in REFERENCE_IONS if s in species}
            rows.append(
                {
                    "T_K": T,
                    "water_factor": water_factor,
                    "converged": bool(e.converged()),
                    "code": int(code),
                    "ionic_strength_mol_per_kg": float(e.ionicStrength()),
                    "pH": float(e.pH()),
                    # The coefficient of each reference ion. With a common ion
                    # size, ions of one charge share one coefficient; a few
                    # species of the export have theirs fixed at 1, which is why
                    # named ions are recorded rather than a class average.
                    "gamma": gamma,
                }
            )

    e.setPT(298.15, 1.0e5)
    g0 = {}
    for s in ("H2O@", "H+", "OH-", "Ca+2", "K+", "Na+", "SO4-2", "AlO2-", "HSiO3-", "SiO2@", "CO3-2",
              "Portlandite", "Cal", "ettringite", "monosulphate12", "monocarbonate", "C3AH6",
              "hydrotalcite", "Gp", "CSHQ-TobD", "CSHQ-JenH", "KSiOH", "NaSiOH"):
        if s in species:
            g0[s] = float(e.standardMolarGibbsEnergy(species.index(s)))

    out = {
        "generator": "test/reference/xgems_cement.py",
        "code": "xgems " + getattr(xgems, "__version__", "unknown"),
        "export": {"repository": REPO, "commit": commit, "path": PREFIX, "sha256": digests,
                   "license": "LGPL-3.0 (xGEMS); the export files are not redistributed here"},
        "database": "CEMDATA18, as exported in the CemGEMS system of xGEMS",
        "aqueous_model": aqueous_parameters(os.path.join(d, "CemHyds-ipm.dat")),
        "rows": rows,
        "standard_gibbs_298K_J_per_mol": g0,
    }
    target = os.path.join(HERE, "xgems_cement.json")
    with open(target, "w") as fh:
        json.dump(out, fh, indent=1)
    print("wrote", target, len(rows), "rows")
    return 0


if __name__ == "__main__":
    sys.exit(main())
