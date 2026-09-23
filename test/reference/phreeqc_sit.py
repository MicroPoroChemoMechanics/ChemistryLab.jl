#!/usr/bin/env python3
"""Activity coefficients under SIT, from PHREEQC, for a cross-code gate.

WHAT THIS COMPARES, AND WHAT IT DELIBERATELY DOES NOT

The fixture carries, for each point, the molality AND the activity coefficient
PHREEQC computed for every reported species. The Julia side then evaluates its
own SIT at *those* molalities. So the comparison is of the activity formula
alone: no speciation, no database of log K, no convergence path enters it.

That matters because the database here is not ours to reproduce. `sit.dat` is
the ANDRA/RWM ThermoChimie-TDB compilation, redistributed with PHREEQC but not
USGS-authored, so this package ships neither the file nor its ε. The ε actually
used are written into the fixture — a handful of coefficients, with the database
they came from named — which is what makes the test runnable in CI without
redistributing a compilation.

Usage:
    python3 test/reference/phreeqc_sit.py --database /path/to/sit.dat
"""

import argparse
import hashlib
import json
import os
import sys

import phreeqpython
from phreeqpython.viphreeqc import VIPhreeqc

# NaCl across four decades of ionic strength, which is where SIT earns its
# pairwise term: at 0.001 m it is Debye-Hückel and at 3 m it is the ε that
# separates it from Davies.
MOLALITIES = [0.001, 0.01, 0.1, 0.5, 1.0, 2.0, 3.0]
REPORTED = ["Na+", "Cl-", "H+", "OH-"]
# The pairs the Julia side needs to evaluate the same formula on the same
# composition. Read from the database, never typed here.
PAIRS = [("Na+", "Cl-"), ("H+", "Cl-"), ("Na+", "OH-")]


def epsilon_from(path, pairs):
    """The ε of `pairs`, read out of the database's own SIT block."""
    want = {tuple(sorted(p)) for p in pairs}
    found, in_block, seen = {}, False, False
    with open(path, encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.split("#")[0].strip()
            if not in_block:
                in_block = line.upper() == "SIT"
                continue
            if not line:
                continue
            if line.startswith("-"):
                seen = line.lower() == "-epsilon"
                continue
            f = line.split()
            if len(f) == 3 and seen:
                key = tuple(sorted(f[:2]))
                if key in want:
                    found[key] = float(f[2])
            else:
                break
    missing = want - set(found)
    if missing:
        raise SystemExit(f"{path} has no ε for {sorted(missing)}")
    return found


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--database",
        default=os.path.join(
            os.path.dirname(phreeqpython.__file__), "database", "sit.dat"
        ),
        help="a PHREEQC database carrying a SIT block (default: phreeqpython's)",
    )
    args = parser.parse_args()
    db = args.database
    if not os.path.exists(db):
        raise SystemExit(
            f"{db} does not exist. `sit.dat` is the ANDRA/RWM ThermoChimie "
            "compilation distributed with PHREEQC; this repository does not "
            "carry it, so point --database at your own copy."
        )

    ip = VIPhreeqc()
    ip.load_database(db)
    if ip.phc_database_error_count:
        raise SystemExit(f"{db} failed to load: {ip.get_error_string()}")

    with open(db, "rb") as handle:
        digest = hashlib.sha256(handle.read()).hexdigest()
    try:
        from importlib.metadata import version

        pp_version = version("phreeqpython")
    except Exception:  # pragma: no cover
        pp_version = "unknown"

    payload = {
        "generator": "test/reference/phreeqc_sit.py",
        "python": sys.version.split()[0],
        "phreeqpython": pp_version,
        "database": os.path.basename(db),
        "database_sha256": digest,
        "database_note": (
            "ANDRA/RWM ThermoChimie-TDB, redistributed with PHREEQC. Not "
            "USGS-authored and not vendored here; only the ε used below travel "
            "with this fixture."
        ),
        "epsilon": [
            {"a": k[0], "b": k[1], "value": v}
            for k, v in sorted(epsilon_from(db, PAIRS).items())
        ],
        "points": [],
    }

    for m in MOLALITIES:
        script = f"""
SOLUTION 1
    units    mol/kgw
    temp     25.0
    water    1.0
    pH       7.0 charge
    Na       {m}
    Cl       {m}
SELECTED_OUTPUT
    -reset      false
    -high_precision true
    -ionic_strength true
    -molalities {' '.join(REPORTED)}
    -activities {' '.join(REPORTED)}
END
"""
        ip.run_string(script)
        rows = ip.get_selected_output_array()
        header, values = rows[0], rows[-1]
        col = lambda name: values[header.index(name)]
        species = {}
        for nm in REPORTED:
            mol = col(f"m_{nm}(mol/kgw)")
            la = col(f"la_{nm}")
            if mol <= 0.0:
                continue
            species[nm] = {"molality": mol, "log_gamma": la - __import__("math").log10(mol)}
        payload["points"].append(
            {"nacl": m, "I": col("mu"), "species": species}
        )

    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(here, "phreeqc_sit.json")
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=False)
        handle.write("\n")
    print(f"wrote {os.path.basename(path)}")


if __name__ == "__main__":
    main()
