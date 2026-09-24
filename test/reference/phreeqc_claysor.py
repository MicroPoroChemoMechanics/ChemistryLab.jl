#!/usr/bin/env python3
"""A ClaySor 2023 subset, run in PHREEQC, for a cross-code gate.

WHAT IS COMPARED, AND WHAT IS NOT CLAIMED

ClaySor 2023 is the 2SPNE SC/CE sorption model of Bradbury and Baeyens for
illite and montmorillonite: two-site protolysis, NON-electrostatic, plus cation
exchange. Its own first lines say it is written against the PSI/Nagra TDB 2020
aqueous database, and its constants are that database's.

This repository does not have TDB 2020 and does not ship it. So what runs here
is the ClaySor SORPTION MODEL over `phreeqc.dat`'s aqueous chemistry, on both
sides of the comparison. That tests the implementation of the model — site
budgets, protolysis, competition, the exchange convention — and is **not** a
reproduction of ClaySor, which would need the database it was fitted with. The
difference is the same one as using a surface constant fitted with a diffuse
layer in a model without one, and it is stated rather than hoped over.

THE SUBSET, NAMED

Na-montmorillonite: the three edge site families (`Mnt_s`, `Mnt_v`, `Mnt_w`),
their protolysis, and Na/Ca exchange on the planar sites. No metals, no illite,
no frayed edges. Site capacities and the specific surface are the ones ClaySor's
own header documents for 1 g of Na-montmorillonite.

Usage:
    python3 test/reference/phreeqc_claysor.py --model /path/to/claysor23_v0.7.dat
"""

import argparse
import hashlib
import json
import os
import re
import sys

from phreeqpython.viphreeqc import VIPhreeqc

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.join(HERE, "phreeqc.dat")

# The subset, by the species each reaction produces.
SURFACE_WANTED = {
    "Mnt_sOH2+", "Mnt_sO-", "Mnt_vOH2+", "Mnt_vO-", "Mnt_wOH2+", "Mnt_wO-",
}
EXCHANGE_WANTED = {"MntxNa", "Mntx2Ca"}
MASTERS = [("Mnt_s", "Mnt_sOH"), ("Mnt_v", "Mnt_vOH"), ("Mnt_w", "Mnt_wOH")]

# Documented by ClaySor's own header, for 1 g of Na-montmorillonite.
SITES = {"Mnt_sOH": 2.0e-6, "Mnt_vOH": 4.0e-5, "Mnt_wOH": 4.0e-5}   # mol/g
CEC = 0.00087            # mol/g
SPECIFIC_AREA = 18.0     # m2/g
GRAMS = 1.0

PH_VALUES = [4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
NA, CA = 0.1, 0.005      # mol/kgw


def parse_model(path):
    """Reactions of the subset, with the log K, reference and error as written."""
    text = open(path, encoding="latin-1").read()
    out = {"surface": [], "exchange": []}
    block = None
    pending = None
    for raw in text.split("\n"):
        line = rstrip = raw.rstrip()
        if re.match(r"^SURFACE_SPECIES", line):
            block = "surface"; pending = None; continue
        if re.match(r"^EXCHANGE_SPECIES", line):
            block = "exchange"; pending = None; continue
        if re.match(r"^[A-Z_]{4,}\s*$", line.strip()) and not line.startswith((" ", "\t")):
            block = None; continue
        if block is None:
            continue
        body, _, comment = line.partition("#")
        s = body.strip()
        if not s:
            continue
        if "=" in s and not s.startswith("-"):
            pending = (s, comment.strip())
            continue
        m = re.match(r"^-log_k\s*=?\s*(\S+)", s, re.I)
        if m and pending:
            equation, comment = pending
            pending = None
            products = equation.split("=", 1)[1]
            names = set(re.findall(r"[A-Za-z_][\w()+\-]*", products))
            wanted = SURFACE_WANTED if block == "surface" else EXCHANGE_WANTED
            if not (names & wanted):
                continue
            ref = re.search(r"ref:\s*(\S+)", comment)
            err = re.search(r"error:\s*([0-9.eE+-]+)", comment)
            out[block].append({
                "equation": equation,
                "log_K": float(m.group(1)),
                "reference": ref.group(1) if ref else None,
                "uncertainty": float(err.group(1)) if err else None,
            })
    missing = SURFACE_WANTED - {
        n for r in out["surface"] for n in re.findall(r"[A-Za-z_][\w()+\-]*", r["equation"])
    }
    if missing:
        raise SystemExit(f"{path} does not define {sorted(missing)}")
    return out


def database_with(model, path_out):
    """`phreeqc.dat` with the subset spliced in BEFORE its terminating END.

    After the `END` a database stops being read, and everything appended there
    loads without error and defines nothing — which shows up only as
    "Master species not in database" when a run uses it.
    """
    base = open(BASE, encoding="latin-1").read()
    lines = ["", "EXCHANGE_MASTER_SPECIES", "    Mntx    Mntx-", "EXCHANGE_SPECIES",
             "    Mntx- = Mntx-", "        -log_k 0.0"]
    for r in model["exchange"]:
        lines += [f"    {r['equation']}", f"        -log_k {r['log_K']}"]
    lines += ["SURFACE_MASTER_SPECIES"]
    lines += [f"    {m}    {ref}" for m, ref in MASTERS]
    lines += ["SURFACE_SPECIES"]
    for _, ref in MASTERS:
        lines += [f"    {ref} = {ref}", "        -log_k 0.0"]
    for r in model["surface"]:
        lines += [f"    {r['equation']}", f"        -log_k {r['log_K']}"]
    i = base.rindex("\nEND")
    open(path_out, "w", encoding="utf-8").write(
        base[:i] + "\n" + "\n".join(lines) + "\n" + base[i:]
    )
    return path_out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", required=True, help="path to claysor23_v0.7.dat")
    args = parser.parse_args()
    if not os.path.exists(args.model):
        raise SystemExit(
            f"{args.model} does not exist. ClaySor 2023 is CC-BY-4.0 and freely "
            "available from doi:10.5281/zenodo.15095062; this repository does "
            "not carry it."
        )

    model = parse_model(args.model)
    spliced = os.path.join(os.path.dirname(args.model), "_claysor_subset.dat")
    database_with(model, spliced)

    ip = VIPhreeqc()
    ip.load_database(spliced)
    if ip.phc_database_error_count:
        raise SystemExit(f"spliced database failed: {ip.get_error_string()}")

    with open(args.model, "rb") as h:
        digest = hashlib.sha256(h.read()).hexdigest()
    reported = ["MntxNa", "Mntx2Ca"] + [
        f"{ref}{suffix}" for _, ref in MASTERS for suffix in ("", "2+", "")
    ]
    reported = ["MntxNa", "Mntx2Ca"]
    for _, ref in MASTERS:
        stem = ref[:-2]          # "Mnt_sOH" -> "Mnt_s"
        reported += [ref, stem + "OH2+", stem + "O-"]

    payload = {
        "generator": "test/reference/phreeqc_claysor.py",
        "python": sys.version.split()[0],
        "model": "ClaySor 2023 (2SPNE SC/CE), Na-montmorillonite subset",
        "model_file": os.path.basename(args.model),
        "model_sha256": digest,
        "model_doi": "10.5281/zenodo.15095062",
        "model_licence": "CC-BY-4.0",
        "aqueous_database": "phreeqc.dat (NOT ClaySor's own PSI/Nagra TDB 2020)",
        "not_a_reproduction": (
            "ClaySor's constants were fitted against PSI/Nagra TDB 2020. Run "
            "here over phreeqc.dat on both sides, this tests the sorption model "
            "and is not a reproduction of ClaySor."
        ),
        "exchange_convention": "Gaines-Thomas (PHREEQC's EXCHANGE default)",
        "grams": GRAMS,
        "specific_area_m2_per_g": SPECIFIC_AREA,
        "cec_mol_per_g": CEC,
        "sites_mol_per_g": SITES,
        "na": NA,
        "ca": CA,
        "reactions": model,
        "points": [],
    }

    for ph in PH_VALUES:
        ip.run_string(f"""
PHASES
Fix_H+
    H+ = H+
    log_k 0.0
HCl
    HCl = H+ + Cl-
    log_k 7.0
END
SOLUTION 1
    units mol/kgw
    temp 25.0
    water 1.0
    pH 7.0
    Na {NA}
    Ca {CA}
    Cl {NA + 2 * CA} charge
EXCHANGE 1
    MntxNa {CEC * GRAMS}
SURFACE 1
    Mnt_sOH {SITES['Mnt_sOH'] * GRAMS} {SPECIFIC_AREA} {GRAMS}
    Mnt_vOH {SITES['Mnt_vOH'] * GRAMS}
    Mnt_wOH {SITES['Mnt_wOH'] * GRAMS}
    -no_edl
EQUILIBRIUM_PHASES 1
    Fix_H+ {-ph} HCl 10.0
SELECTED_OUTPUT
    -reset false
    -high_precision true
    -pH true
    -ionic_strength true
    -activities H+ Na+ Ca+2
    -molalities {' '.join(reported)}
END
""")
        rows = ip.get_selected_output_array()
        header, values = rows[0], rows[-1]
        col = lambda n: values[header.index(n)]
        got = col("pH")
        if abs(got - ph) > 1.0e-6:
            raise SystemExit(f"PHREEQC returned pH {got} for a requested {ph}")
        payload["points"].append({
            "pH": ph,
            "I": col("mu"),
            "la_H": col("la_H+"),
            "la_Na": col("la_Na+"),
            "la_Ca": col("la_Ca+2"),
            "species": {n: col("m_" + n + "(mol/kgw)") for n in reported},
        })

    out = os.path.join(HERE, "phreeqc_claysor.json")
    with open(out, "w", encoding="utf-8") as h:
        json.dump(payload, h, indent=2, sort_keys=False)
        h.write("\n")
    os.remove(spliced)
    print(f"wrote {os.path.basename(out)}")


if __name__ == "__main__":
    main()
