#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)
"""A surface whose site total follows the phase that carries it.

PHREEQC has coupled a `SURFACE` to an `EQUILIBRIUM_PHASES` mineral since v2:

    Hfo_sOH   Fe(OH)3(a)   equilibrium_phase   0.005   53300

reads "0.005 mol of strong sites per mole of Fe(OH)3(a), over 53300 m² per mole".
The engine says so in its own words when the amount changes — its message is
`Resetting number of sites in surface %s to be consistent with moles of phase %s`.

That is exactly what `SITES_FOLLOW_HOST` means in this package, which is why it
is the oracle for it. What this generator captures is the **coupling law**, not
the thermodynamics: a titration that dissolves the sorbent progressively, with
the phase amount and both site totals at every step.

WHAT IS NOT COMPARED, and the fixture says so in a field rather than in a
comment: the two codes are not run over one surface model. PHREEQC scales the
site totals from the phase and applies no back-reaction on the phase's own
stability; this package carries the coupling in the constraint matrix, so the
host's saturation index picks up the site potential. The two therefore agree on
the ratio and on the amounts the acid budget fixes, and are expected to differ
on dissolved iron.

Three traps are avoided here on purpose, all three measured:

1. `Fix_H+ <-pH> HCl 10.0` hands PHREEQC a ten-mole acid reservoir. Harmless
   when the solid cannot dissolve; with a phase-coupled surface it dissolves
   the whole sorbent and every site total comes back zero. The partial
   dissolution path is driven by a FINITE `REACTION` budget instead.
2. Under `-no_edl` the engine stops printing its `[0.005 mol/(mol ...)]`
   annotation and the specific-area line. The coupling is still live. Assert on
   numbers, never on that text.
3. `-sites_units` is inert for a coupled site: the coefficient is always moles
   of sites per mole of phase.

Run it from the repository root, in the `mpcm-oracles` environment:

    python3 test/reference/phreeqc_evolving_surface.py
"""

import ctypes
import hashlib
import importlib.metadata
import json
import os
import sys

from phreeqpython.viphreeqc import VIPhreeqc

HERE = os.path.dirname(os.path.abspath(__file__))

# The acid budget, in moles of HCl. Three moles dissolve one of Fe(OH)3, so this
# walks from the intact sorbent to none of it, with four partial states in
# between — the states that show the site totals tracking rather than switching.
HCL_STEPS = [0.0, 1.0e-3, 2.0e-3, 2.6e-3, 2.9e-3, 2.99e-3, 3.0e-3]

FE_INITIAL = 1.0e-3          # mol of Fe(OH)3(a) present at the start
SITES_STRONG = 0.005         # mol of Hfo_s per mol of Fe(OH)3(a)
SITES_WEAK = 0.2             # mol of Hfo_w per mol of Fe(OH)3(a)
AREA_PER_MOL = 53300.0       # m² per mol; parsed under -no_edl, unused by it


def database_path():
    """The database vendored in this repository, never the wheel's copy.

    `phreeqpython` ships a `phreeqc.dat` that matches no upstream tag — its md5
    is `fca384eb…` against v3.7.3's `4d3f4378…`. A comparison run over a
    database nobody can name is not a comparison.
    """
    path = os.path.join(HERE, "phreeqc.dat")
    if not os.path.exists(path):
        raise SystemExit(f"{path} is missing; see PHREEQC-PROVENANCE.md")
    return path


def engine_version(ip):
    """The PHREEQC engine version, read from the library at run time.

    `phreeqc_hfo_surface.py` records the `phreeqpython` version and the database
    md5 but not this, so the number that decides what the comparison means
    appeared only in prose. Three lines fix that.
    """
    lib = ip.dll
    lib.GetVersionString.restype = ctypes.c_char_p
    return lib.GetVersionString().decode()


def deck():
    steps = " ".join(f"{x:g}" for x in HCL_STEPS)
    return f"""
SOLUTION 1
    units     mol/kgw
    temp      25.0
    water     1.0
    pH        7.0
    Na        0.01
    Cl        0.01 charge
EQUILIBRIUM_PHASES 1
    Fe(OH)3(a)   0.0   {FE_INITIAL}
SURFACE 1
    Hfo_sOH   Fe(OH)3(a)   equilibrium_phase   {SITES_STRONG}   {AREA_PER_MOL}
    Hfo_wOH   Fe(OH)3(a)   equilibrium_phase   {SITES_WEAK}
    -no_edl
REACTION 1
    HCl 1.0
    {steps} moles
SELECTED_OUTPUT
    -reset          false
    -high_precision true
    -pH             true
    -molalities     Hfo_sOH Hfo_sOH2+ Hfo_sO- Hfo_wOH Hfo_wOH2+ Hfo_wO-
USER_PUNCH
    -headings  Fe_phase  kgw
    10 PUNCH EQUI("Fe(OH)3(a)")
    20 PUNCH TOT("water")
END
"""


def run():
    db = database_path()
    ip = VIPhreeqc()
    ip.load_database(db)
    if ip.phc_database_error_count:
        raise SystemExit("PHREEQC refused the database")
    ip.run_string(deck())
    rows = ip.get_selected_output_array()
    # PHREEQC PUNCHES THE INITIAL SOLUTION TOO, before any equilibration, so the
    # output carries one row more than there are reaction steps. Zipping the
    # steps against `rows[1:]` pairs every amount with the state before it and
    # silently drops the last — and the table stays entirely plausible, which is
    # how it went unnoticed until the row count was counted.
    header, values = rows[0], rows[2:]
    if len(values) != len(HCL_STEPS):
        raise SystemExit(
            f"PHREEQC returned {len(values)} states for {len(HCL_STEPS)} reaction "
            "steps; the pairing below would be wrong."
        )

    def col(name):
        return header.index(name)

    strong = [col(f"m_{s}(mol/kgw)") for s in ("Hfo_sOH", "Hfo_sOH2+", "Hfo_sO-")]
    weak = [col(f"m_{s}(mol/kgw)") for s in ("Hfo_wOH", "Hfo_wOH2+", "Hfo_wO-")]
    i_fe, i_ph, i_kgw = col("Fe_phase"), col("pH"), col("kgw")

    points = []
    for hcl, row in zip(HCL_STEPS, values):
        fe = float(row[i_fe])
        # `-molalities` are per KILOGRAM OF WATER, and the water mass is not 1 kg
        # once acid has been added and the solid has dissolved into it. Comparing
        # a molality sum with a mole count without this factor leaves a 7e-6
        # residual that looks like a coupling error and is a units error.
        kgw = float(row[i_kgw])
        points.append(
            {
                "hcl_added": hcl,
                "fe_phase": fe,
                "kgw": kgw,
                "pH": float(row[i_ph]),
                "strong_total": sum(float(row[k]) for k in strong) * kgw,
                "weak_total": sum(float(row[k]) for k in weak) * kgw,
            }
        )

    # The assertion that makes this a fixture rather than a table: PHREEQC's own
    # site totals must be the declared coefficient times the phase amount. If
    # they are not, the coupling was not live and there is nothing to compare.
    for pt in points:
        if pt["fe_phase"] <= 0.0:
            continue
        for key, coef in (("strong_total", SITES_STRONG), ("weak_total", SITES_WEAK)):
            got = pt[key] / pt["fe_phase"]
            if abs(got / coef - 1.0) > 1.0e-6:
                raise SystemExit(
                    f"PHREEQC did not couple {key}: {got} against {coef} "
                    f"at {pt['hcl_added']} mol HCl"
                )

    with open(db, "rb") as handle:
        md5 = hashlib.md5(handle.read()).hexdigest()

    payload = {
        "generator": "test/reference/phreeqc_evolving_surface.py",
        "python": sys.version.split()[0],
        "phreeqpython": importlib.metadata.version("phreeqpython"),
        "phreeqc_engine": engine_version(ip),
        "database": os.path.basename(db),
        "database_md5": md5,
        "phase": "Fe(OH)3(a)",
        "sites_strong_per_mol": SITES_STRONG,
        "sites_weak_per_mol": SITES_WEAK,
        "area_m2_per_mol": AREA_PER_MOL,
        "electrostatics": "-no_edl",
        "not_a_reproduction": (
            "This captures PHREEQC's COUPLING LAW, not a shared surface model. "
            "PHREEQC scales the site totals from the phase and applies no "
            "back-reaction on the phase's own stability, where this package "
            "carries the coupling in the constraint matrix so the host's "
            "saturation index picks up the site potential. The two agree on the "
            "ratio of sites to phase and on the amounts the acid budget fixes; "
            "they are expected to differ on dissolved iron."
        ),
        "points": points,
    }
    out = os.path.join(HERE, "phreeqc_evolving_surface.json")
    with open(out, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=False)
        handle.write("\n")
    print(f"wrote {out} with {len(points)} points, engine {payload['phreeqc_engine']}")


if __name__ == "__main__":
    run()
