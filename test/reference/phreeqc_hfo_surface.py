# Oracle for the surface-complexation tests: hydrous ferric oxide, no electric
# double layer.
#
# The canonical benchmark of the field is Dzombak & Morel's HFO titration, and
# PHREEQC ships its parameters in `phreeqc.dat` as the `Hfo_s` / `Hfo_w`
# surface. This script runs the *non-electrostatic* variant, `-no_edl`, because
# that is the model the first milestone implements: explicit strong and weak
# sites, ideal site mixing, no surface potential. The published D&M calibration
# assumes the diffuse layer and is a separate gate — run this same script with
# `--edl` for it, and do not compare the two.
#
#   conda create -n mpcm-oracles -c conda-forge python=3.12 xgems gems3k
#   /path/to/envs/mpcm-oracles/bin/python -m ensurepip --upgrade
#   /path/to/envs/mpcm-oracles/bin/python -m pip install phreeqpython
#   conda run -n mpcm-oracles python test/reference/phreeqc_hfo_surface.py
#
# Paste the printed block into the fixture of `test/surface_complexation.jl`.
#
# What has to match for the comparison to mean anything, and what does not:
#
#   * the site densities and the reaction set -- taken from `phreeqc.dat`, so
#     they are the reference's own, not ours;
#   * the aqueous activity model -- PHREEQC's own extended Debye-Huckel here,
#     which is *not* one of the four this package ships. A comparison of
#     surface amounts at low ionic strength is still meaningful because the
#     surface reactions involve mostly the same ions on both sides; a
#     comparison of activity coefficients is not, and is not attempted.
#   * `-no_edl` on both sides. Leaving it out changes the model, not the
#     numerics.

import argparse
import json
import hashlib
import os
import re
import sys

import phreeqpython
from phreeqpython.viphreeqc import VIPhreeqc

DATABASE = "phreeqc.dat"

# Dzombak & Morel's ferrihydrite: 0.005 mol strong sites and 0.2 mol weak sites
# per mole of Fe, 600 m2/g, and 89 g/mol of Fe(OH)3. One millimole of Fe here.
N_FE = 1.0e-3
SITES_STRONG = 0.005 * N_FE
SITES_WEAK = 0.2 * N_FE
AREA_PER_GRAM = 600.0
MASS_SOLID = 89.0 * N_FE * 1.0e-3        # kg -> PHREEQC wants grams

ZN_TOTAL = 1.0e-5                        # mol/kgw
PH_VALUES = [4.0, 5.0, 6.0, 7.0, 8.0, 9.0]

REPORTED = [
    "Hfo_sOH", "Hfo_sOH2+", "Hfo_sO-", "Hfo_sOZn+",
    "Hfo_wOH", "Hfo_wOH2+", "Hfo_wO-", "Hfo_wOZn+",
    "Zn+2",
]

# The protolysis-only case: weak sites, no metal. It is the narrowest comparison
# that still exercises everything the first milestone implements — a site
# balance, ideal site mixing, and two states competing for one budget — and it
# is narrow on purpose: with no metal there is nothing in it that depends on the
# aqueous activity model, so the two codes are compared on the surface chemistry
# alone.
PROTOLYSIS_REPORTED = ["Hfo_wOH", "Hfo_wOH2+", "Hfo_wO-"]
SITES_PROTOLYSIS = 2.0e-4


def logk_for(db: str, reaction: str) -> float:
    """The log K of one reaction, read from the database rather than typed.

    A constant copied by hand is a constant that drifts from the file it came
    from; reading it means the comparison is against what PHREEQC actually used.
    The match is on the reaction written with single spaces, so the database's
    mixture of tabs and spaces does not matter.
    """
    def norm(text):
        return " ".join(text.replace("\t", " ").split())

    target = norm(reaction)
    lines = open(db, encoding="utf-8", errors="ignore").read().split("\n")
    for i, line in enumerate(lines):
        if norm(line) != target:
            continue
        for follow in lines[i + 1 : i + 4]:
            if "log_k" in follow:
                return float(follow.split("log_k")[1].split("#")[0])
    raise SystemExit(f"no log K for {reaction!r} in {os.path.basename(db)}")


def emit(name: str, payload: dict) -> None:
    """Write a fixture to `test/reference/<name>.json`.

    JSON rather than a block of Julia pasted into a test file, for three
    reasons that were all learned the hard way. A pasted block is a manual
    transcription, which is the error class these generators exist to remove.
    It has to be written in the *formatter's* dialect — this script used to
    carry a `jl()` that rewrote `6.0e-09` as `6.0e-9` and `5e-06` as `5.0e-6`,
    because Runic gates CI and Python renders neither the way Julia does. And a
    provenance header written as comments is not readable by anything; here it
    is a field.
    """
    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(here, f"{name}.json")
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=False)
        handle.write("\n")
    print(f"wrote {os.path.relpath(path, os.path.dirname(here))}")


def database_path(name: str) -> str:
    """The database this repository ships, not the one the wheel happens to carry.

    `phreeqpython` bundles its own `phreeqc.dat`, and using it made these
    fixtures depend on a file that is not in this repository: the md5 recorded
    beside them pinned something nobody here could check, and that copy matches
    no upstream PHREEQC tag. The vendored one does, and it is the database of
    exactly the engine version bundled here (3.7.3). See PHREEQC-PROVENANCE.md.
    """
    here = os.path.dirname(os.path.abspath(__file__))
    vendored = os.path.join(here, name)
    if os.path.exists(vendored):
        return vendored
    raise SystemExit(
        f"{name} is not in {here}; it is committed there on purpose — "
        "see PHREEQC-PROVENANCE.md"
    )


def provenance(db: str, **extra) -> dict:
    """Everything needed to reproduce this run, read at run time, never assumed."""
    with open(db, "rb") as handle:
        digest = hashlib.md5(handle.read()).hexdigest()
    try:
        from importlib.metadata import version

        pp_version = version("phreeqpython")
    except Exception:  # pragma: no cover - provenance must never fail the run
        pp_version = "unknown"
    return {
        "generator": "test/reference/phreeqc_hfo_surface.py",
        "python": sys.version.split()[0],
        "phreeqpython": pp_version,
        "database": os.path.basename(db),
        "database_md5": digest,
        **extra,
    }


def run_zn_edge(edl: bool):
    """The two-site case: strong and weak sites, protolysis, and a metal.

    Dzombak & Morel's ferrihydrite has two families that differ in exactly two
    ways — the strong sites are forty times scarcer and bind zinc about a
    thousand times more strongly — and share the same protolysis constants. The
    consequence is the shape of a sorption edge: the strong sites take the metal
    first and saturate, the weak ones take over.
    """
    db = database_path(DATABASE)
    ip = VIPhreeqc()
    ip.load_database(db)
    if ip.phc_database_error_count:
        raise SystemExit(f"{DATABASE} failed to load: {ip.get_error_string()}")

    logks = {
        "protonation": logk_for(db, "Hfo_wOH + H+ = Hfo_wOH2+"),
        "deprotonation": logk_for(db, "Hfo_wOH = Hfo_wO- + H+"),
        "zn_strong": logk_for(db, "Hfo_sOH + Zn+2 = Hfo_sOZn+ + H+"),
        "zn_weak": logk_for(db, "Hfo_wOH + Zn+2 = Hfo_wOZn+ + H+"),
    }

    payload = provenance(
        db,
        model="diffuse layer" if edl else "no_edl",
        n_strong=SITES_STRONG,
        n_weak=SITES_WEAK,
        zn_total=ZN_TOTAL,
        area=AREA_PER_GRAM * MASS_SOLID * 1000.0,
    )
    payload.update({f"logK_{k}": v for k, v in logks.items()})

    edl_line = "" if edl else "    -no_edl\n"
    points = []
    for ph in PH_VALUES:
        script = f"""
PHASES
Fix_H+
    H+ = H+
    log_k    0.0
HCl
    HCl = H+ + Cl-
    log_k    7.0
END
SOLUTION 1
    units    mol/kgw
    temp     25.0
    water    1.0
    pH       7.0
    Na       0.01
    Cl       0.01 charge
    Zn       {ZN_TOTAL}
SURFACE 1
    Hfo_sOH  {SITES_STRONG}  {AREA_PER_GRAM}  {MASS_SOLID * 1000.0}
    Hfo_wOH  {SITES_WEAK}
{edl_line}EQUILIBRIUM_PHASES 1
    Fix_H+   {-ph}  HCl  10.0
SELECTED_OUTPUT
    -reset      false
    -high_precision true
    -pH         true
    -ionic_strength true
    -activities H+ Zn+2
    -molalities {' '.join(REPORTED)}
END
"""
        ip.run_string(script)
        rows = ip.get_selected_output_array()
        header, values = rows[0], rows[-1]
        got_ph = values[header.index("pH")]
        if abs(got_ph - ph) > 1.0e-6:
            raise SystemExit(f"PHREEQC returned pH {got_ph} for a requested {ph}")
        col = lambda name: values[header.index(name)]
        m = {nm: col(f"m_{nm}(mol/kgw)") for nm in REPORTED}
        zn_sorbed = m["Hfo_sOZn+"] + m["Hfo_wOZn+"]
        points.append({
            "pH": ph, "la_H": col("la_H+"), "la_Zn": col("la_Zn+2"),
            "I": col("mu"),
            "s_free": m["Hfo_sOH"], "s_prot": m["Hfo_sOH2+"],
            "s_depr": m["Hfo_sO-"], "s_zn": m["Hfo_sOZn+"],
            "w_free": m["Hfo_wOH"], "w_prot": m["Hfo_wOH2+"],
            "w_depr": m["Hfo_wO-"], "w_zn": m["Hfo_wOZn+"],
            "zn_free": m["Zn+2"], "zn_sorbed_fraction": zn_sorbed / ZN_TOTAL,
        })
    payload["points"] = points
    emit("phreeqc_hfo_zn" + ("_ddl" if edl else ""), payload)


def run_protolysis():
    """The acid-base case: one family of weak sites, nothing else."""
    db = database_path(DATABASE)
    ip = VIPhreeqc()
    ip.load_database(db)
    if ip.phc_database_error_count:
        raise SystemExit(f"{DATABASE} failed to load: {ip.get_error_string()}")

    logks = {
        "protonation": logk_for(db, "Hfo_wOH + H+ = Hfo_wOH2+"),
        "deprotonation": logk_for(db, "Hfo_wOH = Hfo_wO- + H+"),
    }
    payload = provenance(
        db,
        model="no_edl, weak sites only, no metal",
        logK_protonation=logks["protonation"],
        logK_deprotonation=logks["deprotonation"],
        n_sites=SITES_PROTOLYSIS,
    )
    points = []
    for ph in PH_VALUES:
        script = f"""
PHASES
Fix_H+
    H+ = H+
    log_k    0.0
HCl
    HCl = H+ + Cl-
    log_k    7.0
END
SOLUTION 1
    units    mol/kgw
    temp     25.0
    water    1.0
    pH       7.0
    Na       0.01
    Cl       0.01 charge
SURFACE 1
    Hfo_wOH  {SITES_PROTOLYSIS}  {AREA_PER_GRAM}  {MASS_SOLID * 1000.0}
    -no_edl
EQUILIBRIUM_PHASES 1
    Fix_H+   {-ph}  HCl  10.0
SELECTED_OUTPUT
    -reset      false
    -high_precision true
    -pH         true
    -activities H+
    -molalities {' '.join(PROTOLYSIS_REPORTED)}
END
"""
        ip.run_string(script)
        rows = ip.get_selected_output_array()
        header, values = rows[0], rows[-1]
        got_ph = values[header.index("pH")]
        if abs(got_ph - ph) > 1.0e-6:
            raise SystemExit(f"PHREEQC returned pH {got_ph} for a requested {ph}")
        la_h = values[header.index("la_H+")]
        amounts = [values[header.index(f"m_{nm}(mol/kgw)")] for nm in PROTOLYSIS_REPORTED]
        total = sum(amounts)
        fractions = [a / total for a in amounts]
        points.append({
            "pH": ph, "la_H": la_h,
            "free": fractions[0], "protonated": fractions[1],
            "deprotonated": fractions[2],
        })

    payload["points"] = points
    emit("phreeqc_protolysis", payload)


IONIC_STRENGTHS = [0.1, 0.01, 0.001]


def run_diffuse_layer():
    """The acid-base case again, with the diffuse layer PHREEQC uses by default.

    Three background electrolyte levels rather than one, because the whole point
    of a diffuse layer is that the screening depends on the ionic strength: a
    model that got the charge right and the screening wrong would reproduce one
    series and miss the other two. The ionic strength PHREEQC computed is
    emitted with each point, so the comparison can be made at a matched one
    instead of assuming the two codes speciate the background identically.
    """
    db = database_path(DATABASE)
    ip = VIPhreeqc()
    ip.load_database(db)
    if ip.phc_database_error_count:
        raise SystemExit(f"{DATABASE} failed to load: {ip.get_error_string()}")

    logks = {
        "protonation": logk_for(db, "Hfo_wOH + H+ = Hfo_wOH2+"),
        "deprotonation": logk_for(db, "Hfo_wOH = Hfo_wO- + H+"),
    }
    area = AREA_PER_GRAM * MASS_SOLID * 1000.0
    payload = provenance(
        db,
        model="diffuse layer (PHREEQC's default SURFACE), weak sites, no metal",
        logK_protonation=logks["protonation"],
        logK_deprotonation=logks["deprotonation"],
        n_sites=SITES_PROTOLYSIS,
        area=area,
        site_density_umol_per_m2=SITES_PROTOLYSIS / area * 1e6,
    )
    series = []
    for ionic in IONIC_STRENGTHS:
        points = []
        for ph in PH_VALUES:
            script = f"""
PHASES
Fix_H+
    H+ = H+
    log_k    0.0
HCl
    HCl = H+ + Cl-
    log_k    7.0
END
SOLUTION 1
    units    mol/kgw
    temp     25.0
    water    1.0
    pH       7.0
    Na       {ionic}
    Cl       {ionic} charge
SURFACE 1
    Hfo_wOH  {SITES_PROTOLYSIS}  {AREA_PER_GRAM}  {MASS_SOLID * 1000.0}
EQUILIBRIUM_PHASES 1
    Fix_H+   {-ph}  HCl  10.0
SELECTED_OUTPUT
    -reset      false
    -high_precision true
    -pH         true
    -ionic_strength true
    -activities H+
    -molalities {' '.join(PROTOLYSIS_REPORTED)}
END
"""
            ip.run_string(script)
            rows = ip.get_selected_output_array()
            header, values = rows[0], rows[-1]
            got_ph = values[header.index("pH")]
            if abs(got_ph - ph) > 1.0e-6:
                raise SystemExit(f"PHREEQC returned pH {got_ph} for a requested {ph}")
            la_h = values[header.index("la_H+")]
            ionic_str = values[header.index("mu")]
            amounts = [values[header.index(f"m_{nm}(mol/kgw)")]
                       for nm in PROTOLYSIS_REPORTED]
            total = sum(amounts)
            fractions = [a / total for a in amounts]
            points.append({
                "pH": ph, "la_H": la_h, "I": ionic_str,
                "free": fractions[0], "protonated": fractions[1],
                "deprotonated": fractions[2],
            })
        series.append({"nacl": ionic, "points": points})

    payload["series"] = series
    emit("phreeqc_diffuse_layer", payload)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case",
        choices=("protolysis", "protolysis-ddl", "zn-edge"),
        default="protolysis",
        help="which comparison to emit (default: the acid-base one)",
    )
    parser.add_argument(
        "--edl",
        action="store_true",
        help="use the diffuse double layer instead of -no_edl (a different model)",
    )
    args = parser.parse_args()
    if args.case == "protolysis":
        run_protolysis()
    elif args.case == "protolysis-ddl":
        run_diffuse_layer()
    else:
        run_zn_edge(args.edl)
