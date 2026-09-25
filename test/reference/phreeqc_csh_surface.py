# Oracle for the C-S-H surface benchmark: silanol sites with a diffuse layer,
# in solutions of NaOH, CaCl2 and NaCl.
#
#   conda run -n mpcm-oracles python test/reference/phreeqc_csh_surface.py
#
# Writes test/reference/phreeqc_csh_surface.json, which test/csh_surface.jl
# reads.
#
# One source of truth for both codes. The surface reactions and their constants
# are read from data/literature/Guo2018.json (the table surface_reactions_phreeqc,
# Guo's Table 1 in PHREEQC syntax), which the Julia side reads too through
# `site_family`; the site density and the specific area come from the same file.
#
# What is held identical, so that the comparison is of the surface and of
# nothing else:
#
#   * the aqueous species -- the free ions and the ion product of water, and no
#     ion pair. The database is written here rather than taken from phreeqc.dat,
#     whose ion pairs and WATEQ activity coefficients have no counterpart in the
#     package's species list;
#   * the activity model -- Davies on both sides: PHREEQC applies it to every
#     species that declares no -gamma, and none here does;
#   * the ion product of water -- its log K is computed from the standard Gibbs
#     energies of data/slop98-inorganic-thermofun.json, the database the Julia
#     side builds its species from;
#   * the diffuse layer -- PHREEQC's default SURFACE, Dzombak and Morel's model,
#     which `DiffuseLayer` reproduces (test/diffuse_layer.jl);
#   * the system -- closed, of the same total composition: the solution, its
#     pH set by its charge balance, and free sites, reacted together in one
#     step. The pH is not prescribed; it comes out of the balance, on both
#     sides.

import hashlib
import json
import math
import os
import sys
import tempfile

import phreeqpython
from phreeqpython.viphreeqc import VIPhreeqc

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(os.path.dirname(_HERE))
R_GAS = 8.31446261815324          # J/(mol K), CODATA 2018, as the package's R_GAS
T_K = 298.15

with open(os.path.join(_ROOT, "data", "literature", "Guo2018.json"), encoding="utf-8") as _f:
    _GUO = json.load(_f)
REACTIONS = [(r[0], float(r[1])) for r in _GUO["tables"]["surface_reactions_phreeqc"]["rows"]]
# The K row is left out: no potassium enters these solutions.
REACTIONS = [(eq, lk) for eq, lk in REACTIONS if "K+" not in eq]
SITE_DENSITY = float(_GUO["quantities"]["csh_site_density"]["value"])          # mol/g
AREA_PER_GRAM = float(_GUO["quantities"]["csh_specific_surface_area"]["value"])  # m2/g

CSH_GRAMS = 1.0                    # per kilogram of water
N_SITES = SITE_DENSITY * CSH_GRAMS

# (NaOH, CaCl2, NaCl), mol per kilogram of water.
GRID = [
    (naoh, cacl2, nacl)
    for naoh in (0.003, 0.03, 0.3)
    for cacl2 in (0.001, 0.01)
    for nacl in (0.01, 0.1, 0.3)
]


def slop98_log_kw():
    """log K of H2O = OH- + H+ from the standard Gibbs energies of slop98."""
    with open(os.path.join(_ROOT, "data", "slop98-inorganic-thermofun.json"), encoding="utf-8") as f:
        subs = json.load(f)["substances"]
    g = {}
    for s in subs:
        v = s.get("sm_gibbs_energy", {}).get("values")
        if v:
            g[s["symbol"]] = float(v[0])
    drg = g["OH-"] + g["H+"] - g["H2O@"]
    return -drg / (R_GAS * T_K * math.log(10)), {k: g[k] for k in ("OH-", "H+", "H2O@")}


def database(log_kw):
    return f"""SOLUTION_MASTER_SPECIES
H        H+     -1.0  H        1.008
H(0)     H2      0.0  H
H(1)     H+     -1.0  0.0
E        e-      0.0  0.0      0.0
O        H2O     0.0  O       16.00
O(0)     O2      0.0  O
O(-2)    H2O     0.0  0.0
Na       Na+     0.0  Na      22.9898
Ca       Ca+2    0.0  Ca      40.078
Cl       Cl-     0.0  Cl      35.453
SOLUTION_SPECIES
H+ = H+
    log_k 0.0
e- = e-
    log_k 0.0
H2O = H2O
    log_k 0.0
Na+ = Na+
    log_k 0.0
Ca+2 = Ca+2
    log_k 0.0
Cl- = Cl-
    log_k 0.0
H2O = OH- + H+
    log_k {log_kw:.10f}
2 H2O = O2 + 4 H+ + 4 e-
    log_k -86.08
2 H+ + 2 e- = H2
    log_k -3.15
SURFACE_MASTER_SPECIES
Csh_w    Csh_wOH
SURFACE_SPECIES
Csh_wOH = Csh_wOH
    log_k 0.0
""" + "".join(f"{eq}\n    log_k {lk}\n" for eq, lk in REACTIONS)


REPORTED = ["Csh_wOH", "Csh_wO-", "Csh_wOCa+", "Csh_wOHCl-", "Csh_wONa"]


def main():
    log_kw, energies = slop98_log_kw()
    text = database(log_kw)
    with tempfile.NamedTemporaryFile("w", suffix=".dat", delete=False) as handle:
        handle.write(text)
        dbfile = handle.name
    ip = VIPhreeqc()
    ip.load_database(dbfile)
    if ip.phc_database_error_count:
        raise SystemExit(f"the generated database failed to load: {ip.get_error_string()}")

    points = []
    for naoh, cacl2, nacl in GRID:
        script = f"""
SOLUTION 1
    units    mol/kgw
    temp     25.0
    water    1.0
    pH       12.0 charge
    Na       {naoh + nacl}
    Ca       {cacl2}
    Cl       {2 * cacl2 + nacl}
SURFACE 1
    Csh_wOH  {N_SITES}  {AREA_PER_GRAM}  {CSH_GRAMS}
USER_PUNCH
    -headings psi_V sigma
    10 PUNCH EDL("psi", "Csh"), EDL("sigma", "Csh")
SELECTED_OUTPUT
    -reset      false
    -high_precision true
    -pH         true
    -ionic_strength true
    -activities H+ Ca+2 Na+ Cl-
    -molalities {' '.join(REPORTED)}
END
"""
        ip.run_string(script)
        if ip.get_error_string():
            raise SystemExit(ip.get_error_string())
        rows = ip.get_selected_output_array()
        header, values = rows[0], rows[-1]
        col = lambda name: values[header.index(name)]
        amounts = {nm: col(f"m_{nm}(mol/kgw)") for nm in REPORTED}
        total = sum(amounts.values())
        points.append({
            "naoh": naoh, "cacl2": cacl2, "nacl": nacl,
            "pH": col("pH"), "I": col("mu"),
            "la_H": col("la_H+"), "la_Ca": col("la_Ca+2"),
            "la_Na": col("la_Na+"), "la_Cl": col("la_Cl-"),
            "psi_V": col("psi_V"), "sigma_C_m2": col("sigma"),
            "sites_total": total,
            "fractions": {nm: amounts[nm] / total for nm in REPORTED},
        })
    os.unlink(dbfile)

    try:
        from importlib.metadata import version
        pp_version = version("phreeqpython")
    except Exception:  # pragma: no cover - provenance must never fail the run
        pp_version = "unknown"
    payload = {
        "generator": "test/reference/phreeqc_csh_surface.py",
        "python": sys.version.split()[0],
        "phreeqpython": pp_version,
        "database": "written by the generator: free ions, the ion product of water, the surface of Guo2018.json",
        "database_md5": hashlib.md5(text.encode()).hexdigest(),
        "model": "PHREEQC default SURFACE (Dzombak-Morel diffuse layer), Davies activity, closed system",
        "log_kw": log_kw,
        "slop98_gibbs_J_per_mol": energies,
        "reactions": [{"equation": eq, "log_K": lk} for eq, lk in REACTIONS],
        "site_density_mol_per_g": SITE_DENSITY,
        "area_m2_per_g": AREA_PER_GRAM,
        "csh_grams_per_kg_water": CSH_GRAMS,
        "n_sites": N_SITES,
        "points": points,
    }
    path = os.path.join(_HERE, "phreeqc_csh_surface.json")
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")
    print(f"wrote {os.path.relpath(path, _ROOT)} ({len(points)} points)")


if __name__ == "__main__":
    main()
