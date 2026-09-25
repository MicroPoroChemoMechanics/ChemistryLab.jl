# Oracle for the C-S-H surface benchmark: silanol sites with a diffuse layer,
# in solutions of NaOH, CaCl2 and NaCl.
#
#   conda run -n mpcm-oracles python test/reference/phreeqc_csh_surface.py
#   conda run -n mpcm-oracles python test/reference/phreeqc_csh_surface.py --case paste
#
# Writes test/reference/phreeqc_csh_surface.json (the surface in solutions of
# NaOH, CaCl2 and NaCl) and, with --case paste, phreeqc_csh_paste.json (the same
# surface in Guo's hydrated paste, with portlandite, the AFm and AFt phases,
# Friedel's and Kuzel's salts, swept in NaCl). test/csh_surface.jl reads both.
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
#   * the standard energies and the molar masses -- asked of ChemistryLab
#     itself (`library_values`): the ΔₐG⁰ of each species at 25 °C and 1 bar as
#     the package evaluates it from data/slop98-inorganic-thermofun.json and
#     data/cemdata18-thermofun.json, and the molar masses it computes from the
#     formulas when it creates the species. No energy, atomic mass or molar mass
#     is read here from a database file or written by hand;
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
import subprocess
import sys
import tempfile

import phreeqpython
from phreeqpython.viphreeqc import VIPhreeqc

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(os.path.dirname(_HERE))
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


# What the package computes, asked of the package. Each line it prints is
# tab-separated: `R`, an `atom` and its molar mass, or a `species` of a database
# with its standard Gibbs energy (J/mol) and its molar mass (g/mol).
_LIBRARY_SCRIPT = r"""
using ChemistryLab, DynamicQuantities
T, P = 298.15u"K", 1.0e5u"Pa"
println("R\t", R_GAS)
for el in split(ARGS[1], ",")
    M = calculate_molar_mass(Dict(Symbol(el) => 1))
    println("atom\t", el, "\t", ustrip(us"g/mol", M))
end
for arg in ARGS[2:end]
    db, list = split(arg, ":"; limit = 2)
    subs = Dict(symbol(s) => s for s in build_species(datapath(db); verbose = false))
    for sym in split(list, ",")
        s = subs[sym]
        G = ustrip(us"J/mol", s[:ΔₐG⁰](T = T, P = P; unit = true))
        M = haskey(s, :M) ? ustrip(us"g/mol", s[:M]) : NaN
        println("species\t", db, "\t", sym, "\t", G, "\t", M)
    end
end
"""

SLOP98 = "slop98-inorganic-thermofun.json"
CEMDATA18 = "cemdata18-thermofun.json"


def library_values(elements, species):
    """R, atomic masses and species data as ChemistryLab computes them.

    `species` maps a database file to the symbols wanted from it. The package
    is the one this file sits in, unless CHEMISTRYLAB_PROJECT names another
    Julia environment that loads it."""
    project = os.environ.get("CHEMISTRYLAB_PROJECT", _ROOT)
    args = [",".join(elements)] + [f"{db}:{','.join(syms)}" for db, syms in species.items()]
    out = subprocess.run(
        ["julia", f"--project={project}", "--startup-file=no", "-e", _LIBRARY_SCRIPT, *args],
        check=True, capture_output=True, text=True,
    ).stdout
    lib = {"R": None, "atoms": {}, "G": {}, "M": {}}
    for line in out.splitlines():
        f = line.split("\t")
        if f[0] == "R":
            lib["R"] = float(f[1])
        elif f[0] == "atom":
            lib["atoms"][f[1]] = float(f[2])
        elif f[0] == "species":
            lib["G"][(f[1], f[2])] = float(f[3])
            lib["M"][(f[1], f[2])] = float(f[4])
    return lib


def _redox_log_k():
    """The two redox couples PHREEQC needs to define O(0) and H(0), from the
    phreeqc.dat committed next to this file. They set the pe of these
    solutions, which nothing compared here depends on."""
    with open(os.path.join(_HERE, "phreeqc.dat"), encoding="utf-8") as f:
        lines = f.read().splitlines()
    out = {}
    for eq in ("2 H2O = O2 + 4 H+ + 4 e-", "2 H+ + 2 e- = H2"):
        k = lines.index(eq)
        out[eq] = float(lines[k + 1].split()[1])
    return out


def log_kw(lib, db):
    g = lambda sp: lib["G"][(db, sp)]
    drg = g("OH-") + g("H+") - g("H2O@")
    return -drg / (lib["R"] * T_K * math.log(10))


def database(kw, lib, extra=()):
    """The minimal database: free ions, water, the surface. `extra` adds
    elements to the master species as (element, master species) pairs; the gram
    formula weights are the library's atomic masses, which PHREEQC uses only to
    convert mass units no input here is given in."""
    a = lib["atoms"]
    master = [("Na", "Na+"), ("Ca", "Ca+2"), ("Cl", "Cl-"), *extra]
    redox = _redox_log_k()
    return f"""SOLUTION_MASTER_SPECIES
H        H+     -1.0  H        {a["H"]}
H(0)     H2      0.0  H
H(1)     H+     -1.0  0.0
E        e-      0.0  0.0      0.0
O        H2O     0.0  O        {a["O"]}
O(0)     O2      0.0  O
O(-2)    H2O     0.0  0.0
""" + "".join(f"{el:<8} {sp:<6}  0.0  {el:<7} {a[el]}\n" for el, sp in master) + \
        "".join(f"S(6)     SO4-2   0.0  S\n" for el, sp in extra if el == "S") + f"""SOLUTION_SPECIES
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
""" + "".join(f"{sp} = {sp}\n    log_k 0.0\n" for el, sp in extra) + f"""H2O = OH- + H+
    log_k {kw:.10f}
2 H2O = O2 + 4 H+ + 4 e-
    log_k {redox["2 H2O = O2 + 4 H+ + 4 e-"]}
2 H+ + 2 e- = H2
    log_k {redox["2 H+ + 2 e- = H2"]}
SURFACE_MASTER_SPECIES
Csh_w    Csh_wOH
SURFACE_SPECIES
Csh_wOH = Csh_wOH
    log_k 0.0
""" + "".join(f"{eq}\n    log_k {lk}\n" for eq, lk in REACTIONS)


REPORTED = ["Csh_wOH", "Csh_wO-", "Csh_wOCa+", "Csh_wOHCl-", "Csh_wONa"]


WATER = ("H2O@", "H+", "OH-")


def main():
    lib = library_values(("H", "O", "Na", "Ca", "Cl"), {SLOP98: WATER})
    kw = log_kw(lib, SLOP98)
    text = database(kw, lib)
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
        "inputs": "standard Gibbs energies and atomic masses from ChemistryLab (library_values)",
        "log_kw": kw,
        "slop98_gibbs_J_per_mol": {sp: lib["G"][(SLOP98, sp)] for sp in WATER},
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


# ── the paste: the surface with the hydrates of Guo's inventory ───────────────

# Dissolution over the free ions, balanced by hand and checked on the Julia
# side; the log K comes from the standard Gibbs energies ChemistryLab computes
# from data/cemdata18-thermofun.json, the database it builds these solids from.
PHASES = {
    "Portlandite": {"Ca+2": 1, "OH-": 2},
    "monosulphate12": {"Ca+2": 4, "AlO2-": 2, "SO4-2": 1, "OH-": 4, "H2O@": 10},
    "monosulphate14": {"Ca+2": 4, "AlO2-": 2, "SO4-2": 1, "OH-": 4, "H2O@": 12},
    "ettringite": {"Ca+2": 6, "AlO2-": 2, "SO4-2": 3, "OH-": 4, "H2O@": 30},
    "C4AClH10": {"Ca+2": 4, "AlO2-": 2, "Cl-": 2, "OH-": 4, "H2O@": 8},
    "C4AsClH12": {"Ca+2": 4, "AlO2-": 2, "Cl-": 1, "SO4-2": 0.5, "OH-": 4, "H2O@": 10},
}
PASTE_IONS = ("H2O@", "H+", "OH-", "Ca+2", "AlO2-", "SO4-2", "Cl-")
PHREEQC_NAME = {"H2O@": "H2O"}
# Guo's inventory per liter of concrete, brought to one kilogram of pore water.
PORE_G = 10 * float(_GUO["quantities"]["porosity_percent"]["value"])
SCALE = 1000.0 / PORE_G


def _grams(q):
    return float(_GUO["quantities"][q]["value"]) * SCALE


def phase_log_k(lib, phase):
    g = lambda sp: lib["G"][(CEMDATA18, sp)]
    drg = sum(nu * g(sp) for sp, nu in PHASES[phase].items()) - g(phase)
    return -drg / (lib["R"] * T_K * math.log(10))


def paste_database(lib):
    kw = log_kw(lib, CEMDATA18)
    extra = (("Al", "AlO2-"), ("S", "SO4-2"))
    full = database(kw, lib, extra)
    base, surface = full.split("SURFACE_MASTER_SPECIES")
    phases = "PHASES\n"
    for ph, prod in PHASES.items():
        # PHREEQC checks each phase's balance against its formula; the formula
        # written is the sum of the products, so the check is of the reaction.
        formula = _formula_of(prod)
        rhs = " + ".join(f"{nu} {PHREEQC_NAME.get(sp, sp)}" for sp, nu in prod.items())
        phases += f"{ph}\n    {formula} = {rhs}\n    log_k {phase_log_k(lib, ph):.10f}\n"
    return base + phases + "SURFACE_MASTER_SPECIES" + surface, kw


_ION_ATOMS = {
    "Ca+2": {"Ca": 1}, "AlO2-": {"Al": 1, "O": 2}, "SO4-2": {"S": 1, "O": 4},
    "OH-": {"O": 1, "H": 1}, "Cl-": {"Cl": 1}, "H2O@": {"H": 2, "O": 1},
}


def _formula_of(prod):
    atoms = {}
    for sp, nu in prod.items():
        for el, k in _ION_ATOMS[sp].items():
            atoms[el] = atoms.get(el, 0) + nu * k
    fmt = lambda x: str(int(x)) if float(x).is_integer() else str(x)
    return "".join(f"{el}{fmt(atoms[el])}" for el in ("Ca", "Al", "S", "Cl", "O", "H") if el in atoms)


NACL = [0.0, 0.02, 0.05, 0.1, 0.2, 0.4]
PASTE_REPORTED = REPORTED


def run_paste():
    lib = library_values(
        ("H", "O", "Na", "Ca", "Cl", "Al", "S"), {CEMDATA18: [*PASTE_IONS, *PHASES]},
    )
    molar_mass = {ph: lib["M"][(CEMDATA18, ph)] for ph in PHASES}
    text, kw = paste_database(lib)
    with tempfile.NamedTemporaryFile("w", suffix=".dat", delete=False) as handle:
        handle.write(text)
        dbfile = handle.name
    ip = VIPhreeqc()
    ip.load_database(dbfile)
    if ip.phc_database_error_count:
        raise SystemExit(f"the generated database failed to load: {ip.get_error_string()}")

    csh_g = _grams("csh_per_liter")
    n_sites = SITE_DENSITY * csh_g
    initial = {
        "Portlandite": _grams("ch_per_liter"),
        "monosulphate12": _grams("afm_per_liter"),
        "ettringite": _grams("aft_per_liter"),
    }
    # Amounts from Guo's grams and the molar masses ChemistryLab gives these
    # solids. The Julia side recomputes them from the same grams and asserts
    # that they are these.
    points = []
    for nacl in NACL:
        eq_phases = "\n".join(
            f"    {ph}  0.0  {initial.get(ph, 0.0) / molar_mass[ph]}" for ph in PHASES
        )
        script = f"""
SOLUTION 1
    units    mol/kgw
    temp     25.0
    water    1.0
    pH       7.0 charge
    Na       {nacl if nacl > 0 else 1.0e-10}
    Cl       {nacl if nacl > 0 else 1.0e-10}
EQUILIBRIUM_PHASES 1
{eq_phases}
SURFACE 1
    Csh_wOH  {n_sites}  {AREA_PER_GRAM}  {csh_g}
USER_PUNCH
    -headings psi_V sigma water_kg
    10 PUNCH EDL("psi", "Csh"), EDL("sigma", "Csh"), TOT("water")
SELECTED_OUTPUT
    -reset      false
    -high_precision true
    -pH         true
    -ionic_strength true
    -totals     Ca Cl Na Al S
    -equilibrium_phases {' '.join(PHASES)}
    -molalities {' '.join(PASTE_REPORTED)}
END
"""
        ip.run_string(script)
        if ip.get_error_string():
            raise SystemExit(ip.get_error_string())
        rows = ip.get_selected_output_array()
        header, values = rows[0], rows[-1]
        col = lambda name: values[header.index(name)]
        # Molalities are per kilogram of the water present, and the hydrates
        # and the surface both trade water with the solution: amounts in moles
        # are the molalities times that mass.
        w = col("water_kg")
        points.append({
            "nacl": nacl, "pH": col("pH"), "I": col("mu"),
            "psi_V": col("psi_V"), "water_kg": w,
            "totals_mol": {el: col(f"{el}(mol/kgw)") * w for el in ("Ca", "Cl", "Na", "Al", "S")},
            "phases": {ph: col(ph) for ph in PHASES},
            "surface": {nm: col(f"m_{nm}(mol/kgw)") * w for nm in PASTE_REPORTED},
        })
    os.unlink(dbfile)
    payload = {
        "generator": "test/reference/phreeqc_csh_surface.py --case paste",
        "python": sys.version.split()[0],
        "database_md5": hashlib.md5(text.encode()).hexdigest(),
        "model": "PHREEQC default SURFACE (Dzombak-Morel diffuse layer), Davies activity, EQUILIBRIUM_PHASES, closed system",
        "inputs": "standard Gibbs energies, molar and atomic masses from ChemistryLab (library_values)",
        "log_kw": kw,
        "phase_log_k": {ph: phase_log_k(lib, ph) for ph in PHASES},
        "phase_reactions": PHASES,
        "grams_per_kg_water": {ph: initial.get(ph, 0.0) for ph in PHASES},
        "molar_mass_g_per_mol": molar_mass,
        "initial_mol": {ph: initial.get(ph, 0.0) / molar_mass[ph] for ph in PHASES},
        "n_sites": n_sites, "area_m2": AREA_PER_GRAM * csh_g,
        "points": points,
    }
    path = os.path.join(_HERE, "phreeqc_csh_paste.json")
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")
    print(f"wrote {os.path.relpath(path, _ROOT)} ({len(points)} points)")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=("surface", "paste"), default="surface")
    args = parser.parse_args()
    run_paste() if args.case == "paste" else main()
