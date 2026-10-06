# Oracle for the surfaces of three charge planes: PHREEQC's CD-MUSIC surface.
#
# Goethite as Hiemstra and Van Riemsdijk (1996) describe it: two kinds of
# surface groups of the 110 face, singly and triply coordinated, each carrying
# the charge -1/2 of the MUSIC model and protonating once, with ion pairs of
# the background electrolyte. Three cases, one per way of placing the charges:
#
#   basic_stern            Fig. 6 of the paper: one capacitance, the ion pairs
#                          where the diffuse layer begins. PHREEQC has three
#                          planes always; a second capacitance of 1e5 F/m2
#                          merges the outer two, as its manual suggests.
#   three_plane_phosphate  Table 2, Case III: the solution-oriented ligands of a
#                          monodentate phosphate complex in the 1-plane, its
#                          charge shared with the 0-plane (f = 0.25), the ion
#                          pairs in the d-plane, a finite second capacitance.
#   triple_layer           the placement of Davis et al. (1978): the ion pairs
#                          in the 1-plane, between their two capacitances.
#
#   conda run -n mpcm-oracles python test/reference/phreeqc_cd_music.py
#
# Writes test/reference/phreeqc_cd_music.json, which test/charge_planes.jl reads.
#
# What has to match for the comparison to mean anything:
#
#   * every constant -- the surface ones from data/literature/Hiemstra1996.json
#     and Davis1978.json, which the Julia side reads too; the aqueous ones
#     (water, phosphate) from phreeqc.dat, recorded in the fixture, from which
#     the Julia side builds its species;
#   * the activity coefficients: one on both sides. Each charged species is
#     redefined here with a Debye-Hueckel ion size of 1e9 angstrom and no
#     linear term, so that log gamma vanishes to 1e-9, and the sodium complexes
#     of phreeqc.dat are switched off. What is compared is the surface model,
#     not two aqueous models;
#   * NaCl as the background electrolyte where the paper has NaNO3, with the
#     same ion-pair constants: PHREEQC would otherwise let nitrate reduce at
#     the default pe.

import hashlib
import json
import os
import sys

from phreeqpython.viphreeqc import VIPhreeqc

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
DATABASE = os.path.join(HERE, "phreeqc.dat")


def literature(key):
    with open(os.path.join(ROOT, "data", "literature", f"{key}.json"), encoding="utf-8") as f:
        return {k: q["value"] for k, q in json.load(f)["quantities"].items()}


HV = literature("Hiemstra1996")
DLJ = literature("Davis1978")
MASS = 1.0                                  # g of goethite in the kilogram of water
PH_VALUES = [4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0]
P_TOTAL = 1.0e-4                            # mol/kgw, a sixth of the singly coordinated groups

UNI = ["Goe_uniOH-0.5", "Goe_uniOH2+0.5", "Goe_uniOHNa+0.5", "Goe_uniOH2Cl-0.5"]
TRI = ["Goe_triO-0.5", "Goe_triOH+0.5", "Goe_triONa+0.5", "Goe_triOHCl-0.5"]
PHOSPHATE = "Goe_uniOPO3-2.5"


def logk_for(reaction):
    """The log K of one reaction of phreeqc.dat, read rather than typed."""
    def norm(text):
        return " ".join(text.replace("\t", " ").split())

    lines = open(DATABASE, encoding="utf-8", errors="ignore").read().split("\n")
    for i, line in enumerate(lines):
        if norm(line) == norm(reaction):
            for follow in lines[i + 1:i + 4]:
                if "log_k" in follow:
                    return float(follow.split("log_k")[1].split("#")[0])
    raise SystemExit(f"no log K for {reaction!r} in phreeqc.dat")


def logk_water():
    """The ionization constant of water at 25 °C, which phreeqc.dat gives by an
    analytic expression in T: asked of PHREEQC rather than evaluated here."""
    ip = VIPhreeqc()
    ip.load_database(DATABASE)
    ip.run_string("""
SOLUTION 1
    temp 25.0
SELECTED_OUTPUT
    -reset false
USER_PUNCH
    -headings lk
    10 PUNCH LK_SPECIES("OH-")
END
""")
    rows = ip.get_selected_output_array()
    return rows[-1][rows[0].index("lk")]


LOGK_W = logk_water()
LOGK_P = {
    "HPO4-2": logk_for("PO4-3 + H+ = HPO4-2"),
    "H2PO4-": logk_for("PO4-3 + 2 H+ = H2PO4-"),
    "H3PO4": logk_for("PO4-3 + 3H+ = H3PO4"),
}

CASES = [
    {
        "name": "basic_stern", "C1": HV["stern_capacitance"], "C2": 1.0e5,
        "log_k_H": HV["log_K_protonation_charging"], "pair_plane": 2,
        "nacl": [0.001, 0.01, 0.1], "phosphate": None,
    },
    {
        "name": "three_plane_phosphate", "C1": HV["stern_capacitance"],
        "C2": HV["outer_capacitance_three_planes"],
        "log_k_H": HV["log_K_protonation_phosphate_sets"], "pair_plane": 2,
        "nacl": [0.01, 0.1],
        "phosphate": {
            "total": P_TOTAL, "log_k": HV["log_K_monodentate_phosphate"],
            "dz": [HV["delta_z0_monodentate_phosphate"], HV["delta_z1_monodentate_phosphate"], 0.0],
        },
    },
    {
        "name": "triple_layer", "C1": DLJ["inner_capacitance"], "C2": DLJ["outer_capacitance"],
        "log_k_H": HV["log_K_protonation_charging"], "pair_plane": 1,
        "nacl": [0.01, 0.1], "phosphate": None,
    },
]


def definitions(case):
    """The species blocks: unit activity coefficients, the goethite surface."""
    pair = [0, 0, 0]
    pair[case["pair_plane"]] = 1
    cation = " ".join(str(x) for x in pair)
    anion = " ".join(str(x) for x in [1 - pair[0], -pair[1], -pair[2]])
    kH, kp = case["log_k_H"], HV["log_K_ion_pair"]
    unit = "    -gamma 1e9 0\n"
    text = "SOLUTION_SPECIES\n"
    for master in ("H+ = H+", "Na+ = Na+", "Cl- = Cl-", "PO4-3 = PO4-3"):
        text += f"{master}\n    log_k 0\n{unit}"
    text += f"H2O = OH- + H+\n    log_k {LOGK_W}\n{unit}"
    for sp, k in LOGK_P.items():
        n = {"HPO4-2": 1, "H2PO4-": 2, "H3PO4": 3}[sp]
        text += f"PO4-3 + {n} H+ = {sp}\n    log_k {k}\n{unit}"
    for off in ("Na+ + OH- = NaOH", "Na+ + HPO4-2 = NaHPO4-"):
        text += f"{off}\n    log_k -50\n"
    text += "SURFACE_MASTER_SPECIES\n    Goe_uni Goe_uniOH-0.5\n    Goe_tri Goe_triO-0.5\n"
    text += "SURFACE_SPECIES\n"
    text += "Goe_uniOH-0.5 = Goe_uniOH-0.5\n    log_k 0\n"
    text += f"Goe_uniOH-0.5 + H+ = Goe_uniOH2+0.5\n    log_k {kH}\n    -cd_music 1 0 0\n"
    text += f"Goe_uniOH-0.5 + Na+ = Goe_uniOHNa+0.5\n    log_k {kp}\n    -cd_music {cation}\n"
    text += f"Goe_uniOH-0.5 + H+ + Cl- = Goe_uniOH2Cl-0.5\n    log_k {kH + kp}\n    -cd_music {anion}\n"
    text += "Goe_triO-0.5 = Goe_triO-0.5\n    log_k 0\n"
    text += f"Goe_triO-0.5 + H+ = Goe_triOH+0.5\n    log_k {kH}\n    -cd_music 1 0 0\n"
    text += f"Goe_triO-0.5 + Na+ = Goe_triONa+0.5\n    log_k {kp}\n    -cd_music {cation}\n"
    text += f"Goe_triO-0.5 + H+ + Cl- = Goe_triOHCl-0.5\n    log_k {kH + kp}\n    -cd_music {anion}\n"
    if case["phosphate"]:
        ph = case["phosphate"]
        dz = " ".join(str(x) for x in ph["dz"])
        text += f"Goe_uniOH-0.5 + H+ + PO4-3 = {PHOSPHATE} + H2O\n    log_k {ph['log_k']}\n    -cd_music {dz}\n"
    return text


def run_point(ip, case, nacl, ph):
    phosphate = case["phosphate"]
    titrant = "HCl" if ph < 7.0 else "NaOH"
    surface = UNI + TRI + ([PHOSPHATE] if phosphate else [])
    script = definitions(case) + f"""
PHASES
Fix_H+
    H+ = H+
    log_k 0.0
END
SOLUTION 1
    units mol/kgw
    temp 25.0
    water 1.0
    pH 7.0
    Na {nacl}
    Cl {nacl} charge
{f"    P {phosphate['total']}" if phosphate else ""}
SURFACE 1
    -sites_units density
    -cd_music
    Goe_uniOH {HV['singly_coordinated_site_density']} {HV['specific_surface_area']} {MASS}
    Goe_triO {HV['triply_coordinated_site_density']}
    -capacitances {case['C1']} {case['C2']}
EQUILIBRIUM_PHASES 1
    Fix_H+ {-ph} {titrant} 10.0
SELECTED_OUTPUT
    -reset false
    -high_precision true
    -pH true
    -ionic_strength true
    -water true
    -activities H+
    -molalities {' '.join(surface)}
    -totals Na Cl P
USER_PUNCH
    -headings psi0 psi1 psi2 sigma0 sigma1 sigma2 gamma_Na
    10 PUNCH EDL("psi", "Goe"), EDL("psi1", "Goe"), EDL("psi2", "Goe")
    20 PUNCH EDL("sigma", "Goe"), EDL("sigma1", "Goe"), EDL("sigma2", "Goe"), GAMMA("Na+")
END
"""
    ip.run_string(script)
    rows = ip.get_selected_output_array()
    header, values = rows[0], rows[-1]
    col = lambda name: values[header.index(name)]
    if abs(col("pH") - ph) > 1.0e-6:
        raise SystemExit(f"PHREEQC returned pH {col('pH')} for a requested {ph}")
    if abs(col("gamma_Na") - 1.0) > 1.0e-6:
        raise SystemExit(f"the activity coefficient of Na+ is {col('gamma_Na')}, not one")
    kgw = col("mass_H2O")
    # Amounts in moles, the kilogram of water having moved by the water the
    # reactions released or took.
    amount = lambda name: col(f"m_{name}(mol/kgw)") * kgw
    return {
        "pH": ph, "la_H": col("la_H+"), "I": col("mu"), "kgw": kgw,
        "Na": col("Na(mol/kgw)") * kgw, "Cl": col("Cl(mol/kgw)") * kgw,
        "P": col("P(mol/kgw)") * kgw,
        "uni": [amount(s) for s in UNI] + ([amount(PHOSPHATE)] if phosphate else []),
        "tri": [amount(s) for s in TRI],
        "psi": [col("psi0"), col("psi1"), col("psi2")],
        "sigma": [col("sigma0"), col("sigma1"), col("sigma2")],
    }


def main():
    ip = VIPhreeqc()
    ip.load_database(DATABASE)
    if ip.phc_database_error_count:
        raise SystemExit(f"phreeqc.dat failed to load: {ip.get_error_string()}")
    with open(DATABASE, "rb") as handle:
        digest = hashlib.md5(handle.read()).hexdigest()
    try:
        from importlib.metadata import version
        pp_version = version("phreeqpython")
    except Exception:  # pragma: no cover - provenance must never fail the run
        pp_version = "unknown"
    payload = {
        "generator": "test/reference/phreeqc_cd_music.py",
        "python": sys.version.split()[0],
        "phreeqpython": pp_version,
        "database": "phreeqc.dat",
        "database_md5": digest,
        "log_k_water": LOGK_W,
        "log_k_phosphate": LOGK_P,
        "log_k_ion_pair": HV["log_K_ion_pair"],
        "area": HV["specific_surface_area"] * MASS,
        "uni_species": UNI + [PHOSPHATE],
        "tri_species": TRI,
        "cases": [],
    }
    for case in CASES:
        out = {k: v for k, v in case.items() if k != "nacl"}
        out["series"] = [
            {"nacl": nacl, "points": [run_point(ip, case, nacl, ph) for ph in PH_VALUES]}
            for nacl in case["nacl"]
        ]
        payload["cases"].append(out)
    path = os.path.join(HERE, "phreeqc_cd_music.json")
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")
    print(f"wrote {os.path.relpath(path, ROOT)}")


if __name__ == "__main__":
    main()
