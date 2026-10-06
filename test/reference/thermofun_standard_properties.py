# Oracle for test/water_eos_reference.jl.
#
# Standard properties of the solvent, of aqueous species (HKF), of minerals and of
# a gas, from 0 to 1000 °C and 1 to 5000 bar, as Reaktoro 2.13 computes them with
# the ThermoFun library it embeds, reading the *same* Cemdata18 ThermoFun file
# this package reads. The solvent is computed by the equation of state of water
# the record declares (`water_eos_hgk84_reaktoro`); Reaktoro's own embedded copy of
# Cemdata18 declares another one and must not be used.
#
#   conda create -n reaktoro-env -c conda-forge reaktoro thermofun
#   conda run -n reaktoro-env python test/reference/thermofun_standard_properties.py <path to cemdata18-thermofun.json>
#
# The library crashes on some states outside the domain of the HKF equations, so
# each species and temperature is computed in a child process, and the aqueous
# species only where the library's own water is dense enough for those
# equations (350 kg/m³, Shock et al. 1992), away from the near-critical states (above 360 °C and below 600
# bar) where the library's HKF departs from Reaktoro's own implementation of the
# same equations (at 400 °C and 500 bar it gives Ca2+ a volume of -706 cm³/mol,
# where Reaktoro's HKF and this package give +79).
#
# Writes `test/reference/thermofun_standard_properties.json`.

import hashlib
import json
import os
import subprocess
import sys

TEMPERATURES_C = (0, 25, 50, 100, 150, 200, 250, 300, 350, 400, 500, 600, 800, 1000)
PRESSURES_BAR = (1, 10, 100, 500, 1000, 2000, 5000)
SOLVENT = "H2O@"
AQUEOUS = ("CO2@", "Ca+2", "OH-", "HCO3-")
OTHERS = ("Portlandite", "Cal", "CO2")

CHILD = r"""
import sys, reaktoro as rkt
db = rkt.ThermoFunDatabase.fromFile(sys.argv[1])
species = db.species().get(sys.argv[2])
T = float(sys.argv[3])
for P in map(float, sys.argv[4:]):
    p = species.standardThermoProps(T, P)
    print("ROW", P, float(p.G0), float(p.H0), float(p.V0), float(p.Cp0), flush=True)
"""


def row(db, name, T, pressures):
    out = subprocess.run(
        [sys.executable, "-c", CHILD, db, name, repr(T)] + [repr(P) for P in pressures],
        capture_output=True, text=True,
    )
    values = {}
    # The library writes its warnings to the same output: only the rows count.
    for line in out.stdout.splitlines():
        if not line.startswith("ROW "):
            continue
        P, G, H, V, Cp = map(float, line.split()[1:])
        values[P] = [G, H, V, Cp]
    return values


def main():
    import reaktoro as rkt

    db = sys.argv[1]
    points = []
    # The density of water at each state, from the library's own solvent: the
    # aqueous species are computed only where it is at least 350 kg/m³.
    density = {}
    for name in (SOLVENT,) + AQUEOUS + OTHERS:
        for TC in TEMPERATURES_C:
            T = TC + 273.15
            pressures = []
            for Pb in PRESSURES_BAR:
                P = Pb * 1e5
                if name in AQUEOUS:
                    if density.get((T, P), 0.0) < 350 or (TC > 360 and Pb < 600):
                        continue
                pressures.append(P)
            for P, (G, H, V, Cp) in row(db, name, T, pressures).items():
                points.append({"species": name, "T": T, "P": P, "G": G, "H": H, "V": V, "Cp": Cp})
                if name == SOLVENT:
                    density[(T, P)] = rkt.waterMolarMass / V
    with open(db, "rb") as handle:
        md5 = hashlib.md5(handle.read()).hexdigest()
    payload = {
        "generator": "test/reference/thermofun_standard_properties.py",
        "python": sys.version.split()[0],
        "reaktoro": rkt.__version__,
        "database": os.path.basename(db),
        "database_md5": md5,
        # The molar mass of water the library's solvent model multiplies the
        # specific properties of the equation of state by, in kg/mol.
        "water_molar_mass": rkt.waterMolarMass,
        "points": points,
    }
    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(here, "thermofun_standard_properties.json")
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=1)
        handle.write("\n")
    print(f"wrote {os.path.basename(path)} ({len(points)} points)")


if __name__ == "__main__":
    main()
