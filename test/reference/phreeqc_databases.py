# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)
"""Reference values for the reader of PHREEQC databases (`read_phreeqc_database`).

Two cases for each database PHREEQC distributes, at 25 and 60 °C:

  - `speciation`: one kilogram of water with Na 0.5, K 0.01, Ca 0.01, Mg 0.02,
    Cl 0.525, S(6) 0.02 and C(4) 0.005 mol/kgw, the pH adjusted to balance the
    charge, which makes it the closed system of these elements whose solutes are
    Na+, K+, Ca+2, Mg+2, Cl-, SO4-2 and HCO3-;
  - `minerals`: one kilogram of pure water at equilibrium with calcite and
    gypsum, 10 mol of each.

For each, the molality of every aqueous species, the ionic strength, the pH, the
activity of water, the mass of water, and the Debye-Hückel A and B PHREEQC used
(`DH_A`, `DH_B`), so that a comparison can separate the gap these make from the
rest. For `minerals`, the moles of each mineral dissolved.

Run with the `mpcm-oracles` environment (phreeqpython, IPhreeqc 3.7.3):

    python test/reference/phreeqc_databases.py > test/reference/phreeqc_databases.json

The databases are those of `usgs_database.py`: no database is stored here.
"""

import json
import sys

import usgs_database
from phreeqpython.viphreeqc import VIPhreeqc

DATABASES = ["phreeqc.dat", "llnl.dat", "minteq.v4.dat", "wateq4f.dat", "sit.dat", "pitzer.dat"]
TEMPERATURES = [25.0, 60.0]

SPECIATION = """
SOLUTION 1
  temp {T}
  units mol/kgw
  pH 8 charge
  Na 0.5
  K 0.01
  Ca 0.01
  Mg 0.02
  Cl 0.525
  S(6) 0.02
  C(4) 0.005
"""

MINERALS = """
SOLUTION 1
  temp {T}
  units mol/kgw
  pH 7 charge
EQUILIBRIUM_PHASES 1
  Calcite 0 10
  Gypsum 0 10
"""

# Every aqueous species, by name and molality, then the scalars.
PUNCH = """
SELECTED_OUTPUT
  -reset false
USER_PUNCH
  -headings name molality
  10 n = SYS("aq", count, name$, type$, moles)
  20 FOR i = 1 TO count
  30   PUNCH name$(i), MOL(name$(i))
  40 NEXT i
  50 PUNCH "$mu", MU, "$pH", -LA("H+"), "$aw", ACT("H2O"), "$mass_water", TOT("water")
  60 PUNCH "$DH_A", DH_A, "$DH_B", DH_B
  70 PUNCH "$Calcite", EQUI_DELTA("Calcite"), "$Gypsum", EQUI_DELTA("Gypsum")
END
"""


def run(database, block, T):
    ip = VIPhreeqc()
    ip.load_database(usgs_database.database_path(database))
    ip.run_string(block.format(T=T) + PUNCH)
    error = ip.get_error_string()
    if error:
        raise SystemExit(f"{database} at {T} °C: {error}")
    row = ip.get_selected_output_array()[-1]
    values = dict(zip(row[0::2], row[1::2]))
    species = {k: v for k, v in values.items() if not k.startswith("$") and k != "H2O"}
    scalars = {k[1:]: v for k, v in values.items() if k.startswith("$")}
    return {"molality": species, **scalars}


def main():
    out = {
        "generator": "test/reference/phreeqc_databases.py",
        "engine": "IPhreeqc 3.7.3 (phreeqpython)",
        "databases": {},
    }
    for database in DATABASES:
        cases = {"speciation": {}, "minerals": {}}
        for T in TEMPERATURES:
            cases["speciation"][str(T)] = run(database, SPECIATION, T)
            minerals = run(database, MINERALS, T)
            # EQUI_DELTA is the change of a phase; a mineral dissolved is a loss.
            minerals["dissolved"] = {"Calcite": -minerals.pop("Calcite"), "Gypsum": -minerals.pop("Gypsum")}
            cases["minerals"][str(T)] = minerals
        out["databases"][database] = cases
    json.dump(out, sys.stdout, indent=1, sort_keys=True)
    print()


if __name__ == "__main__":
    main()
