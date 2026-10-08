# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)
"""Standard properties computed by Reaktoro from a database in its YAML format,
for the reader of such databases (`read_reaktoro_database`, `test/reaktoro_yaml.jl`).

The database is the one `test/reference/reaktoro_yaml_write.jl` writes from the
slop98 and aq17 data of ThermoHub (HKF aqueous species, Maier-Kelley minerals,
minerals of the HollandPowell model, water), so that none is stored here:

    julia --project test/reference/reaktoro_yaml_write.jl /tmp/slop98-subset.yaml
    conda run -n reaktoro-env python test/reference/reaktoro_yaml.py /tmp/slop98-subset.yaml \\
        > test/reference/reaktoro_yaml.json

For every species, G, H, V and Cp at 25, 60, 90 and 150 °C, and at 1 and 100
bar, with the SHA-256 of the file they were computed from.
"""

import hashlib
import json
import sys

import reaktoro as rkt

TEMPERATURES_C = (25.0, 60.0, 90.0, 150.0)
PRESSURES_BAR = (1.0, 100.0)


def main():
    path = sys.argv[1]
    with open(path, "rb") as handle:
        digest = hashlib.sha256(handle.read()).hexdigest()
    db = rkt.Database.fromFile(path)
    out = {
        "generator": "test/reference/reaktoro_yaml.py",
        "reaktoro": rkt.__version__,
        "yaml_sha256": digest,
        "species": {},
    }
    for species in db.species():
        rows = []
        for TC in TEMPERATURES_C:
            for Pb in PRESSURES_BAR:
                p = species.standardThermoProps(TC + 273.15, Pb * 1e5)
                rows.append({"T": TC + 273.15, "P": Pb * 1e5, "G": float(p.G0), "H": float(p.H0), "V": float(p.V0), "Cp": float(p.Cp0)})
        out["species"][species.name()] = rows
    json.dump(out, sys.stdout, indent=1, sort_keys=True)
    print()


if __name__ == "__main__":
    main()
