# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)
"""Standard properties from ThermoFun itself, for the methods of a ThermoFun
database that are not the heat-capacity polynomial and the HKF equations:

  - the substances a database defines by a reaction (440 in PSI/Nagra 12/07),
    whose energy at any temperature is that of the reactants plus the Gibbs
    energy of the reaction, from its log K;
  - the reactions of a database, with their entropy (`drsm_entropy`);
  - the substances whose entropy is obtained by integrating the heat capacity
    (`standard_entropy_cp_integration`);
  - the methods of aq17: Akinfiev and Diamond's equation for the dissolved gases
    (`solute_aknifiev_diamond03`), the Landau transitions and the volumes of
    Holland and Powell (`landau_holland_powell98`, `mv_eos_murnaghan_hp98`), and
    the fluids of Holland and Powell (`fluid_comp_redlich_kwong_hp91`).

For each, G, H, S, Cp and V at 25, 60, 90 and 150 °C and at 1 and 100 bar.
`Pu(OH)+3` of PSI/Nagra, whose reaction gives no coefficients of log K, is not
among them: ThermoFun 0.6 ends on a segmentation fault there.

Run with the environment that carries the ThermoFun Python package:

    conda run -n reaktoro-env python test/reference/thermofun_methods.py \
        <psinagra-12-07-thermofun.json> <cemdata18-thermofun.json> <aq17-thermofun.json> \
        > test/reference/thermofun_methods.json

The files are those `datapath` obtains from ThermoHub: no database is stored
here, only the values computed from them.
"""

import atexit
import hashlib
import importlib
import importlib.metadata
import json
import os
import shutil
import sys
import tempfile

TEMPERATURES_C = (25.0, 60.0, 90.0, 150.0)
PRESSURES_BAR = (1.0, 100.0)

# The substances and reactions to compute, by database.
SUBSTANCES = {
    "psinagra-12-07-thermofun.json": [
        # Defined by a reaction, on species of the database.
        "CaSiO3@", "AlF5-2", "CaSeO4@", "AmCO3+", "Ca2UO2(CO3)3@",
        # Defined by a reaction on a species itself defined by one.
        "Ni(CN)5-3", "BaSeO4(cr)",
        # Entropy by integration of the heat capacity.
        "Np+4", "HCN@",
    ],
    "cemdata18-thermofun.json": [
        "HCN@", "Melanterite",
        # Defined by a reaction, some of them given twice.
        "Fe(OH)3(am)", "zeoliteP_Ca", "chabazite", "MgSiO3@", "S-2", "CN-", "FeCO3(pr)",
    ],
    # Akinfiev and Diamond's dissolved gases; minerals of Holland and Powell,
    # with a Landau transition (Calcite, Quartz) or without, and without a bulk
    # modulus (Gibbsite); fluids of Holland and Powell; the water the dissolved
    # gases are computed with.
    "aq17-thermofun.json": ["CO2@", "CH4@", "H2@", "O2@", "Calcite", "Quartz", "Albite", "Gibbsite", "CO2", "H2O", "CH4", "H2O@"],
}
REACTIONS = {
    "psinagra-12-07-thermofun.json": ["CaSiO3@", "Ni(CN)5-3"],
    "cemdata18-thermofun.json": None,   # all of them
}


def props(p):
    return {
        "G": p.gibbs_energy.val,
        "H": p.enthalpy.val,
        "S": p.entropy.val,
        "Cp": p.heat_capacity_cp.val,
        "V": p.volume.val,
    }


def reaction_props(p):
    return {
        "G": p.reaction_gibbs_energy.val,
        "H": p.reaction_enthalpy.val,
        "S": p.reaction_entropy.val,
        "Cp": p.reaction_heat_capacity_cp.val,
        "V": p.reaction_volume.val,
        "logK": p.log_equilibrium_constant.val,
    }


def main():
    # ThermoFun writes a log, `thermofun.log`, in the working directory: run in a
    # directory of its own, it leaves nothing behind in the repository.
    paths = [os.path.abspath(p) for p in sys.argv[1:]]
    workdir = tempfile.mkdtemp(prefix="thermofun-")
    atexit.register(shutil.rmtree, workdir, True)
    os.chdir(workdir)
    # Imported here, the log it opens on import goes to that directory.
    global fun
    fun = importlib.import_module("thermofun")
    out = {"generator": "test/reference/thermofun_methods.py", "thermofun": importlib.metadata.version("thermofun"), "databases": {}}
    for path in paths:
        name = os.path.basename(path)
        engine = fun.ThermoEngine(fun.Database(path))
        with open(path, "rb") as handle:
            digest = hashlib.sha256(handle.read()).hexdigest()
        entry = {"sha256": digest, "substances": {}, "reactions": {}}
        for symbol in SUBSTANCES.get(name, []):
            entry["substances"][symbol] = [
                {"T": TC + 273.15, "P": Pb * 1e5, **props(engine.thermoPropertiesSubstance(TC + 273.15, Pb * 1e5, symbol))}
                for TC in TEMPERATURES_C for Pb in PRESSURES_BAR
            ]
        if name in REACTIONS:
            symbols = REACTIONS[name]
            if symbols is None:
                with open(path, encoding="utf-8") as handle:
                    symbols = [r["symbol"] for r in json.load(handle)["reactions"]]
            for symbol in symbols:
                entry["reactions"][symbol] = [
                    {"T": TC + 273.15, "P": Pb * 1e5, **reaction_props(engine.thermoPropertiesReaction(TC + 273.15, Pb * 1e5, symbol))}
                    for TC in TEMPERATURES_C for Pb in PRESSURES_BAR
                ]
        out["databases"][name] = entry
    json.dump(out, sys.stdout, indent=1, sort_keys=True)
    print()


if __name__ == "__main__":
    main()
