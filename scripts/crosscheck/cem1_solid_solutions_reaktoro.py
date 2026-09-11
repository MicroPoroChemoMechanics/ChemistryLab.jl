#!/usr/bin/env python3
"""The Reaktoro half of the CEM I solid-solution cross-check.

Reaktoro (https://reaktoro.org) solves the same problem from the same
literature, reads the same ThermoFun file this package ships, and is mature and
widely used -- which is what makes it a useful second opinion.  For a
disagreement to mean anything the three knobs have to match: the same database
file, the same species and phase grouping, the same activity model.  They do
here, and so does a fourth that is specific to a cement: the element vector,
transferred in moles rather than in grams of oxide so that the two codes do not
each convert the datasheet with their own atomic-mass tables.

Run the Julia half first; see README.md.

    conda run -n reaktoro-env python cem1_solid_solutions_reaktoro.py
"""
import json
import os
import sys
from pathlib import Path

import reaktoro as rkt

HERE = Path(__file__).parent
OUT = HERE / "out"
DATABASE = Path(
    os.environ.get("CHEMISTRYLAB_DATA", HERE.parent.parent / "data")
) / "cemdata18-thermofun.json"

ELEMENTS = "H O C Ca Si Al Fe Mg K Na S"

# A full cement assemblage needs several hundred iterations, above the default
# cap; a run stopped at the cap returns an intermediate iterate rather than an
# answer, with an impossible pH.
MAXITERS = 4000

ALL_SOLUTIONS = {
    "CSHQ": ["CSHQ-JenD", "CSHQ-JenH", "CSHQ-TobD", "CSHQ-TobH", "KSiOH", "NaSiOH"],
    "C3(AF)S0.84H": ["C3AFS0.84H4.32", "C3FS0.84H4.32"],
    "AFt_SO4": ["ettringite", "ettringite30"],
    "AFt_SO4_CO3": ["tricarboalu03", "ettringite03_ss"],
    "AFm_SO4_OH": ["C4AH13", "monosulphate12"],
    "straetlingite": ["straetlingite", "straetlingite7"],
    "hydrotalc-pyro": ["Mg3AlC0.5OH", "Mg3FeC0.5OH"],
    "MSH": ["M075SH", "M15SH"],
}
# The three that carry mass on this cement.  Used as a configuration of its own
# because the eight-phase one does not converge here -- see README.md.
THREE = {k: ALL_SOLUTIONS[k] for k in ("CSHQ", "C3(AF)S0.84H", "AFt_SO4")}

PURE = (
    "AlOHmic Kln Gr C12A7 C2S C3A C3S C4AF CA CA2 C2AH7.5 C3AH6 CAH10 "
    "monosulphate10.5 monosulphate12 monosulphate14 monosulphate16 "
    "monosulphate9 chabazite zeoliteP_Ca straetlingite5.5 monocarbonate9 "
    "hemicarbonat10.5 hemicarbonate hemicarbonate9 monocarbonate "
    "ettringite13 ettringite9 Arg Cal C3FH6 C4FH13 C3FS1.34H3.32 "
    "Fe-hemicarbonate Femonocarbonate Dis-Dol Ord-Dol Lim Portlandite "
    "Anh Gp hemihydrate Fe Sd Mag FeOOHmic Py Tro Melanterite K2SO4 "
    "syngenite K2O hydrotalcite Mgs Brc Na2SO4 natrolite zeoliteX "
    "zeoliteY Na2O Sulfur Amor-Sl"
).split()

# `Material.add` reads a formula into element amounts, which is how the charge
# can be entered in the moles the Julia side computed rather than in grams.
FORMULA = {
    "C3S": "Ca3SiO5", "C2S": "Ca2SiO4", "C3A": "Ca3Al2O6",
    "C4AF": "Ca4Al2Fe2O10", "Gp": "CaSO4(H2O)2", "Brc": "Mg(OH)2",
    "Cal": "CaCO3", "K2O": "K2O", "Na2O": "Na2O", "H2O@": "H2O",
}


def build_system(solutions):
    db = rkt.ThermoFunDatabase.fromFile(str(DATABASE))
    aq = rkt.AqueousPhase(rkt.speciate(ELEMENTS))
    # CEMDATA18 carries no ion-size parameter, so the Debye-Huckel limiting law
    # with the non-ideality in the B-dot term and none of it on the neutral
    # species.  Matches HKFActivityModel(a = 0.0, Bdot = 0.097637, Kn = 0.0).
    dhp = rkt.ActivityModelDebyeHuckelParams()
    dhp.aiondefault = 0.0
    dhp.biondefault = 0.097637
    dhp.bneutraldefault = 0.0
    aq.set(rkt.ActivityModelDebyeHuckel(dhp))

    claimed = {m for members in solutions.values() for m in members}
    phases = [aq]
    for name, members in solutions.items():
        phase = rkt.SolidPhase(" ".join(members))
        phase.setName(name)
        phase.set(rkt.ActivityModelIdealSolution(rkt.StateOfMatter.Solid))
        phases.append(phase)
    phases.append(rkt.MineralPhases(" ".join(s for s in PURE if s not in claimed)))
    return rkt.ChemicalSystem(db, *phases)


def options():
    opts = rkt.EquilibriumOptions()
    opts.optima.maxiters = MAXITERS
    return opts


def cold_start(system, payload):
    mat = rkt.Material(system)
    for name, mol in payload["charge"].items():
        if mol > 0.0:
            mat.add(FORMULA.get(name, name), mol, "mol")
    mat.add("H2O", payload["free_water_mol"], "mol")
    state = mat.equilibrate(25.0, "celsius", 1.0, "bar", options())
    return state, mat.result()


def warm_start(system, amounts):
    state = rkt.ChemicalState(system)
    state.temperature(25, "celsius")
    state.pressure(1, "bar")
    known = {sp.name() for sp in system.species()}
    for name, mol in amounts.items():
        if mol > 0.0 and name in known:
            state.set(name, mol, "mol")
    return state, rkt.equilibrate(state, options())


def element_check(system, state, want):
    b = state.componentAmounts()
    got = {el.symbol(): float(b[i]) for i, el in enumerate(system.elements())}
    worst, worst_el = 0.0, ""
    for el, target in want.items():
        if el == "Zz" or abs(target) < 1e-8:   # the charge row, and trace carbon
            continue
        rel = abs(got.get(el, 0.0) - target) / abs(target)
        if rel > worst:
            worst, worst_el = rel, el
    print(f"    element vectors agree to {worst:.1e}"
          + (f" (worst: {worst_el})" if worst_el else ""))


def report(label, system, state, result, reference=None):
    ok = result.succeeded()
    print(f"\n  {label}: converged={ok} ({result.iterations()} iterations)")
    if not ok:
        print("    no answer to compare -- see README.md")
        return
    ph = float(rkt.AqueousProps(state).pH())
    vol = float(rkt.ChemicalProps(state).volume()) * 1e6
    print(f"    pH = {ph:.4f}     total volume = {vol:.4f} cm3")
    if reference is None:
        return
    print(f"    {'phase':18s} {'Reaktoro':>12s} {'ChemistryLab':>14s}  {'dev':>8s}")
    known = {sp.name() for sp in system.species()}
    worst = 0.0
    for name, want in sorted(reference.items(), key=lambda t: -t[1]):
        got = float(state.speciesAmount(name)) if name in known else 0.0
        dev = 100 * (got / want - 1)
        worst = max(worst, abs(dev))
        print(f"    {name:18s} {got:12.6g} {want:14.6g}  {dev:+7.2f} %")
    print(f"    worst deviation: {worst:.2f} %")


def main():
    charge_file, solution_file = OUT / "charge.json", OUT / "solution.json"
    for f in (charge_file, solution_file):
        if not f.exists():
            sys.exit(f"missing {f}: run cem1_solid_solutions.jl first")
    payload = json.loads(charge_file.read_text())
    converged = json.loads(solution_file.read_text())

    system = build_system(THREE)
    reference = {
        k: v for k, v in converged.items()
        if v > 1e-8 and any(
            sp.name() == k and sp.aggregateState() != rkt.AggregateState.Aqueous
            for sp in system.species()
        )
    }

    print(f"  three solid solutions declared, {len(system.species())} species")
    state, result = cold_start(system, payload)
    element_check(system, state, payload["elements"])
    report("cold start, from the anhydrous charge", system, state, result, reference)

    # The eight-phase configuration, for the record: it is the one the Julia side
    # certifies, and the one that does not converge here.  A solid solution whose
    # end-members are all at zero has no mole fractions, so its ideal-mixing term
    # is undefined there -- a difficulty of the formulation that every code has
    # to regularize somehow.
    system8 = build_system(ALL_SOLUTIONS)
    print(f"\n  eight solid solutions declared, {len(system8.species())} species")
    state8, result8 = cold_start(system8, payload)
    report("cold start", system8, state8, result8)
    state8, result8 = warm_start(system8, converged)
    report("warm start, from the certified composition",
           system8, state8, result8)


if __name__ == "__main__":
    main()
