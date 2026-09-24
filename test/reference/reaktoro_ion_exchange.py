# Oracle for the ion-exchange tests: the same exchange, in both conventions.
#
# Reaktoro carries **both** — `ActivityModelIonExchangeGainesThomas` and
# `...Vanselow` — which makes it the one reference that can settle the question
# this comparison is about: the two conventions are not interchangeable, and a
# heterovalent exchange is where they part company.
#
#   conda create -n reaktoro-env -c conda-forge reaktoro thermofun
#   conda run -n reaktoro-env python test/reference/reaktoro_ion_exchange.py
#
# Writes `test/reference/reaktoro_ion_exchange.json`, which the test reads.
#
# What matches, and what is deliberately not matched:
#
#   * the exchange constants — read from PHREEQC's `phreeqc.dat` through
#     Reaktoro, and reported below so the Julia side declares the same ones
#     rather than a transcription of them;
#   * the aqueous activity model — **ideal on both sides**. Not because ideal is
#     right, but because an exchange constant and an aqueous activity
#     coefficient are two different things and this comparison is about the
#     first. The surface work makes the same choice for the same reason;
#   * the convention, which is the whole subject: each block below is one, and
#     they are reported separately because no conversion factor relates them.

import json
import os
import sys

import reaktoro as rkt

DATABASE = "phreeqc.dat"
EXCHANGERS = ["NaX", "KX", "CaX2"]
CHARGE = {"NaX": 1, "KX": 1, "CaX2": 2}

# One millimole of exchange capacity, and a solution holding all three cations
# in comparable amounts — the composition where the two conventions are furthest
# apart, rather than a trace limit where they agree by construction.
CAPACITY_MMOL = 1.0
AQUEOUS = {"Na+": 1.0, "K+": 1.0, "Ca+2": 0.5, "Cl-": 3.0}   # mmol


def formation_lgk(db, reaction: str) -> float:
    """The formation constant of one exchanger species, as Reaktoro read it.

    Note the syntax: Reaktoro wants `2*X-`, not `2X-`, and the second form
    raises rather than being parsed as one.
    """
    return float(db.reaction(reaction).props(298.15, "K", 1.0, "bar").lgK)


def exchange_constants(db):
    """The two exchange constants **relative to the sodium form**.

    The database gives formation constants from the free exchanger `X-`; the
    reaction a Julia `SiteFamily` writes is an exchange against its reference
    member. Subtracting is the conversion, and doing it here rather than by hand
    is what keeps the two sides declaring the same number.
    """
    lg_na = formation_lgk(db, "Na+ + X- = NaX")
    return {
        "K": formation_lgk(db, "K+ + X- = KX") - lg_na,
        "Ca": formation_lgk(db, "Ca+2 + 2*X- = CaX2") - 2 * lg_na,
    }


def run(convention: str):
    db = rkt.PhreeqcDatabase(DATABASE)
    solution = rkt.AqueousPhase(" ".join(["H2O", "H+", "OH-"] + list(AQUEOUS)))
    solution.set(rkt.ActivityModelIdealAqueous())
    exchange = rkt.IonExchangePhase(" ".join(EXCHANGERS))
    exchange.set(
        rkt.ActivityModelIonExchangeGainesThomas()
        if convention == "GainesThomas"
        else rkt.ActivityModelIonExchangeVanselow()
    )
    system = rkt.ChemicalSystem(db, solution, exchange)

    state = rkt.ChemicalState(system)
    state.temperature(25, "celsius")
    state.pressure(1, "bar")
    state.set("H2O", 1.0, "kg")
    for name, mmol in AQUEOUS.items():
        state.set(name, mmol, "mmol")
    state.set("NaX", CAPACITY_MMOL, "mmol")

    result = rkt.equilibrate(state)
    if not result.succeeded():
        raise SystemExit(f"Reaktoro failed to equilibrate ({convention})")

    amounts = {n: float(state.speciesAmount(n)) for n in EXCHANGERS}
    equivalents = sum(CHARGE[n] * a for n, a in amounts.items())
    return amounts, equivalents


if __name__ == "__main__":
    db = rkt.PhreeqcDatabase(DATABASE)
    constants = exchange_constants(db)
    payload = {
        "generator": "test/reference/reaktoro_ion_exchange.py",
        "python": sys.version.split()[0],
        "reaktoro": rkt.__version__,
        "database": f"{DATABASE} (via Reaktoro's PhreeqcDatabase)",
        "aqueous_model": "ideal on both sides, deliberately",
        "capacity_mmol_of_charge": CAPACITY_MMOL,
        "solution_mmol": AQUEOUS,
        # exchange constants relative to the sodium form, read from the database
        "logK_K": constants["K"],     # K+ + NaX = KX + Na+
        "logK_Ca": constants["Ca"],   # Ca+2 + 2 NaX = CaX2 + 2 Na+
        "capacity": CAPACITY_MMOL * 1.0e-3,
    }
    for convention in ("Vanselow", "GainesThomas"):
        amounts, equivalents = run(convention)
        payload[convention.lower()] = {**amounts, "equivalents": equivalents}

    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(here, "reaktoro_ion_exchange.json")
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=False)
        handle.write("\n")
    print(f"wrote {os.path.basename(path)}")
