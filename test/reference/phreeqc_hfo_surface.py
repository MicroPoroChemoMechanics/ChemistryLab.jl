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


def jl(value: float) -> str:
    """A float as Julia's formatter writes it.

    Two differences, both of which Runic — the formatter this repository's CI
    gates on — would otherwise flag. Python renders an exponent with a leading
    zero, `6.0e-09` against `6.0e-9`, and it drops the fractional part of a
    round mantissa, `5e-06` against `5.0e-6`.

    The fix belongs here rather than in the file this emits: correcting the
    output alone means the next regeneration puts both back, and a formatting
    job then fails on a file nobody edited.
    """
    text = repr(value)
    if "e" not in text:
        return text
    mantissa, exponent = text.split("e")
    if "." not in mantissa:
        mantissa += ".0"
    exponent = re.sub(r"^([+-])0(\d)$", r"\1\2", exponent)
    return f"{mantissa}e{exponent}"


def database_path(name: str) -> str:
    return os.path.join(os.path.dirname(phreeqpython.__file__), "database", name)


def provenance(db: str) -> None:
    """Everything needed to reproduce this run, read at run time, never assumed."""
    with open(db, "rb") as handle:
        digest = hashlib.md5(handle.read()).hexdigest()
    try:
        from importlib.metadata import version

        pp_version = version("phreeqpython")
    except Exception:  # pragma: no cover - provenance must never fail the run
        pp_version = "unknown"
    print("# generated by test/reference/phreeqc_hfo_surface.py")
    print(f"# python       {sys.version.split()[0]}")
    print(f"# phreeqpython {pp_version}")
    print(f"# database     {os.path.basename(db)}  md5 {digest}")


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

    provenance(db)
    print(f"# model        {'diffuse double layer' if edl else 'no_edl'}")
    print(f"# Fe           {N_FE} mol -> strong {SITES_STRONG}, weak {SITES_WEAK} mol")
    print(f"# Zn total     {ZN_TOTAL} mol/kgw")
    print()
    print("const PHREEQC_HFO_ZN = (")
    for key, value in logks.items():
        print(f"    logK_{key} = {jl(value)},")
    print(f"    n_strong = {jl(SITES_STRONG)},")
    print(f"    n_weak = {jl(SITES_WEAK)},")
    print(f"    zn_total = {jl(ZN_TOTAL)},")
    print("    points = [")

    edl_line = "" if edl else "    -no_edl\n"
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
        print(
            f"        (pH = {jl(ph)}, la_H = {jl(col('la_H+'))}, "
            f"la_Zn = {jl(col('la_Zn+2'))}, "
            f"s_free = {jl(m['Hfo_sOH'])}, s_prot = {jl(m['Hfo_sOH2+'])}, "
            f"s_depr = {jl(m['Hfo_sO-'])}, s_zn = {jl(m['Hfo_sOZn+'])}, "
            f"w_free = {jl(m['Hfo_wOH'])}, w_prot = {jl(m['Hfo_wOH2+'])}, "
            f"w_depr = {jl(m['Hfo_wO-'])}, w_zn = {jl(m['Hfo_wOZn+'])}, "
            f"zn_free = {jl(m['Zn+2'])}, zn_sorbed_fraction = {jl(zn_sorbed / ZN_TOTAL)}),"
        )

    print("    ],")
    print(")")


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
    provenance(db)
    print("# model        no_edl, weak sites only, no metal")
    print(f"# sites        {SITES_PROTOLYSIS} mol")
    print(f"# log K        protonation {logks['protonation']}, "
          f"deprotonation {logks['deprotonation']}")
    print()
    print("const PHREEQC_PROTOLYSIS = (")
    print(f"    logK_protonation = {jl(logks['protonation'])},")
    print(f"    logK_deprotonation = {jl(logks['deprotonation'])},")
    print(f"    n_sites = {jl(SITES_PROTOLYSIS)},")
    print("    points = [")

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
        print(
            f"        (pH = {jl(ph)}, la_H = {jl(la_h)}, "
            f"free = {jl(fractions[0])}, protonated = {jl(fractions[1])}, "
            f"deprotonated = {jl(fractions[2])}),"
        )

    print("    ],")
    print(")")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case",
        choices=("protolysis", "zn-edge"),
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
    else:
        run_zn_edge(args.edl)
