# Oracle for test/pitzer.jl: Pitzer's J(x), as Reaktoro tabulates it.
#
# Reaktoro's ActivityModelPitzer.cpp carries J(x) as tables of values computed
# on grids of x (`J0region1` to `J0region3`), independently of the Chebyshev
# series of Harvie that this package evaluates. The test compares the two, which
# is what says the 42 coefficients of the series were transcribed right.
#
#   python3 test/reference/reaktoro_pitzer_j0.py
#
# Reads the source at a pinned commit and writes
# `test/reference/reaktoro_pitzer_j0.json`. No Reaktoro installation is needed.

import json
import re
import urllib.request

COMMIT = "f587235e692b7168faa272c8892a4c079a1f7dd9"
URL = ("https://raw.githubusercontent.com/reaktoro/reaktoro/" + COMMIT +
       "/Reaktoro/Models/ActivityModels/ActivityModelPitzer.cpp")

src = urllib.request.urlopen(URL).read().decode()


def region(name):
    m = re.search(name + r"\s*=\s*\{(.*?)\};", src, re.S)
    return [float(v) for v in re.findall(r"[-+]?\d\.\d+e[-+]\d+", m.group(1))]


r1, r2, r3 = region("J0region1"), region("J0region2"), region("J0region3")
# Region 1: x = 0.01 i on [0, 1]; region 2: x = 1 + 0.1 i on [1, 10];
# region 3: x = 10 + i on [10, 100].
points = [(0.01, r1[1]), (0.1, r1[10]), (0.5, r1[50]), (1.0, r1[100]),
          (2.0, r2[10]), (5.0, r2[40]), (10.0, r3[0]), (50.0, r3[40]),
          (100.0, r3[90])]

out = {
    "source": "Reaktoro, Reaktoro/Models/ActivityModels/ActivityModelPitzer.cpp",
    "commit": COMMIT,
    "license": "LGPL-2.1-or-later",
    "x": [x for x, _ in points],
    "J": [j for _, j in points],
}
with open("test/reference/reaktoro_pitzer_j0.json", "w") as f:
    json.dump(out, f, indent=2)
    f.write("\n")
