# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)
"""The PHREEQC databases the oracle generators load, obtained as ChemistryLab
obtains them.

No database is stored in the repository. The generators load the very files the
package downloads: `database/<name>` of the PHREEQC repository at tag v3.7.3,
the release of the IPhreeqc engine bundled with phreeqpython (a database of a
later release uses keywords that engine refuses), checked against the same
SHA-256 as `THIRD_PARTY_DATABASES` in `src/databases/remote.jl`.

A copy is looked for in `$CHEMISTRYLAB_DATABASE_DIR` first, then in a cache of
this user's (`~/.cache/chemistrylab-oracles`), then downloaded. The files come
under the User Rights Notice of PHREEQC: Parkhurst & Appelo (2013), U.S.
Geological Survey Techniques and Methods 6-A43, doi:10.3133/tm6A43.
"""

import hashlib
import os
import urllib.request

TAG = "v3.7.3"

SHA256 = {
    "phreeqc.dat": "3e819f36a78a134b9e53557e9fd9e640d3616b84282cf7095c60337c49c6357b",
    "llnl.dat": "24d9266ff5c02aab84b2ba295567c8a5328d26f0a00fc69a43b6bb33b3efade2",
    "minteq.v4.dat": "1d4dd3f14932ccc4236f862f56689be6b1d1bbaa14067b23277bfcfa7c9bea3b",
    "wateq4f.dat": "b12c4a9818c946a882c675c458499d76a71fa0ec0becc8e5b9e175391ffaeb39",
    "sit.dat": "427d6114ed3f3135054683882319a0852b593d2685f0ccc80855f10ae1c4b840",
    "pitzer.dat": "3895bc5caf3f843abbb6deade72c9515c2be1c2efa3dd57052d0bb10cba68996",
}

URL = "https://raw.githubusercontent.com/phreeqc-dev/phreeqc3/{tag}/database/{name}"


def _sha256(path):
    with open(path, "rb") as handle:
        return hashlib.sha256(handle.read()).hexdigest()


def database_path(name="phreeqc.dat"):
    """The path of the PHREEQC database `name`, downloaded once if needed."""
    if name not in SHA256:
        raise SystemExit(f"{name} is not a PHREEQC database the oracles know: {sorted(SHA256)}")
    expected = SHA256[name]
    own = os.environ.get("CHEMISTRYLAB_DATABASE_DIR")
    if own:
        for root in [own] + [os.path.join(own, d) for d in sorted(os.listdir(own))]:
            candidate = os.path.join(root, name)
            if os.path.isfile(candidate) and _sha256(candidate) == expected:
                return candidate
    cache = os.path.join(os.path.expanduser("~"), ".cache", "chemistrylab-oracles")
    target = os.path.join(cache, name)
    if os.path.isfile(target) and _sha256(target) == expected:
        return target
    os.makedirs(cache, exist_ok=True)
    partial = target + ".part"
    urllib.request.urlretrieve(URL.format(tag=TAG, name=name), partial)
    found = _sha256(partial)
    if found != expected:
        os.remove(partial)
        raise SystemExit(f"{name}: downloaded SHA-256 {found}, expected {expected}")
    os.replace(partial, target)
    return target
