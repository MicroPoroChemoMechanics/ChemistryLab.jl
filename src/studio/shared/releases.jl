# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2026 Jean-François Barthélémy (Cerema, UMR MCD)

# The key of a release of the solids in a thermogram, "Portlandite (water)" or
# "Cal (CO2)": a carbonated hydrate releases both. Included by the server, which
# reads the windows the page sets, and by the compute process, which reports the
# windows it used, so that both spell a key the same way.

const RELEASE_KEY = r"^(.+) \((water|CO2)\)$"

"The key of a release: its phase, and what it releases (`:water` or `:carbon_dioxide`)."
release_key(phase::AbstractString, releases::Symbol) = "$phase ($(releases === :water ? "water" : "CO2"))"

"The phase and the product a key names, `nothing` for a key that names no release."
function release_of(key::AbstractString)
    m = match(RELEASE_KEY, key)
    m === nothing && return nothing
    return String(m[1]) => (m[2] == "water" ? :water : :carbon_dioxide)
end
