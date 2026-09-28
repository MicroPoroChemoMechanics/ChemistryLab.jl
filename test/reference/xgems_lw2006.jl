# The Julia half of the GEMS3K replay of Lothenbach & Winnefeld (2006).
#
#   julia --project=docs test/reference/xgems_lw2006.jl
#   conda run -n mpcm-oracles python test/reference/xgems_replay.py lw2006
#
# Writes, at each age of their Table 3, the element budget the recipe of
# `scripts/lothenbach_winnefeld_2006.jl` puts into the equilibrium, and the list
# of species of its system, into the oracle cache (outside this repository). The
# Python half reads them, runs GEMS3K on the same budgets with the same phases,
# and writes `xgems_lw2006.json` beside itself.

include(joinpath(@__DIR__, "..", "..", "scripts", "lothenbach_winnefeld_2006.jl"))
using JSON

cs = lw06_system()
recipe = lw06_recipe()
rows = [
    Dict("time_h" => h, "elements" => lw06_elements(cs, budget(recipe, cs; t = h / 24).b))
        for h in lw06_hours()
]
cache = get(ENV, "CHEMISTRYLAB_ORACLE_CACHE", joinpath(homedir(), ".cache", "chemistrylab-oracles"))
mkpath(joinpath(cache, "lw2006"))
path = joinpath(cache, "lw2006", "budgets.json")
open(path, "w") do io
    JSON.print(io, Dict("species" => [String(symbol(s)) for s in cs.species], "temperature_K" => recipe.T, "rows" => rows), 1)
end
println("wrote ", path, ": ", length(rows), " budgets")
