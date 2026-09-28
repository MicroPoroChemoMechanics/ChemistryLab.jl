# The Julia half of the GEMS3K replay of De Weerdt et al. (2011).
#
#   julia --project=docs test/reference/xgems_deweerdt2011.jl
#   conda run -n mpcm-oracles python test/reference/xgems_replay.py deweerdt2011
#
# Writes, for each of the four mixes and each age of their Table 8, the element
# budget the recipe of `scripts/de_weerdt_2011.jl` puts into the equilibrium, and
# the list of species of its system, into the oracle cache (outside this
# repository). The Python half reads them, runs GEMS3K on the same budgets with
# the same phases, and writes `xgems_deweerdt2011.json` beside itself.

include(joinpath(@__DIR__, "..", "..", "scripts", "de_weerdt_2011.jl"))
using JSON

cs = dw11_system()
rows = Dict{String, Any}[]
for mix in DW11_MIXES
    recipe = dw11_recipe(mix)
    for d in dw11_days()
        push!(rows, Dict("mix" => mix, "time_d" => d, "elements" => budget_elements(cs, budget(recipe, cs; t = d).b)))
    end
end
cache = get(ENV, "CHEMISTRYLAB_ORACLE_CACHE", joinpath(homedir(), ".cache", "chemistrylab-oracles"))
mkpath(joinpath(cache, "deweerdt2011"))
path = joinpath(cache, "deweerdt2011", "budgets.json")
open(path, "w") do io
    JSON.print(
        io, Dict(
            "species" => [String(symbol(s)) for s in cs.species],
            "temperature_K" => dw11_value("curing_temperature"), "rows" => rows,
        ), 1,
    )
end
println("wrote ", path, ": ", length(rows), " budgets")
