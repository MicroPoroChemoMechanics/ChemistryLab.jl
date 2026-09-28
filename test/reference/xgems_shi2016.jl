# The Julia half of the GEMS3K replay of Shi et al. (2016).
#
#   julia --project=docs test/reference/xgems_shi2016.jl
#   conda run -n mpcm-oracles python test/reference/xgems_replay.py shi2016
#
# Writes, for each of the four mortars, the element budget of its paste at 91
# days and of the same paste with each amount of carbon dioxide of
# `shi16_co2_grams()` added, and the list of species of its system, into the
# oracle cache (outside this repository). The Python half reads them, runs GEMS3K
# on the same budgets with the same phases, and writes `xgems_shi2016.json`
# beside itself.

include(joinpath(@__DIR__, "..", "..", "scripts", "shi_2016.jl"))
using JSON

cs = shi16_system()
rows = Dict{String, Any}[]
for mix in SHI16_MIXES
    paste = budget(shi16_recipe(mix), cs; t = SHI16_AGE)
    for g in shi16_co2_grams()
        b = shi16_carbonated_budget(cs, paste.b, g)
        push!(rows, Dict("mix" => mix, "co2_g" => g, "elements" => budget_elements(cs, b)))
    end
end
cache = get(ENV, "CHEMISTRYLAB_ORACLE_CACHE", joinpath(homedir(), ".cache", "chemistrylab-oracles"))
mkpath(joinpath(cache, "shi2016"))
path = joinpath(cache, "shi2016", "budgets.json")
open(path, "w") do io
    JSON.print(
        io, Dict(
            "species" => [String(symbol(s)) for s in cs.species],
            "temperature_K" => shi16_value("curing_temperature"), "rows" => rows,
        ), 1,
    )
end
println("wrote ", path, ": ", length(rows), " budgets")
