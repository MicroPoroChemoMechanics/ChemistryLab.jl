# Build `cemdata18-zeolites.json`, CEMDATA18 plus the Empa zeolite series, and
# print the checks the build makes.
#
#   julia --project=docs data/zeolites/regenerate.jl
#
# ChemistryLab builds this database itself the first time it is asked for
# (`datapath("cemdata18-zeolites.json")`, src/databases/derived.jl), from the
# Cemdata18 file it obtains from its publisher and the published values of
# data/literature/MaLothenbach2020.json and MaLothenbach2021.json. This script
# runs the same build and shows what it checked, so that the result can be
# audited rather than trusted.
#
# What the build refuses to do: overwrite a CEMDATA18 substance (the 28 zeolites
# take new symbols, so `natrolite` and `NAT-Na` coexist); emit a dissolution that
# does not balance in atoms and in charge; emit a phase whose log Ksp, recomputed
# from its published ΔfG⁰ through CEMDATA18's own aqueous Gibbs energies, differs
# from the published one by more than 0.05 log units. That round trip can only
# close if the two datasets share a reference state and the transcription is
# exact, which is what makes the merge defensible.

import ChemistryLab
using Printf

function main()
    base = ChemistryLab.datapath("cemdata18-thermofun.json")
    println(rpad("zeolite", 14), lpad("log Ksp", 9), lpad("recomputed", 12), lpad("Δ", 8))
    for c in ChemistryLab.zeolite_logK_check(base)
        @printf("%-14s%9.2f%12.2f%8.3f\n", c.symbol, c.published, c.recomputed, c.difference)
    end
    path = ChemistryLab.datapath("cemdata18-zeolites.json")
    db = ChemistryLab.JSON.parsefile(path)
    v = db["zeolite_extension"]["verification"]
    @printf("\nlargest discrepancy: %.4f log units (tolerance %.2f)\n", v["worst_discrepancy_log_units"], v["tolerance_log_units"])
    println("built $(length(db["substances"])) substances: ", path)
    return nothing
end

abspath(PROGRAM_FILE) == (@__FILE__) && main()
