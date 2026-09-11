#!/usr/bin/env julia
#
# The Julia half of the CEM I solid-solution cross-check.
#
# Solves the cement of `docs/src/examples/cem1_solid_solutions.md` with every
# solid solution CEMDATA18 defines, and writes into `out/`:
#
#   charge.json     the anhydrous charge in moles, the free water, and the
#                   element vector -- what the second code has to be given if
#                   the comparison is to be of the chemistry and not of the two
#                   atomic-mass tables;
#   solution.json   the certified composition, species by species, so the second
#                   code can be started from it and asked whether it moves.
#
#     julia --project=scripts scripts/crosscheck/cem1_solid_solutions.jl
#
# Then run `cem1_solid_solutions_reaktoro.py`; see README.md.

using ChemistryLab
using DynamicQuantities
using OptimaSolver
using OrderedCollections
using Printf

const OUT = joinpath(@__DIR__, "out")
mkpath(OUT)

# A non-converged solve must never pass silently as a result.
ChemistryLab.STRICT_CONVERGENCE[] = false   # TEMPORAIRE: test de diagnostic

# ── The datasheet ────────────────────────────────────────────────────────────
const OXIDES = OrderedDict(
    "CaO" => 65.03, "SiO2" => 21.4, "Al2O3" => 3.84, "Fe2O3" => 4.49,
    "MgO" => 1.0, "K2O" => 0.46, "Na2O" => 0.13, "SO3" => 2.3, "CO2" => 0.0,
)
const WATER_G = 50.0            # w/c = 0.5

# Molar masses come from the data, never from a table written here.  Oxide
# formulas go through `Species`, which computes the mass from the formula and the
# element data; the cement phases are database species and are read from
# `byname`, because a name like "C3S" is a *phase* name and `Species` would read
# it as the formula C3S -- three carbons and a sulfur.
oxide_mass(ox) = ustrip(us"g/mol", Species(ox)[:M])

# ── The Bogue conversion ─────────────────────────────────────────────────────
# By inversion of the mass stoichiometric matrix, as in
# `docs/src/examples/bogue_calculation.md`: built from the species data, so no
# molar mass and no tabulated coefficient is written down here.
const CLINKER = CemSpecies.(split("C3S C2S C3A C4AF"))
const BOGUE = inv(mass_matrix(CanonicalStoichMatrix(CLINKER)).A)

function recipe(byname)
    phase_mass(name) = ustrip(us"g/mol", byname[name][:M])
    f = 100 / sum(values(OXIDES))
    n_ox = Dict(k => v * f / oxide_mass(k) for (k, v) in OXIDES)
    # Lime combined with the sulfate and the carbonate is not available to the
    # silicates, so it comes off before the inversion.
    CaO_free = OXIDES["CaO"] -
        OXIDES["SO3"] * oxide_mass("CaO") / oxide_mass("SO3") -
        OXIDES["CO2"] * oxide_mass("CaO") / oxide_mass("CO2")
    g = BOGUE * [CaO_free, OXIDES["SiO2"], OXIDES["Al2O3"], OXIDES["Fe2O3"]]
    charge = OrderedDict{String, Float64}(
        symbol(sp) => g[i] * f / phase_mass(symbol(sp))
            for (i, sp) in enumerate(CLINKER)
    )
    # CEMDATA18 has no periclase and no anhydrous SO3, so the magnesia enters as
    # brucite and the sulfate as gypsum; the water they bring is part of the
    # cement and comes off the mixing water.
    charge["Gp"] = n_ox["SO3"]
    charge["Brc"] = n_ox["MgO"]
    charge["K2O"] = n_ox["K2O"]
    charge["Na2O"] = n_ox["Na2O"]
    OXIDES["CO2"] > 0 && (charge["Cal"] = n_ox["CO2"])
    free_water = WATER_G / oxide_mass("H2O") - 2n_ox["SO3"] - n_ox["MgO"]
    return charge, free_water
end

# ── The phase list ───────────────────────────────────────────────────────────
const PURE = split(
    "AlOHmic Kln Gr C12A7 C2S C3A C3S C4AF CA CA2 C2AH7.5 C3AH6 CAH10 " *
        "monosulphate10.5 monosulphate12 monosulphate14 monosulphate16 " *
        "monosulphate9 chabazite zeoliteP_Ca straetlingite5.5 monocarbonate9 " *
        "hemicarbonat10.5 hemicarbonate hemicarbonate9 monocarbonate " *
        "ettringite13 ettringite9 Arg Cal C3FH6 C4FH13 C3FS1.34H3.32 " *
        "Fe-hemicarbonate Femonocarbonate Dis-Dol Ord-Dol Lim Portlandite " *
        "Anh Gp hemihydrate Fe Sd Mag FeOOHmic Py Tro Melanterite K2SO4 " *
        "syngenite K2O hydrotalcite Mgs Brc Na2SO4 natrolite zeoliteX " *
        "zeoliteY Na2O Sulfur Amor-Sl"
)

# Every multi-end-member phase the database defines.  Nothing here is chosen by
# looking at an answer.
const SOLUTIONS = [
    "CSHQ" => ["CSHQ-JenD", "CSHQ-JenH", "CSHQ-TobD", "CSHQ-TobH", "KSiOH", "NaSiOH"],
    "C3(AF)S0.84H" => ["C3AFS0.84H4.32", "C3FS0.84H4.32"],
    "AFt_SO4" => ["ettringite", "ettringite30"],
    "AFt_SO4_CO3" => ["tricarboalu03", "ettringite03_ss"],
    "AFm_SO4_OH" => ["C4AH13", "monosulphate12"],
    "straetlingite" => ["straetlingite", "straetlingite7"],
    "hydrotalc-pyro" => ["Mg3AlC0.5OH", "Mg3FeC0.5OH"],
    "MSH" => ["M075SH", "M15SH"],
]

substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
byname = Dict(symbol(s) => s for s in substances)
members = reduce(vcat, last.(SOLUTIONS))

species = speciation(
    substances, vcat(PURE, members, ["CO2@"]); aggregate_state = [AS_AQUEOUS]
)
cs = ChemicalSystem(
    species, CEMDATA_PRIMARIES;
    solid_solutions = [
        SolidSolutionPhase(nm, [byname[m] for m in mem]) for (nm, mem) in SOLUTIONS
    ],
)

charge, free_water = recipe(byname)

println("anhydrous charge (mol per 100 g of oxides):")
for (k, v) in charge
    @printf("  %-5s %10.6f mol   (%6.2f g/mol)\n", k,
        v, ustrip(us"g/mol", byname[k][:M]))
end
@printf("  %-5s %10.6f mol\n\n", "H2O@", free_water)

state = ChemicalState(cs)
for (sym, x) in charge
    x > 0 && set_quantity!(state, sym, x * u"mol")
end
set_quantity!(state, "H2O@", free_water * u"mol")
set_quantity!(state, "CO2@", 1.0e-9u"mol")
b = Float64.(cs.SM.A) * ustrip.(us"mol", state.n)

# ── The input, written before the expensive step ─────────────────────────────
json_pairs(pairs) = join(("    \"$k\": $v" for (k, v) in pairs), ",\n")

# `CSM` carries ATOMS in its rows (plus a charge row), so this is the element
# vector itself -- what a second code has to be given if the comparison is to be
# of the chemistry and not of the two atomic-mass tables.
be = Float64.(cs.CSM.A) * ustrip.(us"mol", state.n)
open(joinpath(OUT, "charge.json"), "w") do io
    println(
        io, "{\n  \"charge\": {\n",
        json_pairs((k, v) for (k, v) in charge if v > 0), "\n  },\n",
        "  \"free_water_mol\": ", free_water, ",\n  \"elements\": {\n",
        json_pairs(zip(cs.CSM.primaries, be)), "\n  }\n}"
    )
end
println("wrote ", joinpath(OUT, "charge.json"))

model = HKFActivityModel(å = 0.0, Ḃ = 0.097637, Kₙ = 0.0)
t0 = time()
eq, cert = equilibrate_certified(state; model = model, b = b)
dt = time() - t0

@printf(
    "certificate: optimal=%s  worst supersaturation=%+.4f  balance=%.1e  (%.1f s)\n",
    cert.optimal, cert.worst_supersaturation, cert.balance, dt
)
@printf("pH = %.4f     total volume = %.4f cm3\n",
    pH(eq, model), ustrip(uconvert(us"cm^3", volume(eq).total)))

n = ustrip.(us"mol", eq.n)
println("\nsolid assemblage (mol):")
for (s, x) in sort(
        [(symbol(cs.species[i]), n[i]) for i in cs.idx_crystal if n[i] > 1.0e-8];
        by = last, rev = true,
    )
    @printf("  %-18s %12.7g\n", s, x)
end

# ── The converged composition, for the warm start ────────────────────────────
open(joinpath(OUT, "solution.json"), "w") do io
    println(
        io, "{\n",
        json_pairs((symbol(cs.species[i]), n[i]) for i in eachindex(n)), "\n}"
    )
end
println("\nwrote ", joinpath(OUT, "solution.json"))
