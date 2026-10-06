# =============================================================================
#  csh_aluminum.jl — two ways of giving the C-S-H gel aluminum, tried and tested
#
#  None of the gels the package ships holds both the calcium of a gel beside
#  portlandite and aluminum. This file tries two extensions, each with ONE energy
#  fitted by ChemistryLab, and tests them on data they were not fitted on. Nothing
#  here is part of a database; the page docs/src/examples/csh_aluminum.md says
#  what was found, and why no extension is shipped.
#
#    way A: CASH+ (Kulik et al. 2022, Miron et al. 2022a) with the aluminate
#           AlO(OH)2⁻ its authors describe for the bridging site: one compound per
#           occupant of the interlayer site (TAvh, TACh, TANh, TAKh), written as
#           the silicate compound with its bridging silicate exchanged for an
#           aluminate, + gibbsite − amorphous silica − water, plus δ.
#    way B: CSHQ (Kulik 2011) with one more ideal member, four formula units of
#           5CA of CNASH_ss (Myers et al. 2014), so that it carries one aluminum,
#           its Gibbs energy four times that of 5CA plus δ.
#
#  The fit and the tests use the uptake isotherm: for each synthesis, the gel at
#  the Al/Si measured by mass balance and the solution in equilibrium with it,
#  with the aluminum-bearing pure phases left out. L'Hôpital et al. (2015,
#  Appendix C) compute their solutions undersaturated with respect to the
#  strätlingite and the katoite they hold, so these phases are not at equilibrium
#  with them, and the Al/Si the gel keeps beside them is not an equilibrium datum.
#
#  ASSUMED, each where it is made: the 90 mL of solution are 90 g of water; the
#  syntheses of Yan et al. (2022), whose masses are not printed, are made like
#  those of L'Hôpital et al. (2 g of solids per 90 g of solution, Ca/Si 1.0 with
#  the calcium of the aluminate counted); the entropy, heat capacity and volume of
#  an aluminate compound are those of its silicate compound changed by the
#  exchange reaction; the aluminate has no site interaction of its own; the
#  same δ for every occupant of the interlayer site.
# =============================================================================

isdefined(@__MODULE__, :lh16a_table) || include(joinpath(@__DIR__, "lhopital2016_aluminum.jl"))

const _CSH_AL_JSON = ChemistryLab.JSON
const _CSH_AL_C18 = Ref{Any}(nothing)
const _CSH_AL_CP = Ref{Any}(nothing)
const _CSH_AL_CPDB = Ref{Any}(nothing)
# Read once, on first use.
_cached!(r, f) = (r[] === nothing && (r[] = f()); r[])
_csh_al_c18() = _cached!(_CSH_AL_C18, () -> _CSH_AL_JSON.parsefile(datapath("cemdata18-thermofun.json")))
_csh_al_cp() = _cached!(_CSH_AL_CP, () -> _CSH_AL_JSON.parsefile(datapath("cemdata18-cashplus.json")))
"""The species of `cemdata18-cashplus.json`, by symbol."""
csh_al_cashplus() = _cached!(_CSH_AL_CPDB, () -> Dict(symbol(s) => s for s in build_species(datapath("cemdata18-cashplus.json"); verbose = false)))

_csh_al_record(db, s) = only(r for r in db["substances"] if r["symbol"] == s)
_csh_al_value(r, k) = Float64(r[k]["values"][1])

# A ThermoFun record built into a species, as `build_species` reads a file; a
# record with a constant heat capacity has no temperature methods.
function _csh_al_species(e)
    haskey(e, "TPMethods") || (e["TPMethods"] = missing)
    haskey(e, "mass_per_mole") && delete!(e, "mass_per_mole")
    return only(ChemistryLab.build_species(ChemistryLab.DataFrame(ChemistryLab.Tables.dictrowtable([e]))))
end

"""The aluminum-bearing pure phases of the syntheses, left out of the isotherm."""
const CSH_AL_PHASES = ("C3AH6", "C3AS0.84H4.32", "straetlingite", "C4AH13", "CAH10", "AlOHmic")

# ── way A ────────────────────────────────────────────────────────────────────

"""The silicate compound each aluminate compound is written from, and its formula, by occupant of the interlayer site."""
const CSH_AL_A_MEMBERS = Dict(
    "v" => ("TSvh", "TAvh", "Ca2Si2AlO11H7"), "C" => ("TSCh", "TACh", "Ca3Si2AlO13H9"),
    "N" => ("TSNh", "TANh", "Ca2Si2AlNaO12H8"), "K" => ("TSKh", "TAKh", "Ca2Si2AlKO12H8"),
)

"""
    csh_al_aluminate(occupant, δ) -> Species

The compound of the aluminate on the bridging site of CASH+ with `occupant` on
the interlayer site (`"v"`, `"C"`, `"N"` or `"K"`): its silicate compound plus
gibbsite, minus amorphous silica and water, the Gibbs energy and the enthalpy
raised by `δ` (J/mol); entropy, heat capacity and volume follow the exchange.
"""
function csh_al_aluminate(occupant, δ)
    base, sym, formula = CSH_AL_A_MEMBERS[occupant]
    b = _csh_al_record(_csh_al_cp(), base)
    c18 = _csh_al_c18()
    gbs, sil, w = (_csh_al_record(c18, s) for s in ("Gbs", "Amor-Sl", "H2O@"))
    exchange(k) = _csh_al_value(gbs, k) - _csh_al_value(sil, k) - _csh_al_value(w, k)
    e = _CSH_AL_JSON.parse(_CSH_AL_JSON.json(b))
    e["symbol"], e["formula"] = sym, formula
    e["name"] = "$sym, aluminate on the bridging site of CASH+ (fitted by ChemistryLab, not published)"
    for k in ("sm_entropy_abs", "sm_heat_capacity_p", "sm_volume")
        e[k]["values"] = [_csh_al_value(b, k) + exchange(k)]
    end
    for k in ("sm_gibbs_energy", "sm_enthalpy")
        e[k]["values"] = [_csh_al_value(b, k) + exchange(k) + δ]
    end
    return _csh_al_species(e)
end

"""
    csh_al_gel_A(δ, occupants) -> SolidSolutionPhase

The gel of way A, `"C-S-H"`: the compounds of CASH+ (with sodium and potassium
from Miron et al. 2022a when `occupants` holds `"N"` or `"K"`) whose interlayer
occupant is among `occupants`, and the aluminate compound of each, mixed by
the compound energy formalism with the published site interactions.
"""
function csh_al_gel_A(δ, occupants)
    alkali = any(in(("N", "K")), occupants)
    ref, table = alkali ? ("Miron2022a", "cashplus_nk") : ("Kulik2022", "cashplus")
    occ = literature_table(ref, "$(table)_occupancy")
    on(m, site) = only(occ.species[i] for i in eachindex(occ.site) if occ.end_member[i] == m && occ.site[i] == site)
    core = [m for m in unique(occ.end_member) if on(m, "IC") in occupants]
    db = csh_al_cashplus()
    members = vcat([db[m] for m in core], [csh_al_aluminate(c, δ) for c in occupants])
    labels = permutedims(
        hcat(
            vcat([on(m, "BT") for m in core], fill("A", length(occupants))),
            vcat([on(m, "IC") for m in core], collect(occupants))
        )
    )
    lattice = SublatticeModel([1.0, 1.0], labels; sites = ["BT", "IC"])
    t = literature_table(ref, "$(table)_interactions")
    W = ustrip.(us"J/mol", t.W)
    keep(r) = t.site[r] == "BT" || (t.species_1[r] in occupants && t.species_2[r] in occupants)
    model = CompoundEnergyModel(lattice; interactions = [(t.site[r], t.species_1[r], t.species_2[r], W[r]) for r in eachindex(W) if keep(r)])
    return SolidSolutionPhase("C-S-H", members; model)
end

# ── way B ────────────────────────────────────────────────────────────────────

"""
    csh_al_member_B(δ) -> Species

Four formula units of 5CA, (CaO)5(SiO2)4(Al2O3)0.5(H2O)6.5, one aluminum per
formula unit: its Gibbs energy and enthalpy four times those of 5CA plus `δ`
(J/mol), its entropy, heat capacity and volume four times those of 5CA.
"""
function csh_al_member_B(δ)
    e = _CSH_AL_JSON.parse(_CSH_AL_JSON.json(_csh_al_record(_csh_al_c18(), "5CA")))
    for k in ("sm_gibbs_energy", "sm_enthalpy", "sm_entropy_abs", "sm_heat_capacity_p", "sm_volume")
        e[k]["values"] = [4 * _csh_al_value(e, k)]
    end
    for m in something(get(e, "TPMethods", nothing), [])
        haskey(m, "m_heat_capacity_ft_coeffs") && (m["m_heat_capacity_ft_coeffs"]["values"] = 4 .* m["m_heat_capacity_ft_coeffs"]["values"])
    end
    for k in ("sm_gibbs_energy", "sm_enthalpy")
        e[k]["values"] = [_csh_al_value(e, k) + δ]
    end
    e["symbol"] = "CSHQ-Al"
    e["formula"] = "(CaO)5(SiO2)4(Al2O3)0.5(H2O)6.5"
    e["name"] = "aluminum end member of CSHQ, four formula units of 5CA (fitted by ChemistryLab, not published)"
    return _csh_al_species(e)
end

# ── the isotherm ─────────────────────────────────────────────────────────────

"""
    csh_al_system(way, δ, alkali) -> ChemicalSystem

The system of a synthesis in `alkali` (`"none"`, `"NaOH"` or `"KOH"`) with the
gel of `way` (`:A` or `:B`) at `δ` (J/mol), and the pure phases of the
syntheses without the aluminum-bearing ones.
"""
function csh_al_system(way, δ, alkali)
    pure = [p for p in LH16_PURE if !(p in CSH_AL_PHASES)]
    if way === :A
        occupants = alkali == "none" ? ["v", "C"] : alkali == "NaOH" ? ["v", "C", "N"] : ["v", "C", "K"]
        gel = csh_al_gel_A(δ, occupants)
        excluded = vcat(split("H2@ O2@ CH4@"), alkali == "none" ? ["Ca(OH)2@"] : LH16_CASHPLUS_EXCLUDED)
        ions = alkali == "none" ? String[] : [first(LH16_ALKALI[alkali])]
        sp = speciation(collect(values(csh_al_cashplus())), vcat(pure, ions); aggregate_state = [AS_AQUEOUS], exclude_species = excluded)
    else
        names = alkali == "none" ? LH16_CSHQ : vcat(LH16_CSHQ, [last(LH16_ALKALI[alkali])])
        gel = SolidSolutionPhase("C-S-H", vcat([LH16_DB[m] for m in names], [csh_al_member_B(δ)]))
        sp = speciation(collect(values(LH16_DB)), vcat(pure, names); aggregate_state = [AS_AQUEOUS], exclude_species = split("H2@ O2@ CH4@"))
    end
    gel_symbols = Set(symbol.(end_members(gel)))
    return ChemicalSystem(vcat([s for s in sp if !(symbol(s) in gel_symbols)], end_members(gel)), CEMDATA_PRIMARIES; solid_solutions = [gel])
end

"""
    csh_al_synthesis(cs, Ca_Si, Al_Si; alkali = "none", c = 0.0) -> ChemicalState

2 g of lime, silica fume and monocalcium aluminate at `Ca_Si` (the calcium of
the aluminate counted) and `Al_Si`, in 90 g of water with `c` mol/L of `alkali`
hydroxide, at 20 °C: the recipe of L'Hôpital et al. (2015, 2016a), whose
Appendix A it reproduces to the printed digit.
"""
function csh_al_synthesis(cs, Ca_Si, Al_Si; alkali = "none", c = 0.0)
    M(s) = ustrip(us"g/mol", LH16_DB[s][:M])
    n_Si = 2.0 / (M("Amor-Sl") + Ca_Si * M("Lim") + Al_Si / 2 * (M("CA") - M("Lim")))
    st = ChemicalState(cs; T = 293.15u"K")
    set_quantity!(st, "H2O@", 90.0u"g")
    set_quantity!(st, "Lim", (Ca_Si - Al_Si / 2) * n_Si * u"mol")
    set_quantity!(st, "Amor-Sl", n_Si * u"mol")
    Al_Si > 0 && set_quantity!(st, "CA", Al_Si / 2 * n_Si * u"mol")
    if alkali != "none"
        set_quantity!(st, first(LH16_ALKALI[alkali]), c * 0.09u"mol")
        set_quantity!(st, "OH-", c * 0.09u"mol")
    end
    return st
end

"""
    csh_al_samples(set) -> Vector{NamedTuple}

The syntheses whose gel Al/Si and dissolved aluminum were both measured:
`"fit"`, those of L'Hôpital et al. (2016a) without alkali, which the energies
are fitted on; and, not used in the fit, `"lh16b"` (L'Hôpital et al. 2016b, in
KOH and NaOH), `"lh15"` (L'Hôpital et al. 2015, in 0.5 M KOH) and `"yan"` (Yan
et al. 2022, in NaOH and KOH). Each: a label, the target Ca/Si, the measured gel
Al/Si, the hydroxide and its concentration (mol/L), the measured aluminum
(mmol/L) and its qualifier.
"""
function csh_al_samples(set)
    out = NamedTuple[]
    add(label, Ca_Si, Al_Si, alkali, c, Al, q) = push!(out, (; label, Ca_Si, Al_Si, alkali, c, Al, qualifier = q))
    if set == "fit"
        g, p = lh16a_table("gel_composition"), lh16a_table("pore_solution")
        for i in eachindex(g.time_d)
            ismissing(g.Al_Si[i]) && continue
            same(k) = p.Ca_Si_target[k] == g.Ca_Si_target[i] && p.Al_Si_target[k] == g.Al_Si_target[i] && p.time_d[k] == g.time_d[i] && p.element[k] == "Al"
            j = findfirst(same, eachindex(p.element))
            j === nothing || add("$(g.Ca_Si_target[i])/$(g.Al_Si_target[i])", g.Ca_Si_target[i], g.Al_Si[i], "none", 0.0, ustrip(us"mmol/L", p.concentration[j]), p.qualifier[j])
        end
    elseif set == "lh16b"
        g, p = literature_table("LHopital2016b", "solid_composition"), literature_table("LHopital2016b", "pore_solution")
        # Two local functions of one name would be one function: each has its own.
        batch(t, k) = (t.Al_Si_target[k], t.Ca_Si_target[k], t.alkali[k], t.alkali_concentration[k], t.time_d[k])
        for i in eachindex(g.time_d)
            (g.alkali[i] == "none" || ismissing(g.Al_Si[i]) || iszero(g.Al_Si_target[i])) && continue
            j = findfirst(k -> batch(p, k) == batch(g, i) && p.element[k] == "Al", eachindex(p.element))
            j === nothing || add(
                "$(g.Ca_Si_target[i]) $(g.alkali[i])", g.Ca_Si_target[i], g.Al_Si[i], g.alkali[i],
                ustrip(us"mol/L", g.alkali_concentration[i]), ustrip(us"mmol/L", p.concentration[j]), p.qualifier[j]
            )
        end
    elseif set == "lh15"
        g, p = literature_table("LHopital2015", "gel_composition_KOH"), literature_table("LHopital2015", "pore_solution_KOH")
        c = ustrip(u"mol/L", literature_value("LHopital2015", "koh_concentration"))
        for i in eachindex(g.time_d)
            ismissing(g.Al_Si[i]) && continue
            j = findfirst(k -> p.Al_Si_target[k] == g.Al_Si_target[i] && p.time_d[k] == g.time_d[i] && p.element[k] == "Al", eachindex(p.element))
            j === nothing || add("1.0 KOH", 1.0, g.Al_Si[i], "KOH", c, ustrip(us"mmol/L", p.concentration[j]), p.qualifier[j])
        end
    elseif set == "yan"
        g, p = literature_table("Yan2022", "gel_composition"), literature_table("Yan2022", "pore_solution")
        sample(t, k) = (t.hydroxide[k], t.hydroxide_concentration[k], t.time_months[k], t.Al_Si_target[k])
        for i in eachindex(g.Al_Si)
            iszero(g.Al_Si[i]) && continue
            j = findfirst(k -> sample(p, k) == sample(g, i) && p.element[k] == "Al", eachindex(p.element))
            j === nothing || add(
                "1.0 $(g.hydroxide[i])", 1.0, g.Al_Si[i], g.hydroxide[i],
                ustrip(u"mol/L", g.hydroxide_concentration[i]), ustrip(us"mmol/L", p.concentration[j]), p.qualifier[j]
            )
        end
    else
        throw(ArgumentError("csh_al_samples: the set is \"fit\", \"lh16b\", \"lh15\" or \"yan\"; got \"$set\"."))
    end
    return out
end

"""
    csh_al_isotherm(way, δ, samples) -> Vector{NamedTuple}

Each of `samples` with its gel, at the measured Al/Si, in equilibrium with its
solution under the gel of `way` at `δ` (J/mol): the sample, the dissolved
aluminum the model gives (mmol/L), and whether the equilibrium is certified.
"""
function csh_al_isotherm(way, δ, samples)
    systems = Dict{String, Any}()
    return map(samples) do s
        cs = get!(() -> csh_al_system(way, δ, s.alkali), systems, s.alkali)
        eq, cert = equilibrate_certified(csh_al_synthesis(cs, s.Ca_Si, s.Al_Si; alkali = s.alkali, c = s.c); model = lh16_model(s.alkali))
        n = ustrip.(us"mol", eq.n)
        V = ustrip(uconvert(us"L", volume(eq).liquid))
        Al = 1000 * sum(n[k] * Float64(get(atoms(cs.species[k]), :Al, 0)) for k in cs.idx_aqueous) / V
        (; s..., model = Al, certified = cert.optimal)
    end
end

"""
    csh_al_misfit(results) -> NamedTuple

The decimal logarithm of the computed over the measured aluminum, over the
certified results whose aluminum was measured: its root mean square `rms`, its
mean `bias`, their number `n`; and `excess`, the largest amount by which a
computed aluminum exceeds a detection limit, in decimal logarithm. A result whose
equilibrium did not certify is left out.
"""
function csh_al_misfit(results)
    r = [log10(x.model / x.Al) for x in results if x.certified && x.qualifier == "measured"]
    excess = maximum((max(0.0, log10(x.model / x.Al)) for x in results if x.certified && x.qualifier != "measured"); init = 0.0)
    return (; rms = sqrt(sum(abs2, r) / length(r)), bias = sum(r) / length(r), n = length(r), excess)
end
