# The chloride end member of CSHQ, and the sorption test it is fitted on.
#
# Included by data/chloride/regenerate.jl, which fits the end member on the test
# and writes data/cemdata18-chloride.json, and by test/cshq_chloride.jl, which
# checks what that file ships. The test is that of Hirao et al. (2005): 1 g of
# hydrated alite, C-S-H with the portlandite it formed, in 10 cm³ of NaCl
# solution at 20 °C, the bound chloride computed from what the solution lost.
# This file models that measurement, not the binding alone.
#
# Two candidates are written here, and both are fitted. The first follows the
# sodium end member of CSHQ; the second carries calcium chloride and nothing
# else. data/chloride/README.md says why the second is the one shipped.

module ChlorideMember

using ChemistryLab
using DynamicQuantities

"""
A candidate end member: its symbol, its formula, the species its composition is
written over (`composition(c)`: pairs `symbol => ν`, the member being `Σ ν s`),
and the rule for its entropy, heat capacity and volume.
"""
abstract type Candidate end

"""
Two formula units of NaSiOH, ((NaOH)2.5SiO2H2O)0.2, with their NaOH exchanged
for half as much CaCl2. Entropy, heat capacity and volume: those of the two
NaSiOH units.
"""
struct AfterNaSiOH <: Candidate end

"""
Half a CaCl2 in the gel, with no silica. Entropy and volume: those of the ions
it is made of, so that its dissolution has neither; heat capacity zero, the
negative heat capacities of the ions having no counterpart in a solid.
"""
struct CalciumChloride <: Candidate end

member_symbol(::AfterNaSiOH) = "CaClSiOH"
member_symbol(::CalciumChloride) = "CSHQ-Cl"
member_formula(::AfterNaSiOH) = "((CaCl2)1.25SiO2H2O)0.4"
member_formula(::CalciumChloride) = "(CaCl2)0.5"
member_name(::AfterNaSiOH) = "(CaCl2)1.25SiO2H2O, chloride end member of CSHQ written after NaSiOH, fitted by ChemistryLab"
member_name(::CalciumChloride) = "(CaCl2)0.5, chloride end member of CSHQ fitted by ChemistryLab"
composition(::AfterNaSiOH) = ("NaSiOH" => 2, "Na+" => -1, "OH-" => -1, "Ca+2" => 1 // 2, "Cl-" => 1)
composition(::CalciumChloride) = ("Ca+2" => 1 // 2, "Cl-" => 1)

"The candidate the database ships."
const SHIPPED = CalciumChloride()

# NaSiOH's reference temperature, and that of the test.
const T_REF = 293.15

# Above 1 mol/L the solution leaves the range of the B-dot model the equilibria
# are computed with.
const MAX_CHLORIDE = 1.0

at_ref(s, key, unit) = ustrip(unit, s[key](T = T_REF * u"K", P = 1.0e5u"Pa"; unit = true))
along(c, byname, key, unit) = sum(ν * at_ref(byname[s], key, unit) for (s, ν) in composition(c))

entropy(::AfterNaSiOH, byname) = 2 * at_ref(byname["NaSiOH"], :S⁰, us"J/(mol*K)")
entropy(c::CalciumChloride, byname) = along(c, byname, :S⁰, us"J/(mol*K)")

"""
    standard_properties(c, byname, δ) -> (; G, H, S, ΔrS)

The Gibbs energy, enthalpy and entropy of candidate `c` at 293.15 K and 1 bar
(J/mol, J/(mol K)), for `δ` (J/mol) the Gibbs energy of forming it from the
species of `composition(c)`: the one fitted number. The enthalpy follows from
the Gibbs energy and the entropy through that reaction, whose entropy `ΔrS` the
rule for the entropy implies.
"""
function standard_properties(c::Candidate, byname, δ)
    G = along(c, byname, :ΔₐG⁰, us"J/mol") + δ
    S = entropy(c, byname)
    ΔrS = S - along(c, byname, :S⁰, us"J/(mol*K)")
    H = along(c, byname, :ΔₐH⁰, us"J/mol") + δ + T_REF * ΔrS
    return (; G, H, S, ΔrS)
end

# The heat capacity and the volume of the record, from NaSiOH's.
function thermal!(e, ::AfterNaSiOH, base, byname)
    e["sm_heat_capacity_p"]["values"] = 2 .* base["sm_heat_capacity_p"]["values"]
    e["sm_volume"]["values"] = 2 .* base["sm_volume"]["values"]
    for m in e["TPMethods"]
        haskey(m, "m_heat_capacity_ft_coeffs") || continue
        m["m_heat_capacity_ft_coeffs"]["values"] = 2 .* m["m_heat_capacity_ft_coeffs"]["values"]
    end
    return e
end
function thermal!(e, c::CalciumChloride, base, byname)
    e["sm_heat_capacity_p"]["values"] = [0.0]
    # The partial molar volumes of the ions, in the record's J/bar.
    e["sm_volume"]["values"] = [along(c, byname, :V⁰, us"J/(bar*mol)")]
    for m in e["TPMethods"]
        haskey(m, "m_heat_capacity_ft_coeffs") || continue
        m["m_heat_capacity_ft_coeffs"]["values"] = zero.(m["m_heat_capacity_ft_coeffs"]["values"])
    end
    return e
end

"""
    member_entry(c, db, byname, δ; provenance = nothing) -> ThermoFun substance

Candidate `c` as a ThermoFun substance, on the pattern of the NaSiOH record of
`db` (the parsed CEMDATA18 file), with `byname` the species built from it.
"""
function member_entry(c::Candidate, db, byname, δ; provenance = nothing)
    base = only(s for s in db["substances"] if s["symbol"] == "NaSiOH")
    base["Tst"] == T_REF || error("NaSiOH is no longer referred to $T_REF K; the end member must be redone.")
    JSON = ChemistryLab.JSON
    e = JSON.parse(JSON.json(base))
    p = standard_properties(c, byname, δ)
    e["name"] = member_name(c)
    e["symbol"] = member_symbol(c)
    e["formula"] = member_formula(c)
    e["sm_gibbs_energy"]["values"] = [p.G]
    e["sm_enthalpy"]["values"] = [p.H]
    e["sm_entropy_abs"]["values"] = [p.S]
    thermal!(e, c, base, byname)
    e["datasources"] = ["Hirao2005: fitted by ChemistryLab, data/chloride/regenerate.jl"]
    provenance === nothing || (e["chloride_provenance"] = provenance)
    return e
end

"""
    member_species(entry) -> Species

The species a ThermoFun substance record gives, built by the reader
`build_species` uses on a database file.
"""
member_species(entry) = only(
    ChemistryLab.build_species(ChemistryLab.DataFrame(ChemistryLab.Tables.dictrowtable([entry])))
)

"""
    gel_system(subs, member) -> ChemicalSystem

The aqueous species of Ca, Si, Na and Cl, portlandite, and CSHQ with its sodium
end member and `member`, as `CSHQ_Cl`. Without potassium in the budget KSiOH
would have nothing to hold, and it is left out.
"""
function gel_system(subs, member)
    byname = Dict(symbol(s) => s for s in subs)
    gel = vcat([byname[m] for m in ("CSHQ-TobD", "CSHQ-TobH", "CSHQ-JenH", "CSHQ-JenD", "NaSiOH")], [member])
    sp = speciation(
        subs, ["Portlandite", "Amor-Sl", "Lim", "Na+", "Cl-"];
        aggregate_state = [AS_AQUEOUS], exclude_species = split("H2@ O2@ CH4@"),
    )
    return ChemicalSystem(
        vcat(sp, gel), CEMDATA_PRIMARIES; solid_solutions = [SolidSolutionPhase("CSHQ_Cl", gel)],
    )
end

"""
    calibration_points() -> (; chloride, bound)

The points of Hirao et al. (2005), Fig. 5, up to 1 mol/L, the origin left out:
chloride in solution after the test (mol/L) and bound chloride per gram of
hydrate (mmol/g).
"""
function calibration_points()
    t = literature_table("Hirao2005", "binding_hydrate")
    c = ustrip.(u"mol/L", t.chloride)
    b = ustrip.(u"mmol/g", t.bound)
    keep = [i for i in eachindex(c) if 0 < c[i] <= MAX_CHLORIDE]
    return (; chloride = c[keep], bound = b[keep])
end

"""
    forward(c, subs, db; model, water_fraction) -> δ -> NamedTuple

The test with candidate `c` as a function of δ (J/mol), at each calibration
point: `bound`, the bound chloride as the depletion measures it (mmol per gram
of hydrate); `held`, the chloride the end member holds (mmol/g); the states.

The hydrate is 1 g: C-S-H and portlandite in the proportions the two capacities
printed by Hirao et al. give, 0.61602 per gram of C-S-H against 0.47001 per gram
of hydrate. The C-S-H enters as lime, silica and water, at the Ca/Si the model
gives the gel at portlandite saturation, with `water_fraction` of its mass as
water. The solution is 10 cm³ at the initial concentration, `c₀ = c + b m/V`,
its water chosen so that the model's own liquid volume is 10 cm³. The bound
chloride is then computed as the measurement computes it, `(c₀ − c) V/m`, from
the concentration the model leaves in its liquid: the water the gel takes up or
gives back moves that concentration, in the model as in the test.
"""
function forward(c::Candidate, subs, db; model, water_fraction)
    byname = Dict(symbol(s) => s for s in subs)
    hq(q) = literature_value("Hirao2005", q)
    m = ustrip(us"g", hq("solid_mass"))
    V = ustrip(u"L", hq("solution_volume"))
    T = (273.15 + ustrip(hq("temperature_C"))) * u"K"
    csh_fraction = ustrip(hq("langmuir_capacity_hydrate")) / ustrip(hq("langmuir_capacity_csh"))
    pts = calibration_points()
    c0 = pts.chloride .+ pts.bound .* 1.0e-3 .* m ./ V         # mol/L
    M(s) = ustrip(us"g/mol", byname[s][:M])

    function charged(cs, r, n_Si, conc)
        idx = Dict(symbol(s) => i for (i, s) in enumerate(cs.species))
        n = zeros(length(cs.species))
        n[idx["Lim"]] = r * n_Si
        n[idx["Amor-Sl"]] = n_Si
        n[idx["Portlandite"]] = (1 - csh_fraction) * m / M("Portlandite")
        gel_water = water_fraction * csh_fraction * m / M("H2O@")
        # The solution's water, so that the model's liquid is V liters.
        W = V * 1000 / M("H2O@")
        for _ in 1:4
            s = zeros(length(cs.species))
            s[idx["H2O@"]] = W
            s[idx["Na+"]] = s[idx["Cl-"]] = max(conc, 1.0e-12) * V
            W *= V / ustrip(u"L", volume(ChemicalState(cs, s .* u"mol"; T)).liquid)
        end
        n[idx["H2O@"]] = W + gel_water
        n[idx["Na+"]] = n[idx["Cl-"]] = max(conc, 1.0e-12) * V
        return ChemicalState(cs, n .* u"mol"; T)
    end
    # The Ca/Si of the gel at portlandite saturation, which does not depend on
    # how much gel there is: from one equilibrium without chloride.
    n_Si(r) = (1 - water_fraction) * csh_fraction * m / (r * M("Lim") + M("Amor-Sl"))
    cs0 = gel_system(subs, member_species(member_entry(c, db, byname, 0.0)))
    eq0, cert0 = equilibrate_certified(charged(cs0, 1.7, n_Si(1.7), 0.0); model)
    cert0.optimal || error("the test without chloride did not certify.")
    t0 = solid_solution_totals(eq0, "CSHQ_Cl")
    r = t0.elements[:Ca] / t0.elements[:Si]

    return function (δ)
        cs = gel_system(subs, member_species(member_entry(c, db, byname, δ)))
        k_member = findfirst(s -> symbol(s) == member_symbol(c), cs.species)
        per_member = Float64(atoms(cs.species[k_member])[:Cl])
        order = sortperm(c0; rev = true)
        states0 = [charged(cs, r, n_Si(r), c0[k]) for k in order]
        A = Float64.(cs.SM.A)
        eqs, certs = equilibrate_path(first(states0), [A * ustrip.(us"mol", s.n) for s in states0]; model)
        all(x -> x.optimal, certs) || error("a point of the test did not certify at δ = $δ J/mol.")
        states = similar(eqs)
        states[order] = eqs
        chlorine = [Float64(get(atoms(s), :Cl, 0)) for s in cs.species]
        bound = map(eachindex(states)) do k
            n = ustrip.(us"mol", states[k].n)
            conc = sum(chlorine[i] * n[i] for i in cs.idx_aqueous) / ustrip(u"L", volume(states[k]).liquid)
            (c0[k] - conc) * V / m * 1.0e3
        end
        held = [per_member * ustrip(us"mol", eq.n[k_member]) * 1.0e3 / m for eq in states]
        return (; bound, held, states, c0, r)
    end
end

"""
    calcium_trend(subs, member; model) -> (; low, high, ratio_low, ratio_high)

Chloride held per silicon by a gel below portlandite saturation, one mole of
silicon made at Ca/Si 1.0 and at 1.4, in 0.3 mol of NaCl per kilogram of water
at 20 °C; and the Ca/Si the gel ends at. Measured, the binding grows with the
Ca/Si of the C-S-H.
"""
function calcium_trend(subs, member; model)
    cs = gel_system(subs, member)
    k = findfirst(s -> symbol(s) == symbol(member), cs.species)
    per_member = Float64(atoms(member)[:Cl])
    function held(r)
        st = ChemicalState(cs; T = 293.15u"K")
        set_quantity!(st, "Lim", r * u"mol")
        set_quantity!(st, "Amor-Sl", 1.0u"mol")
        set_quantity!(st, "H2O@", (1000 / ustrip(us"g/mol", Species("H2O")[:M])) * u"mol")
        set_quantity!(st, "Na+", 0.3u"mol")
        set_quantity!(st, "Cl-", 0.3u"mol")
        eq, cert = equilibrate_certified(st; model)
        cert.optimal || error("the gel at Ca/Si $r did not certify.")
        t = solid_solution_totals(eq, "CSHQ_Cl")
        return per_member * ustrip(us"mol", eq.n[k]) / t.elements[:Si], t.elements[:Ca] / t.elements[:Si]
    end
    (low, ratio_low), (high, ratio_high) = held(1.0), held(1.4)
    return (; low, high, ratio_low, ratio_high)
end

end # module
