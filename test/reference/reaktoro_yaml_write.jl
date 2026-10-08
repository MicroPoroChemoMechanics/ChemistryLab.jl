# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)
#
# Writes a database of Reaktoro, in its YAML format, from the slop98 and aq17 data
# that `datapath` obtains from ThermoHub: aqueous species with their HKF
# parameters, minerals with a Maier-Kelley heat capacity, minerals of the
# HollandPowell model, and water. It is written during
# `test/reaktoro_yaml.jl`, and written once by hand for the oracle,
#
#     julia --project test/reference/reaktoro_yaml_write.jl /tmp/slop98-subset.yaml
#     conda run -n reaktoro-env python test/reference/reaktoro_yaml.py /tmp/slop98-subset.yaml
#
# so that no database of Reaktoro is stored or downloaded; the oracle records the
# SHA-256 of the file it was computed from, which the test checks.

using ChemistryLab
using DynamicQuantities

# ThermoFun's units as Reaktoro's: J/bar and kbar in SI.
const J_PER_BAR = ustrip(us"m^3/mol", 1.0u"J/(bar*mol)")
const KBAR = 1000 * ustrip(us"Pa", 1.0u"bar")

const REAKTORO_AQUEOUS = ["Ca+2", "Na+", "Cl-", "OH-", "HCO3-", "CO3-2", "CO2@", "SiO2@", "SO4-2"]
const REAKTORO_MINERALS = ["Cal", "Lim", "Per"]
# Minerals of the HollandPowell model from aq17: its heat capacity, and the
# thermal expansion and the bulk modulus of Holland and Powell (1998) as `alpha0`
# and `kappa0` of the modified Tait equation, with `kappa0p` = 4 and `kappa0pp`
# the article's -kappa0p/kappa0, which Reaktoro requires, except for Forsterite.
# They test the reading and the equations, not a data
# set. Gibbsite has no bulk modulus in aq17, and Albite keeps its Landau
# transition.
const REAKTORO_HOLLAND_POWELL = ["Periclase", "Brucite", "Forsterite", "Albite", "Gibbsite"]

# A ThermoFun value at the reference, as a number in SI units.
_val(row, key) = Float64(only(row[key]["values"]))

function write_reaktoro_yaml(path)
    _, subs, _ = read_thermofun_database(datapath("slop98-inorganic-thermofun.json"))
    by = Dict(String(r.symbol) => r for r in eachrow(subs))
    composition(sp) = sort!([(e, n) for (e, n) in ChemistryLab.atoms(sp) if e != :Zz]; by = first)
    elements(sp) = join(("$(n):$(e)" for (e, n) in composition(sp)), " ")
    formula(sp) = join(("$(e)$(isone(n) ? "" : n)" for (e, n) in composition(sp)), "")
    open(path, "w") do io
        println(io, "Species:")
        sp = only(build_species(subs, ["H2O@"]))
        println(io, "  H2O(aq):\n    Name: H2O(aq)\n    Formula: H2O\n    Elements: 2:H 1:O\n    AggregateState: Aqueous")
        println(io, "    StandardThermoModel:\n      WaterHKF:\n        Ttr: 273.16\n        Str: 63.312288\n        Gtr: -235517.36\n        Htr: -287721.13")
        for s in REAKTORO_AQUEOUS
            row = by[s]
            hkf = only(m for m in row.TPMethods if only(values(m["method"])) == "solute_hkf88_reaktoro")
            c = Float64.(hkf["eos_hkf_coeffs"]["values"])
            f = ChemistryLab.HKF_SI_CONVERSIONS
            name = endswith(s, "@") ? replace(s, "@" => "(aq)") : s
            z = Float64(something(row.formula_charge, 0))
            sp = only(build_species(subs, [s]))
            println(io, "  $name:\n    Name: $name\n    Formula: $(replace(s, "@" => ""))\n    Elements: $(elements(sp))")
            iszero(z) || println(io, "    Charge: $z")
            println(io, "    AggregateState: Aqueous\n    StandardThermoModel:\n      HKF:")
            println(io, "        Gf: $(_val(row, "sm_gibbs_energy"))\n        Hf: $(_val(row, "sm_enthalpy"))\n        Sr: $(_val(row, "sm_entropy_abs"))")
            for (k, key) in enumerate((:a1, :a2, :a3, :a4, :c1, :c2, :wref))
                println(io, "        $(key): $(c[k] * f[key])")
            end
            println(io, "        charge: $z")
        end
        for s in REAKTORO_MINERALS
            row = by[s]
            cp = only(m for m in row.TPMethods if only(values(m["method"])) == "cp_ft_equation")
            a = Float64.(cp["m_heat_capacity_ft_coeffs"]["values"])
            sp = only(build_species(subs, [s]))
            println(io, "  $s:\n    Name: $s\n    Formula: $(formula(sp))\n    Elements: $(elements(sp))\n    AggregateState: Solid")
            println(io, "    StandardThermoModel:\n      MaierKelley:")
            println(io, "        Gf: $(_val(row, "sm_gibbs_energy"))\n        Hf: $(_val(row, "sm_enthalpy"))\n        Sr: $(_val(row, "sm_entropy_abs"))")
            println(io, "        Vr: $(J_PER_BAR * _val(row, "sm_volume"))\n        a: $(a[1])\n        b: $(a[2])\n        c: $(a[3])\n        Tmax: 1000.0")
        end
        _, aq17, _ = read_thermofun_database(datapath("aq17-thermofun.json"))
        hp = Dict(String(r.symbol) => r for r in eachrow(aq17))
        for s in REAKTORO_HOLLAND_POWELL
            row = hp[s]
            methods = Dict(only(values(m["method"])) => m for m in row.TPMethods)
            a = Float64.(methods["cp_ft_equation"]["m_heat_capacity_ft_coeffs"]["values"])
            sp = Species(String(row.formula))
            name = s * "-HP"
            println(io, "  $name:\n    Name: $name\n    Formula: $(formula(sp))\n    Elements: $(elements(sp))\n    AggregateState: Solid")
            println(io, "    StandardThermoModel:\n      HollandPowell:")
            println(io, "        Gf: $(_val(row, "sm_gibbs_energy"))\n        Hf: $(_val(row, "sm_enthalpy"))\n        Sr: $(_val(row, "sm_entropy_abs"))")
            println(io, "        Vr: $(J_PER_BAR * _val(row, "sm_volume"))\n        a: $(a[1])\n        b: $(a[2])\n        c: $(a[3])\n        d: $(a[4])")
            k = row.m_compressibility
            if !ismissing(k) && k !== nothing
                κ = KBAR * only(k["values"])
                println(io, "        alpha0: $(only(row.m_expansivity["values"]))\n        kappa0: $κ\n        kappa0p: 4.0")
                println(io, "        kappa0pp: $((s == "Forsterite" ? -4.5 : -4.0) / κ)")
                println(io, "        numatoms: $(sum(last, composition(sp)))")
            end
            if haskey(methods, "landau_holland_powell98")
                Tc, Smax, Vmax = Float64.(methods["landau_holland_powell98"]["m_landau_phase_trans_props"]["values"])
                println(io, "        Tcr: $(Tc + T_ZERO_CELSIUS)\n        Smax: $Smax\n        Vmax: $(J_PER_BAR * Vmax)")
            end
            println(io, "        Tmax: 2000.0")
        end
    end
    return path
end

if abspath(PROGRAM_FILE) == @__FILE__
    write_reaktoro_yaml(ARGS[1])
end
