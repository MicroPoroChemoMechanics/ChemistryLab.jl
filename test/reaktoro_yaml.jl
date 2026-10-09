using JSON
using ChemistryLab.DataFrames: metadata

isdefined(@__MODULE__, :write_reaktoro_yaml) || include(joinpath(@__DIR__, "reference", "reaktoro_yaml_write.jl"))

# A database of Reaktoro, read and computed against Reaktoro itself
# (`test/reference/reaktoro_yaml.py`). The file is written during the test from
# ThermoHub data (`test/reference/reaktoro_yaml_write.jl`): no database of
# Reaktoro is stored or downloaded, and the oracle records the SHA-256 of the
# file it was computed from.
#
# Reaktoro's water has the densities of IAPWS-95, 7e-5 from those of the
# equation of Haar, Gallagher and Kell the package takes at 150 °C, as
# ThermoFun's do (`test/thermofun_methods.jl`); the HKF equations, which take
# the density and the permittivity of water, carry that difference, up to
# 2.5 J/mol in G and 3.4 J/(mol K) in Cp at 150 °C. Reaktoro integrates the
# volume of `HollandPowell` from zero pressure, as Eq. 13 of Holland and Powell
# (2011) writes it: its energies carry the integral up to the reference
# pressure P°, which the package's leave out; and away from P° its enthalpy and
# heat capacity leave out the temperature derivative of ∫ V dP, 9.4 J/mol in H
# for brucite at 100 bar, its Gibbs energies and volumes agreeing to 1e-10.
@testsection "Reaktoro YAML against Reaktoro" begin
    oracle = JSON.parsefile(joinpath(@__DIR__, "reference", "reaktoro_yaml.json"))
    path = write_reaktoro_yaml(joinpath(mktempdir(), "subset.yaml"))
    @test bytes2hex(open(ChemistryLab.sha256, path)) == oracle["yaml_sha256"]
    _, subs, _ = read_reaktoro_database(path)
    @test metadata(subs, "format") == "reaktoro"
    @test isequal(import_database(path)[2], subs)
    built = Dict(s.name => s for s in build_species(subs))
    @test Set(keys(built)) == Set(keys(oracle["species"]))
    at(f, row) = f(T = row["T"], P = row["P"])
    Pr = 1.0e5
    for (name, rows) in sort!(collect(oracle["species"]); by = first)
        s = built[name]
        model = only(subs[subs.name .== name, :standard_model])
        tait = model == "HollandPowell" && haskey(only(subs[subs.name .== name, :parameters]), "kappa0")
        @testset "$name" begin
            if haskey(ChemistryLab.properties(s), :refused_method)
                # A Landau transition of HollandPowell, at its reference state only.
                @test s[:refused_method] == "HollandPowell with a Landau transition"
                @test_throws ArgumentError at(s[:ΔₐG⁰], last(rows))
                continue
            end
            for row in rows
                if tait
                    # The integral up to P° taken as P° V, the compression over
                    # one bar leaving 3e-6 J/mol.
                    V(t) = s[:V⁰](T = t, P = Pr)
                    @test at(s[:ΔₐG⁰], row) + Pr * V(row["T"]) ≈ row["G"] atol = 1.0e-5
                    @test at(s[:V⁰], row) ≈ row["V"] rtol = 1.0e-10
                    if row["P"] == Pr
                        @test at(s[:ΔₐH⁰], row) + Pr * V(row["T"]) ≈ row["H"] atol = 1.0e-5
                        @test at(s[:Cp⁰], row) ≈ row["Cp"] atol = 1.0e-9
                    end
                elseif model in ("MaierKelley", "HollandPowell")
                    @test at(s[:ΔₐG⁰], row) ≈ row["G"] atol = 1.0e-6
                    @test at(s[:ΔₐH⁰], row) ≈ row["H"] atol = 1.0e-6
                    @test at(s[:Cp⁰], row) ≈ row["Cp"] atol = 1.0e-9
                    @test at(s[:V⁰], row) ≈ row["V"] rtol = 1.0e-12
                elseif model == "HKF"
                    @test at(s[:ΔₐG⁰], row) ≈ row["G"] atol = 3.0
                    @test at(s[:ΔₐH⁰], row) ≈ row["H"] atol = 15.0
                    @test at(s[:Cp⁰], row) ≈ row["Cp"] atol = 4.0
                    @test at(s[:V⁰], row) ≈ row["V"] atol = 6.0e-8
                else
                    # Water.
                    @test at(s[:ΔₐG⁰], row) ≈ row["G"] atol = 0.5
                    @test at(s[:ΔₐH⁰], row) ≈ row["H"] atol = 3.0
                    @test at(s[:Cp⁰], row) ≈ row["Cp"] atol = 0.1
                    @test at(s[:V⁰], row) ≈ row["V"] rtol = 1.0e-4
                end
            end
        end
    end
end

@testsection "Reaktoro YAML: the equations and the file" begin
    # The modified Tait equation without `kappa0pp` takes the article's
    # -kappa0p/kappa0; its volume is the derivative of its Gibbs energy.
    p = Dict(
        "Gf" => -569196.0, "Hf" => -601550.0, "Sr" => 26.94, "Vr" => 1.125e-5, "a" => 60.5, "b" => 3.62e-4,
        "c" => -535800.0, "d" => -299.2, "alpha0" => 3.11e-5, "kappa0" => 1.616e11, "kappa0p" => 3.95, "numatoms" => 2.0,
    )
    q = ChemistryLab._tait_parameters(p, 298.15)
    q2 = ChemistryLab._tait_parameters(merge(p, Dict("kappa0pp" => -3.95 / 1.616e11)), 298.15)
    @test q == q2
    for (T, P) in ((298.15, 1.0e5), (500.0, 2.0e8))
        dIdP = ForwardDiff.derivative(x -> ChemistryLab._tait_integral(q, T, x), P)
        Pth = q.ακ * q.θ / q.ξ0 * (1 / expm1(q.θ / T) - 1 / expm1(q.u0))
        @test dIdP ≈ q.V0 * (1 - q.a * (1 - (1 + q.b * (P - Pth))^(-q.c))) rtol = 1.0e-12
    end
    # A species of the model with `kappa0` and without its other parameters is
    # refused, naming them.
    dir = mktempdir()
    bad = joinpath(dir, "bad.yaml")
    write(
        bad, """
        Species:
          X:
            Name: X
            Formula: MgO
            Elements: 1:Mg 1:O
            AggregateState: Solid
            StandardThermoModel:
              HollandPowell:
                Gf: -569196.0
                Hf: -601550.0
                Sr: 26.94
                Vr: 1.125e-5
                a: 60.5
                kappa0: 1.616e11
        """,
    )
    _, subs, _ = read_reaktoro_database(bad)
    @test_throws "needs alpha0, kappa0p, numatoms" build_species(subs)
    # A species without a single standard model is not read, and one of a model
    # the package does not compute is left out with a warning.
    write(
        bad, """
        Species:
          Y:
            Name: Y
            Elements: 1:O
            AggregateState: Aqueous
          Z:
            Name: Z
            Elements: 1:Mg 1:O
            AggregateState: Solid
            StandardThermoModel:
              MineralHKF:
                Gf: -1.0
        """,
    )
    _, subs, _ = read_reaktoro_database(bad)
    @test any(occursin("Y has no single standard model", n) for n in metadata(subs, "notes"))
    @test isempty(@test_logs (:warn, r"left out.*Z \(`MineralHKF`\)") build_species(subs))
    # The excerpts the manual shows: a record of the zeolites ChemistryLab adds
    # to Cemdata18, and the same zeolite as a database of Reaktoro.
    page = read(joinpath(pkgdir(ChemistryLab), "docs", "src", "manual", "importing_databases.md"), String)
    shown = JSON.parse(only(m.captures[1] for m in eachmatch(r"```json\n(.*?)```"s, page)))
    record = only(r for r in JSON.parsefile(datapath("cemdata18-zeolites.json"))["substances"] if r["symbol"] == shown["symbol"])
    for (k, v) in shown
        @test v isa AbstractDict && haskey(v, "values") ? only(v["values"]) ≈ only(record[k]["values"]) : v == record[k]
    end
    write(bad, only(m.captures[1] for m in eachmatch(r"```yaml\n(.*?)```"s, page)))
    _, subs, _ = read_reaktoro_database(bad)
    z = only(build_species(subs))
    @test z[:ΔₐG⁰](T = 298.15, P = 1.0e5) ≈ only(record["sm_gibbs_energy"]["values"])
    @test z[:V⁰](T = 298.15, P = 1.0e5) ≈ 1.0e-5 * only(record["sm_volume"]["values"])
    @test z[:Cp⁰](T = 350.0, P = 1.0e5) ≈ only(record["sm_heat_capacity_p"]["values"])
    # A list is not read; a negative number is a scalar.
    write(bad, "Species:\n  X:\n    Name: X\n    Elements: [1:Mg, 1:O]\n")
    @test_throws "holds a list" read_reaktoro_database(bad)
    write(bad, "Species:\n  X:\n    Name: X\n    Charge: -2\n    Elements: 1:O\n    AggregateState: Aqueous\n    StandardThermoModel:\n      HKF:\n        Gf: -1.0\n")
    _, subs, _ = read_reaktoro_database(bad)
    @test only(subs.charge) == -2.0
end
