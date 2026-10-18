# The model: defaults, the JSON round trip, a model of the first studio read as
# the same calculation, the refusal of malformed fields, and each problem a
# model can have, raised by a model built for it.

m0 = CLS.default_model()
@test CLS.to_dict(CLS.model_from_dict(CLS.to_dict(m0))) == CLS.to_dict(m0)
d = JSON.parse(JSON.json(CLS.to_dict(m0)); dicttype = Dict{String, Any})
@test CLS.to_dict(CLS.model_from_dict(d)) == CLS.to_dict(m0)

# omitted parts take the defaults
m = CLS.model_from_dict(Dict{String, Any}("mix" => Dict("kind" => "mix", "w_b" => 0.4)))
@test m.mix.w_b == 0.4 && m.mix.temperature == 20.0
@test m.system == m0.system && length(m.materials) == 1

# a model of the first studio
v1 = Dict{String, Any}(
    "version" => 1, "title" => "old", "database" => "cemdata18-thermofun.json",
    "phase_list" => "Portland paste (Lothenbach and Winnefeld 2006)", "activity_model" => "cemdata18_NaOH",
    "materials" => [
        Dict(
            "template" => "PC (Durdzinski 2017), Bogue", "fraction" => 0.7,
            "extents" => [Dict("constituent" => "C3S", "extent" => 0.8)], "material_extent" => nothing
        ),
        Dict("template" => "S1 slag (Durdzinski 2017)", "fraction" => 0.3, "extents" => [], "material_extent" => 0.4),
    ],
    "w_b" => 0.45, "temperature_c" => 25.0, "mode" => "equilibrium",
)
old = CLS.model_from_dict(v1)
@test old.title == "old" && old.activity == CLS.Cemdata18Activity("NaOH", true)
@test old.mix.w_b == 0.45 && old.mix.temperature == 25.0
@test old.materials[1].extents == ["C3S" => CLS.ConstantExtentChoice(0.8)]
@test old.materials[2].material_extent == CLS.ConstantExtentChoice(0.4)

# malformed fields are refused by name
err(d) = try
    CLS.model_from_dict(d)
    nothing
catch e
    e
end
@test occursin("\"w_b\"", err(Dict{String, Any}("mix" => Dict("kind" => "mix", "w_b" => "x"))).msg)
@test occursin("\"title\"", err(Dict{String, Any}("title" => 3)).msg)
@test occursin("materials", err(Dict{String, Any}("materials" => "x")).msg)
@test occursin("version", err(Dict{String, Any}("version" => 99)).msg)
@test occursin("unknown calculation", err(Dict{String, Any}("calculation" => Dict("kind" => "magic"))).msg)

# problems: none for the default, each code for the model built to raise it
codes(m) = [p.code for p in CLS.problems(m)]
@test isempty(CLS.problems(m0))
with(f) = (m = CLS.model_from_dict(CLS.to_dict(m0)); f(m); m)
template(m) = m.materials[1]
@test codes(with(m -> empty!(m.materials))) == ["no_material"]
@test "fraction_not_positive" in codes(with(m -> (m.materials[1] = CLS.TemplateMaterial(template(m).template, 0.0, template(m).extents, nothing))))
@test codes(with(m -> (m.materials[1] = CLS.TemplateMaterial(template(m).template, 0.7, template(m).extents, nothing)))) == ["fractions_sum"]
@test codes(with(m -> (m.mix = CLS.Mix(0.5, 120.0, 100.0, 1.0)))) == ["out_of_range"]
@test codes(with(m -> (m.materials[1] = CLS.TemplateMaterial(template(m).template, 1.0, ["C3S" => CLS.TabulatedExtentChoice([1.0], [0.5])], nothing)))) == ["tabulated_extent"]
@test codes(with(m -> (m.materials[1] = CLS.TemplateMaterial(template(m).template, 1.0, ["Gp" => CLS.ParrottKillohExtentChoice(385.0, true)], nothing)))) == ["parrott_killoh_phase"]
@test codes(with(m -> (m.materials[1] = CLS.TemplateMaterial(template(m).template, 1.0, Pair{String, CLS.ExtentChoice}[], CLS.CappedExtentChoice(CLS.ParrottKillohExtentChoice(385.0, true), 0.5))))) == ["parrott_killoh_phase"]
oxide(rows) = CLS.OxideMaterial("mine", "cement", 1.0, rows, "bogue", "Gp", 0.0, Pair{String, CLS.ExtentChoice}[], nothing)
mineral(rows) = CLS.MineralMaterial("mine", "cement", 1.0, rows, Pair{String, CLS.ExtentChoice}[], nothing)
@test codes(with(m -> (m.materials[1] = oxide(CLS.AmountRow[])))) == ["empty_table"]
@test codes(with(m -> (m.materials[1] = oxide([CLS.AmountRow("CaO", 60, "g")])))) == ["oxide_unit"]
@test codes(with(m -> (m.materials[1] = oxide([CLS.AmountRow("CaO", 80, "%"), CLS.AmountRow("SiO2", 30, "%")])))) == ["oxide_total"]
@test isempty(codes(with(m -> (m.materials[1] = oxide([CLS.AmountRow("CaO", 64, "%"), CLS.AmountRow("SiO2", 21, "%")])))))
@test codes(with(m -> (m.materials[1] = mineral(CLS.AmountRow[])))) == ["empty_table"]
@test codes(with(m -> (m.materials[1] = mineral([CLS.AmountRow("C3S", 1, "mol")])))) == ["mineral_unit"]
@test codes(with(m -> (m.materials[1] = mineral([CLS.AmountRow("C3S", 80, "%"), CLS.AmountRow("C2S", 30, "%")])))) == ["mineral_total"]
@test isempty(codes(with(m -> (m.materials[1] = mineral([CLS.AmountRow("C3S", 80, "%")])))))
@test codes(with(m -> (m.database = CLS.ImportedDatabase("  ", "auto")))) == ["no_database_file"]
@test isempty(codes(with(m -> (m.database = CLS.ImportedDatabase("/some/file.dat", "auto")))))
@test codes(with(m -> (m.system = CLS.CustomSystem(String[], String[], String[])))) == ["empty_system"]
@test codes(with(m -> (m.system = CLS.CustomSystem(["C3S"], ["CSHQ", "CNASH_ss"], String[])))) == ["two_gels"]
# every path with an empty list of values
for T in CLS.members(CLS.CalculationChoice)
    list = findfirst(f -> f.kind === :numbers, CLS.fields(T))
    list === nothing && continue
    spec = Dict{String, Any}("kind" => CLS.kind(T), CLS.fields(T)[list].name => Float64[])
    @test "empty_path" in codes(with(m -> (m.calculation = CLS.choice_from_dict(CLS.CalculationChoice, spec))))
end
@test codes(with(m -> (m.calculation = CLS.AttackPath(CLS.AmountRow[], 500.0, 3, 28.0)))) == ["empty_solution"]
@test codes(with(m -> (m.calculation = CLS.AttackPath([CLS.AmountRow("Na2SO4", 1, "g")], 500.0, 3, 28.0)))) == ["solution_unit"]
@test codes(with(m -> (m.calculation = CLS.HydrationPath([1.0, 7.0])))) == ["no_law_of_age"]
@test isempty(codes(with(m -> (m.calculation = CLS.HydrationPath([1.0, 7.0]); m.materials[1] = CLS.TemplateMaterial(template(m).template, 1.0, ["C3S" => CLS.LogisticExtentChoice(0.8, 3.0, 1.5, 0.0)], nothing)))))
@test "activity_not_positive" in codes(with(m -> (m.constraint = CLS.FixedActivityChoice("CO2@", 0.0))))
@test "no_species" in codes(with(m -> (m.constraint = CLS.FixedActivityChoice(" ", 1.0))))
p = CLS.problems(with(m -> (m.materials[1] = CLS.TemplateMaterial(template(m).template, 0.7, template(m).extents, nothing))))[1]
@test p.values["total"] == 0.7
@test CLS.to_dict(p) == Dict{String, Any}("code" => p.code, "text" => p.text, "values" => p.values)

# laws that depend on the age, or not
@test !CLS.depends_on_age(m0)
@test CLS.depends_on_age(CLS.CappedExtentChoice(CLS.LogisticExtentChoice(0.8, 3.0, 1.5, 0.0), 0.5))
@test !CLS.depends_on_age(CLS.CappedExtentChoice(CLS.ConstantExtentChoice(0.5), 0.5))
@test CLS.depends_on_age(CLS.TemplateMaterial("x", 1.0, Pair{String, CLS.ExtentChoice}[], CLS.LogisticExtentChoice(0.8, 3.0, 1.5, 0.0)))
@test CLS.label(CLS.default_model().materials[1]) == "CEM I 52.5 N (Gruyaert 2010), Bogue"

# choices compare and hash by their fields
a, b = CLS.TabulatedExtentChoice([1.0, 7.0], [0.2, 0.5]), CLS.TabulatedExtentChoice([1.0, 7.0], [0.2, 0.5])
@test a == b && hash(a) == hash(b) && a !== b
@test a != CLS.TabulatedExtentChoice([1.0, 7.0], [0.2, 0.6])
@test length(Set([a, b])) == 1
