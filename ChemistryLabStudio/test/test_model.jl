# The model: defaults, the JSON round trip, the refusal of malformed fields,
# and each problem a model can have, triggered by a model built for it.

using JSON

m0 = CLS.default_model()

@test CLS.to_dict(CLS.model_from_dict(CLS.to_dict(m0))) == CLS.to_dict(m0)
# through the JSON text, as the page and the embedded block carry it
d = JSON.parse(JSON.json(CLS.to_dict(m0)); dicttype = Dict{String, Any})
@test CLS.to_dict(CLS.model_from_dict(d)) == CLS.to_dict(m0)

# omitted fields take the defaults
m = CLS.model_from_dict(Dict{String, Any}("w_b" => 0.4))
@test m.w_b == 0.4
@test m.phase_list == m0.phase_list
@test length(m.materials) == 1

# an integer from JSON is a number like any other
@test CLS.model_from_dict(Dict{String, Any}("temperature_c" => 25)).temperature_c === 25.0

# malformed fields are refused by name
err(d) = try
    CLS.model_from_dict(d)
    nothing
catch e
    e
end
@test occursin("\"w_b\"", err(Dict{String, Any}("w_b" => "x")).msg)
@test occursin("\"title\"", err(Dict{String, Any}("title" => 3)).msg)
@test err(Dict{String, Any}("w_b" => true)) isa ArgumentError
@test occursin("materials", err(Dict{String, Any}("materials" => "x")).msg)
@test occursin("object", err(Dict{String, Any}("materials" => [1])).msg)
@test occursin("extents", err(Dict{String, Any}("materials" => [Dict("template" => "x", "fraction" => 1, "extents" => 3)])).msg)
@test occursin("version", err(Dict{String, Any}("version" => 99)).msg)

# material_extent: omitted, null or a number
mat = Dict{String, Any}("template" => "t", "fraction" => 1.0, "extents" => Any[])
@test CLS.model_from_dict(Dict{String, Any}("materials" => [mat])).materials[1].material_extent === nothing
mat["material_extent"] = 0.3
@test CLS.model_from_dict(Dict{String, Any}("materials" => [mat])).materials[1].material_extent == 0.3

# problems: none for the default, each code for the model built to raise it
codes(m) = [p.code for p in CLS.problems(m)]
@test isempty(CLS.problems(m0))
with(f) = (m = CLS.model_from_dict(CLS.to_dict(m0)); f(m); m)
@test codes(with(m -> empty!(m.materials))) == ["no_material"]
@test "fraction_not_positive" in codes(with(m -> (m.materials[1].fraction = 0.0)))
@test codes(with(m -> (m.materials[1].extents[1] = "C3S" => 1.5))) == ["extent_out_of_range"]
@test codes(with(m -> (m.materials[1].material_extent = -0.1))) == ["extent_out_of_range"]
@test codes(with(m -> (m.materials[1].fraction = 0.7))) == ["fractions_sum"]
@test CLS.problems(with(m -> (m.materials[1].fraction = 0.7)))[1].values["total"] == 0.7
@test codes(with(m -> (m.w_b = 0.0))) == ["w_b_not_positive"]
@test codes(with(m -> (m.temperature_c = 120.0))) == ["temperature_out_of_range"]
@test codes(with(m -> (m.activity_model = "pitzer"))) == ["unknown_activity_model"]
@test codes(with(m -> (m.mode = "kinetics"))) == ["unknown_mode"]
p = CLS.problems(with(m -> (m.w_b = 0.0)))[1]
@test CLS.to_dict(p) == Dict{String, Any}("code" => p.code, "text" => p.text, "values" => p.values)
