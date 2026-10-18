# The declarations of the choices: for every member of every family, the fields
# declared are the fields of the struct, in its order; the defaults build; an
# object round-trips; a model holding it writes a script that parses. Then the
# reading of each kind of field refuses what it must, by name.

using JSON

parses(src) = !any(a -> a isa Expr && a.head in (:error, :incomplete), Meta.parseall(src).args)

"A default model with `x` put where its family belongs."
function model_with(x)
    m = CLS.default_model()
    x isa CLS.DatabaseChoice && (m.database = x)
    x isa CLS.SystemChoice && (m.system = x)
    x isa CLS.ActivityChoice && (m.activity = x)
    x isa CLS.MaterialChoice && (m.materials = CLS.MaterialChoice[x])
    x isa CLS.MixChoice && (m.mix = x)
    x isa CLS.CalculationChoice && (m.calculation = x)
    x isa CLS.ConstraintChoice && (m.constraint = x)
    if x isa CLS.ExtentChoice
        m.materials = CLS.MaterialChoice[CLS.TemplateMaterial(m.materials[1].template, 1.0, ["C3S" => x], nothing)]
    end
    return m
end

for (name, F) in CLS.FAMILIES, T in CLS.members(F)
    @testset "$name / $(CLS.kind(T))" begin
        @test [f.name for f in CLS.fields(T)] == collect(string.(fieldnames(T)))
        d = Dict{String, Any}("kind" => CLS.kind(T))
        x = CLS.choice_from_dict(F, d)
        @test x isa T
        @test CLS.to_dict(CLS.choice_from_dict(F, CLS.to_dict(x))) == CLS.to_dict(x)
        @test CLS.to_dict(CLS.choice_from_dict(F, JSON.parse(JSON.json(CLS.to_dict(x)); dicttype = Dict{String, Any}))) == CLS.to_dict(x)
        m = model_with(x)
        g = CLS.generate(m)
        @test parses(g.source)
        @test startswith(g.source, g.compute)
        @test parses(CLS.diagnostic_source(m))
        @test CLS.to_dict(CLS.read_script(g.source).model) == CLS.to_dict(m)
    end
end

s = CLS.schema()
@test sort(collect(keys(s))) == sort(collect(keys(CLS.FAMILIES)))
# a kind names one choice: the page's labels are keyed by it
kinds = [CLS.kind(T) for (_, F) in CLS.FAMILIES for T in CLS.members(F)]
@test allunique(kinds)
@test all(m -> haskey(m, "kind") && haskey(m, "fields"), s["calculation"])
# the page's buttons: one equilibrium, or a path of equilibria
modes = Dict(m["kind"] => m["mode"] for m in s["calculation"])
@test modes["equilibrium"] == "equilibrium" && modes["hydration"] == "extents" && modes["water_sweep"] == "extents"
@test !any(m -> haskey(m, "mode"), s["material"])

# a family added the way a new feature adds one, then taken away
abstract type ToyChoice <: CLS.Choice end
struct ToyMember <: ToyChoice
    level::Float64
end
CLS.kind(::Type{ToyMember}) = "toy"
CLS.fields(::Type{ToyMember}) = [CLS.Field("level", :number; default = 1.0, min = 0)]
CLS.members(::Type{ToyChoice}) = (ToyMember,)
CLS.register_family!("toy", ToyChoice)
try
    @test CLS.family_name(ToyChoice) == "toy"
    @test only(CLS.schema()["toy"])["fields"][1]["name"] == "level"
    @test CLS.choice_from_dict(ToyChoice, Dict("kind" => "toy", "level" => 2)) == ToyMember(2.0)
    @test only(CLS.bound_problems(ToyMember(-1.0), "the toy")).code == "out_of_range"
finally
    delete!(CLS.FAMILIES, "toy")
end
@test !haskey(CLS.schema(), "toy")

# the reading of each kind of field
f(kind; kw...) = CLS.Field("x", kind; kw...)
@test_throws ArgumentError CLS.Field("x", :furlong)
@test_throws ArgumentError CLS.Field("x", :number; level = :guru)
@test CLS.field_value(f(:number), 2) === 2.0
@test_throws ArgumentError CLS.field_value(f(:number), "2")
@test_throws ArgumentError CLS.field_value(f(:number), true)
@test CLS.field_value(f(:integer), 3.0) === 3
@test_throws ArgumentError CLS.field_value(f(:integer), 3.5)
@test_throws ArgumentError CLS.field_value(f(:integer), true)
@test CLS.field_value(f(:choice; options = ["a", "b"]), "b") == "b"
@test_throws ArgumentError CLS.field_value(f(:choice; options = ["a", "b"]), "c")
@test_throws ArgumentError CLS.field_value(f(:text), 1)
@test CLS.field_value(f(:bool), false) === false
@test_throws ArgumentError CLS.field_value(f(:bool), 0)
@test CLS.field_value(f(:numbers), [1, 2.5]) == [1.0, 2.5]
@test_throws ArgumentError CLS.field_value(f(:numbers), [1, "x"])
@test CLS.field_value(f(:texts), ["a"]) == ["a"]
@test_throws ArgumentError CLS.field_value(f(:texts), [1])
@test CLS.field_value(f(:rows), [Dict("formula" => "CaO", "amount" => 1)]) == [CLS.AmountRow("CaO", 1, "%")]
@test_throws ArgumentError CLS.field_value(f(:rows), 3)
@test CLS.field_value(f(:nested; source = "extent"), nothing) === nothing
@test CLS.field_value(f(:nested; source = "extent"), Dict("kind" => "constant", "value" => 0.5)) == CLS.ConstantExtentChoice(0.5)
keyed = CLS.field_value(f(:keyed; source = "extent"), [Dict("key" => "C3S", "value" => Dict("kind" => "constant", "value" => 0.5))])
@test keyed == ["C3S" => CLS.ConstantExtentChoice(0.5)]
@test CLS.plain(first(keyed)) == Dict{String, Any}("key" => "C3S", "value" => Dict{String, Any}("kind" => "constant", "value" => 0.5))
@test_throws ArgumentError CLS.field_value(f(:keyed; source = "extent"), 3)
@test_throws ArgumentError CLS.field_value(f(:keyed; source = "extent"), [Dict("value" => 1)])
@test_throws ArgumentError CLS.choice_from_dict(CLS.ExtentChoice, 3)
e = try
    CLS.choice_from_dict(CLS.ExtentChoice, Dict("kind" => "magic"))
catch err
    err
end
@test e isa ArgumentError && occursin("constant", e.msg)
@test CLS.family_name(CLS.ExtentChoice) == "extent"
@test CLS.to_dict(CLS.Field("x", :number; default = 1.0, unit = "d", min = 0))["unit"] == "d"

# bounds of numbers, lists, nested and keyed choices
@test isempty(CLS.bound_problems(CLS.ConstantExtentChoice(0.5), "here"))
p = only(CLS.bound_problems(CLS.ConstantExtentChoice(1.5), "here"))
@test p.code == "out_of_range" && p.values["field"] == "value" && occursin("up to 1.0", p.text)
@test only(CLS.bound_problems(CLS.TabulatedExtentChoice([1.0, 2.0], [0.2, 1.4]), "here")).values["value"] == 1.4
@test length(CLS.bound_problems(CLS.CappedExtentChoice(CLS.ConstantExtentChoice(-1.0), 0.5), "here")) == 1
t = CLS.TemplateMaterial("x", 1.0, ["C3S" => CLS.ConstantExtentChoice(2.0)], nothing)
@test occursin("C3S", only(CLS.bound_problems(t, "x")).values["where"])
@test occursin("from 0.01", only(CLS.bound_problems(CLS.Mix(0.0, 20.0, 100.0, 1.0), "the mix")).text)
