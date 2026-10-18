# The generated script: it parses, its compute part is sections 1 to 5, and the
# values the model holds are the ones the code passes, read back from the
# parsed expressions rather than from the text the generator wrote.

using Dates

m0 = CLS.default_model()
g = CLS.generate(m0; version = v"0.1.0", date = Date(2026, 10, 18))

"Every call to `f` in an expression tree."
function calls(ex, f::Symbol, found = Expr[])
    ex isa Expr || return found
    ex.head === :call && ex.args[1] === f && push!(found, ex)
    foreach(a -> calls(a, f, found), ex.args)
    return found
end
"The value of keyword `name` in a call, written after `;` or not."
function kwarg(call::Expr, name::Symbol)
    params = [a for a in call.args if a isa Expr && a.head === :parameters]
    kws = vcat((p.args for p in params)..., call.args)
    return only(a.args[2] for a in kws if a isa Expr && a.head === :kw && a.args[1] === name)
end

parsed = Meta.parseall(g.source)
@test !any(a -> a isa Expr && a.head in (:error, :incomplete), parsed.args)
@test startswith(g.source, g.compute)
@test occursin("## 5. Calculation", g.compute)
@test !occursin("## 6.", g.compute)
@test all(occursin("## $k. ", g.source) for k in 1:7)
@test occursin("Written by ChemistryLab Studio 0.1.0 on 2026-10-18.", g.source)
@test occursin("julia --project=@chemistrylab-studio cem_i_52_5_n_paste.jl", g.source)

compute = Meta.parseall(g.compute)
recipe = only(calls(compute, :Recipe))
@test kwarg(recipe, :w_b) == m0.w_b
@test recipe.args[end] == :(cem_i_52_5_n => 1.0)
extents = only(calls(compute, :with_extents))
@test extents.args[3].args[1] == :(Dict{String, Any})
@test Dict(e.args[2] => e.args[3] for e in extents.args[3].args[2:end]) == Dict(m0.materials[1].extents)
@test only(calls(compute, :phase_list_system)).args[2] == m0.phase_list
@test only(calls(compute, :cemdata18_activity_model)).args[2] == QuoteNode(:KOH)

# two materials, one with a degree for the whole material and none per constituent
m = CLS.model_from_dict(CLS.to_dict(m0))
m.materials[1].fraction = 0.7
push!(m.materials, CLS.MaterialChoice("S1 slag (Durdzinski 2017)", 0.3, Pair{String, Float64}[], 0.4))
m.activity_model = "cemdata18_NaOH"
m.temperature_c = 5
c = Meta.parseall(CLS.generate(m).compute)
@test only(calls(c, :Recipe)).args[(end - 1):end] == [:(cem_i_52_5_n => 0.7), :(s1_slag => 0.3)]
slag = only(e for e in calls(c, :with_extents) if e.args[3] == :s1_slag)
@test only(a.args[2] for a in slag.args[2].args if a.head === :kw) == 0.4
@test only(calls(c, :cemdata18_activity_model)).args[2] == QuoteNode(:NaOH)
# no extents and no material degree: the template is used as it is
m.materials[2].material_extent = nothing
@test length(calls(Meta.parseall(CLS.generate(m).compute), :with_extents)) == 1

@test CLS.material_name("S1 slag (Durdzinski 2017)", String[]) == "s1_slag"
@test CLS.material_name("PC (Durdzinski 2017), Bogue", ["pc"]) == "pc_2"
@test CLS.material_name("PC (x)", ["pc", "pc_2"]) == "pc_3"
@test CLS.material_name("42 (x)", String[]) == "material_42"
@test CLS.material_name("Recipe", String[]) == "recipe_material"
@test CLS.material_name("(x)", String[]) == "material_"

@test CLS.julia_number(1) == "1.0"
@test CLS.julia_number(0.45) == "0.45"
@test CLS.script_stem("CEM I 52.5 N, at 28 d!") == "cem_i_52_5_n_at_28_d"
@test CLS.script_stem("…") == "chemistrylab_study"
@test textwidth(CLS.heading(3, "Chemical system")) == CLS.RULE_WIDTH
@test_throws ArgumentError CLS.activity_model_code("pitzer")
