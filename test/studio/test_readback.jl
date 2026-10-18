# Reading a script back: the model returns exactly; an edit since the studio
# wrote it is detected; a file that is not a studio script is refused; a title
# holding the delimiters of a Julia block comment survives.

m0 = CLS.default_model()
src = CLS.generate(m0).source

rb = CLS.read_script(src)
@test CLS.to_dict(rb.model) == CLS.to_dict(m0)
@test !rb.edited
@test isempty(rb.message)

edited = replace(src, "w_b = 0.5 " => "w_b = 0.4 ")
rb = CLS.read_script(edited)
@test rb.edited
@test rb.model.mix.w_b == 0.5       # the model of the block, not the hand edit
@test occursin("w_b = 0.4 ", edited)
@test occursin("edited by hand", rb.message)

rb = CLS.read_script("x = 1\n")
@test rb.model === nothing
@test occursin("not written by ChemistryLab Studio", rb.message)

rb = CLS.read_script("x = 1\n" * CLS.MODEL_TAG)
@test rb.model === nothing
@test occursin("cut short", rb.message)

@test_throws Exception CLS.read_script("x = 1\n" * CLS.MODEL_TAG * " sha256=0\n{not json\n=#\n")

# delimiters of block comments and non-ASCII text in the title
m = CLS.model_from_dict(CLS.to_dict(m0))
m.title = "Pâte #= piège =# à 20 °C"
s = CLS.generate(m).source
@test !any(a -> a isa Expr && a.head in (:error, :incomplete), Meta.parseall(s).args)
rb = CLS.read_script(s)
@test rb.model.title == m.title
@test !rb.edited
