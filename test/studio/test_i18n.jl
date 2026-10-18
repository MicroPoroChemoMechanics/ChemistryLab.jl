# The words of the page: English and French hold the same keys; every kind and
# field of the schema, every problem the model can raise, every key the page
# asks for has its words in both.

words(lang) = JSON.parse(read(joinpath(CLS.WEB_DIR, "i18n", "$lang.json"), String); dicttype = Dict{String, Any})
en, fr = words("en"), words("fr")

@test sort(collect(keys(en))) == sort(collect(keys(fr)))
@test all(v -> v isa String && !isempty(v), values(en)) && all(v -> v isa String && !isempty(v), values(fr))

needed = Set{String}()
for (_, F) in CLS.FAMILIES, T in CLS.members(F)
    push!(needed, "k_" * CLS.kind(T))
    for f in CLS.fields(T)
        push!(needed, "f_" * f.name)
    end
end
for dir in (CLS.STUDIO_ROOT, joinpath(CLS.STUDIO_ROOT, "features")), file in readdir(dir; join = true)
    endswith(file, ".jl") || continue
    for m in eachmatch(r"Problem\(\s*\"([a-z_]+)\"", read(file, String))
        push!(needed, "problem_" * m[1])
    end
end
page = join((read(joinpath(CLS.WEB_DIR, f), String) for f in ("index.html", "app.js", "forms.js", "picker.js", "api.js")), '\n')
for m in eachmatch(r"data-i18n(?:-title)?=\"([a-z_0-9]+)\"", page)
    push!(needed, m[1])
end
for m in eachmatch(r"\bt\(\"([a-z_0-9]+)\"[,)]", page)
    push!(needed, m[1])
end
missing_words = sort!([k for k in needed if !haskey(en, k)])
@test isempty(missing_words)
