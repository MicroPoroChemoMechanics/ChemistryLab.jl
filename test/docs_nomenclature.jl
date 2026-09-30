# The nomenclature of the documentation, docs/nomenclature.toml: the page
# `nomenclature.md` renders it and the hints over the equations read it. Every
# page an entry is restricted to exists, a symbol has at most one meaning that
# holds everywhere, and every constant an entry names is one that
# docs/nomenclature.jl reads from the library.

using TOML

@testsection "the nomenclature of the documentation" begin
    docs = joinpath(pkgdir(ChemistryLab), "docs")
    entries = TOML.parsefile(joinpath(docs, "nomenclature.toml"))["symbol"]
    @test length(entries) > 50
    for e in entries
        @test all(k -> haskey(e, k), ("group", "tex", "label", "name"))
        for p in get(e, "pages", String[])
            @test endswith(p, "/") ? isdir(joinpath(docs, "src", p)) : isfile(joinpath(docs, "src", p))
        end
    end
    # A subscript and a symbol written alike are two forms (the `r` of Δ_r, the
    # rate `r`), each with at most one meaning that holds everywhere.
    form(e) = (e["tex"], get(e, "script", false))
    for f in unique(form.(entries))
        @test count(e -> form(e) == f && isempty(get(e, "pages", String[])), entries) <= 1
    end
    # The label, the name and the unit are shown as HTML over the equations
    # (theme/symbol-hints.ts), so they hold no tag but <sub>, <sup> and <b>, and
    # every one is closed.
    for e in entries, field in ("label", "name", "unit")
        text = get(e, field, "")
        tags = [m.match for m in eachmatch(r"</?[^>]*>", text)]
        @test all(in(("<sub>", "</sub>", "<sup>", "</sup>", "<b>", "</b>")), tags)
        for t in ("sub", "sup", "b")
            @test count("<$t>", text) == count("</$t>", text)
        end
    end
    # One definition per concept: two entries never carry the same name, so a
    # quantity is written one way throughout (C_p, not C_p here and C_P there).
    names = [lowercase(strip(e["name"])) for e in entries]
    @test allunique(names)

    # The typography of the formulas (the note at the top of nomenclature.md):
    # a vector or a matrix in bold upright, its components in italic, and the
    # transpose written ^\mathsf{T}. What would bring back a second convention
    # is refused: an underlined vector, `\top` or `^{\mathsf T}` for the
    # transpose, and the heat capacity written C_P.
    src = joinpath(docs, "src")
    pages = [joinpath(r, f) for (r, _, fs) in walkdir(src) for f in fs if endswith(f, ".md") && !occursin("/.", r)]
    formulas(text) = vcat(
        [m.captures[1] for m in eachmatch(r"```math\n(.*?)```"s, text)],
        [m.captures[1] for m in eachmatch(r"``([^`]+)``", replace(text, r"```.*?```"s => " "))],
    )
    for p in pages, f in formulas(read(p, String))
        bad = [x for x in (r"\\underline", r"\\top\b", r"\^\{\\mathsf ?T\}", r"C_P\b") if occursin(x, f)]
        isempty(bad) || @error "a formula of $(relpath(p, src)) breaks the notation" f bad
        @test isempty(bad)
    end

    # The constants docs/nomenclature.jl maps to the library's values.
    @test issubset(Set(e["constant"] for e in entries if haskey(e, "constant")), Set(["R", "F", "N_A", "k_B", "e", "epsilon_0"]))
end
