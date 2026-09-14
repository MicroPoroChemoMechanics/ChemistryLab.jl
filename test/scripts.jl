@testsection "the generated scripts still match their pages" begin
    # The binder-family pages carry their calculation in a dozen small blocks
    # separated by prose, each showing an intermediate result. That is good for a
    # reader and useless for someone who wants to RUN the calculation, so
    # `scripts/_page_to_script.jl` concatenates those blocks into one file per
    # page.
    #
    # A generated file that nobody checks is a file that rots. This regenerates
    # every one of them and compares.
    #
    # BY SYNTAX TREE, not by text. The scripts are formatted with Runic like the
    # rest of the repository, so their bytes differ from the page's by whitespace
    # alone; comparing text would turn this into a formatting check that fails on
    # the next Runic release. Comparing parsed expressions says the thing that is
    # actually meant -- it is the same program -- and a real divergence still
    # fails.
    root = dirname(@__DIR__)
    include(joinpath(root, "scripts", "_page_to_script.jl"))

    # Two normalizations, and both are needed rather than defensive.
    #
    # Line numbers differ between a page and the file extracted from it and carry
    # no meaning here. And Runic inserts an EXPLICIT `return` on the last
    # expression of a function, which is a real difference in the tree — a
    # `println()` becomes `return println()` — and the only one it makes that a
    # parser can see. Checked: with those two removed all seven scripts match
    # their pages exactly, so this normalization is `Runic`'s return rule and
    # nothing wider.
    normalize(x) = x
    function normalize(ex::Expr)
        ex.head === :return && length(ex.args) == 1 && return normalize(ex.args[1])
        args = Any[normalize(a) for a in ex.args if !(a isa LineNumberNode)]
        return Expr(ex.head, args...)
    end
    program(code) = normalize(Meta.parseall(code))

    @test !isempty(PAGES)
    stale = String[]
    for rel in PAGES
        page = joinpath(root, "docs", "src", rel)
        script = joinpath(root, "scripts", replace(basename(rel), ".md" => ".jl"))
        @test isfile(page)
        @test isfile(script)
        program(read(script, String)) == program(_script_text(page, rel)) ||
            push!(stale, basename(script))
    end
    isempty(stale) || @info "regenerate with: julia --project=. scripts/_page_to_script.jl"
    @test isempty(stale)

    # And the property that makes the concatenation sound in the first place:
    # every block of such a page shares one name, so Documenter runs them in ONE
    # module and the script is the same program. A page that mixed names would
    # produce a script that runs differently from the page it came from, and the
    # extractor refuses it rather than emitting one.
    for rel in PAGES
        name, codes = _page_blocks(joinpath(root, "docs", "src", rel))
        @test !isempty(name)
        @test !isempty(codes)
    end

    # The comparison must be able to FAIL. A tree comparison that accepted
    # anything would be worse than no test, so a deliberate change is shown to
    # break it while a pure reformatting does not.
    @test program("a = 1 + 2") == program("a  =  1 +\n    2")   # whitespace only
    @test program("f() = g()") == program("f() = return g()")   # Runic's return
    @test program("a = 1 + 2") != program("a = 1 + 3")          # a real change
    @test program("f(x) = x") != program("f(x) = -x")           # and another
end
