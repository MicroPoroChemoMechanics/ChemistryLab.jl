# ── Partial builds: `CHEMLAB_DOCS_ONLY` ──────────────────────────────────────
#
# The full site costs over two hours, nearly all of it inside `@example` blocks,
# so a change touching one page pays for every other. This builds only the pages
# whose path contains one of a comma-separated list of patterns:
#
#     CHEMLAB_DOCS_ONLY=examples/cem3_slag julia --project=docs docs/make.jl
#     CHEMLAB_DOCS_ONLY=manual/,theory/    julia --project=docs docs/make.jl
#
# FILTERING `pages` IS NOT ENOUGH, and this is the whole reason the file exists.
# Documenter's `SetupBuildDirectory` walks the **source directory** and adds
# every `.md` it finds; `pages` only decides the navigation. So pruning the page
# tree alone leaves every `@example` on the site still running, and a "partial"
# build costs exactly what the full one costs. Measured: filtering `pages` down
# to one manual page still had the build inside `examples/` fifty minutes later.
#
# What does work is giving Documenter a different source tree. It walks with
# `follow_symlinks = true`, so a kept page is linked rather than copied and
# cannot go stale against its original -- falling back to a copy on a platform
# that refuses to create a symlink, which Windows does without Developer Mode.
# The page tree is filtered to match, so the navigation shows what was built.
#
# WHAT A PARTIAL BUILD PROVES: that the `@example` blocks of the pages it kept
# run, and that those pages render.
#
# WHAT IT CANNOT PROVE, and says so out loud: anything about links or about
# docstring coverage. Every `@ref` into a pruned page points at nothing, and
# every docstring whose API page was pruned is "missing", so `cross_references`
# and `missing_docs` are demoted to warnings -- for that build only. It also
# refuses to deploy. On a pull request the gate is therefore two jobs together:
# the partial builds of `docs/shards.jl`, which between them run every page, and
# the draft pre-flight of `docs/make.jl`, which checks links and docstring
# coverage on the whole tree. The full build, with `warnonly` kept down to
# `[:docs_block]`, runs on `main` and is the only one that deploys.
#
# `index.md` and `references.md` are always kept: the site root the nav is
# rendered around, and the page every `[@cite]` on a kept page resolves into.
#
# THE THIRD PIECE IS IN `docs/src/.vitepress/config.mts`. Demoting
# `cross_references` makes Documenter emit the unresolved references literally,
# as `./@ref` links, and VitePress then refuses to build over dead links -- which
# it is right to do for the full build. Its `ignoreDeadLinks` is therefore driven
# by this same environment variable, so the Julia half and the Node half of one
# decision cannot disagree. Measured before that was added: the Documenter half
# succeeded and `vitepress build` failed with "3 dead link(s) found".

const ALWAYS_KEPT = ("index.md", "references.md")

"""
    page_leaves(node) -> Vector{String}

Every page path under `node`, in build order. A leaf is a bare path or the value
of a `"Title" => "path"` pair; a `"Title" => [...]` pair is a section.
"""
page_leaves(node::AbstractString) = [node]
page_leaves(node::Pair) = page_leaves(node.second)
page_leaves(node::AbstractVector) = reduce(vcat, page_leaves.(node); init = String[])

"""
    is_kept(page, patterns) -> Bool

Whether a source-relative page path survives the filter.
"""
is_kept(page, pats) = page in ALWAYS_KEPT || any(p -> occursin(p, page), pats)

"""
    keep_pages(node, patterns) -> node or `nothing`

The page tree with only the leaves that [`is_kept`](@ref) admits. Returns
`nothing` when nothing under `node` survives, which prunes a section header
whose whole contents were filtered out.
"""
keep_pages(node::AbstractString, pats) = is_kept(node, pats) ? node : nothing
function keep_pages(node::Pair, pats)
    kept = keep_pages(node.second, pats)
    return kept === nothing ? nothing : (node.first => kept)
end
function keep_pages(node::AbstractVector, pats)
    kept = Any[]
    for child in node
        k = keep_pages(child, pats)
        k === nothing || push!(kept, k)
    end
    return isempty(kept) ? nothing : kept
end

"""
    pruned_source(srcdir, patterns) -> String

A temporary directory mirroring `srcdir`, with the Markdown pages the filter
rejects left out. Everything that is not a page -- images, `refs.bib`,
`.vitepress` -- is carried through unconditionally, because a kept page may
reference any of it. The directory is removed when the process exits.

# Pages are linked, everything else is COPIED, and the asymmetry is not a detail

Documenter walks with `follow_symlinks = true`, so a symlinked page is read
through and cannot drift from the file it stands for. Nothing writes back to a
page's source, so linking them is safe.

Everything else must be copied, and the reason is a defect this file caused once
already. `SetupBuildDirectory` copies each non-Markdown file with
`cp(src, dst; force = true)`, and Julia's `cp` reproduces a **symlink as a
symlink**. So `build/.vitepress/config.mts` became a link pointing back at
`docs/src/.vitepress/config.mts`, and DocumenterVitepress -- which fills that
file's `REPLACE_ME_DOCUMENTER_VITEPRESS` markers and writes the result to what it
believes is the build copy -- wrote straight through the link and **overwrote the
template in the repository**, baking in the navigation of the three-page build
that produced it. Measured: a committed `config.mts` with its markers gone and a
sidebar listing one manual page.

Copying is 0.6 MB across 28 files.
"""
function pruned_source(srcdir, pats)
    tmp = mktempdir(; prefix = "chemlab_docs_")
    for (root, _, files) in walkdir(srcdir)
        rel = relpath(root, srcdir)
        mkpath(normpath(joinpath(tmp, rel)))
        for f in files
            src = abspath(joinpath(root, f))
            dst = normpath(joinpath(tmp, rel, f))
            if endswith(f, ".md")
                is_kept(replace(normpath(joinpath(rel, f)), '\\' => '/'), pats) || continue
                # A symlink where the platform allows one, a copy where it does
                # not. Creating a symlink on Windows needs Developer Mode or an
                # elevated process, and a documentation build must not depend on
                # either; the copy is equivalent for a page, which is only ever
                # read. What it loses is the guarantee that a pruned page cannot
                # drift from its source -- and it cannot drift far, since the
                # tree lives only for the length of the build.
                try
                    symlink(src, dst)
                catch
                    cp(src, dst; force = true)
                end
            else
                cp(src, dst; force = true)
            end
        end
    end
    return tmp
end

const DOCS_ONLY = strip(get(ENV, "CHEMLAB_DOCS_ONLY", ""))
const PARTIAL_BUILD = !isempty(DOCS_ONLY)

const DOCS_SOURCE = if PARTIAL_BUILD
    patterns = filter(!isempty, strip.(split(DOCS_ONLY, ',')))

    # Tested against the leaves BEFORE filtering, not against the result. The
    # two pages above are kept unconditionally, so the filtered tree is never
    # empty and a typo in the variable would otherwise build a two-page site and
    # report success -- the failure this mechanism exists to prevent.
    matched = filter(pg -> any(p -> occursin(p, pg), patterns), page_leaves(pages))
    isempty(matched) && error(
        "CHEMLAB_DOCS_ONLY=\"$DOCS_ONLY\" matches no page. A pattern is matched " *
            "against the path as it appears in docs/pages.jl, for example " *
            "`examples/cem3_slag` or `manual/`.",
    )

    pages = keep_pages(pages, patterns)
    @warn """
    PARTIAL DOCUMENTATION BUILD -- not deployable, and it checks neither links
    nor docstring coverage. `cross_references` and `missing_docs` are demoted to
    warnings, because a reference into a pruned page, and a docstring whose API
    page was pruned, would otherwise be reported as broken. Run the full build
    before merging.
    Patterns: $(join(patterns, ", "))
    Pages   : $(join(matched, ", "))"""

    pruned_source(joinpath(@__DIR__, "src"), patterns)
else
    "src"
end
