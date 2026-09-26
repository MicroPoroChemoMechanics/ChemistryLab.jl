# ── Documentation shards: the pages of a pull request, built in parallel ─────
#
# The full site is built in one process, page after page, and a pull request
# waited for all of it before learning whether one `@example` block failed. A
# shard is a group of pages built on its own runner with `CHEMLAB_DOCS_ONLY`
# (see `partial.jl`), so the pull-request jobs run side by side and the feedback
# arrives after the slowest shard rather than after the sum.
#
# What the shards prove together is what a partial build proves, for every page
# at once: that every executed block runs and every page renders. What none of
# them proves is that the cross references resolve, since each prunes the pages
# of the others; the draft pre-flight job checks those on the whole tree, and the
# full build on `main` remains the one that deploys.
#
# The groups are balanced on measured cost (`scripts/docs_timing.jl`), and the
# pages no group names fall into the last one, `rest`, so a new page is built
# by a pull request without anyone having to assign it. A pattern that matches
# no page is refused rather than ignored: a renamed page would otherwise leave
# its group silently empty.

const DOC_SHARDS = [
    # The heaviest page alone: a sweep of certified solves on a 109-species
    # cement, several of which are refusals that pay the whole cascade.
    "cem4" => ["examples/cem4_pozzolanic.md"],
    # The other binders with solid solutions declared, each a cold certified
    # solve on a large system.
    "binders" => [
        "examples/cem1_solid_solutions.md",
        "examples/cem3_slag.md",
        "examples/cem5_composite.md",
        "examples/chloride_binding_blended.md",
    ],
    # Kinetic trajectories, coupled or not. The coupled ones are integrated by
    # `scripts/precomputed.jl`, on the runner's threads.
    "hydration" => [
        "tutorials/self_desiccation_kinetics.md",
        "examples/coupled_hydration.md",
        "examples/ionic_hydration.md",
        "examples/cement_clinker_kinetics.md",
    ],
    # The calibration on its own: three sets of coupled trajectories, distinct
    # from those of the pages above, so separating them duplicates nothing.
    "calibration" => [
        "examples/hydration_calibration.md",
    ],
    # The outputs of a calculation, each a handful of certified solves on the
    # cement system of `scripts/gruyaert2010.jl`, and the semi-adiabatic run,
    # a coupled trajectory with the temperature among its unknowns.
    "calorimetry" => [
        "examples/isothermal_calorimetry.md",
        "examples/thermogravimetry.md",
        "examples/semiadiabatic_calorimetry.md",
    ],
    # Everything else: theory, manual, the other tutorials and applications.
    "rest" => String[],
]

"""
    doc_page_leaves(node) -> Vector{String}

Every page path of a `pages` tree, in order. Duplicated from `partial.jl` so that
this file can be read without triggering a partial build.
"""
doc_page_leaves(node::AbstractString) = [node]
doc_page_leaves(node::Pair) = doc_page_leaves(node.second)
doc_page_leaves(node::AbstractVector) =
    reduce(vcat, doc_page_leaves.(node); init = String[])

"""
    shard_pages(name, leaves) -> Vector{String}

The pages of shard `name` among `leaves`. The shard `rest` holds every page no
other shard claims. Throws on an unknown name and on a pattern that matches no
page.
"""
function shard_pages(name::AbstractString, leaves::Vector{String})
    named = [s for s in DOC_SHARDS if !isempty(s.second)]
    for (shard, pats) in named, p in pats
        any(pg -> occursin(p, pg), leaves) || error(
            "documentation shard `$shard`: pattern `$p` matches no page of " *
                "docs/pages.jl; rename it with the page it stood for.",
        )
    end
    claimed(pg) = any(s -> any(p -> occursin(p, pg), s.second), named)
    i = findfirst(s -> s.first == name, DOC_SHARDS)
    i === nothing && error(
        "unknown documentation shard `$name`; known: " *
            join(first.(DOC_SHARDS), ", "),
    )
    pats = DOC_SHARDS[i].second
    return isempty(pats) ? filter(!claimed, leaves) :
        filter(pg -> any(p -> occursin(p, pg), pats), leaves)
end
