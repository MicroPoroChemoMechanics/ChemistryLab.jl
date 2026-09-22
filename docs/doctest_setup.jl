# The doctest preamble, in ONE place.
#
# `docs/make.jl` needs it so that a `jldoctest` block runs with the same modules
# in scope as the docstring's reader would have; the `Doctests` job of
# `.github/workflows/Documentation.yml` needs the identical one, because it runs
# the same blocks in a process that never loads `make.jl`.
#
# They were written out twice, and they drifted: the workflow's copy named
# `ModelingToolkit`, which this package does not depend on and which no doctest
# uses, while `Symbolics` — which several of them do use — was missing from it.
# The copy was commented out, so nothing ever caught the difference. Sharing the
# definition is what makes that class of drift impossible rather than unlikely.
#
# Every module listed here is a direct dependency of ChemistryLab, so this file
# loads in the package's own environment and needs nothing from `docs/`.

using Documenter: DocMeta
using ChemistryLab

DocMeta.setdocmeta!(
    ChemistryLab,
    :DocTestSetup,
    :(using ChemistryLab, DynamicQuantities, OrderedCollections, Symbolics);
    recursive = true,
)
