# The page tree, grouped by what a chapter is *about* rather than by the order the
# pages were written in.
#
# Four chapters carry the documentation, and the difference between them is the
# question each answers.
#
#   **Theory**    — why is this the right calculation? Readable without running
#                   anything; every equation is the one the code evaluates.
#   **Manual**    — how is this written? One page per kind of object, syntax
#                   first, on the smallest example that shows it.
#   **Tutorials** — how do I drive a calculation from end to end? Narrative, and
#                   each one goes somewhere.
#   **Applications** — what does a real case look like, and what do the choices
#                   cost in numbers? Everything executed lives here.
#   **API**       — the docstrings, generated.
#
# A page that answers two of those questions is worth splitting; a page in the
# wrong chapter is worth moving. `Manual` exists because nine of the fourteen
# pages once filed under `Tutorials` were object syntax rather than a
# calculation, which is what made the chapter hard to navigate.
#
# **Cementitious media is a subsection of four chapters, not a chapter of its
# own.** The material is the subject of the package, so it appears wherever its
# question is being asked: the theory of its water budget, the syntax of its
# species, the trajectory of its hydration, the worked cases. Grouping it into
# one chapter would have separated each of those from the general treatment it
# specializes.
#
# Three paths are load-bearing and must not move: `tutorials/self_desiccation.md`,
# `tutorials/equilibrium.md` and `examples/hydration_calibration.md` are linked by
# URL from released CHANGELOG sections, which are not retro-edited.

pages = [
    "Home" => "index.md",
    "Getting Started" => "quickstart.md",
    "Theory" => [
        "theory/index.md",
        # The definitions and identities the rest is written in. Read first: the
        # remaining pages use its notation, which is the code's.
        "Foundations" => [
            "theory/thermodynamics.md",
            "theory/equilibrium.md",
        ],
        # The two places a mixture stops being ideal, and the only two where a
        # standard state has to be argued about rather than looked up.
        "Non-ideal mixtures" => [
            "theory/activity_models.md",
            "theory/solid_solutions.md",
        ],
        # The conservation law that is not an element, and the potential
        # conjugate to it. Needed by any binder whose sulfur is not all
        # sulfate -- which is every binder containing slag.
        "Oxidation state" => [
            "theory/redox.md",
        ],
        # The second conserved quantity that is not an element, and the third
        # place a mixture stops being ideal. It reads after the two above
        # because it uses both: a site balance is written like the charge row,
        # and site mixing is written like a solid solution's.
        "Surfaces and interfaces" => [
            "theory/surface_complexation.md",
        ],
        # When, rather than what: the rate laws, their parameters and the
        # provenance of every number in them. Equilibrium says nothing about
        # time, and a binder's engineering behavior is entirely about time.
        "Kinetics" => [
            "theory/kinetics.md",
        ],
        # The material this package exists for. What a Gibbs minimization can
        # predict about a drying paste, and what is not a thermodynamic
        # quantity at all.
        "Cementitious media" => [
            "theory/cement_water_budget.md",
        ],
    ],
    "Manual" => [
        # What a formula, a species and a reaction *are* here. Everything else
        # consumes these.
        "Chemical description" => [
            "manual/where_the_numbers_come_from.md",
            "manual/formula_manipulation.md",
            "manual/species.md",
            "manual/reactions.md",
            "manual/stoich_matrices.md",
        ],
        "Databases and thermodynamic data" => [
            "manual/databases.md",
            "manual/thermodynamic_data.md",
        ],
        "Systems and states" => [
            "manual/chemical_system_state.md",
        ],
        # How much area a solid offers, and to what. One page because the same
        # objects serve a dissolution rate law and, later, the sites a surface
        # binds with -- the numbers differ, the abstraction does not.
        "Surfaces" => [
            "manual/surfaces.md",
        ],
        # The cement-specific syntax: phase names, the Bogue notation, the
        # shorthand the literature uses.
        "Cementitious media" => [
            # The shorthand first: every page below writes phases in it, and a
            # reader who has not met `C3S` cannot follow them.
            # The check-list that turns a missing phase from a silent wrong
            # answer into a decision made on purpose. Every trap in it was met
            # while writing this manual.
            "manual/choosing_species.md",
            "manual/cement_notation.md",
            "manual/cement_species.md",
            # Then the map: which binder is which, what each constituent brings,
            # and which of the package's models a given family needs.
            "manual/binder_families.md",
        ],
        "Appendices" => [
            "manual/advanced.md",
        ],
    ],
    "Tutorials" => [
        "Equilibrium" => [
            "tutorials/equilibrium.md",
        ],
        # The kinetics calls the equilibrium solver, so it reads after it.
        "Kinetics and coupling" => [
            "tutorials/kinetics.md",
            "tutorials/coupling.md",
        ],
        # Uses the whole chain, which is why it comes last.
        "Cementitious media" => [
            "tutorials/self_desiccation.md",
            "tutorials/self_desiccation_kinetics.md",
        ],
        "Validation against other codes" => [
            "tutorials/reaktoro_comparison.md",
        ],
    ],
    "Applications" => [
        # From an oxide analysis to a species list — the entry point for someone
        # holding a cement datasheet rather than a database.
        "From a cement analysis to a chemical system" => [
            "examples/bogue_calculation.md",
            "examples/example_stoich_matrix.md",
            "examples/from_scratch.md",
        ],
        # The executed counterpart of the Theory chapter: what the modeling
        # choices cost, in numbers, on a composition small enough to check.
        "Non-ideal mixtures, measured" => [
            "examples/activity_models_compared.md",
            "examples/solid_solution_models.md",
            # The case a non-ideal model is written for and a single-composition
            # formulation cannot hold: the published AFm binary, run three ways.
            "examples/miscibility_gap.md",
            "examples/pitzer_model.md",
        ],
        # The smallest complete surface calculation, against the closed form it
        # is supposed to reproduce -- and it reproduces it rather than assuming
        # it, which is the whole claim.
        "Surfaces" => [
            "examples/surface_langmuir.md",
            # Two families on one support, which is what makes a sorption edge
            # bend, and the cross-code comparison that says where the remaining
            # difference lives.
            "examples/hfo_titration.md",
            "examples/hfo_diffuse_layer.md",
        ],
        # Small, checkable aqueous cases with an analytical answer to compare to.
        "Aqueous equilibria" => [
            "examples/titration_acetic_acid.md",
            "examples/titration_malonic_acid.md",
            "examples/co2_carbonate_system.md",
        ],
        # The entry point to the cement material: clinker in, everything else
        # computed.
        "Cementitious media from the clinker up" => [
            "examples/cem1_from_clinker.md",
            # The phase list as a modeling decision: declare every solid
            # solution the database defines and let the minimization choose,
            # rather than choosing for it.
            "examples/cem1_solid_solutions.md",
        ],
        "Cementitious media at equilibrium" => [
            "examples/simplified_clinker_dissolution.md",
            "examples/cement_wc_ratio.md",
            "examples/cement_carbonation.md",
        ],
        # Beyond the Portland cement. One page per EN 197-1 family, in order of
        # how much of the clinker is replaced and of what the replacement asks
        # of the models: a carbonate and a first slag (CEM II), then the slag
        # binder where the oxidation state stops being ignorable (CEM III), the
        # pozzolanic binder where the C-S-H must carry the aluminum (CEM IV),
        # and the composite where all of it holds at once (CEM V). The map of
        # the families themselves is in the manual.
        "Blended binders" => [
            "examples/cem2_blended.md",
            "examples/cem3_slag.md",
            "examples/cem4_pozzolanic.md",
            "examples/cem5_composite.md",
        ],
        "Cementitious media in time" => [
            "examples/cement_clinker_kinetics.md",
            "examples/coupled_hydration.md",
            "examples/ionic_hydration.md",
            "examples/hydration_calibration.md",
        ],
    ],
    "API" => Any[
        "Chemical description" => [
            "Formulas" => "api/formulas.md",
            "Species" => "api/species.md",
            "Parsing tools" => "api/parsing_tools.md",
            "Element order" => "api/element_order.md",
            "Stoichiometric Matrix" => "api/stoich_matrices.md",
            "Reactions" => "api/reactions.md",
            "Databases" => "api/databases.md",
        ],
        "Thermodynamics" => [
            "Thermodynamical functions" => "api/thermo_functions.md",
            "Thermodynamical models" => "api/thermo_models.md",
            "Water properties" => "api/water_properties.md",
        ],
        "Systems, states and equilibrium" => [
            "Chemical systems and states" => "api/chemical_systems.md",
            "Equilibrium" => "api/equilibrium.md",
        ],
        "Kinetics" => [
            "Kinetics" => "api/kinetics.md",
            "Surfaces" => "api/surfaces.md",
        ],
        "Utilities" => [
            "Utilities" => "api/utils.md",
        ],
    ],
    "References" => "references.md",
]
