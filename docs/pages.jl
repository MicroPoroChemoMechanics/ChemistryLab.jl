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
            # Three pages in the order a first reader's questions arise. The two
            # laws written as balances on a boundary: the enthalpy a calorimeter
            # measures, the Gibbs energy the solver minimizes, and what each
            # one is not.
            "theory/the_two_laws.md",
            # The composition allowed to change: the chemical potential, Euler and
            # Gibbs-Duhem, the Gibbs energy of a reaction, the equilibrium
            # constant, and why a metastable state is not a local minimum.
            "theory/energies_and_potentials.md",
            # What a database stores and how it is measured: formation
            # quantities, entropies, primaries, temperature and pressure, the
            # apparent convention. Reference material; a first reading can
            # come back to it.
            "theory/formation_quantities.md",
            "theory/thermodynamics.md",
            # What an activity is measured from. Before the certificate, because
            # every potential the certificate compares rests on these conventions.
            "theory/standard_states.md",
            "theory/equilibrium.md",
        ],
        # The places a mixture stops being ideal: the solution and the solids,
        # where a standard state has to be argued about rather than looked up,
        # and the gas under pressure, which keeps its standard state and gains a
        # fugacity coefficient.
        "Non-ideal mixtures" => [
            "theory/activity_models.md",
            "theory/solid_solutions.md",
            "theory/real_gases.md",
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
            # How a run under partial equilibrium is integrated: the partition and
            # why the state carries element amounts, the right-hand side and its
            # exact Jacobian, the implicit step, and the energy balance of a
            # calorimeter. The tutorials and applications run it and point here.
            "theory/partial_equilibrium_kinetics.md",
        ],
        # The material this package exists for. What a Gibbs minimization can
        # predict about a drying paste, and what is not a thermodynamic
        # quantity at all.
        "Cementitious media" => [
            "theory/cement_water_budget.md",
            # From materials, masses and degrees of reaction to the budget the
            # minimization conserves, and what is kept out of it.
            "theory/recipe_bookkeeping.md",
            # The enthalpy of a glass that has no formula, from measured
            # silicate glasses, for the heat of a blend.
            "theory/glass_thermochemistry.md",
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
        # How a database of any format is read, then what is built on it: the
        # extensions, the filters and the solid solutions, and the functions of
        # temperature its species carry.
        "Databases and thermodynamic data" => [
            "manual/importing_databases.md",
            "manual/databases.md",
            "manual/thermodynamic_data.md",
        ],
        "Systems and states" => [
            "manual/chemical_system_state.md",
        ],
        # The options of a calculation, split out of the two tutorials so that
        # a tutorial reads as a calculation and the reference is found in one
        # place. Both assume the corresponding tutorial has been read.
        "Solving and integrating" => [
            "manual/solving.md",
            "manual/kinetics_syntax.md",
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
            # From materials to a budget: extents, the unreacted residue, and the
            # processes built on the certified equilibrium.
            "manual/recipes.md",
        ],
        "Appendices" => [
            "manual/advanced.md",
        ],
    ],
    "Tutorials" => [
        "Equilibrium" => [
            "tutorials/equilibrium.md",
            # Solid solutions as candidates: which are present, with which
            # composition, what the user chooses, and how to read the certificate.
            "tutorials/solid_solutions.md",
        ],
        # The kinetics calls the equilibrium solver, so it reads after it.
        "Kinetics and coupling" => [
            "tutorials/kinetics.md",
            "tutorials/coupling.md",
            # A site family at a rate: needs the coupling and a surface.
            "tutorials/slow_surface.md",
        ],
        # Uses the whole chain, which is why it comes last.
        "Cementitious media" => [
            "tutorials/self_desiccation.md",
            "tutorials/self_desiccation_kinetics.md",
        ],
        # Two halves of the same question, and neither substitutes for the
        # other: a second code reading the same database cannot see an error in
        # the database, and a paper's own tables cannot see an error in the
        # solver.
        "Validation" => [
            "tutorials/reaktoro_comparison.md",
            "tutorials/published_data_validation.md",
            # The temperature dependence of the database's constants, against
            # the fits to measured constants that PHREEQC ships.
            "tutorials/validation_logk_temperature.md",
            # The pressure: carbon dioxide dissolved up to 500 atm, the gas
            # ideal and following an equation of state.
            "tutorials/validation_co2_solubility.md",
            # The range: the ionization constant of water and the solubility of
            # quartz to 1000 °C and 5 kbar, against fits to the measurements.
            "tutorials/validation_high_temperature.md",
            # A measured paste through its first year, and the same budgets
            # through GEMS3K: the model against the paste, the code against the code.
            "tutorials/validation_measured_pastes.md",
            # Alkali uptake by C-A-S-H batches: CSHQ predicted on data it was not
            # fitted on, CASH+NK checked on data it was.
            "tutorials/validation_alkali_uptake.md",
            # And by C-S-H of the aluminum, the companion syntheses: CNASH_ss, the
            # gel that takes it, against CSHQ, which does not.
            "tutorials/validation_aluminum_uptake.md",
            # Four blended pastes, their degrees of reaction from the measurement,
            # the chemistry against the measurement and against GEMS3K.
            "tutorials/validation_blended_pastes.md",
            # A low-pH shotcrete, 40 % silica fume, and the formate of its set
            # accelerator: the alkali end-members of CASH+NK on a gel of low Ca/Si.
            "tutorials/validation_silica_fume_paste.md",
            # Forty-eight pore solutions of the first six hours, speciated at
            # their measured pH: the aqueous model against the paper's indices.
            "tutorials/validation_early_pore_solutions.md",
            "tutorials/validation_fly_ash_pore_solutions.md",
            # The same pastes cured from 7 to 80 °C: the temperature dependence
            # of the activity model and of the solubility products.
            "tutorials/validation_temperature_pore_solutions.md",
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
            # Two gels mixed on their sites: CSH3T both ways, and the chain
            # length of the CNASH gel.
            "examples/sublattice_csh.md",
            # The CASH+ gel in the compound energy formalism: the invariant points
            # of the paper, the alkalis, and the authors' 110 gel compositions.
            "examples/cashplus_csh.md",
            # A gel that holds aluminum at a high Ca/Si: two extensions fitted
            # here, tested on syntheses and pastes, and why neither is shipped.
            "examples/csh_aluminum.md",
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
            # And the case a fixed budget cannot describe: the solid carrying
            # the sites is itself dissolving, so the budget has to follow it.
            "examples/evolving_sorbent.md",
            "examples/claysor_clay.md",
            # A surface in a cement paste, competing with the AFm salts for the
            # same chloride, against PHREEQC on the same model.
            "examples/csh_chloride_binding.md",
            # The same surface on a C-S-H frozen after a first equilibrium with
            # CSHQ, and the chloride end member of CSHQ, in slag cements.
            "examples/chloride_binding_blended.md",
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
            # The four processes of the recipe layer on one measured paste.
            "examples/cement_processes.md",
            # The hydrates of two cements from 0 to 60 °C, against the same
            # calculation on cemdata2007.
            "examples/hydrates_temperature.md",
            "examples/chloride_temperature.md",
        ],
        # What the surroundings do to a hydrated paste, at equilibrium and in
        # zero dimensions: a sequence of equilibria in which what enters or
        # leaves the paste is the parameter, each against measurements.
        "Durability in zero dimensions" => [
            "examples/leaching.md",
            "examples/sulfate_attack.md",
            "examples/seawater.md",
            "examples/seawater_flushing.md",
            "examples/delayed_ettringite.md",
            "examples/hemicarbonate.md",
            "examples/chloride_afm.md",
            # The products of the alkali-silica reaction synthesized at 80 °C:
            # two sets of constants for them, and an equilibrium on the middle
            # root of the ionic strength.
            "examples/asr_products.md",
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
        # The semi-adiabatic calorimeter belongs here rather than with the
        # outputs below: its temperature is an unknown of the kinetics, raised
        # by the heat and raising the rates in turn.
        "Cementitious media in time" => [
            "examples/cement_clinker_kinetics.md",
            "examples/coupled_hydration.md",
            "examples/ionic_hydration.md",
            "examples/hydration_calibration.md",
            # The comparisons of Lavergne et al. (2018), the semi-adiabatic
            # calorimeter among them.
            "examples/lavergne2018.md",
            # The first blended cement integrated in time: the slag pastes the
            # two output pages read at measured degrees of hydration.
            "examples/blended_slag_kinetics.md",
            # Four materials, ten blends, two glasses under their own laws and
            # the limestone at equilibrium, against the thermogravimetry.
            "examples/quaternary_kinetics.md",
            "examples/ternary_kinetics.md",
            # The clinker and the slag laws at 5, 20 and 40 °C, against the
            # degrees of reaction measured at each temperature.
            "examples/slag_temperature.md",
            "examples/slag_temperature_pastes.md",
            # The glasses of the supplementary materials dissolving at pH 13,
            # and what dissolved calcium and aluminum do to them.
            "examples/glass_dissolution.md",
        ],
        # What a laboratory measures on a paste, read off computed states: the
        # heat an isothermal calorimeter records and the mass a thermobalance
        # loses. Neither feeds back on the calculation.
        "Outputs of a calculation" => [
            "examples/isothermal_calorimetry.md",
            # The heat of a blend at its measured degrees of reaction, the slag
            # glass at the enthalpy of the measured glasses.
            "examples/slag_heat.md",
            "examples/thermogravimetry.md",
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
    "Nomenclature" => "nomenclature.md",
    "References" => "references.md",
]

"""
    page_leaves(node) -> Vector{String}

Every page path under `node`, in build order. A leaf is a bare path or the value
of a `"Title" => "path"` pair; a `"Title" => [...]` pair is a section. The full
build, the partial builds, the shards of the CI and `scripts/docs_timing.jl` all
walk the tree with it.
"""
page_leaves(node::AbstractString) = [node]
page_leaves(node::Pair) = page_leaves(node.second)
page_leaves(node::AbstractVector) = reduce(vcat, page_leaves.(node); init = String[])
