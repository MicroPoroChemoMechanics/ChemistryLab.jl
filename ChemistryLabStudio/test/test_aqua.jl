using Aqua

# `persistent_tasks` precompiles the package in an environment resolved from the
# registry, which cannot hold a package before its first registration: it is
# enabled once ChemistryLabStudio is in General.
Aqua.test_all(ChemistryLabStudio; persistent_tasks = false)
