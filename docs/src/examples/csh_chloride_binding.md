# [Chloride binding by C-S-H and Friedel's salt](@id sec-example-csh-chloride)

!!! info "Before this page"
    [A charged surface, screened](@ref sec-example-diffuse-layer), and
    [Chemistry that happens on a surface](@ref sec-theory-surface) §8 and §9.

A hydrated paste takes chloride out of its pore solution in two ways. The AFm
hydrates exchange their sulfate for it and become Kuzel's salt, then Friedel's
salt; the C-S-H adsorbs it on its silanol sites. The two compete for the same
chloride and neither can be read off the other, so a calculation that keeps only
the salts assigns all of the bound chloride to them.

The surface model is that of [Elakneswaran2010](@cite), as [Guo2018](@cite) use
it: silanol sites `≡SiOH` on the C-S-H, which deprotonate and bind calcium,
sodium and chloride, with the potential set by a diffuse layer. The proton,
calcium and chloride constants were fitted to zeta potentials by
[Elakneswaran2009](@cite), and the first two lie close to the values
[Pointeau2006](@cite) obtained by titration, −12.0 and −9.2 against −12.7 and
−9.4.

!!! warning "What this calculation is, and what it is not"
      - The C-S-H has a **fixed composition** and only carries the sites: it
        neither dissolves nor precipitates, and it is not the CSHQ solid
        solution of the cement pages. Putting the sites on CSHQ would count its
        calcium and its alkalis twice, once in the solid and once on the surface.
      - The diffuse layer sets the potential, following Dzombak and Morel as
        PHREEQC does by default. The ions it holds are counted only in
        [the last section](@ref "The chloride of the diffuse layer"), and there
        with a thickness that nothing published fixes for this paste.
      - The hydrates are those of CEMDATA18, and Guo used Cemdata07. This is
        their model on the package's database, not a reproduction of their
        figure; [the validation chapter](@ref "Validation against published data") compares
        the salts with that figure on its own terms.

## The paste

Guo state their inventory per liter of concrete (C-S-H, portlandite, AFm, AFt,
and the porosity), and every amount below is brought to one kilogram of pore
water. The site density and the specific area of the C-S-H are theirs, from
the same file.

```@example cshcl
using ChemistryLab, DynamicQuantities, OptimaSolver, Printf

const DB = Dict(symbol(s) => s for s in
                build_species(datapath("cemdata18-thermofun.json"); verbose = false))
ions = [DB[s] for s in ("H2O@", "H+", "OH-", "Na+", "Cl-", "Ca+2", "AlO2-", "SO4-2")]
hydrates = [DB[s] for s in ("Portlandite", "monosulphate12", "monosulphate14",
                            "ettringite", "C4AClH10", "C4AsClH12")]

guo(q) = literature_value("Guo2018", q)
# A liter of concrete holds 10 × porosity grams of pore water.
const SCALE = 1000 / (10 * ustrip(guo("porosity_percent")))
grams = Dict("Portlandite" => guo("ch_per_liter"),
             "monosulphate12" => guo("afm_per_liter"),
             "ettringite" => guo("aft_per_liter"))
csh_g = ustrip(us"g", guo("csh_per_liter")) * SCALE
n_sites = ustrip(u"mol/g", guo("csh_site_density")) * csh_g
area = ustrip(u"m^2/g", guo("csh_specific_surface_area")) * csh_g
@printf("C-S-H %.0f g, %.2f mol of sites on %.3g m², per kg of pore water\n",
        csh_g, n_sites, area)
```

The surface reactions are read from the table of the same file, written in
PHREEQC's syntax, and [`site_family`](@ref) turns each constant into the
standard energy of its complex. The potassium row is left out, since no
potassium enters here.

```@example cshcl
t = literature_table("Guo2018", "surface_reactions_phreeqc")
reactions = [eq => lk for (eq, lk) in zip(t.reaction, t.log_K) if !occursin("K+", eq)]
for (eq, lk) in reactions
    @printf("  %-34s log K = %6.2f\n", eq, lk)
end

dl = DiffuseLayer(; area)
family = site_family(
    "Csh_w", reactions, ions; master = "Csh_w", site = "Xw",
    capacity = TotalSiteAmount(n_sites * u"mol"),
    support = SurfaceSupport("C-S-H", nothing, FixedSurfaceArea(area)), model = dl,
)
members = vcat([family.free_site], family.complexes)
cs = ChemicalSystem(
    vcat(ions, hydrates, members),
    [DB["H2O@"], DB["H+"], DB["Na+"], DB["Ca+2"], DB["Cl-"], DB["AlO2-"],
     DB["SO4-2"], family.free_site];
    site_families = [family],
)
idx = Dict(symbol(sp) => k for (k, sp) in enumerate(cs.species))
nothing # hide
```

## A sweep in sodium chloride

The paste is charged with sodium chloride up to 0.4 mol per kilogram of water
and solved at each step, from the saltiest point down, each solve starting from
the previous answer. Every point is certified.

```@example cshcl
function paste(nacl)
    n = Any[fill(1.0e-14u"mol", length(cs.species))...]
    n[idx["H2O@"]] = 1.0u"kg" / DB["H2O@"][:M]
    for (ph, g) in grams
        n[idx[ph]] = g * SCALE / DB[ph][:M]
    end
    n[idx["XwOH"]] = n_sites * u"mol"
    n[idx["Na+"]] = n[idx["Cl-"]] = max(nacl, 1.0e-10) * u"mol"
    return ChemicalState(cs, n)
end

model = DaviesActivityModel()
NACL = [0.4, 0.3, 0.2, 0.15, 0.1, 0.07, 0.05, 0.03, 0.02, 0.01, 0.0]
A = Float64.(cs.SM.A)
budgets = [A * ustrip.(us"mol", paste(c).n) for c in NACL]
states, certs = equilibrate_path(paste(first(NACL)), budgets; model)
println("certified: ", count(c -> c.optimal, certs), " of ", length(certs))
```

The chloride out of solution sits in three places: two per Friedel's salt, one
per Kuzel's salt, and one per `≡SiOHCl⁻` complex.

```@example cshcl
const MV = 1000 * R_GAS * 298.15 / FARADAY          # mV per unit of FΨ/RT
mol(eq, s) = ustrip(us"mol", eq.n[idx[s]])
z = Float64[charge(sp) for sp in members]
function partition(eq)
    friedel = 2mol(eq, "C4AClH10")
    kuzel = mol(eq, "C4AsClH12")
    surface = mol(eq, "XwOHCl-")
    psi = diffuse_layer_potential(dl, z, [mol(eq, symbol(sp)) for sp in members],
                                  ionic_strength(eq), 298.15)
    return (; friedel, kuzel, surface, bound = friedel + kuzel + surface,
            pH = pH(eq, model), psi_mV = psi * MV)
end
rows = [partition(eq) for eq in states]

println(" NaCl     pH     Ψ/mV   Cl in Friedel  Kuzel   surface   bound   surface share")
for (c, r) in zip(NACL, rows)
    share = r.bound > 1.0e-9 ? @sprintf("%5.1f %%", 100r.surface / r.bound) : "   –"
    @printf("%5.2f  %6.3f  %6.1f  %9.4f  %9.4f  %8.4f  %7.4f   %s\n",
            c, r.pH, r.psi_mV, r.friedel, r.kuzel, r.surface, r.bound, share)
end
```

## Against PHREEQC

The same paste is computed by PHREEQC on the same model: the same reactions
and constants, the hydrates with the log K the package gives them, the Davies
equation on both sides and Dzombak and Morel's diffuse layer. The generator,
`test/reference/phreeqc_csh_surface.py`, asks the package for every energy and
molar mass it uses, so the comparison is of the two solvers and of nothing else.

```@example cshcl
using JSON
ORACLE = JSON.parsefile(joinpath(pkgdir(ChemistryLab), "test", "reference",
                                 "phreeqc_csh_paste.json"))
worst = (phase = 0.0, surface = 0.0, pH = 0.0, psi_mV = 0.0)
for pt in ORACLE["points"]
    k = findfirst(≈(pt["nacl"]), NACL)
    eq, r = states[k], rows[k]
    for (ph, m) in pt["phases"]
        global worst = merge(worst, (phase = max(worst.phase, abs(mol(eq, ph) - m)),))
    end
    for (sp, m) in pt["surface"]
        s = replace(sp, "Csh_w" => "Xw"; count = 1)
        global worst = merge(worst, (surface = max(worst.surface, abs(mol(eq, s) - m)),))
    end
    global worst = merge(worst, (pH = max(worst.pH, abs(r.pH - pt["pH"])),
                                 psi_mV = max(worst.psi_mV, abs(r.psi_mV - 1000pt["psi_V"]))))
end
@printf("over %d points: %.1e mol on a phase, %.1e mol on a surface species,\n",
        length(ORACLE["points"]), worst.phase, worst.surface)
@printf("%.1e on the pH, %.3f mV on the potential\n", worst.pH, worst.psi_mV)
```

## Reading the partition

The surface is positive throughout. Calcium outnumbers the deprotonated sites
and reverses the sign that a silanol surface has in a calcium-free solution,
which is what lets chloride bind to it at all. Salt then lowers the potential by
screening it and by adding negative complexes.

The surface takes the first chloride. Up to 0.03 mol/kg no salt has formed and
all of the bound chloride is on the C-S-H. The monosulfate then turns into
Kuzel's salt, which holds one chloride per formula unit, and between 0.3 and
0.4 mol/kg Kuzel's salt gives way to Friedel's salt, which holds two.

From 0.05 to 0.1 mol/kg the chloride on the surface does not move, while
Kuzel's salt grows fivefold. The monosulfate and Kuzel's salt coexist there,
with ettringite and portlandite, and coexisting hydrates fix the activity of the
chloride they exchange: the salt added goes into Kuzel's salt, and the surface,
which responds to activities alone, sees no change. Once the monosulfate is
gone the activity rises again, and so does the surface's share, to 44 % of the
bound chloride at 0.3 mol/kg. A model that kept only the salts would miss a
third of the bound chloride at the top of the sweep and all of it at the
bottom.

The chloride retained in the diffuse layer is the term missing from these
totals. With a positive surface every anion is in excess in the layer, so the
next section counts it.

## The chloride of the diffuse layer

A [`DonnanLayer`](@ref) makes that layer explicit, as PHREEQC's `SURFACE -Donnan`
does, after [AppeloWersin2007](@cite): a layer of water of fixed thickness on the
surface, holding each solute at
the average Boltzmann enrichment whose charge balances the surface's, while the
surface keeps its Gouy-Chapman potential. [`equilibrate_donnan`](@ref) withdraws
what the layer holds from the solution and solves again until the two agree.
Against PHREEQC, on the surface alone in the eighteen solutions of the test
suite, the layer's chloride agrees to 3 × 10⁻⁴ relative.

The thickness is the parameter of the approach, and nothing published fixes it
here. One Debye length of the solution, computed from its ionic strength, is
taken below; the layer's water is taken from the pore solution, as it must be
in a closed paste.

```@example cshcl
debye_length(I) = sqrt(water_relative_permittivity(298.15, 1.0e5) * VACUUM_PERMITTIVITY *
                       R_GAS * 298.15 / (2 * FARADAY^2 * I * 1000))
println(" NaCl  thickness  layer water  Cl excess  surface  Kuzel  Friedel   bound  (no layer)")
for c in (0.4, 0.2, 0.1)
    k = findfirst(==(c), NACL)
    t = debye_length(ionic_strength(states[k]))
    res = equilibrate_donnan(states[k], DonnanLayer(thickness = t); model, water = :taken)
    r = partition(res.state)
    ex = res.layer.excess[idx["Cl-"]]
    @printf("%5.2f  %6.2f nm  %8.3f kg  %9.4f  %7.4f  %6.4f  %7.4f  %6.4f  (%6.4f)\n",
            c, 1.0e9t, res.layer.water, ex, r.surface, r.kuzel, r.friedel, r.bound + ex,
            rows[k].bound)
end
```

The layer adds 2 to 8 % to the bound chloride. At 0.1 mol/kg its excess is two
thirds of what the surface complexes hold, and part of it is chloride the salts
would otherwise have taken: Kuzel's salt falls from 0.0555 to 0.0503 mol, since
the layer lowers the concentration the AFm phases see. A layer one Debye length
thick holds 41 to 70 % of the pore water here, because the pores of this paste
are barely wider than the layer: this is where a Donnan average is a coarse
description, and where the number depends on the thickness chosen more than on
anything measured.

## See also

  - [A charged surface, screened](@ref sec-example-diffuse-layer) — the diffuse
    layer checked alone, on ferrihydrite.
  - [Chemistry that happens on a surface](@ref sec-theory-surface) — the site
    families and the potential.
  - [Validation against published data](@ref) — Guo's
    salts against their figure, without the surface.
  - `test/csh_surface.jl` — the Donnan layer against PHREEQC's `-Donnan`.
