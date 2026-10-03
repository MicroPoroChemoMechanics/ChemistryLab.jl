# [A CEM I paste replaced by fly ash, carbonated, salted and leached](@id ex-cement-processes)

!!! info "Before this page"
    [Recipes: materials, extents and what has not reacted](@ref man-recipes), and
    the materials of
    [Validation against measured blended pastes](@ref ex-validation-blended).

A **process** is a sequence of equilibria in which one quantity changes from each
equilibrium to the next: the share of the binder replaced by a fly ash, the
carbon dioxide taken up, the salt added, the water that has flowed through. Each
equilibrium starts from the answer before it, which is faster than starting cold
and keeps the sequence on one branch. This page runs the four processes of the
recipe layer on one paste.

The paste is the CEM I of [DeWeerdt2011](@citet), a clinker
interground with 3.7 % gypsum and hydrated at a water/binder ratio of 0.5 and
20 °C, at 90 days. Its clinker phases have reacted as the paper measured them by
XRD.

```@example processes
using ChemistryLab
using DynamicQuantities
using Printf
using Plots
default(framestyle = :box, grid = false)

# The materials and recipes of De Weerdt et al. (2011), shared with the
# validation page.
include(joinpath(pkgdir(ChemistryLab), "scripts", "de_weerdt_2011.jl"))
model = cemdata18_activity_model(:KOH)
nothing # hide
```

## 1. The paste and its phases

The phases a paste may form are a choice, and `data/phase_lists.toml` holds that
choice for published pastes, each list written after the paper it names. This
list is written for a Portland cement paste:

```@example processes
list = "Portland paste (Lothenbach and Winnefeld 2006)"
pl = phase_list(list)
println("reactants: ", join(pl.reactants, ", "))
println("products:  ", join(pl.products, ", "))
println("solid solutions: ", join(pl.solid_solutions, ", "))
```

The list has no chloride phase. Section 4 adds salt, so this page departs from
the list, and says so: `add` declares Friedel's salt (`C4AClH10`) and Kuzel's
salt (`C4AsClH12`), the AFm phases that bind chloride. Chlorine then enters the
system, and so do its aqueous species.

```@example processes
cs = phase_list_system(list, DW11_SUBSTANCES; add = ["C4AClH10", "C4AsClH12"])
recipe = dw11_recipe("OPC")
paste, cert = equilibrate_certified(recipe, cs; t = 90, model)
cert.optimal
```

Two helpers read what the sections below compare: the mass of a few solids, in
grams per 100 g of binder, and the C-S-H, its mass and its calcium-to-silicon
ratio.

```@example processes
solids(rs, names) = [get(phase_masses(rs), n, 0.0) for n in names]
function csh(rs)
    gel = solid_solution_totals(rs.state, "CSHQ")
    si = get(gel.elements, :Si, 0.0)
    return (g = 1000 * ustrip(gel.mass), ca_si = si > 0 ? gel.elements[:Ca] / si : NaN)
end
ratio(x) = isnan(x) ? "    –" : @sprintf("%5.2f", x)
@printf("portlandite %.1f g, C-S-H %.1f g at Ca/Si %.2f, pH %.2f\n",
        solids(paste, ["Portlandite"])[1], csh(paste).g, csh(paste).ca_si, pH(paste.state, model))
```

## 2. Replacing part of the cement by fly ash

[`blend`](@ref) replaces the fraction `f` of the binder by another material and
keeps the binder mass and the water/binder ratio. Here the material is the
siliceous fly ash of the same paper. Its glass reacts at the rate the paper
fitted, 30 % of the fly ash at 90 days, and its crystals stay inert. The
clinker keeps the degrees of reaction measured in the CEM I, although in the
blended pastes of the paper the clinker reacted somewhat faster.

```@example processes
fractions = 0.0:0.1:0.5
bl = blend(recipe, dw11_fly_ash(), fractions, cs; t = 90, model)
println("  f     portlandite   C-S-H   Ca/Si    pH")
for (f, rs) in zip(fractions, bl)
    @printf("%4.1f  %9.1f g  %6.1f g  %5.2f  %5.2f\n", f, solids(rs, ["Portlandite"])[1], csh(rs).g, csh(rs).ca_si,
            pH(rs.state, model))
end
```

Portlandite falls faster than the cement is removed. With 40 % of fly ash, the
60 % of cement left would make 16.7 g of it, and 2.3 g remain: the silica of the
glass has taken the rest to form C-S-H, which is the pozzolanic reaction. The
C-S-H therefore falls much less than the cement, from 57.1 g to 54.0 g where the
cement alone would make 34 g. Its Ca/Si stays at 1.58, the value in equilibrium
with portlandite, for as long as portlandite remains; at 50 % portlandite is gone
and the gel starts to lose calcium. [The validation page](@ref ex-validation-blended)
measured how far this goes against the real pastes: at 35 % the model consumed
portlandite faster than the pastes did, so this sweep overstates it too.

## 3. Carbonation

[`carbonate`](@ref) adds carbon dioxide to the budget of the paste, here up to
1 mol per 100 g of binder. The unreacted clinker stays aside.

```@example processes
co2 = 0.0:0.1:1.0
co = carbonate(paste, co2)
println("CO2 (mol)  portlandite  calcite   C-S-H   Ca/Si    pH")
for (x, rs) in zip(co2, co)
    p, c = solids(rs, ["Portlandite", "Cal"])
    @printf("%6.1f  %10.1f g  %6.1f g  %6.1f g  %s  %5.2f\n", x, p, c, csh(rs).g, ratio(csh(rs).ca_si),
            pH(rs.state, model))
end
```

The carbon dioxide first turns portlandite into calcite, which is gone by
0.4 mol. The pH hardly moves meanwhile, from 13.64 to 13.62: it is set by the
sodium and potassium of the pore solution, not by portlandite. Calcite then takes
its calcium from the C-S-H. The Ca/Si of the gel falls from 1.58 to 0.67 by
0.8 mol, and the pH with it, to 9.90. At a Ca/Si of 0.67 the gel is its most
calcium-poor end member, and the pH holds at that value while it is carbonated in
turn. At 1.0 mol no C-S-H is left, and the pH is 7.99. The
[carbonation page](@ref sec-cement-carbonation) compares this sequence with the
calculation and the measurements of Shi et al. on four mortars.

```@example processes
names = ["Portlandite", "Cal", "ettringite", "monocarbonate", "Gp", "Amor-Sl"]
p1 = plot(; xlabel = "CO₂ added (mol per 100 g of binder)", ylabel = "g per 100 g of binder", legend = :topleft)
for n in names
    plot!(p1, co2, [solids(rs, [n])[1] for rs in co]; label = n, linewidth = 2)
end
plot!(p1, co2, [csh(rs).g for rs in co]; label = "C-S-H", linewidth = 2, color = :black)
p2 = plot(co2, [pH(rs.state, model) for rs in co]; xlabel = "CO₂ added (mol per 100 g of binder)",
          ylabel = "pH", legend = false, linewidth = 2, color = :black)
fig = plot(p1, p2; layout = (1, 2), size = (900, 380), left_margin = 6Plots.mm, bottom_margin = 6Plots.mm)
savefig(fig, "processes-carbonation.svg"); nothing # hide
```

![](processes-carbonation.svg)

## 4. Salt

[`add_salt`](@ref) adds a salt by its formula. Sodium chloride is not a phase of
this system, and its ions enter together, so the budget stays neutral. Section 1
declared Friedel's and Kuzel's salts, so the paste can bind chloride. Every
equilibrium's certificate also reports its ionic strength against the range
the Cemdata18 activity model is stated for, about 1 mol/kg.

```@example processes
nacl = 0.0:0.01:0.06
sa = add_salt(paste, "NaCl", nacl)
water_kg(rs) = ustrip(us"kg", mass(rs.state, cs.species[only(cs.idx_solvent)]))
afm(rs) = sum(m for (k, m) in phase_masses(rs) if startswith(k, "monosulphate12") || startswith(k, "C4AH13"); init = 0.0)
println("NaCl (mol)  free Cl (mol/kg)  bound Cl   Friedel  Kuzel  monosulfate   I (mol/kg)  in range")
for (x, rs) in zip(nacl, sa)
    free = get(pore_solution(rs).elements, :Cl, 0.0)
    bound = x - free * water_kg(rs)
    f, k = solids(rs, ["C4AClH10", "C4AsClH12"])
    m = afm(rs)
    @printf("%7.2f  %12.3f  %11.1f %%  %6.1f g  %5.1f g  %8.1f g  %9.2f   %s\n", x, free,
            x > 0 ? 100 * bound / x : 0.0, f, k, m, rs.certificate.ionic_strength, rs.certificate.within_activity_range)
end
```

This paste holds no calcite, so its AFm is the monosulfate, and the first
chloride goes into Kuzel's salt, which is formed from it: 81 to 89 % of the
chloride added is bound. Friedel's salt appears from 0.03 mol, and the Kuzel's
salt then gives way to it. From 0.05 mol the Kuzel's salt is gone and Friedel's
salt stops growing, at 10.3 g, so what is added stays in solution: 61 % of the
chloride is bound at 0.06 mol. From
0.04 mol the ionic strength is past the range of the model, and the certificate
says so. The answer is certified, but for an activity model applied outside its
stated range.

## 5. Leaching

[`leach`](@ref) removes the whole pore solution and replaces it with fresh water,
step after step, as water flowing through a paste would. Each step here uses 5 kg
of water, a hundred times the water of the paste, so that a few steps show the
whole sequence.

```@example processes
le = leach(paste, 12; renewal = 5000u"g")
println("step   portlandite   C-S-H   Ca/Si    pH")
for (k, rs) in enumerate(le)
    @printf("%4d  %9.1f g  %6.1f g  %5.2f  %5.2f\n", k, solids(rs, ["Portlandite"])[1], csh(rs).g, csh(rs).ca_si,
            pH(rs.state, model))
end
```

The first renewal takes the sodium and potassium away, and the pH falls from
13.64 to 12.68, the value set by portlandite. The Ca/Si of the gel rises at the
same time, from 1.58 to 1.63, because its alkali end members, which hold silicon
and no calcium, leave with the alkalis. Portlandite then dissolves step by step
at an almost constant pH, and is gone at the fourth renewal. The C-S-H then loses
its calcium in turn: its Ca/Si falls to 0.97 by the twelfth renewal, and the pH
to 11.71. A paste in contact with flowing soft water degrades in this order.

## Where to go next

The materials, and how the model compares with the measured pastes, are on
[Validation against measured blended pastes](@ref ex-validation-blended). The
layer itself is described on [Recipes](@ref man-recipes).
