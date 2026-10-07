# [The enthalpy of a glass](@id sec-theory-glass)

The heat a paste releases is the fall of its enthalpy, and the enthalpy of a
paste that still holds unreacted slag or fly ash includes that of their glass.
A clinker phase has a formula and a database record; a glass has neither. This
page builds the enthalpy of a glass from measurements on silicate glasses, and
says at each step what is measured and what is assumed. The function is
[`glass_enthalpy`](@ref); its use in a recipe is in
[Recipes: materials, extents and what has not reacted](@ref man-recipes).

## 1. What is measured

Two properties of silicate glasses have been measured on enough compositions
to be used here.

**The heat capacity.** Below the glass transition, the heat capacity of a
silicate glass is the sum of contributions of its oxides, each independent of
the composition [Richet1987](@cite):

```math
C_p = \sum_i x_i\, C_{p,i}(T),
\qquad
C_{p,i} = a_i + b_i T + \frac{c_i}{T^2} + \frac{d_i}{\sqrt{T}} ,
```

``x_i`` the mole fraction of oxide ``i``, with coefficients for eleven oxides
fitted on 36 glasses, to about 1 % between 270 K and the glass transition.

**The enthalpy of vitrification**, ``\Delta H_v``, the enthalpy by which a glass
exceeds its crystal of the same composition, from the difference of their
enthalpies of solution in a calorimeter. At 298 K it is measured for six
glasses near a slag's composition, compiled by [RichetBottinga1986](@citet),
[RichetBottinga1984](@citet) and [Navrotsky1980](@citet):

```@example glass
using ChemistryLab, DynamicQuantities, Printf
for g in ChemistryLab._vitrification_enthalpies()
    @printf("%-19s %-11s ΔHv = %5.1f kJ/mol%s  (%s)\n", g.name, g.formula, g.dH_v / 1000,
            isnan(g.uncertainty) ? "       " : @sprintf(" ± %3.1f", g.uncertainty / 1000), g.key)
end
```

The crystal of CaSiO₃ is pseudowollastonite, the one the calorimetry used.

## 2. The construction

Write the analysis as ``n_{ox}`` moles of each oxide per gram. The CaO, MgO,
Al₂O₃ and SiO₂ part is written as ``\nu_k`` moles of the six measured glasses
plus ``s_{ox}`` moles of oxide that no measured glass takes,

```math
\sum_k a_{k,ox}\, \nu_k + s_{ox} = n_{ox},
\qquad \nu_k \ge 0,\quad s_{ox} \ge 0 ,
```

``a_{k,ox}`` the moles of oxide ``ox`` in one formula of glass ``k``. Two
assumptions are made. The glasses are mixed ideally, with no enthalpy of
mixing. And the oxide left over is counted as its crystal, with an enthalpy of
vitrification of zero. The ``\nu_k`` are chosen to leave as little as possible
to the second, that is to minimize ``\sum_{ox} s_{ox}``. This is a linear
program with four equations; its solutions are among the vertices of the
polytope, the choices of four columns whose square system has a non-negative
solution, and with ten columns there are 210 to try.

The enthalpy of formation of the glass is then

```math
\Delta_f H_{\text{glass}} = \sum_{ox} n_{ox}\, \Delta_f H_{ox}
 + \sum_k \nu_k \Bigl[\Delta_f H_k - \sum_{ox} a_{k,ox} \Delta_f H_{ox}\Bigr]
 + \sum_k \nu_k\, \Delta H_{v,k} ,
```

the oxides, the formation of each crystal from its oxides, and its
vitrification. The bracket is taken within one database, aq17, which holds the
six crystals and the four oxides: the enthalpy of forming a crystal from its
oxides is consistent only when both come from the same assessment. The
databases do not agree on the crystals themselves:

```@example glass
aq = ChemistryLab._aq17_crystals()
slop = Dict(symbol(s) => s for s in build_species(datapath("slop98-inorganic-thermofun.json"); verbose = false))
H(sp) = ustrip(us"kJ/mol", sp[:ΔₐH⁰](T = 298.15u"K", P = 1.0e5u"Pa"; unit = true))
for (a, s) in (("Gehlenite", "Gh"), ("Akermanite", "Ak"), ("Anorthite", "An"), ("Diopside", "Di"), ("Corundum", "Crn"))
    @printf("%-11s aq17 %9.2f   slop98 %9.2f kJ/mol   difference %6.2f\n", a, H(aq[a]), H(slop[s]), H(aq[a]) - H(slop[s]))
end
```

so the choice of database moves a glass by up to twenty kilojoules per mole of
crystal, several times the uncertainty of its enthalpy of vitrification. The
enthalpies of formation of the oxides themselves, ``\Delta_f H_{ox}``, come
from a reference the caller chooses, so that the glass can be put on the scale
of the database of the calculation it enters; the oxides of the analysis other
than the four, Fe₂O₃, Na₂O, K₂O, TiO₂, MnO, are counted as their crystals, with
no enthalpy of vitrification (assumed). At a temperature other than 298.15 K,
the heat capacity of section 1 is integrated, an oxide it does not cover taking
that of its crystal, as [Richet1987](@citet) prescribes.

## 3. What ideal mixing leaves open

When several combinations leave the same amount of oxide over, ideal mixing
should give them the same enthalpy, and the measured glasses do not. The
composition Ca₃MgAl₂Si₄O₁₅ is anorthite plus akermanite, and also gehlenite
plus diopside plus silica, and also half gehlenite, half anorthite, diopside
and half pseudowollastonite; every combination with nothing left over lies
between these three. Their enthalpies of vitrification are 104, 151 and
168 kJ per formula. The function returns the midpoint of the extremes and their
half-difference as `span`:

```@example glass
M(f) = ustrip(us"g/mol", Species(f)[:M])
mol = Dict("CaO" => 3, "MgO" => 1, "Al2O3" => 1, "SiO2" => 4)
Mf = sum(n * M(ox) for (ox, n) in mol)
g = glass_enthalpy(Dict(ox => n * M(ox) / Mf for (ox, n) in mol))
@printf("vitrification %.1f ± %.1f kJ per Ca3MgAl2Si4O15\n",
        ustrip(us"J/g", g.vitrification) * Mf / 1000, ustrip(us"J/g", g.span) * Mf / 1000)
```

A span of 32 kJ per formula is large against the enthalpies of mixing measured
between silicate glasses, which [RichetBottinga1986](@citet) find below about
5 kJ/mol near 1000 K. It is larger too than what they attribute to the fictive
temperature of a glass, the temperature its structure was frozen at: an error
of 100 K on it moves the enthalpy of a diopside or anorthite glass by about
10 kJ/mol. The six measurements were
made on glasses of different thermal histories, by several laboratories, over
half a century: the span is the measure of how far they can be combined, and
it is reported rather than hidden in the midpoint.

## 4. What a slag gets

A blast-furnace slag holds more lime than the measured glasses can take: its
alumina goes to gehlenite, its magnesia to akermanite, its remaining silica to
pseudowollastonite, and some lime is left over, which no combination can avoid,
so the combination is unique and the span zero. For the slag of
[Snellings2022](@citet):

```@example glass
rows = literature_table("Snellings2022", "chemical_composition"; material = "GGBFS")
slag = Dict(String(o) => ustrip(p) / 100 for (o, p) in zip(rows.oxide, rows.percent) if o != "Total")
s = glass_enthalpy(slag; ignore = ("P2O5",))
for k in (:formation_298, :from_oxides, :vitrification, :uncertainty)
    @printf("%-14s %9.1f J/g\n", k, ustrip(us"J/g", getfield(s, k)))
end
@printf("lime left over: %.2f %% of the slag\n", 100 * s.unassigned["CaO"] * M("CaO"))
```

The enthalpy of vitrification, about 160 J/g, is the part of the heat a slag
can release that a crystalline assemblage of the same composition would not;
its uncertainty, from the printed uncertainties of the three glasses, is a few
joules per gram. The phosphorus of the analysis has no crystal in the
reference databases and is left out by name (`ignore`), which the result
reports.

## 5. What it does not give

The construction gives an enthalpy, which is what a heat needs. It does not
give an entropy or a Gibbs energy: the residual entropy of these glasses is not
measured, so the glass cannot be a phase whose stability an equilibrium
decides, and it enters a calculation as the residue of a constituent that
reacts at a prescribed extent. The fictive temperature of a granulated slag,
quenched in water, is not known, and section 3 gives its order of magnitude.
And the glass of a siliceous fly ash, rich in silica and alumina, with
potassium and iron, lies outside the measured glasses: its alumina has no
measured glass to go to without lime, and it is left over as corundum, with the
assumption that carries.
