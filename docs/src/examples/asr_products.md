# [The products of the alkali-silica reaction, synthesized at 80 °C](@id ex-asr-products)

!!! info "Before this page"
    [The products of the alkali-silica reaction](@ref sec-asr-extension), for
    the database extension this page uses, and [the certifying solver](@ref
    sec-theory-certificate), for the ionic strength of section 4.

The alkali-silica reaction dissolves reactive silica of the aggregate in the
alkaline pore solution of a concrete, and what precipitates in its place swells.
[ShiLothenbach2019](@citet) synthesized the crystalline products at 80 °C: 4 g of
silica fume, lime at a Ca/Si from 0 to 0.5 and KOH or NaOH at an alkali/Si of
0.5, in 30 to 100 g of water. After 90 days they identified the solids by X-ray
diffraction (XRD) and analyzed the solutions: amorphous silica without lime,
K- or Na-shlykovite at a low Ca/Si, and above it, in KOH, a potassium calcium
silicate hydrate they name ASR-P1. From the solutions they derived the
solubility products of the three solids at 80 °C. This page computes the 17
syntheses at equilibrium with Cemdata18 [Lothenbach2019](@cite), its C-S-H (CSHQ
with the alkali
end-member), and the three products.

```@example asr
using ChemistryLab, DynamicQuantities, Logging, OptimaSolver, Printf, Plots
default(framestyle = :box, grid = false)
include(joinpath(pkgdir(ChemistryLab), "scripts", "asr_products.jl"))
nothing # hide
```

## 1. Two sets of constants for the products

The constants of the products come in two sets, and they disagree.
[Jin2023](@citet) gave the two shlykovites standard properties at 25 °C, which
`cemdata18-asr.json` carries, by refining the measurements of
[ShiLothenbach2019](@citet) with the pH recalculated at 80 °C, where the
original authors had corrected the pH measured at 23 °C by subtracting 1.47.
ASR-P1 has only the solubility product of [ShiLothenbach2019](@citet), at 80 °C.
The set `:jin` takes the shlykovites of the database and ASR-P1 of
[ShiLothenbach2019](@citet); the set `:shi` takes all three from
[ShiLothenbach2019](@citet), consistent with one another and valid at 80 °C
only. Their solubility products at 80 °C, computed through the aqueous species
of Cemdata18:

```@example asr
g(sp) = ustrip(us"J/mol", sp[:ΔₐG⁰](T = SL19_T * u"K", P = 1.0e5u"Pa"; unit = true))
# The dissolution of each product as Table 5 writes it (sl19_reaction).
logK(sp) = -(sum(ν * g(SL19_DB[s]) for (s, ν) in sl19_reaction(symbol(sp))) - g(sp)) / (R_GAS * SL19_T * log(10))
println("product         set :jin   set :shi")
jin, shi = Dict(symbol(p) => p for p in sl19_products(:jin)), Dict(symbol(p) => p for p in sl19_products(:shi))
for name in ("K-shlykovite", "Na-shlykovite", "ASR-P1")
    @printf("%-14s %9.2f %10.2f\n", name, logK(jin[name]), logK(shi[name]))
end
```

The shlykovites of the first set are two to four orders of magnitude less
soluble. Within the second, the three products share one treatment of the pH;
the first mixes two.

## 2. The 17 syntheses

```@example asr
# Each synthesis at equilibrium at 80 °C, in each set. The warnings of the solver,
# if any, are kept out of the page; the certificates are printed.
runs = with_logger(NullLogger()) do # hide
runs = Dict(set => sl19_syntheses(set) for set in (:jin, :shi))
end # hide
println("sample   series    XRD                        set :jin                      set :shi")
for (a, b) in zip(runs[:jin], runs[:shi])
    @printf("%-8s %-8s  %-26s %-29s %s\n", a.sample, a.series, isempty(a.xrd) ? "(no sample)" : a.xrd, a.solids, b.solids)
end
for set in (:jin, :shi)
    println("set $set: ", count(r -> r.certified, runs[set]), " of ", length(runs[set]), " certified")
end
```

Every synthesis certifies in both sets. Without lime both give amorphous silica,
as observed. In the sodium series both give Na-shlykovite from a Ca/Si of 0.1,
as observed, with a C-S-H beside it earlier than the XRD finds one: from 0.3 in
the set `:jin` and 0.2 in the set `:shi`, where the XRD shows it at 0.4. At 0.5,
where the XRD finds the C-S-H alone, the set `:shi` gives it alone and the set
`:jin` keeps the shlykovite. The potassium series separates them. The XRD shows K-shlykovite at a Ca/Si of 0.1 and 0.2 and ASR-P1
from 0.3 up. The set `:jin` gives K-shlykovite at every Ca/Si, the set `:shi`
ASR-P1 at every Ca/Si below 0.5: each reproduces one half of the series. The
two solids have nearly the same composition, a Ca/Si of 0.25 and 0.29 and a
K/Si of 0.25 and 0.13, so that which one forms is decided by a few kilojoules,
less than the gap between the two treatments of the pH. Neither set is
consistent with both halves of the measurements, and the first is not
consistent with itself, its two potassium products having been derived from the
same solutions with two pH treatments.

## 3. The solutions

```@example asr
fig = plot([1.0, 2000.0], [1.0, 2000.0]; xscale = :log10, yscale = :log10, color = :gray, linestyle = :dash,
           label = "", xlabel = "silicon measured (mmol/L)", ylabel = "silicon computed (mmol/L)",
           title = "syntheses at 80 °C, 90 days", legend = :topleft, size = (680, 460),
           left_margin = 6Plots.mm, bottom_margin = 6Plots.mm)
for (set, color, shape) in ((:jin, :steelblue, :circle), (:shi, :darkorange, :diamond))
    rows = [r for r in runs[set] if !isnan(r.Si_meas)]
    scatter!(fig, [r.Si_meas for r in rows], [r.Si for r in rows]; color, markershape = shape, label = "set :$set")
end
savefig(fig, "asr-silicon.svg"); nothing # hide
```

![The silicon computed in the solutions of the syntheses against the silicon measured, in the two sets of constants, with the line of equality.](asr-silicon.svg)

Where the solid is amorphous silica, or a shlykovite at a Ca/Si up to 0.2, the
computed silicon is within a factor of about two of the measured one, between
400 and 1400 mmol/L. Beyond, the two sets part. In KOH from a Ca/Si of 0.3, where
the measured silicon falls to 10 to 400 mmol/L, both keep several hundred: the
products that hold the silicon there are not the ones the model forms, or not
with the energies it gives them. In NaOH from 0.3, the set `:jin` goes the other
way, its Na-shlykovite, four orders of magnitude less soluble than that of the
set `:shi`, taking the silicon down to 7 to 19 mmol/L where 145 to 442 were
measured; the set `:shi` keeps it within a factor of about two.
[ShiLothenbach2019; Table 6, footnote c](@citet) warn that at such
concentrations polynuclear silicate species dominate the solution, of which
Cemdata18 carries one, the tetramer; the solutions are the less telling
comparison here, and they still side with the set whose constants were derived
from them, consistently.

## 4. An equilibrium on the middle root of the ionic strength

The solution of SKC0, silica fume and KOH without lime, holds most of its
dissolved silicon as the tetramer ``\mathrm{Si_4O_{10}^{4-}}``. At the potentials
of its equilibrium the equation of the ionic strength the solver inverts at
every step, ``F(\ln I) = \ln\big(\tfrac12 \sum_i z_i^2 m_i(I)\big) - \ln I = 0``,
crosses zero three times, and the equilibrium is the middle root
([the certifying solver](@ref sec-theory-certificate) explains why, and which
root the solver takes). The function, drawn at the equilibrium of the set
`:jin`:

```@example asr
skc0 = runs[:jin][1]
F = sl19_ionic_strength_equation(skc0.state, skc0.alkali)
@printf("%s: certified %s, ionic strength %.3f mol/kg, F there %.1e\n", skc0.sample, skc0.certified, skc0.I, F(log(skc0.I)))
s = range(log(0.05), log(8.0); length = 400)
fig = plot(exp.(s), F.(s); xscale = :log10, color = :steelblue, linewidth = 2, label = "F",
           xlabel = "ionic strength I (mol/kg)", ylabel = "F", title = "SKC0: SiO₂ and KOH at 80 °C",
           legend = :topright, size = (680, 420), left_margin = 6Plots.mm, bottom_margin = 6Plots.mm)
hline!(fig, [0.0]; color = :gray, label = "")
vline!(fig, [skc0.I]; color = :darkorange, linestyle = :dash, label = "the equilibrium")
savefig(fig, "asr-ionic-roots.svg"); nothing # hide
```

![The function whose zero is the ionic strength, at the potentials of the equilibrium of SKC0, crossing zero three times, the equilibrium on the middle root.](asr-ionic-roots.svg)

Its ionic strength lies above the range the extended Debye-Hückel model of
Cemdata18 states for itself, which the certificate reports: the answer is that
of the model past its range, and the model, not the solver, is what limits it.

## Where to go next

[The products of the alkali-silica reaction](@ref sec-asr-extension) describes
the database extension and its checks, and `data/literature/ShiLothenbach2019.json`
and `data/literature/Jin2023.json` record the published values with where each
was read.
