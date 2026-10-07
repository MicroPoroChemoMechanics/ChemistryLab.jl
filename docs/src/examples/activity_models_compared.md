# [What the choice of activity model costs](@id sec-app-activity-models)

!!! info "Before this page"
    [Activity models](@ref sec-theory-activity).

[Activity models](@ref sec-theory-activity) sets out what the three built-in
models are and where each comes from. This page puts numbers on the difference:
the models disagree about the activity coefficients from a tenth molal, barely
about the water activity below a molal, and, in a single electrolyte, not at all
about whether the solvent and the solutes form one thermodynamic system; that
last question is decided in a mixed solution with a neutral species.

A fourth column, `Cemdata18`, is the B-dot model [Helgeson1969](@cite) again,
with the parameters Cemdata18 prescribes for the pore solution of a cement
([`cemdata18_activity_model`](@ref)): one ion size, 3.67 Å, for every ion, and
the same linear coefficient on the ions and on the neutral species.

Everything below is evaluated on an imposed NaCl composition. Nothing is solved,
so the whole page costs a few milliseconds — which is also the point: comparing
models does not require a converged equilibrium.

```@example am
using ChemistryLab
using DynamicQuantities
using Printf
using LinearAlgebra
using ForwardDiff

substances = build_species(datapath("slop98-inorganic-thermofun.json"); verbose = false)
dict = Dict(symbol(s) => s for s in substances)
cs = ChemicalSystem([dict[s] for s in split("H2O@ H+ OH- Na+ Cl-")],
                    ["H2O@", "H+", "Na+", "Cl-", "Zz"])
sym_w = symbol(cs.species[only(cs.idx_solvent)])

function nacl(m)                       # m mol NaCl per kg of water, imposed
    st = ChemicalState(cs)
    set_quantity!(st, "H2O@", 1.0u"kg")
    set_quantity!(st, "Na+", m * u"mol")
    set_quantity!(st, "Cl-", m * u"mol")
    set_quantity!(st, "H+", 1.0e-7u"mol")
    set_quantity!(st, "OH-", 1.0e-7u"mol")
    return st
end

# The fourth is the B-dot model again, with the parameters Cemdata18 prescribes for
# a KOH pore solution: one ion size for every ion, and the same b on the neutral
# species (Lothenbach et al. 2019, Eq. C.1).
models = ["dilute" => DiluteSolutionModel(),
          "Davies" => DaviesActivityModel(),
          "B-dot" => HKFActivityModel(),
          "Cemdata18" => cemdata18_activity_model(:KOH)]
nothing # hide
```

## 1. ``A`` and ``B`` are properties of water

The Debye-Hückel coefficients are not fitting constants: they follow from the
density and the dielectric constant of water, and
[`hkf_debye_huckel_params`](@ref) evaluates them from this package's own equation
of state. The defaults the models carry, ``A = 0.5114`` and ``B = 0.3288`` at
25 °C, are therefore *derived* — and they agree with the values the LLNL aqueous
model tabulates [ParkhurstAppelo2013; p. 118](@cite), printed beside them:

```@example am
llnl = literature_table("ParkhurstAppelo2013", "llnl_debye_huckel")
for θ in (25.0, 60.0, 100.0)
    T = 273.15 + θ
    p = hkf_debye_huckel_params(T, 1.0e5)
    w = water_thermo_props(T, 1.0e5)
    e = water_electro_props_jn(T, 1.0e5, w)
    i = findfirst(==(θ), llnl.temperature_C)
    @printf("T = %6.2f K   ρ = %.4f g/cm³   ε = %6.2f   A = %.4f (%.4f)   B = %.4f (%.4f)\n",
            T, w.D / 1000, e.epsilon, p.A, llnl.A[i], p.B, llnl.B[i])
end
```

[Helgeson1981](@citet), Table 1, computed from the water properties of the time,
gives slightly lower values, 0.5091 and 0.3283 at 25 °C.

Both rise with temperature, because water's dielectric constant falls faster
than ``T`` rises: hot water screens worse, so the same ionic strength costs more.

## 2. The screening length, in nanometers

``\kappa^{-1} = 1/(B\sqrt{I})`` is in ångström when ``B`` is in
``\text{Å}^{-1}\,(\text{kg/mol})^{1/2}``; the table gives it in nanometers:

```@example am
B25 = hkf_debye_huckel_params(298.15, 1.0e5).B
println("     I (mol/kg)    Debye length (nm)")
for I in (0.001, 0.01, 0.1, 0.3, 1.0, 3.0)
    @printf("   %10.3f    %14.3f\n", I, 1 / (B25 * sqrt(I)) / 10)
end
```

A cement pore solution sits near ``I \approx 0.1``–``0.5`` mol/kg, so its
screening length is **a few ångström** — two or three water molecules, and the
width of the gel pores where a hydrating paste keeps its last water. That
coincidence is discussed in
[the water budget](@ref sec-theory-water-budget); it is the reason to treat an
extended Debye-Hückel model as a correlation in bulk solution rather than a
theory of confined water.

## 3. The activity coefficients, and the water activity

```@example am
println("                    γ(Na⁺)                                a_w")
println("  m      dilute   Davies    B-dot  Cemdata18    dilute   Davies    B-dot  Cemdata18")
for m in (0.001, 0.01, 0.1, 0.5, 1.0, 3.0)
    st = nacl(m)
    γ = [activity_coefficients(st, mod)["Na+"] for (_, mod) in models]
    aw = [exp(log_activities(st, mod)[sym_w]) for (_, mod) in models]
    @printf("%6.3f  %7.4f  %7.4f  %7.4f  %7.4f      %7.5f  %7.5f  %7.5f  %7.5f\n", m, γ..., aw...)
end
```

Read the ``\gamma`` columns first. The ideal model is already several percent off
at a **millimolal**, which is worth knowing before treating ideality as a safe
default. The two corrections agree to about 1 % up to a tenth molal and part
company above it — 8 % apart at 0.5 mol/kg — and by 3 mol/kg Davies
[Davies1962](@cite) has returned ``\gamma > 1`` while the B-dot model is still
below 1: the ``bI`` term has taken over, which is the ceiling of the B-dot
construction arriving.

The `Cemdata18` column is the same formula with other constants, and it stays
within 2 % of the B-dot column up to a tenth molal. Above, its larger linear term
takes over sooner: 0.78 against 0.65 at 1 mol/kg, and above 1 at 3 mol/kg, past
the ionic strength of about 1 mol/kg up to which Cemdata18 states it.

Now the ``a_w`` columns: **they barely separate below a molal**, within 0.4 %
of one another at 1 mol/kg. At 3 mol/kg the ideal 0.898 and the B-dot 0.895
still sit together, while Davies and the Cemdata18 constants, whose linear terms
weigh more, take the water activity to 0.860 and 0.868. Each of these is the
solvent row the Gibbs-Duhem relation pairs with the model's own activity
coefficients, so a difference in ``a_w`` is a difference in ``\gamma`` seen from
the solvent.

Seen as curves rather than as a table, the separation is a matter of where each
model leaves the limiting law:

```@example am
using Plots

ms = exp10.(range(-3, log10(3.0); length = 120))
γ(mod, m) = activity_coefficients(nacl(m), mod)["Na+"]
aw(mod, m) = exp(log_activities(nacl(m), mod)[sym_w])
A25 = hkf_debye_huckel_params(298.15, 1.0e5).A

p1 = plot(; xscale = :log10, xlabel = "molality m (mol/kg)", ylabel = "γ(Na⁺)",
    title = "Four models, one electrolyte", legend = :bottomleft)
for ((name, mod), col) in zip(models, (:gray, :firebrick, :steelblue, :darkorange))
    plot!(p1, ms, [γ(mod, m) for m in ms]; label = name, linewidth = 2, color = col)
end
plot!(p1, ms, [exp(-A25 * sqrt(m) * log(10)) for m in ms];
    label = "Debye-Hückel limiting law", linestyle = :dot, color = :black, linewidth = 2)
plot(p1; size = (720, 430), left_margin = 8Plots.mm, bottom_margin = 8Plots.mm)
```

```@example am
p2 = plot(; xscale = :log10, xlabel = "molality m (mol/kg)", ylabel = "water activity a_w",
    title = "The water activity barely separates", legend = :bottomleft)
for ((name, mod), col) in zip(models, (:gray, :firebrick, :steelblue, :darkorange))
    plot!(p2, ms, [aw(mod, m) for m in ms]; label = name, linewidth = 2, color = col,
        linestyle = name == "Davies" ? :dash : :solid)
end
plot(p2; size = (720, 430), left_margin = 10Plots.mm, bottom_margin = 8Plots.mm)
```

The four curves lie together up to a molal, and the two with the larger linear
term leave the others above it.

## 4. Gibbs-Duhem: one system, or two halves

Equilibrium is set by chemical potentials, that is by derivatives of the
activities with respect to composition. So the test that matters is whether
``\sum_i n_i\,\mathrm{d}\mu_i = 0`` holds along a composition change, measured
here along three directions, with the derivative taken by automatic
differentiation so that nothing but the algebra is measured.

```@example am
cs3 = ChemicalSystem([dict[s] for s in split("H2O@ Na+ Cl-")], ["H2O@", "Na+", "Cl-"])
M_W = ustrip(us"kg/mol", cs3.species[only(cs3.idx_solvent)][:M])
n_w = 1.0 / M_W

function gd_residual(mod, m, dn)
    μ = build_potentials(cs3, mod)
    p = (ΔₐG⁰overRT = zeros(3), T = 298.15, P = 1.0e5, ϵ = 1.0e-30)
    n0 = [n_w, m, m]
    dμ = ForwardDiff.jacobian(n -> μ(n, p), n0) * dn
    return abs(sum(n0 .* dμ)) / max(norm(n0 .* abs.(dμ)), 1.0)
end

for (name, dn) in ("dissolution   dn = (0, +1, +1)" => [0.0, 1.0, 1.0],
                   "ion exchange  dn = (0, +1, -1)" => [0.0, 1.0, -1.0],
                   "water removal dn = (-1, 0, 0)" => [-1.0, 0.0, 0.0])
    println("\n── ", name)
    println("   m         dilute       Davies        B-dot    Cemdata18")
    for m in (0.1, 0.3, 1.0, 3.0)
        r = [gd_residual(mod, m, dn) for (_, mod) in models]
        @printf("%6.2f    %10.3e   %10.3e   %10.3e   %10.3e\n", m, r...)
    end
end
```

In sodium chloride every model is consistent to rounding, along every
direction: each builds its solvent row from its own solutes' terms, the ideal
model as ``-M_w \sum m``, Davies by a closed form, the B-dot model by integrating
its ``A``, ``B`` and ``\dot{B}``, and with two ions of one charge and one size
nothing is left to break the relation. What breaks it is what a single salt
cannot show: ions of different sizes or charges under a linear term, and neutral
species carrying a salting-out coefficient that the ions' coefficients do not
return. A mixed solution, sodium and calcium chlorides and sulfates with
dissolved carbon dioxide, shows it, through the symmetry of the Jacobian of the
log activities and the Gibbs-Duhem relation over all its columns:

```@example am
csm = ChemicalSystem([dict[s] for s in split("H2O@ Na+ Ca+2 Cl- SO4-2 CO2@")],
                     ["H2O@", "Na+", "Ca+2", "Cl-", "SO4-2", "CO2@"])
amt = Dict("H2O@" => n_w, "Na+" => 0.2, "Ca+2" => 0.02, "Cl-" => 0.2,
           "SO4-2" => 0.02, "CO2@" => 0.01)
nm = [amt[symbol(s)] for s in csm.species]
pm = ChemistryLab._build_params(ChemicalState(csm, nm .* u"mol"))
println("             symmetry   Gibbs-Duhem")
for (name, mod) in models
    J = ForwardDiff.jacobian(n -> activity_model(csm, mod)(n, pm), nm)
    sym = ChemistryLab._jacobian_asymmetry(J)[1]
    gd = ChemistryLab._gibbs_duhem_defect(J, nm)[1]
    @printf("%-10s  %10.3e  %10.3e\n", name, sym, gd)
end
```

The ideal model stays exact. The other three are not the gradient of a Gibbs
energy here: Davies through the salting-out term of the dissolved carbon dioxide
alone (with `bₙ = 0` it is exact), the B-dot model through its ion sizes, its
linear term and that same salting-out, and the Cemdata18 constants, which give
the neutral species the ions' linear coefficient, most of all. Their certified
equilibria are then compositions consistent with their own activities rather
than minima of an energy, which is what their certificate's `scope` says, and an
optimizer minimizing ``\mathbf{n}\cdot\boldsymbol{\mu}(\mathbf{n})`` would not find
them: see [One equilibrium, whatever the back end](@ref).

## What to take from this

  - below ``I \approx 0.01`` mol/kg the choice hardly matters, and
    [`DiluteSolutionModel`](@ref) is the best-conditioned objective;
  - between there and about a molal, use [`HKFActivityModel`](@ref), whose ion
    sizes distinguish the ions Davies treats alike, and which follows the
    measured coefficients further than Davies once their linear terms take
    over. For the pore
    solution of a cement, [`cemdata18_activity_model`](@ref) gives it the
    parameters Cemdata18 prescribes, which [Lothenbach2019](@citet) state valid
    to about 1 mol/kg, and a certified answer reports its ionic strength against
    that;
  - above a few molal, none of these is defensible, and no warning is issued
    because none of them knows. [`solvent_fraction`](@ref) is the guard that
    does.

## Where to go next

The expressions compared here are derived in
[Activity models](@ref sec-theory-activity), and the same exercise on the
mole-fraction side is [Solid solution models, in numbers](@ref sec-app-solid-solutions).
