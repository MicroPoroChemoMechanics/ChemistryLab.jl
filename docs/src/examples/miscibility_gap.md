# [The AFm of a CEM I 52.5 N and its miscibility gap: three ways to declare a binary that can unmix](@id ex-miscibility-gap)

!!! info "Before this page"
    [Solid solutions](@ref sec-theory-solid-solutions) §6 and [Solid solution
    models, in numbers](@ref sec-app-solid-solutions).

The AFm and AFt phases of a Portland cement are binaries — sulfate against
hydroxide, sulfate against carbonate — and the Redlich-Kister parameters
published for them in Cemdata18 [Lothenbach2019](@cite) are strong enough that the mixing energy is **concave over an
interval**. Where an energy is concave the Gibbs minimum is not one composition
but two: the phase unmixes, and the equilibrium is a **miscibility gap**.

This page runs the same cement three ways on the AFm sulfate/hydroxide binary,
because the three are genuinely different claims:

1. **ideal mixing** — an answer, certified, to a question that was changed;
2. **the published parameters with one composition** — the question kept; the
   answer is a minimum only if the composition falls outside the gap, and the
   certificate tests it;
3. **the published parameters with two compositions** — room for the pair,
   which only a composition inside the gap needs.

The theory is in [the solid-solution chapter](@ref sec-theory-solid-solutions);
this page is the measurement.

```@example gap
using ChemistryLab
using DynamicQuantities
using Logging
using OptimaSolver
using OrderedCollections
using Printf
using Plots
default(framestyle = :box, grid = false)

substances = build_species(datapath("cemdata18-thermofun.json"); verbose = false)
byname = Dict(symbol(s) => s for s in substances)
molar_mass(n) = ustrip(us"g/mol", byname[n][:M])
nothing # hide
```

## 1. Where the gap is, before any cement is involved

`spinodal_interval` scans the second derivative of the molar mixing energy and
returns the interval over which it is negative, or `nothing` when the model is
convex. It needs no solve and no system:

```@example gap
RT = R_GAS * 298.15

# Cemdata18's dimensionless Guggenheim parameters, from
# data/literature/Lothenbach2019.json: a₀ = α₀RT and a₁ = α₁RT.
function published_rk(pair)
    p = literature_row("Lothenbach2019", "guggenheim_parameters", pair)
    return RedlichKisterModel(a0 = p.alpha0 * RT, a1 = p.alpha1 * RT)
end

models = OrderedDict(
    "ideal" => IdealSolidSolutionModel(),
    "AFm SO4/OH (published)" => published_rk("AFm SO4/OH"),
    "AFt SO4/CO3 (published)" => published_rk("AFt SO4/CO3"),
    "regular, W = 1.9 RT" => RegularSolutionModel([0.0 1.9RT; 1.9RT 0.0]),
    "regular, W = 2.1 RT" => RegularSolutionModel([0.0 2.1RT; 2.1RT 0.0]),
)

for (label, m) in models
    gap = spinodal_interval(m, 2)
    if gap === nothing
        @printf("  %-26s convex everywhere\n", label)
    else
        @printf("  %-26s concave for x in [%.3f, %.3f]\n", label, gap[1], gap[2])
    end
end
```

The last two are the classical check: a symmetric regular solution unmixes above
``W = 2RT`` exactly, because ``d^2 g/dx^2`` at ``x = 1/2`` is ``4 - 2W/RT``. The two published
cement sets are well past it.

The curve itself makes the shape of the problem visible. A convex energy has one
minimum; a concave stretch means the straight line between two compositions lies
*below* the curve, and the system takes the line — which is to say, it takes the
two compositions and not any single one between them.

```@example gap
xs = range(0.001, 0.999; length = 400)
g(m, x) = begin
    # molar Gibbs energy of mixing, in units of RT
    ideal = x * log(x) + (1 - x) * log(1 - x)
    excess = if m isa RedlichKisterModel
        x * (1 - x) * (m.a0 + m.a1 * (2x - 1) + m.a2 * (2x - 1)^2) / RT
    else
        0.0
    end
    ideal + excess
end

afm = models["AFm SO4/OH (published)"]
gap = spinodal_interval(afm, 2)

fig = plot(xs, [g(afm, x) for x in xs]; label = "published Redlich-Kister",
           color = :firebrick, linewidth = 2,
           xlabel = "x (sulfate end-member)", ylabel = "g / RT",
           title = "The mixing energy of the AFm sulfate/hydroxide binary",
           size = (760, 420), bottom_margin = 8Plots.mm, left_margin = 8Plots.mm)
plot!(fig, xs, [g(IdealSolidSolutionModel(), x) for x in xs];
      label = "ideal mixing", color = :steelblue, linewidth = 2, linestyle = :dash)
vspan!(fig, [gap[1], gap[2]]; color = :firebrick, alpha = 0.12, label = "spinodal")
savefig(fig, "gap-energy.svg"); nothing # hide
```

![](gap-energy.svg)

## 2. The cement, declared three ways

A CEM I paste, with the AFm binary as the only phase whose model changes between
the three runs:

```@example gap
# The Bogue composition of the CEM I 52.5 N of [Lavergne2018](@cite), Table 9,
# and its gypsum.
bogue = literature_table("Lavergne2018", "cement_bogue")
CLINKER = OrderedDict(zip(bogue.phase, bogue.percent ./ 100))
GYPSUM = literature_value("Lavergne2018", "gypsum_percent") / 100
WB = 0.50
BINDER_G = 100.0

pure = split(
    "C3S C2S C3A C4AF Gp Anh Cal Portlandite ettringite monocarbonate " *
    "hemicarbonate C3AH6 C3FH6 straetlingite FeOOHmic AlOHmic Amor-Sl Brc"
)
CSHQ = ["CSHQ-JenD", "CSHQ-JenH", "CSHQ-TobD", "CSHQ-TobH", "KSiOH", "NaSiOH"]
AFM = ["C4AH13", "monosulphate12"]
aqueous = ["SO4-2", "CO2@"]

# The activity model Cemdata18 prescribes (its Eq. C.1): extended Debye-Hückel,
# with the common ion size and B-dot the paper gives for KOH solutions (it also
# gives them for NaOH). This paste carries no alkali, so neither set describes
# its calcium hydroxide and sulfate solution more closely; the KOH set is used,
# as on the other cement pages.
model = cemdata18_activity_model(:KOH)

"""
One system, differing only in how the AFm binary is declared.

`autostart` is a keyword because two of the three cases below are *meant* not to
certify. For those the multi-start cascade -- every backend, then the ideal
pre-solve, then the homotopy continuation -- cannot help: there is no admissible
single-composition minimum for it to find, so it spends minutes arriving at the
answer the first route already gave. The certificate reported is the same one
either way; declining the cascade declines only the search for a better start.
"""
function run_case(afm_phase; autostart = true)
    sp = speciation(substances, vcat(pure, CSHQ, AFM, aqueous);
                    aggregate_state = [AS_AQUEOUS])
    ss = [SolidSolutionPhase("CSHQ", [byname[m] for m in CSHQ]), afm_phase]
    cs = ChemicalSystem(sp, CEMDATA_PRIMARIES; solid_solutions = ss)

    st = ChemicalState(cs)
    for (phase, frac) in CLINKER
        set_quantity!(st, phase,
            BINDER_G * (1 - GYPSUM) * frac / molar_mass(phase) * u"mol")
    end
    set_quantity!(st, "Gp", BINDER_G * GYPSUM / molar_mass("Gp") * u"mol")
    set_quantity!(st, "H2O@", BINDER_G * WB / molar_mass("H2O@") * u"mol")
    b = Float64.(cs.SM.A) * ustrip.(us"mol", st.n)

    # Every case prints its certificate; the warning of a refusal would only repeat it.
    eq, cert = with_logger(NullLogger()) do
        equilibrate_certified(st; model = model, b = b, autostart = autostart)
    end
    return cs, eq, cert
end

afm_em = [byname[m] for m in AFM]
published = published_rk("AFm SO4/OH")
nothing # hide
```

### Case 1 — ideal mixing

```@example gap
cs1, eq1, c1 = run_case(SolidSolutionPhase("AFm_SO4_OH", afm_em))
@printf("optimal=%-5s  worst SI=%+.2e  balance=%.1e  pH=%.4f\n",
        c1.optimal, c1.worst_supersaturation, c1.balance, pH(eq1, model))
```

### Case 2 — the published parameters, one composition

`SolidSolutionPhase` refuses this declaration outright, and names the interval:

```@example gap
err = try
    SolidSolutionPhase("AFm_SO4_OH", afm_em; model = published)
    nothing
catch e
    e
end
println(err.msg)
```

Waiving the refusal is what a caller does who believes the answer stays outside
the gap. The certificate decides it, testing the present phase against splitting
as well as for stationarity:

```@example gap
cs2, eq2, c2 = run_case(SolidSolutionPhase("AFm_SO4_OH", afm_em;
                                           model = published, check_convexity = false);
                        autostart = false)
@printf("optimal=%-5s  worst SI=%+.2e  balance=%.1e  pH=%.4f\n",
        c2.optimal, c2.worst_supersaturation, c2.balance, pH(eq2, model))
if hasproperty(c2, :worst_violation_split)
    @printf("phases reported as wanting to split: %d   worst split measure %+.2e\n",
            length(c2.split_phases), c2.worst_violation_split)
end
```

Here the belief is right: the AFm of this paste has a C4AH13 fraction of 0.28,
outside the pair section 3 computes, and the certificate finds no phase that
wants to split. The answer is one composition, and a minimum. Until a phase with
an excess term was inverted by Newton's method, this solve did not converge, and
the tangent-plane test, applied to the point it stopped at, reported the phase as
wanting to split: a failure of the search read as one of the chemistry.

### Case 3 — the published parameters, two compositions

```@example gap
cs3, eq3, c3 = run_case(SolidSolutionPhase("AFm_SO4_OH", afm_em;
                                           model = published, instances = 2);
                        autostart = false)
@printf("optimal=%-5s  worst SI=%+.2e  balance=%.1e  pH=%.4f\n",
        c3.optimal, c3.worst_supersaturation, c3.balance, pH(eq3, model))

n3 = ustrip.(us"mol", eq3.n)
for (grp, ph) in zip(cs3.ss_groups, cs3.solid_solutions)
    startswith(name(ph), "AFm") || continue
    tot = sum(n3[i] for i in grp)
    x = tot > 1.0e-12 ? n3[grp[2]] / tot : NaN
    @printf("  %-14s total %8.5f mol   x(sulfate) = %s\n",
            name(ph), tot, isnan(x) ? "—" : @sprintf("%.4f", x))
end
```

With a second instance offered, the paste splits nothing: both instances hold
the composition of case 2, at the same pH, and the answer is certified. How the
AFm is shared between them is not a result. Two instances of a phase outside its
gap are degenerate: any division of one composition between them has the same
Gibbs energy, so only their sum is determined, and the division printed is
wherever the search stopped. That is why `instances = :auto` exists: it gives a
phase its second instance only when the certificate of the one-composition answer
asks for it.

## 3. The pair itself is computable, from the model alone

The two compositions are not unknown. They follow from the mixing model alone,
by the construction Glynn and Reardon set out [GlynnReardon1990](@cite) and that
PHREEQC uses for a binary solid solution: the pair at which a single straight
line is tangent to ``g`` twice, equivalently at which both end-members have equal
chemical potentials in the two phases.

```@example gap
ct = common_tangent(published, 2)
sp = spinodal_interval(published, 2)
@printf("common tangent : x = %.4f and %.4f\n", ct[1], ct[2])
@printf("spinodal       : x = %.4f to %.4f  (contained in it, as it must be)\n",
        sp[1], sp[2])
```

[`common_tangent`](@ref) solves two equations in two unknowns by Newton and costs
microseconds — it involves no chemical system at all. Its answer is checkable
without trusting the code: for a **symmetric** model ``g(1-x) = g(x)``, so
``g'(1-x) = -g'(x)``, and the condition collapses to ``g'(x) = 0``:

```@example gap
A = 3.0                                  # W/RT, well past the threshold of 2
sym = RegularSolutionModel([0.0 A*RT; A*RT 0.0])
a, b = common_tangent(sym, 2)
@printf("symmetric model: x = %.6f and %.6f, and b = 1 - a to %.1e\n",
        a, b, abs(b - (1 - a)))
@printf("residual of  ln(x/(1-x)) + A(1-2x) = 0  :  %+.2e and %+.2e\n",
        log(a / (1 - a)) + A * (1 - 2a), log(b / (1 - b)) + A * (1 - 2b))
```

So the endpoints of the gap are available to a user who needs them. What does
**not** work is handing them to the minimization as a starting point and
expecting it to keep the two lobes apart:

### And the split follows, by mass balance

Inside a gap **only the proportions move**: the two compositions are the same
whatever the overall composition is. So once the pair is known, how a given
overall composition ``\bar{x}`` separates is the lever rule and nothing more —
which [`miscibility_split`](@ref) returns, together with the Gibbs energy the
separation releases:

```@example gap
# The composition of the AFm of this paste, as the single-composition answer of
# case 2 returns it: x is the fraction of C4AH13, the first end-member.
n2 = ustrip.(us"mol", eq2.n)
grp2 = only(g for (g, ph) in zip(cs2.ss_groups, cs2.solid_solutions) if startswith(name(ph), "AFm"))
x̄_paste = n2[grp2[1]] / sum(n2[i] for i in grp2)

for x̄ in (0.20, 0.40, x̄_paste, 0.80, 0.96)
    r = miscibility_split(published, x̄, 2)
    if r.f_beta == 0
        @printf("x̄ = %.4f : homogeneous (outside the pair)\n", x̄)
    else
        # Two calls rather than one: `@printf` takes a LITERAL format string,
        # and `"a" * "b"` is an expression -- `ArgumentError: First argument to
        # @printf after io must be a format string`, raised while the block is
        # lowered, so it kills the build rather than the line.
        @printf("x̄ = %.4f : %.1f %% at x=%.4f and %.1f %% at x=%.4f",
                x̄, 100r.f_alpha, r.x_alpha, 100r.f_beta, r.x_beta)
        @printf("   —  releases %6.1f J/mol   (mass balance %+.0e)\n",
                r.Δg, r.f_alpha * r.x_alpha + r.f_beta * r.x_beta - x̄)
    end
end
```

``\Delta g`` is the distance from the curve down to the common tangent: it says, in
joules per mole of binary, **how much a single-composition answer overstates the
Gibbs energy**. The third line is the AFm of this paste, at the composition case
2 certifies: outside the pair, and homogeneous, which is what the certificate
found. The two lines inside the pair show what a split looks like and what it is
worth.

### Inside the gap

Everything above costs microseconds and needs no solver. What it does not give is
``\bar{x}`` itself, the overall composition the paste puts in the phase, which
the aqueous equilibrium decides. This paste's lies outside the gap. Where one
falls inside, `instances = :auto` gives the phase its second instance and the
split passes look for the pair, with the aqueous solution iterated along; the
theory chapter runs that on a binary its element budget holds in the gap
([Solid solutions](@ref sec-theory-solid-solutions), section 6). PHREEQC, for its
part, treats a binary solid solution with a dedicated construction rather than
the global minimization.

!!! info "Where this leaves the three questions"
    | | |
    |:--|:--|
    | **detect** a gap | yes — the certificate refuses a single composition inside one, and names the phase |
    | **locate** it | yes — `spinodal_interval` and `common_tangent`, from the model alone |
    | **represent** it | yes — `instances = 2`, or `instances = :auto`, which adds the second when the certificate asks for it; the species exist, the groups stay disjoint, conservation is untouched |
    | **split** a given overall composition | yes — `miscibility_split`, exact, with the energy it releases |
    | iterate that back through the aqueous equilibrium | yes, with `instances = :auto` where the budget holds the composition inside the gap ([Solid solutions](@ref sec-theory-solid-solutions), section 6); this paste does not need it |

    The first four are what a user needs to answer "is this phase homogeneous,
    and if not, into what?". The last is the coupling with the rest of the
    paste.

    In GEM-Selektor and Reaktoro the second declaration is the user's, and
    PHREEQC draws the same line: it treats a binary solid solution with a
    dedicated construction rather than handing it to the global minimization.

    Stated plainly because a reader deciding whether to trust a number here
    deserves to know which of these five they are relying on.

!!! note "Why the second instance is refused for a convex model"
    Two instances of a convex phase are **degenerate**: every way of splitting
    the amount between them has the same energy, so the minimum is a flat
    manifold and the optimizer is asked to choose a point on it for no reason.
    Inside a spinodal the common-tangent pair is unique, and the degeneracy does
    not arise. `SolidSolutionPhase` therefore admits `instances > 1` exactly
    where `spinodal_interval` reports a gap.

## 4. The three answers side by side

```@example gap
labels = ["ideal\nmixing", "published,\none composition", "published,\ntwo compositions"]
certs = [c1, c2, c3]
eqs = [eq1, eq2, eq3]
vols = [ustrip(uconvert(us"cm^3", volume(e).total)) for e in eqs]
phs = [pH(e, model) for e in eqs]
ok = [c.optimal for c in certs]

p1 = bar(labels, vols; legend = false, ylabel = "total volume (cm³)",
         color = [o ? :seagreen : :firebrick for o in ok], title = "Volume")
p2 = bar(labels, phs; legend = false, ylabel = "pH", title = "Pore solution pH",
         color = [o ? :seagreen : :firebrick for o in ok],
         ylims = (12.0, 13.5))
fig3 = plot(p1, p2; layout = (1, 2), size = (900, 420),
            bottom_margin = 16Plots.mm, left_margin = 8Plots.mm)
savefig(fig3, "gap-summary.svg"); nothing # hide
```

![](gap-summary.svg)

Green is a proof and red is its absence — not a claim that a red bar is far from
the truth, which nothing here establishes. The three give the same pH to four
figures. The middle one is proved with the published parameters kept, because
this paste's AFm lies outside the gap; the third, which offers a second instance
that nothing needs, is the one the search cannot close.

!!! danger "What this does and does not settle"
    `optimal = true` in case 2 is a proof of a **KKT point at which no present
    phase wants to split and no absent phase is supersaturated**. For a convex
    problem those conditions are sufficient for a global minimum, and that is why
    the certificate means what it means on every other page. With a concave
    mixing model the problem is no longer convex, so global optimality is not
    implied by stationarity alone — what the certificate adds here is the
    tangent-plane test, which is exactly the condition that separates a
    one-composition minimum from a spurious stationary point; its scope does not
    claim a global minimum.

## 5. How other codes represent the same thing

GEM-Selektor computes the same criterion — its phase stability index
``Λ_k = \log_{10} Ω_k`` is, term for term, the quantity this package's `Ω` and
`phase_split_measure` compute, derived independently from the same KKT conditions
[Kulik2013](@cite). Where the two differ is not detection but **declaration**:
CEMDATA18 ships the AFm and AFt binaries under two names each, so a GEMS user
represents a gap by declaring the binary twice, in the database. `instances = 2`
is the same representation asked for by a keyword instead.

A GEMS user gets the right answer because the database ships the binary twice
and the solver is handed two declarations to populate. Here `instances = :auto`
lets the certificate decide: a phase is given its second instance when the
stability test finds it wanting to split, and the split passes then look for the
pair ([Solid solutions](@ref sec-theory-solid-solutions), section 6).

## Where to go next

A cement declared with every solid solution of its database, and the phase list
treated as a modeling decision, is
[A CEM I at equilibrium, with every solid solution declared](@ref). The aqueous
counterpart of this chapter continues with
[The Pitzer model, and what it can be used with](@ref sec-app-pitzer).
