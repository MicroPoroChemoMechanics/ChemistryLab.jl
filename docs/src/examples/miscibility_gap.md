# [A miscibility gap, and the three answers a formulation can give](@id ex-miscibility-gap)

The AFm and AFt phases of a Portland cement are binaries — sulfate against
hydroxide, sulfate against carbonate — and the Redlich-Kister parameters
published for them are strong enough that the mixing energy is **concave over an
interval**. Where an energy is concave the Gibbs minimum is not one composition
but two: the phase unmixes, and the equilibrium is a **miscibility gap**.

This page runs the same cement three ways on the AFm sulfate/hydroxide binary,
because the three are genuinely different claims:

1. **ideal mixing** — an answer, certified, to a question that was changed;
2. **the published parameters with one composition** — the question kept, and
   the answer is not a minimum. The certificate now says so;
3. **the published parameters with two compositions** — the question kept and
   answered, at the common tangent.

The theory is in [the solid-solution chapter](@ref sec-theory-solid-solutions);
this page is the measurement.

```@example gap
using ChemistryLab
using DynamicQuantities
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
RT = 8.31446261815324 * 298.15

models = OrderedDict(
    "ideal" => IdealSolidSolutionModel(),
    "AFm SO4/OH (published)" => RedlichKisterModel(a0 = 0.188RT, a1 = 2.49RT),
    "AFt SO4/CO3 (published)" => RedlichKisterModel(a0 = 1.67RT, a1 = 0.946RT),
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
`W = 2RT` exactly, because `d²g/dx²` at `x = ½` is `4 − 2W/RT`. The two published
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
CLINKER = OrderedDict("C3S" => 0.65, "C2S" => 0.11, "C3A" => 0.11, "C4AF" => 0.08)
GYPSUM = 0.046
WB = 0.50
BINDER_G = 100.0

pure = split(
    "C3S C2S C3A C4AF Gp Anh Cal Portlandite ettringite monocarbonate " *
    "hemicarbonate C3AH6 C3FH6 straetlingite FeOOHmic AlOHmic Amor-Sl Brc"
)
CSHQ = ["CSHQ-JenD", "CSHQ-JenH", "CSHQ-TobD", "CSHQ-TobH", "KSiOH", "NaSiOH"]
AFM = ["C4AH13", "monosulphate12"]
aqueous = ["SO4-2", "CO2@"]

model = HKFActivityModel(å = 0.0, Ḃ = 0.097637, Kₙ = 0.0)

"""One system, differing only in how the AFm binary is declared."""
function run_case(afm_phase)
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

    eq, cert = equilibrate_certified(st; model = model, b = b)
    return cs, eq, cert
end

afm_em = [byname[m] for m in AFM]
published = RedlichKisterModel(a0 = 0.188RT, a1 = 2.49RT)
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
the gap. Here it does not, and the certificate is what says so:

```@example gap
cs2, eq2, c2 = run_case(SolidSolutionPhase("AFm_SO4_OH", afm_em;
                                           model = published, check_convexity = false))
@printf("optimal=%-5s  worst SI=%+.2e  balance=%.1e  pH=%.4f\n",
        c2.optimal, c2.worst_supersaturation, c2.balance, pH(eq2, model))
if hasproperty(c2, :worst_violation_split)
    @printf("phases reported as wanting to split: %d   worst split measure %+.2e\n",
            length(c2.split_phases), c2.worst_violation_split)
end
```

Before the tangent-plane test was applied to phases that are **present**, this
case certified: every member was stationary, nothing absent was supersaturated,
and the element balance closed. Stationarity is blind to the one failure that
matters for a non-ideal phase — that the minimum is two compositions rather than
the one reported — so the answer was a KKT point and not a minimum, and nothing
in the output said which.

### Case 3 — the published parameters, two compositions

```@example gap
cs3, eq3, c3 = run_case(SolidSolutionPhase("AFm_SO4_OH", afm_em;
                                           model = published, instances = 2))
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

The two instances are the **common-tangent pair**: one sits on each side of the
spinodal, and together they hold the amount that a single composition could not.
That is the whole content of `instances = 2` — the substance is present twice,
under derived symbols sharing one thermodynamic record, so the composition vector
has somewhere to put each lobe.

!!! note "Why the second instance is refused for a convex model"
    Two instances of a convex phase are **degenerate**: every way of splitting
    the amount between them has the same energy, so the minimum is a flat
    manifold and the optimizer is asked to choose a point on it for no reason.
    Inside a spinodal the common-tangent pair is unique, and the degeneracy does
    not arise. `SolidSolutionPhase` therefore admits `instances > 1` exactly
    where `spinodal_interval` reports a gap.

## 3. The three answers side by side

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
         ylims = (12.5, 13.5))
fig3 = plot(p1, p2; layout = (1, 2), size = (900, 420),
            bottom_margin = 16Plots.mm, left_margin = 8Plots.mm)
savefig(fig3, "gap-summary.svg"); nothing # hide
```

![](gap-summary.svg)

Green is a proof and red is its absence — not a claim that the red bar is far
from the truth, which nothing here establishes. The point of the middle case is
precisely that it *looks* like the others.

!!! danger "What this does and does not settle"
    `optimal = true` in case 3 is a proof of a **KKT point at which no present
    phase wants to split and no absent phase is supersaturated**. For a convex
    problem those conditions are sufficient for a global minimum, and that is why
    the certificate means what it means on every other page. With a concave
    mixing model the problem is no longer convex, so global optimality is not
    implied by stationarity alone — what the certificate adds here is the
    tangent-plane test, which is exactly the condition that separates a
    common-tangent pair from a spurious stationary point.

    Read case 3 as: the answer satisfies every first-order condition *and* the
    stability test that case 2 fails. That is a stronger statement than case 2's,
    and a weaker one than the proof the convex pages carry.

## 4. How other codes represent the same thing

GEM-Selektor computes the same criterion — its phase stability index
``Λ_k = \log_{10} Ω_k`` is, term for term, the quantity this package's `Ω` and
`phase_split_measure` compute, derived independently from the same KKT conditions
[Kulik2013](@cite). Where the two differ is not detection but **declaration**:
CEMDATA18 ships the AFm and AFt binaries under two names each, so a GEMS user
represents a gap by declaring the binary twice, in the database. `instances = 2`
is the same representation asked for by a keyword instead.

Neither code splits a phase by itself. That is the honest statement of where
things stand, and it is a point of agreement rather than of difference.
