# ── An oxide analysis as an element budget ───────────────────────────────────
#
# A clinker phase has a formula, so it can be named as a species and its amount
# set. A **glass** has none: ground granulated blast-furnace slag, a fly ash and
# a natural pozzolana are reported by their oxide analysis and by nothing else,
# because there is no phase to report. Bogue does not help — it inverts a
# decomposition over phases that exist, and here they do not.
#
# What a glass can still do is contribute its elements, and that is all an
# equilibrium calculation needs: the constraint is `A n = b`, and `b` is a
# vector of component totals. So a glass enters as its oxide analysis, converted
# to moles and decomposed over the system's own components.

"""
    primary_decomposition(species, primaries) -> Vector{Float64}

The coefficients writing `species` as a linear combination of `primaries`.

Solved by column-pivoted QR rather than `\\`, because the matrix is singular
whenever the primaries are not all exercised, and **refused** above a residual
of `1e-8`: a species outside the span of the primaries has no decomposition, and
returning a least-squares approximation of one would put elements into the
budget that the species does not contain.
"""
function primary_decomposition(species::AbstractSpecies, primaries)
    els = Symbol[]
    for p in Iterators.flatten((primaries, (species,))), k in keys(atoms(p))
        k in els || push!(els, k)
    end
    E = zeros(length(els), length(primaries))
    for (j, p) in enumerate(primaries), (k, v) in atoms(p)
        E[findfirst(==(k), els), j] = Float64(v)
    end
    t = zeros(length(els))
    for (k, v) in atoms(species)
        t[findfirst(==(k), els)] = Float64(v)
    end
    x = qr(E, ColumnNorm()) \ t
    r = norm(E * x - t)
    r > 1.0e-8 && throw(
        ArgumentError(
            "`$(symbol(species))` is not in the span of the primaries " *
                "(residual $r). Its composition contains something they cannot " *
                "express, so no decomposition exists — adding one anyway would " *
                "put elements into the budget that this species does not carry.",
        ),
    )
    return x
end

"""
    oxide_budget(oxides, primaries; mass = 100.0u"g") -> Vector{Float64}

The component totals `b` contributed by a material reported as an **oxide
analysis**, for use as the right-hand side of `A n = b`.

`oxides` maps an oxide formula to its mass fraction — the shape of a cement or
slag datasheet, `"CaO" => 0.41` and so on. `mass` is how much of the material
the budget is for. The analysis is **not** renormalized: if the fractions do not
sum to one, what is missing is loss on ignition and minor oxides the sheet does
not report, and silently scaling them up would invent material.

Each oxide's molar mass is computed from its formula by [`Species`](@ref), never
written down, and each is decomposed over `primaries` by
[`primary_decomposition`](@ref), which refuses an oxide the primaries cannot
express.

```julia
# A ground granulated blast-furnace slag, as a datasheet reports it
slag = Dict("CaO" => 0.41, "SiO2" => 0.36, "Al2O3" => 0.11, "MgO" => 0.08)
b = oxide_budget(slag, cs.SM.primaries; mass = 60.0u"g")
```

!!! warning "This says what a glass contains, not what it does"
    An oxide budget is a statement of composition and carries no information
    about **reactivity**. A slag and a quartz sand of the same analysis give the
    same `b`, and an equilibrium calculation will dissolve both completely. How
    much of the glass has actually reacted is a kinetic quantity, supplied from
    outside — by a degree of reaction, or by a rate law such as
    [`waller`](@ref) — and multiplying the budget by it is the caller's
    responsibility, not this function's.
"""
function oxide_budget(
        oxides::AbstractDict{<:AbstractString, <:Real}, primaries;
        mass = 100.0u"g",
    )
    m_g = ustrip(us"g", mass)
    prim = collect(primaries)
    b = zeros(length(prim))
    for (formula, frac) in oxides
        frac == 0 && continue
        frac < 0 && throw(
            ArgumentError("the mass fraction of `$formula` is negative ($frac)"),
        )
        ox = Species(formula)
        n_ox = m_g * float(frac) / ustrip(us"g/mol", ox[:M])
        b .+= n_ox .* primary_decomposition(ox, prim)
    end
    return b
end
