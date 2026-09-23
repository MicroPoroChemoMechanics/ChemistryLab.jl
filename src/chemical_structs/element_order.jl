# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2025-2026 Jean-François Barthélémy and Anthony Soive (Cerema, UMR MCD)

using OrderedCollections

"""
    SITE_SYMBOLS :: NTuple{24, Symbol}

The pseudo-element symbols reserved for **surface site families**.

A site is a conserved quantity that is not a chemical element: a species
occupying one carries the family's symbol in its formula, exactly as a charged
species carries `:Zz`, and the conservation row then falls out of the ordinary
matrix assembly instead of being bolted on at the solve.

# Why a fixed list rather than a registry

A registry would be mutable global state shared between calculations, which is
the defect release 0.19.0 removed from the solver. A fixed list is a pure
predicate, [`is_site_symbol`](@ref), and it cannot drift between two systems
built in the same session.

# Why these twenty-five

The formula parser accepts one uppercase letter and at most one lowercase one,
so a site symbol has two characters; `Xs1` would parse as `Xs` followed by the
stoichiometric coefficient 1. That leaves `X` plus a lowercase letter, minus two:

  - **`:Xe` is xenon**, a real element;
  - **`:Xx` is the conventional placeholder for "not an element"**, and this
    package already uses it as one. Making it a site would mean a typo could
    quietly become a valid site family.

Twenty-four families is a ceiling, and it is a real one — say so rather than
work around it if a system ever needs a twenty-fifth.

By convention `:Xs`, `:Xw` and `:Xv` read as strong, weak and second-weak sites,
which is the naming the clay and oxide literature uses, but nothing enforces it.

See also: [`is_site_symbol`](@ref), [`ATOMIC_ORDER`](@ref).
"""
const SITE_SYMBOLS = (
    :Xa, :Xb, :Xc, :Xd, :Xf, :Xg, :Xh, :Xi, :Xj, :Xk, :Xl, :Xm,
    :Xn, :Xo, :Xp, :Xq, :Xr, :Xs, :Xt, :Xu, :Xv, :Xw, :Xy, :Xz,
)

"""
    is_site_symbol(s::Symbol) -> Bool

Whether `s` is one of the [`SITE_SYMBOLS`](@ref) reserved for surface sites.

Pure, and false for every real element — `:Xe` included.

# Examples

```jldoctest
julia> is_site_symbol(:Xs), is_site_symbol(:Xe), is_site_symbol(:Xx)
(true, false, false)
```
"""
is_site_symbol(s::Symbol) = s in SITE_SYMBOLS

"""
    ATOMIC_ORDER :: Vector{Symbol}

Canonical atomic ordering used across the package for serialization,
stringification and deterministic ordering of formula fields.

This vector lists element symbols in the preferred display/serialization
order and includes the special placeholder :Zz which represents a unit
positive charge in compositions.

# Examples

```julia
julia> print(ATOMIC_ORDER)
[:Ca, :Na, :K, :Mg, :Sr, :Ba, :Al, :Fe, :Ti, :Mn, :Cr, :Si, :C, :H, :N, :S, :O, :P, :B, :F, :Cl, :Br, :I, :Zz]
```
"""
const ATOMIC_ORDER = [
    :Ca,
    :Na,
    :K,
    :Mg,
    :Sr,
    :Ba,
    :Al,
    :Fe,
    :Ti,
    :Mn,
    :Cr,
    :Si,
    :C,
    :H,
    :N,
    :S,
    :O,
    :P,
    :B,
    :F,
    :Cl,
    :Br,
    :I,
    :U,
    # The site families, between the elements and the charge. Their position
    # matters twice: `Formula` sorts a composition by this list and has no
    # fallback for a symbol absent from it, and `speciation` relies on `:Zz`
    # staying last.
    SITE_SYMBOLS...,
    :Zz,
]

"""
    CEMENT_TO_MENDELEEV :: Vector{Pair{Symbol,OrderedDict{Symbol,Int}}}

Mapping from cement shorthand symbols to their corresponding oxide or
elemental compositions. Each Pair maps a cement shorthand Symbol (key)
to an OrderedDict (value) describing composition in terms of element
symbols and integer stoichiometric coefficients.

This mapping is used to translate cement shorthand notation into full
elemental compositions for formula construction and serialization.

# Examples

```julia
julia> haskey(Dict(CEMENT_TO_MENDELEEV), :C)
true

julia> for (k,v) in CEMENT_TO_MENDELEEV println(k, " ≡ ", unicode(Species(v))) end
C ≡ CaO
M ≡ MgO
S ≡ SiO₂
A ≡ Al₂O₃
F ≡ Fe₂O₃
K ≡ K₂O
N ≡ Na₂O
P ≡ O₅P₂
T ≡ TiO₂
C̄ ≡ CO₂
S̄ ≡ SO₃
N̄ ≡ NO₃
H ≡ H₂O
```
"""
const CEMENT_TO_MENDELEEV = [
    :C => OrderedDict(:Ca => 1, :O => 1),
    :M => OrderedDict(:Mg => 1, :O => 1),
    :S => OrderedDict(:Si => 1, :O => 2),
    :A => OrderedDict(:Al => 2, :O => 3),
    :F => OrderedDict(:Fe => 2, :O => 3),
    :K => OrderedDict(:K => 2, :O => 1),
    :N => OrderedDict(:Na => 2, :O => 1),
    :P => OrderedDict(:P => 2, :O => 5),
    :T => OrderedDict(:Ti => 1, :O => 2),
    :C̄ => OrderedDict(:C => 1, :O => 2),
    :S̄ => OrderedDict(:S => 1, :O => 3),
    :N̄ => OrderedDict(:N => 1, :O => 3),
    :H => OrderedDict(:H => 2, :O => 1),
]

"""
    OXIDE_ORDER :: Vector{Symbol}

Derived ordered list of cement oxide shorthand symbols, extracted from
`CEMENT_TO_MENDELEEV` while preserving the original sequence. Useful for
deterministic iteration over oxide types.

# Examples

```julia
julia> print(OXIDE_ORDER)
[:C, :M, :S, :A, :F, :K, :N, :P, :T, :C̄, :S̄, :N̄, :H]
```
"""
const OXIDE_ORDER = collect(first.(CEMENT_TO_MENDELEEV))

"""
    CEMDATA_PRIMARIES :: Vector{String}

List of primaries chosen in CEMDATA.

# Examples

```julia
julia> print(CEMDATA_PRIMARIES)
["AlO2-", "Ca+2", "Cl-", "CO3-2", "FeO2-", "H2O@", "H+", "K+", "Mg+2", "Na+", "NO3-", "SiO2@", "SO4-2", "Sr+2", "Zz"]
```
"""
const CEMDATA_PRIMARIES = [
    "AlO2-",
    "Ca+2",
    "Cl-",
    "CO3-2",
    "FeO2-",
    "H2O@",
    "H+",
    "K+",
    "Mg+2",
    "Na+",
    "NO3-",
    "SiO2@",
    "SO4-2",
    "Sr+2",
    "Zz",
]
