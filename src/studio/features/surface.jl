# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright © 2026 Jean-François Barthélémy (Cerema, UMR MCD)

# The surface of the C-S-H: the sites of the published model of Guo et al.
# (2018) on the gel of the paste, frozen at its composition, and the paste
# solved again with them (ChemistryLab's `with_csh_surface`), with a Donnan
# layer on the diffuse layer if asked. It follows an equilibrium, or each point
# of a path: an addition of salt gives the binding of chloride.

"""
    SurfaceChoice

Whether, and how, the surface of the C-S-H is computed after the equilibrium.
"""
abstract type SurfaceChoice <: Choice end

"No surface: the paste as its phases make it."
struct NoSurfaceChoice <: SurfaceChoice end
kind(::Type{NoSurfaceChoice}) = "no_surface"
fields(::Type{NoSurfaceChoice}) = Field[]

"""
The surface of the C-S-H after Guo et al. (2018), with a Donnan layer of
`thickness` nm on its diffuse layer if `donnan` (0 for one Debye length).
"""
struct CshSurfaceChoice <: SurfaceChoice
    donnan::Bool
    thickness::Float64
end
kind(::Type{CshSurfaceChoice}) = "csh_surface"
fields(::Type{CshSurfaceChoice}) = [
    Field("donnan", :bool; default = false),
    Field("thickness", :number; default = 0.0, unit = "nm", min = 0, max = 1000, level = :expert),
]

members(::Type{SurfaceChoice}) = (NoSurfaceChoice, CshSurfaceChoice)

"What keeps the surface from following the calculation `c`."
surface_problems(::NoSurfaceChoice, c) = Problem[]
function surface_problems(::CshSurfaceChoice, c)
    mode(typeof(c)) == "kinetics" || return Problem[]
    return [Problem("surface_kinetics", "The surface of the C-S-H follows an equilibrium or a path of equilibria, not a kinetic run.")]
end

"""
    surface_lines(s, c) -> Vector{String}

The lines that follow the calculation `c` to compute the surface of its paste:
none without a surface.
"""
surface_lines(::NoSurfaceChoice, c) = String[]
function surface_lines(s::CshSurfaceChoice, c)
    kw = !s.donnan ? "" : s.thickness > 0 ? "; donnan = $(julia_number(s.thickness))e-9" : "; donnan = :debye"
    target = mode(typeof(c)) == "equilibrium" ? "with_csh_surface(state$kw)" : "[with_csh_surface(rs$kw) for rs in path]"
    return vcat(
        [
            "# The surface of the C-S-H, after Guo et al. (2018): the sites of their model on the",
            "# gel, frozen at its composition (portlandite holding the activities it shares with",
            "# it), and the paste solved again with them.",
        ],
        s.donnan ? [
                "# A Donnan layer of water on the diffuse layer, $(s.thickness > 0 ? "$(julia_number(s.thickness)) nm thick" : "one Debye length thick"), its water",
                "# taken from the pore solution.",
            ] : String[],
        [
            "surface = $target",
            "study = merge(study, (; surface))",
        ],
    )
end
