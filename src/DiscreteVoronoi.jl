module DiscreteVoronoi

export find_closest_site # Helper functions
export naive_voronoi!, jfa_voronoi!, dac_voronoi! # Traditional approaches

# Novel work (Joerg Walter and Jacobus Smit):
export exact_condition, exact_elimination, centre_anchor_elimination # Elimination methods
export redac_voronoi! # Reduce-Divide-and-Conquer

using LinearAlgebra: norm
using StaticArrays

include("helper_functions.jl")
include("traditional_approaches.jl")
include("elimination_methods.jl")

# Reduce, Divide, and Conquer
function redac_voronoi!(grid, sites; p=2, elimination=exact_elimination)
    TL = 1, 1
    BR = size(grid)
    _redac_voronoi!(grid, TL, BR, sites, p, elimination)
    return nothing
end

function _redac_voronoi!(grid, TL, BR, sites, p, elimination)
    # Again, first, if the grid is a single cell, we are done
    if all(BR .== 1)
        grid[TL] = find_closest_site(TL, sites, p)
    else
        # Otherwise we check if all corners have the same closest site
        corners = (TL, (TL[1], BR[2]), (BR[1], TL[2]), BR)
        closest_corners = (find_closest_site(corner, sites, p) for corner in corners)
        if allequal(closest_corners)
            grid[TL[1]:BR[1], TL[2]:BR[2]] .= first(closest_corners)
        else
            # And if not we eliminate faraway seeds from subsequent steps
            # (`elimination` returns a view into sites, keeping indexing the same)
            local_sites = elimination(sites, TL, BR)
            # then divide the grid into quadrants and "conquer" each one
            for (quadrant_TL, quadrant_BR) in get_quadrants(TL, BR)
                _dac_voronoi!(grid, quadrant_TL, quadrant_BR, local_sites, p)
            end
        end
    end
    return nothing
end
end # module