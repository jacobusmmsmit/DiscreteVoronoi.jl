function distance(x, y, p=2)
    return norm(x .- y, p)
end

"""
    get_corners(TL, BR)

Returns the corners of the rectangle defined by its top-left (TL) and bottom-right (BR) corners
"""
function get_corners(TL, BR)
    return (TL, (TL[1], BR[2]), (BR[1], TL[2]), BR)
end

"""
    swap!(v, i::Int, j::Int)

Swap (in-place) the elements of `v` indexed by `i` and `j`. Does nothing if `i == j`.
"""
function swap!(v, i::Int, j::Int)
    if i != j
        v[i], v[j] = v[j], v[i]
    end
    return nothing
end

"""
    get_quadrants(TL, BR)

Returns a tuple containing top-left (TL) and bottom-right (BR) corners for each quadrant of input TL-BR
"""
function get_quadrants(TL, BR)
    MID_TL = (TL .+ BR) .÷ 2
    MID_BR = TL .+ BR .- MID_TL
    return (
        (TL, MID_TL),
        ((MID_BR[1], TL[2]), (BR[1], MID_TL[2])),
        ((TL[1], MID_BR[2]), (MID_TL[1], BR[2])),
        (MID_BR, BR),
    )
end

"""
    find_closest_site(cell, sites, p=2)

See name. `p` is which Lᵖ norm to use.
"""
function find_closest_site(cell, sites, p=2)
    first_site, rest_sites = Iterators.peel(sites)
    isnothing(rest_sites) && return first_site
    min_dist = distance(cell, first_site, p)
    min_site = first_site
    for site in rest_sites
        cur_dist = distance(cell, site, p)
        if cur_dist < min_dist
            min_dist = cur_dist
            min_site = site
        end
    end
    return min_site
end

"""
    label_voronoi_grid(grid; shuffle_cells=true)

Given a grid of un-plottable things, it assigns each unique value in the grid a number to plot.
`shuffle_cells` determines whether the assigned values are randomly shuffled before assignment.
If they're not you might get a pretty gradient or other pattern :)
"""
function label_voronoi_grid(grid; shuffle_cells=true)
    if shuffle_cells
        maybe_shuffle = shuffle
    else
        maybe_shuffle = identity
    end
    labelled_grid = similar(grid, Int)
    for (i, val) in (enumerate ∘ maybe_shuffle ∘ unique)(grid)
        labelled_grid[grid.==Ref(val)] .= i
    end
    return labelled_grid
end

"""
    voronoi_equality(grids...)

Checks equality of Voronoi diagrams accounting for the fact that there may be
multiple correct/equivalent solutions as some sites may be the same distance
from some cells.
"""
function voronoi_equality(grid1, grid2)
    size(grid1) == size(grid2) || throw(
        ArgumentError(
            "Grids should have the same shape, got shapes $(size(grid1)) and $(size(grid2))",
        ),
    )
    for I in CartesianIndices(grid1)
        cell = SVector{2,Int}(Tuple(I))
        grid1[I] == grid2[I] && continue
        distance(cell, grid1[I]) == distance(cell, grid2[I]) || return false
    end
    return true
end

voronoi_equality(grids...) = reduce(voronoi_equality, grids)
