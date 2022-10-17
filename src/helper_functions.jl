"""
    get_corners(TL, BR)

Returns the corners of the rectangle defined by its top-left (TL) and
bottom-right (BR) corners. The order of returning is TL, TR, BR, BL such that
the pairs of corners can be looped over to get edges.
"""
function get_corners(TL, BR)
    return (TL, (TL[1], BR[2]), BR, (BR[1], TL[2]))
end


"""
    get_edge(TL, BR)

Returns the edge (cells) between two vertically or horizontally aligned points.
"""
function get_edge(a, b)
    all(a .!= b) && throw(ArgumentError("Points $a and $b must be vertically or horizontally aligned"))
    if a[1] == b[1]
        x, y = extrema((a[2], b[2]))
        return (SVector{2,Int}(a[1], j) for j in x:1:y)
    elseif a[2] == b[2]
        x, y = extrema((a[1], b[1]))
        return (SVector{2,Int}(i, a[2]) for i in x:1:y)
    end
end

"""
    get_edges(TL, BR)

Returns a generator containing the cells on the border edges of the rectangle
defined by its top-left (TL) and bottom-right (BR) corners.
"""
function get_edges(TL, BR)
    corners = get_corners(TL, BR)
    edge_cells = map(1:4) do i
        j = mod(i, 4) + 1
        return get_edge(corners[i], corners[j])
    end
    return Iterators.flatten(edge_cells)
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
    MID_BR = MID_TL .+ 1
    return (
        (TL, MID_TL),
        ((MID_BR[1], TL[2]), (BR[1], MID_TL[2])),
        ((TL[1], MID_BR[2]), (MID_TL[1], BR[2])),
        (MID_BR, BR),
    )
end

"""
    find_closest_site(cell, sites; distance=euclidean)

Returns the closest site to `cell` in `sites` determined by `distance`.
"""
function find_closest_site(cell, sites; distance=euclidean)
    first_site, rest_sites = Iterators.peel(sites)
    isnothing(rest_sites) && return first_site
    min_dist = distance(cell, first_site)
    min_site = first_site
    for site in rest_sites
        cur_dist = distance(cell, site)
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
    voronoi_equality(grid1, grid2; distance=euclidean)

Checks equality of Voronoi diagrams accounting for the fact that there may be
multiple correct/equivalent solutions as some sites may be the same distance
from some cells.
"""
function voronoi_equality(grid1, grid2; distance=euclidean)
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
