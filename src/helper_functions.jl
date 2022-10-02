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

Swap (in-place) the elements of `v` indexed by `i` and `j`.
"""
function swap!(v, i::Int, j::Int) 
    v[i], v[j] = v[j], v[i]
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
        (MID_BR, BR)
    )
end

"""
    find_closest_site(cell, sites, p=2)

See name. `p` is which Lᵖ norm to use.
"""
function find_closest_site(cell, sites, p=2)
    first_site, rest_sites = Iterators.peel(sites)
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
    label(grid; shuffle_cells = true)

Given a grid of un-plottable things, it assigns each unique value in the grid a number to plot.
`shuffle_cells` determines whether the assigned values are randomly shuffled before assignment.
If they're not you might get a pretty gradient or other pattern :)
"""
function label(grid; shuffle_cells = true)
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