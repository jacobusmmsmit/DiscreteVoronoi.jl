using Random: shuffle

function distance(x, y, p=2)
    return norm(x .- y, p)
end


function get_corners(TL, BR)
    return (TL, (TL[1], BR[2]), (BR[1], TL[2]), BR)
end


"""
    get_quadrants(TL, BR)

Returns a tuple containing TL and BR corners for each quadrant of input TL-BR
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

function label(grid)
    labelled_grid = similar(grid, Int)
    # Remove shuffle if you want a pretty gradient
    for (i, val) in (enumerate ∘ shuffle ∘ unique)(grid)
        labelled_grid[grid.==Ref(val)] .= i
    end
    return labelled_grid
end

# # Testing allocations of static bit vector
# sites = [@SVector([1, 1]), @SVector([2, 2]), @SVector([3, 3])]
# new_sites = exact_elimination(sites, (1, 1), (2, 2))
# new_sites2 = exact_elimination(new_sites, (1, 1), (2, 2))
# @ballocated exact_elimination($sites, (1, 1), (1, 1))
# @ballocated find_closest_site((1, 1), $sites)