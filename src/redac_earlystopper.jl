# Reduce, Divide, and Conquer
"""
    site_sort!(sites::ES, keep_predicate, TL, BR) where {ES<:EarlyStopper}

Produce an anonymous function to pass to `early_stop_sort!`.
"""
function site_sort!(sites::ES, keep_predicate, TL, BR) where {ES<:EarlyStopper}
    predicate(site) = keep_predicate(site, sites, TL, BR)
    return early_stop_sort!(sites, predicate)
end

site_sort!(sites, keep_predicate) = site_sort!(EarlyStopper(sites), keep_predicate)

"""
    redac_voronoi!(grid, sites; p=2, predicate=exact_condition)

Performs a divide-and-conquer method similar to `dac_voronoi!` but has an additional site-elimination
step which aims to reduce the work of subsequent steps.
"""
function redac_voronoi!(grid, sites; p=2, predicate=exact_condition)
    TL = 1, 1
    BR = size(grid)
    _redac_voronoi!(grid, TL, BR, EarlyStopper(sites), p, predicate)
    return nothing
end

function _redac_voronoi!(grid, TL, BR, sites::ES, p, predicate) where {ES<:EarlyStopper}
    # Again, first, if the grid is a single cell, we are done
    if all((BR .- TL) .== 0)
        grid[TL...] = find_closest_site(TL, sites, p)
    else
        # Otherwise we check if all corners have the same closest site
        corners = (TL, (TL[1], BR[2]), (BR[1], TL[2]), BR)
        closest_corners = (find_closest_site(corner, sites, p) for corner in corners)
        if allequal(closest_corners)
            view(grid, TL[1]:BR[1], TL[2]:BR[2]) .= Ref(first(closest_corners))
        else
            # And if not we eliminate faraway seeds from subsequent steps
            # `site_sort` sorts sites by whether the predicate is true and stores how many are true.
            local_sites = site_sort!(sites, predicate, TL, BR)
            # then divide the grid into quadrants and "conquer" each one
            for (quadrant_TL, quadrant_BR) in get_quadrants(TL, BR)
                _redac_voronoi!(grid, quadrant_TL, quadrant_BR, local_sites, p, predicate)
            end
        end
    end
    return nothing
end