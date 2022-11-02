# This is non-exported, non-userfacing code.
# This file contains the elimination methods used in `redac_voronoi!`.

function exact_condition(site, sites, TL, BR; distance=euclidean)
    for other_site in sites
        site == other_site && continue
        if all(
            distance(site, corner) > distance(other_site, corner) for
            corner in get_corners(TL, BR)
        )
            return false
        end
    end
    return true
end

function centre_anchor_condition(site, sites, TL, BR; distance=euclidean)
    centre = @. TL + BR / 2
    anchor = find_closest_site(centre, sites)
    corners = get_corners(TL, BR)
    return any(distance(site, corner) <= distance(anchor, corner) for corner in corners)
end

function centre_anchor_aux(grid, sites, TL, BR; distance=euclidean)
    centre = @. TL + BR / 2
    anchor = find_closest_site(centre, sites)
    corners = get_corners(TL, BR)
    function predicate(site)
        return all(TL .<= site .<= BR) ||
               any(distance(corner, site) <= distance(corner, anchor) for corner in corners)
    end
    return early_stop_sort!(sites, predicate)
end

function exact_aux(grid, sites, TL, BR; distance=euclidean)
    predicate(site) = all(TL .<= site .<= BR) || exact_condition(site, sites, TL, BR; distance=distance)
    return early_stop_sort!(sites, predicate)
end

function naive_edge_aux!(grid, sites, TL, BR; distance=euclidean)
    corners = get_corners(TL, BR)
    edges = (CartesianIndex(corners[a]...):CartesianIndex(corners[b]...) for (a, b) in zip((1, 1, 2, 4), (2, 4, 3, 3)))
    sites_to_keep = (find_closest_site!(grid, Coord(Tuple(I)), sites; distance=distance) for edge in edges for I in edge)
    predicate(site) = all(TL .<= site .<= BR) || site in sites_to_keep
    return early_stop_sort!(sites, predicate)
end

function dac_edge_aux!(grid, sites, TL, BR; distance=euclidean)
    corners = get_corners(TL, BR)
    edges = ((corners[a], corners[b]) for (a, b) in zip((1, 1, 2, 4), (2, 4, 3, 3)))
    sites_to_keep = (edge_dac!(grid, sites, edge[1], edge[2], distance=distance) for edge in edges)
    predicate(site) = all(TL .<= site .<= BR) || site in sites_to_keep
    return early_stop_sort!(sites, predicate)
end

function edge_dac!(grid, sites, a, b; distance=euclidean)
    all(a .!= b) && throw(ArgumentError("Points $a and $b must be vertically or horizontally aligned"))
    closest = find_closest_site!(grid, a, sites; distance=distance), find_closest_site!(grid, b, sites; distance=distance)
    sum(abs, a .- b) == 1 && return closest
    if closest[1] == closest[2]
        return Iterators.repeated(closest[1], sum(abs, a .- b))
    else
        c = (a .+ b) .÷ 2
        d = c .+ (a[1] == b[1] ? (0, 1) : (1, 0))
        return (edge_dac!(grid, sites, a, c), edge_dac!(grid, sites, b, d))
    end
end