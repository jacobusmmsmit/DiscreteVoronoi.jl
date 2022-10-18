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

function centre_anchor_aux(sites, TL, BR; distance=euclidean)
    centre = @. TL + BR / 2
    anchor = find_closest_site(centre, sites)
    corners = get_corners(TL, BR)
    function predicate(site)
        return all(TL .<= site .<= BR) ||
               any(distance(corner, site) <= distance(corner, anchor) for corner in corners)
    end
    return early_stop_sort!(sites, predicate)
end

function exact_aux(sites, TL, BR; distance=euclidean)
    predicate(site) = all(TL .<= site .<= BR) || exact_condition(site, sites, TL, BR; distance=distance)
    return early_stop_sort!(sites, predicate)
end

function edge_aux(sites, TL, BR; distance=euclidean)
    corners = get_corners(TL, BR)
    keep_per_edge = Iterators.map(1:4) do i
        j = mod(i, 4) + 1
        edge_cells = get_edge(corners[i], corners[j])
        return Iterators.map(cell -> find_closest_site(cell, sites; distance=distance), edge_cells)
    end
    sites_to_keep = Iterators.flatten(keep_per_edge)
    predicate(site) = all(TL .<= site .<= BR) || site in sites_to_keep
    return early_stop_sort!(sites, predicate)
end

function dac_edge_aux(sites, TL, BR; distance=euclidean)
    corners = get_corners(TL, BR)
    keep_per_edge = Iterators.map(1:4) do i
        j = mod(i, 4) + 1
        return edge_dac(sites, corners[i], corners[j]; distance=distance)
    end
    sites_to_keep = Iterators.flatten(keep_per_edge)
    predicate(site) = all(TL .<= site .<= BR) || site in sites_to_keep
    return early_stop_sort!(sites, predicate)
end

function edge_dac(sites, a, b; distance=euclidean)
    all(a .!= b) && throw(ArgumentError("Points $a and $b must be vertically or horizontally aligned"))
    closest = find_closest_site.((a, b), Ref(sites); distance=distance)
    sum(abs, a .- b) == 1 && return closest
    if closest[1] == closest[2]
        return Iterators.repeated(closest, sum(abs, a .- b))
    else
        c = (a .+ b) .÷ 2
        d = c .+ (a[1] == b[1] ? (0, 1) : (1, 0))
        return Iterators.flatten((edge_dac(sites, a, c), edge_dac(sites, b, d)))
    end
end