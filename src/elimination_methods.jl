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
    anchor = find_closest_site(centre, sites; distance=distance)
    corners = get_corners(TL, BR)
    return any(distance(site, corner) <= distance(anchor, corner) for corner in corners)
end

function centre_anchor_aux(grid, sites, TL, BR; distance=euclidean)
    centre = @. TL + BR / 2
    anchor = find_closest_site(centre, sites; distance=distance)
    corners = get_corners(TL, BR)
    function predicate(site)
        return all(TL .<= site .<= BR) ||
               any(distance(corner, site) <= distance(corner, anchor) for corner in corners)
    end
    return unstable_partition!(predicate, sites)[1]
end

function exact_aux(grid, sites, TL, BR; distance=euclidean)
    predicate(site) = all(TL .<= site .<= BR) || exact_condition(site, sites, TL, BR; distance=distance)
    return unstable_partition!(predicate, sites)[1]
end
