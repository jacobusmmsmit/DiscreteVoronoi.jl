function exact_condition(site, sites, TL, BR)
    for other_site in sites
        site == other_site && continue
        if all(distance(site, corner) > distance(other_site, corner) for corner in get_corners(TL, BR))
            return false
        end
    end
    return true
end

function exact_elimination(sites, TL, BR)
    mask = SVector{length(sites),Bool}(exact_condition(site, sites, TL, BR) for site in sites)
    return view(sites, mask)
end

function centre_anchor_elimination(sites, TL, BR)
    centre = @. TL + BR / 2
    anchor = find_closest_site(centre, sites)
    corners = get_corners(TL, BR)

    function cond(site)
        for corner in corners
            distance(site, corner) < distance(anchor, corner) && return true
        end
        return false
    end

    mask = SVector{length(sites),Bool}([cond(site) for site in sites])
    return view(sites, mask)
end