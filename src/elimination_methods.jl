function exact_condition(site, sites, TL, BR)
    for other_site in sites
        if site == other_site
            continue
        end
        if all(distance(site, corner) > distance(other_site, corner) for corner in get_corners(TL, BR))
            return false
        end
    end
    return true
end

function centre_anchor_condition(site, sites, TL, BR)
    centre = @. TL + BR / 2
    anchor = find_closest_site(centre, sites)
    corners = get_corners(TL, BR)
    return any(distance(site, corner) <= distance(anchor, corner) for corner in corners)
end