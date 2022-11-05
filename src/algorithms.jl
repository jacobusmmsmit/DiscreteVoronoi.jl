using Distances
using Polyester

export preset_voronoi!, preset_voronoi_rounded!
export naive_voronoi!
export jfa_voronoi!, jfa_voronoi_parallel!
export original_site_find, center_site_find
export dac_voronoi!
export original_site_filter, center_site_filter, anchor_site_filter, corner_site_filter, edge_site_filter
export redac_voronoi!

@inbounds function preset_voronoi!(grid, sites)
    for (color, (x, y)) in sites
        if checkbounds(Bool, grid, x, y)
            grid[x, y] = convert(eltype(grid), color)
        end
    end
    return grid
end

@inbounds function preset_voronoi_rounded!(grid, sites)
    for (color, (x, y)) in sites
        xp, yp = round_tuple((x, y))
        if checkbounds(Bool, grid, xp, yp)
            grid[xp, yp] = convert(eltype(grid), color)
        end
    end
    return grid
end

@inbounds function naive_voronoi!(grid, points, distance=euclidean)
    for y in axes(grid, 2), x in axes(grid, 1)
        _, min_color = findmin(points) do point
            distance((x, y), point)
        end
        grid[x, y] = convert(eltype(grid), min_color)
    end
    return grid
end

@inbounds function jfa_voronoi!(grid, points, distance=euclidean)
    for (color, (x, y)) in enumerate(points)
        grid[x, y] = convert(eltype(grid), color)
    end
    k = max(size(grid)...)
    while k > 1
        k = k ÷ 2 + k % 2
        for y in axes(grid, 2), x in axes(grid, 1)
            for j in (-k, 0, k), i in (-k, 0, k)
                ((i !== 0 || j !== 0) && checkbounds(Bool, grid, x + i, y + j) && (colorq = grid[x + i, y + j]) !== 0) || continue
                if (colorp = grid[x, y]) === 0 || distance(points[colorp], (x, y)) > distance(points[colorq], (x, y))
                    grid[x, y] = colorq
                end
            end
        end
    end
    return grid
end

@inbounds function jfa_voronoi_parallel!(grid, points, distance=euclidean)
    for (color, (x, y)) in enumerate(points)
        grid[x, y] = convert(eltype(grid), color)
    end
    gridp = copy(grid)
    k = max(size(grid)...)
    while k > 1
        k = k ÷ 2 + k % 2
        @batch for y in axes(grid, 2), x in axes(grid, 1)
            for j in (-k, 0, k), i in (-k, 0, k)
                ((i !== 0 || j !== 0) && checkbounds(Bool, grid, x + i, y + j) && (colorq = grid[x + i, y + j]) !== 0) || continue
                if (colorp = grid[x, y]) === 0 || distance(points[colorp], (x, y)) > distance(points[colorq], (x, y))
                    gridp[x, y] = colorq
                end
            end
        end
        grid .= gridp
    end
    return grid
end

@inbounds function base_cases!(grid, sites, rect, distance)
    (t, l), (b, r) = rect
    if length(sites) == 1
        grid[t:b, l:r] .= convert(eltype(grid), sites[1][1])
        return true
    end
    if t == b && l == r 
        if grid[t, l] == zero(eltype(grid))
            _, min_color = findmin(sites) do site
                distance((t, l), site[2])
            end
            grid[t, l] = convert(eltype(grid), sites[min_color][1])
        end
        return true
    end
    return all(!=(zero(eltype(grid))), @view grid[t:b, l:r])
end

@inbounds function original_site_find(grid, sites, rect, distance)
    corners = get_corners(rect)
    mins = ((findmin(sites) do site
        distance(corner, site[2])
    end for corner in corners)...,)
    allequal(color for (dist, color) in mins) || return false
    min_color = minimum(color for (dist, color) in mins)
    dists = ((minimum(point for (color, point) in sites if color != min_color) do point
        distance(corner, point)
    end for corner in corners)...,)
    any(zip(mins, dists)) do ((min_dist, _), dist)
        dist <= min_dist
    end && return false
    (t, l), (b, r) = rect
    grid[t:b, l:r] .= convert(eltype(grid), min_color)
    return true
end

@inbounds function center_site_find(grid, sites, rect, distance)
    center = get_center(rect)
    min_dist, min_color = findmin(sites) do site
        distance(center, site[2])
    end
    dist = minimum(point for (color, point) in sites if color != min_color) do point
        distance(center, point)
    end
    (t, l), (b, r) = rect
    max_dist = min_dist + distance((b, r), (t, l)) + 1
    if dist > max_dist
        grid[t:b, l:r] .= convert(eltype(grid), sites[min_color][1])
        return true
    end
    return false
end

@inbounds function dac_voronoi!(grid, sites, site_find, distance=euclidean, depth::Int=1, rect=((1,1), size(grid)))
    base_cases!(grid, sites, rect, distance) && return grid
    site_find(grid, sites, rect, distance) && return grid
    quadrants = get_quadrants(rect)
    if depth > 0
        Threads.@threads for quadrant in quadrants
            dac_voronoi!(grid, sites, site_find, distance, depth - 1, quadrant)
        end
    else
        for quadrant in quadrants
            dac_voronoi!(grid, sites, site_find, distance, depth - 1, quadrant)
        end
    end
    return grid
end

@inbounds function original_site_filter(::Val{:filter}, grid, sites, rect, distance)
    corners = get_corners(rect)
    return collect(filter(sites) do site1
        in_rectangle(rect, site1[2]) ||
        all(sites) do site2
            any(corners) do corner
                distance(site1[2], corner) <= distance(site2[2], corner)
            end
        end
    end)
end
@inbounds function original_site_filter(::Val{:partition}, grid, sites, rect, distance)
    corners = get_corners(rect)
    return unstable_partition!(sites) do site1
        in_rectangle(rect, site1[2]) ||
        all(sites) do site2
            any(corners) do corner
                distance(site1[2], corner) <= distance(site2[2], corner)
            end
        end
    end[1]
end

@inbounds function center_site_filter(::Val{:filter}, grid, sites, rect, distance)
    center = get_center(rect)
    min_dist = minimum(sites) do site
        distance(center, site[2])
    end
    (t, l), (b, r) = rect
    max_dist = min_dist + distance((b, r), (t, l)) + 1
    return collect(filter(sites) do site
        in_rectangle(rect, site[2]) ||
        distance(center, site[2]) <= max_dist
    end)
end
@inbounds function center_site_filter(::Val{:partition}, grid, sites, rect, distance)
    center = get_center(rect)
    min_dist = minimum(sites) do site
        distance(center, site[2])
    end
    (t, l), (b, r) = rect
    max_dist = min_dist + distance((b, r), (t, l)) + 1
    return unstable_partition!(sites) do site
        in_rectangle(rect, site[2]) ||
        distance(center, site[2]) <= max_dist
    end[1]
end

@inbounds function anchor_site_filter(::Val{:filter}, grid, sites, rect, distance)
    center = get_center(rect)
    _, min_color = findmin(sites) do site
        distance(center, site[2])
    end
    corners = get_corners(rect)
    return collect(filter(sites) do site
        in_rectangle(rect, site[2]) ||
        any(corners) do corner
            distance(site[2], corner) <= distance(sites[min_color][2], corner)
        end
    end)
end
@inbounds function anchor_site_filter(::Val{:partition}, grid, sites, rect, distance)
    center = get_center(rect)
    _, min_color = findmin(sites) do site
        distance(center, site[2])
    end
    corners = get_corners(rect)
    return unstable_partition!(sites) do site
        in_rectangle(rect, site[2]) ||
        any(corners) do corner
            distance(site[2], corner) <= distance(sites[min_color][2], corner)
        end
    end[1]
end

@inbounds function corner_site_filter(::Val{:filter}, grid, sites, rect, distance)
    corners = get_corners(rect)
    mins = ((findmin(sites) do site
        distance(corner, site[2])
    end for corner in corners)...,)
    return collect(filter(sites) do site
        in_rectangle(rect, site[2]) ||
        all(mins) do (_, min_color)
            any(corners) do corner
                distance(site[2], corner) <= distance(sites[min_color][2], corner)
            end
        end
    end)
end
@inbounds function corner_site_filter(::Val{:partition}, grid, sites, rect, distance)
    corners = get_corners(rect)
    mins = ((findmin(sites) do site
        distance(corner, site[2])
    end for corner in corners)...,)
    return unstable_partition!(sites) do site
        in_rectangle(rect, site[2]) ||
        all(mins) do (_, min_color)
            any(corners) do corner
                distance(site[2], corner) <= distance(sites[min_color][2], corner)
            end
        end
    end[1]
end

@inbounds function edge_site_filter(::Val{:filter}, grid, sites, rect, distance)
    edges = get_edges(rect)
    mins = unique(((findmin(sites) do site
        distance(edge, site[2])
    end for edge in edges)...,))
    corners = get_corners(rect)
    return collect(filter(sites) do site
        in_rectangle(rect, site[2]) ||
        all(mins) do (_, min_color)
            any(corners) do corner
                distance(site[2], corner) <= distance(sites[min_color][2], corner)
            end
        end
    end)
end
@inbounds function edge_site_filter(::Val{:partition}, grid, sites, rect, distance)
    edges = get_edges(rect)
    mins = unique(((findmin(sites) do site
        distance(edge, site[2])
    end for edge in edges)...,))
    corners = get_corners(rect)
    return unstable_partition!(sites) do site
        in_rectangle(rect, site[2]) ||
        all(mins) do (_, min_color)
            any(corners) do corner
                distance(site[2], corner) <= distance(sites[min_color][2], corner)
            end
        end
    end[1]
end

@inbounds function redac_voronoi!(method::Val{:filter}, grid, sites, site_filter, distance=euclidean, depth::Int=1, rect=((1, 1), size(grid)))
    local_sites = site_filter(method, grid, sites, rect, distance)
    base_cases!(grid, local_sites, rect, distance) && return grid
    quadrants = get_quadrants(rect)
    if depth > 0
        Threads.@threads for quadrant in quadrants
            redac_voronoi!(method, grid, local_sites, site_filter, distance, depth - 1, quadrant)
        end
    else
        for quadrant in quadrants
            redac_voronoi!(method, grid, local_sites, site_filter, distance, depth - 1, quadrant)
        end
    end
    grid
end
@inbounds function redac_voronoi!(method::Val{:partition}, grid, sites, site_filter, distance=euclidean, depth::Int=1, rect=((1, 1), size(grid)))
    local_sites = site_filter(method, grid, sites, rect, distance)
    base_cases!(grid, local_sites, rect, distance) && return grid
    quadrants = get_quadrants(rect)
    if depth > 0
        Threads.@threads for quadrant in quadrants
            redac_voronoi!(method, grid, copy(local_sites), site_filter, distance, depth - 1, quadrant)
        end
    else
        for quadrant in quadrants
            redac_voronoi!(method, grid, local_sites, site_filter, distance, depth - 1, quadrant)
        end
    end
    grid
end
