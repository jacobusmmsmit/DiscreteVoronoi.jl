using Distances
using Polyester

export naive_voronoi
export preset_voronoi!, preset_voronoi_rounded!
export jfa_voronoi!, jfa_voronoi_parallel!
export original_site_find
export center_site_find, center_site_filter, center_site_filter!
export anchor_site_filter, anchor_site_filter!
export dac_voronoi!, redac_voronoi!, redac_voronoi_optimized!

naive_voronoi(grid, sites, distance=euclidean) = map(cell -> findmin(site -> distance(cell, site), sites)[2], grid)

function preset_voronoi!(grid, sites)
    for (color, point) in sites
        if 1 <= point[1] <= size(grid, 1) && 1 <= point[2] <= size(grid, 2)
            grid[point...] = convert(eltype(grid), color)
        end
    end
    grid
end

function preset_voronoi_rounded!(grid, sites)
    for (color, point) in sites
        closest_point = round_tuple(point)
        if 1 <= closest_point[1] <= size(grid, 1) && 1 <= closest_point[2] <= size(grid, 2)
            grid[closest_point...] = convert(eltype(grid), color)
        end
    end
    grid
end

function jfa_voronoi!(grid, points, distance=euclidean)
    for (color, (x, y)) in enumerate(points)
        grid[x, y] = convert(eltype(grid), color)
    end
    k = max(size(grid)...)
    while k > 1
        k = k ÷ 2 + k % 2
        @inbounds for y in axes(grid, 2), x in axes(grid, 1)
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

function jfa_voronoi_parallel!(grid, points, distance=euclidean)
    for (color, (x, y)) in enumerate(points)
        grid[x, y] = convert(eltype(grid), color)
    end
    gridp = copy(grid)
    k = max(size(grid)...)
    while k > 1
        k = k ÷ 2 + k % 2
        @inbounds @batch for y in axes(grid, 2), x in axes(grid, 1)
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

@inbounds function conquer_base_cases!(grid, sites, rect, distance)
    (t, l), (b, r) = rect
    all(!=(0), @view grid[t:b, l:r]) && return true
    if length(sites) == 1
        grid[t:b, l:r] .= convert(eltype(grid), sites[1][1])
        return true
    end
    if t == b && l == r
        _, min_color = _findmin(sites) do site
            distance((t, l), site[2])
        end
        grid[t, l] = convert(eltype(grid), sites[min_color][1])
        return true
    end
    return false
end

@inbounds function original_site_find(grid, sites, rect, distance)
    conquer_base_cases!(grid, sites, rect, distance) && return true
    color = 0
    corners = get_corners(rect)
    for corner in corners
        min_dist, min_color = _findmin(sites) do site
            distance(corner, site[2])
        end
        dist = _minimum(point for (color, point) in sites if color != min_color) do point
            distance(corner, point)
        end
        dist > min_dist || return false
        color == 0 && (color = min_color)
        min_color == color || return false
    end
    (t, l), (b, r) = rect
    grid[t:b, l:r] .= convert(eltype(grid), sites[color][1])
    return true
end

@inbounds function center_site_find(grid, sites, rect, distance)
    conquer_base_cases!(grid, sites, rect, distance) && return true
    center = get_center(rect)
    min_dist, min_color = _findmin(sites) do site
        distance(center, site[2])
    end
    dist = _minimum(point for (color, point) in sites if color != min_color) do point
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

@inbounds function dac_voronoi!(grid, sites, conquer, distance=euclidean, depth::Int=1, rect=((1,1), size(grid)))
    if !conquer(grid, sites, rect, distance)
        quadrants = get_quadrants(rect)
        if depth > 0
            Threads.@threads for quadrant in quadrants
                dac_voronoi!(grid, sites, conquer, distance, depth - 1, quadrant)
            end
        else
            for quadrant in quadrants
                dac_voronoi!(grid, sites, conquer, distance, depth - 1, quadrant)
            end
        end
    end
    return grid
end

@inbounds function center_site_filter(grid, sites, rect, distance)
    conquer_base_cases!(grid, sites, rect, distance) && return true, ()
    center = get_center(rect)
    min_dist = _minimum(sites) do site
        distance(center, site[2])
    end
    (t, l), (b, r) = rect
    max_dist = min_dist + distance((b, r), (t, l)) + 1
    return false, collect(_filter(sites) do site
        distance(center, site[2]) <= max_dist
    end)
end

@inbounds function anchor_site_filter(grid, sites, rect, distance)
    conquer_base_cases!(grid, sites, rect, distance) && return true, ()
    center = get_center(rect)
    _, min_color = _findmin(sites) do site
        distance(center, site[2])
    end
    corners = get_corners(rect)
    return false, collect(_filter(sites) do site
        _any(corners) do corner
            distance(corner, site[2]) <= distance(corner, sites[min_color][2])
        end
    end)
end

@inbounds function redac_voronoi!(grid, sites, conquer, distance=euclidean, depth::Int=1, rect=((1, 1), size(grid)))
    conquered, sites = conquer(grid, sites, rect, distance)
    if !conquered
        quadrants = get_quadrants(rect)
        if depth > 0
            Threads.@threads for quadrant in quadrants
                redac_voronoi!(grid, sites, conquer, distance, depth - 1, quadrant)
            end
        else
            for quadrant in quadrants
                redac_voronoi!(grid, sites, conquer, distance, depth - 1, quadrant)
            end
        end
    end
    grid
end

@inbounds function center_site_filter!(grid, sites, rect, distance, stack)
    conquer_base_cases!(grid, sites, rect, distance) && return true
    center = get_center(rect)
    min_dist = _minimum(sites) do site
        distance(center, site[2])
    end
    (t, l), (b, r) = rect
    max_dist = min_dist + distance((b, r), (t, l)) + 1
    set_sites!(stack, _filter(sites) do site
        distance(center, site[2]) <= max_dist
    end)
    return false
end

@inbounds function anchor_site_filter!(grid, sites, rect, distance, stack)
    conquer_base_cases!(grid, sites, rect, distance) && return true
    center = get_center(rect)
    _, min_color = _findmin(sites) do site
        distance(center, site[2])
    end
    corners = get_corners(rect)
    set_sites!(stack, _filter(sites) do site
        _any(corners) do corner
            distance(corner, site[2]) <= distance(corner, sites[min_color][2])
        end
    end)
    return false
end

@inbounds function redac_voronoi_optimized!(grid, sites, conquer!, distance=euclidean, depth::Int=1, rect=((1, 1), size(grid)), stack=SiteStack{eltype(sites)}())
    push_empty!(stack)
    if !conquer!(grid, sites, rect, distance, stack)
        sites = get_sites(stack)
        quadrants = get_quadrants(rect)
        if depth > 0
            Threads.@threads for quadrant in quadrants
                redac_voronoi_optimized!(grid, sites, conquer!, distance, depth - 1, quadrant, SiteStack{eltype(sites)}())
            end
        else
            for quadrant in quadrants
                redac_voronoi_optimized!(grid, sites, conquer!, distance, depth - 1, quadrant, stack)
            end
        end
    end
    pop!(stack)
    grid
end
