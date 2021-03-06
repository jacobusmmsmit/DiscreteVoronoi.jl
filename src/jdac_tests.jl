export jdac!, jdacx!, jdac_aux0!
export jdac_aux1a!, jdac_aux1b!, jdac_aux1c!
export jdac_aux2a!, jdac_aux2b!, jdac_aux2c!
export jdac_aux3a!, jdac_aux3b!, jdac_aux3c!
export jdac_aux4a!, jdac_aux4c!
export jdac_aux5a!, jdac_aux5c!


"""
    exact_site_filter(sites, corners, p)
Returns {s ā š® : ā s' ā š® s.t. ā C ā {Cā, Cā, Cā, Cā}, d(s', C) < d(s, C)}
Returns sites such that there does not exist a site such that all corners
are closer to the site than the candidate.
"""
function exact_site_filter(sites, corners, p)
    filter(sites) do i_s
        s = i_s[2]
        include = true
        for (_, t) in sites
            nc = 0
            for c in corners
                if distance(s, c, p) <= distance(t, c, p)
                    break
                end
                nc += 1
            end
            if nc == 4
                include = false
            end
        end
        include
    end
end

@inbounds function jdac_aux0!(grid, sites, p, depth, rect, stack)
    (t, l), (N, M) = rect
    if all(>(0), (N, M)) && any(==(0), @view grid[t:t+N-1, l:l+M-1])
        corners = (t, l), (t + N - 1, l), (t, l + M - 1), (t + N - 1, l + M - 1)
        stack_sites = exact_site_filter(sites, corners, p)
        jdac!(grid, stack_sites, jdac_aux0!, p, depth - 1, rect, stack)
    end
    grid
end

@inbounds function jdac_aux1a!(grid, sites, p, depth, rect, stack)
    (t, l), (N, M) = rect
    if all(>(0), (N, M)) && any(==(0), @view grid[t:t+N-1, l:l+M-1])
        center = (t+(N-1)/2, l+(M-1)/2) 
        max_dist = findmin(site -> distance(center, site[2], p), sites)[1] + distance((1, 1), (N, M), p) + 1
        stack_sites = filter(site -> distance(center, site[2], p) <= max_dist, sites)
        jdac!(grid, stack_sites, jdac_aux1a!, p, depth - 1, rect, stack)
    end
    grid
end

@inbounds function jdac_aux1b!(grid, sites, p, depth, rect, stack)
    (t, l), (N, M) = rect
    if all(>(0), (N, M)) && any(==(0), @view grid[t:t+N-1, l:l+M-1])
        center = (t+(N-1)/2, l+(M-1)/2) 
        dists = map(site -> distance(center, site[2], p), sites)
        max_dist = findmin(dists)[1] + distance((1, 1), (N, M), p) + 1
        stack_sites = [site for (dist, site) in zip(dists, sites) if dist <= max_dist]
        jdac!(grid, stack_sites, jdac_aux1b!, p, depth - 1, rect, stack)
    end
    grid
end

@inbounds function jdac_aux1c!(grid, sites, p, depth, rect, stack)
    (t, l), (N, M) = rect
    if all(>(0), (N, M)) && any(==(0), @view grid[t:t+N-1, l:l+M-1])
        center = (t+(N-1)/2, l+(M-1)/2) 
        push_empty!(stack)
        fill_dists!(stack, center, sites, p)
        max_dist = findmin(get_dists(stack))[1] + distance((1, 1), (N, M), p) + 1
        fill_sites!(stack, max_dist, sites)
        jdac!(grid, get_sites(stack), jdac_aux1c!, p, depth - 1, rect, stack)
        pop!(stack)
    end
    grid
end

@inbounds function jdac_aux2a!(grid, sites, p, depth, rect, stack)
    (t, l), (N, M) = rect
    if all(>(0), (N, M)) && any(==(0), @view grid[t:t+N-1, l:l+M-1])
        center = (t+(N-1)/2, l+(M-1)/2) 
        dist, index = findmin(site -> distance(center, site[2], p), sites)
        max_dist = dist + distance((1, 1), (N, M), p) + 1
        min_site = sites[index][2]
        corners = (t, l), (t + N - 1, l), (t + N - 1, l + M - 1), (t, l + M - 1)
        dist = maximum(corner -> distance(corner, min_site, p), corners) + distance((1, 1), center, p) + 1
        max_dist = min(max_dist, dist)
        stack_sites = filter(site -> distance(center, site[2], p) <= max_dist, sites)
        jdac!(grid, stack_sites, jdac_aux2a!, p, depth - 1, rect, stack)
    end
    grid
end

@inbounds function jdac_aux2b!(grid, sites, p, depth, rect, stack)
    (t, l), (N, M) = rect
    if all(>(0), (N, M)) && any(==(0), @view grid[t:t+N-1, l:l+M-1])
        center = (t+(N-1)/2, l+(M-1)/2) 
        dists = map(site -> distance(center, site[2], p), sites)
        dist, index = findmin(dists)
        max_dist = dist + distance((1, 1), (N, M), p) + 1
        min_site = sites[index][2]
        corners = (t, l), (t + N - 1, l), (t + N - 1, l + M - 1), (t, l + M - 1)
        dist = maximum(corner -> distance(corner, min_site, p), corners) + distance((1, 1), center, p) + 1
        max_dist = min(max_dist, dist)
        stack_sites = [site for (dist, site) in zip(dists, sites) if dist <= max_dist]
        jdac!(grid, stack_sites, jdac_aux2b!, p, depth - 1, rect, stack)
    end
    grid
end

@inbounds function jdac_aux2c!(grid, sites, p, depth, rect, stack)
    (t, l), (N, M) = rect
    if all(>(0), (N, M)) && any(==(0), @view grid[t:t+N-1, l:l+M-1])
        center = (t+(N-1)/2, l+(M-1)/2) 
        push_empty!(stack)
        fill_dists!(stack, center, sites, p)
        dist, index = findmin(get_dists(stack))
        max_dist = dist + distance((1, 1), (N, M), p) + 1
        min_site = sites[index][2]
        corners = (t, l), (t + N - 1, l), (t + N - 1, l + M - 1), (t, l + M - 1)
        dist = maximum(corner -> distance(corner, min_site, p), corners) + distance((1, 1), center, p) + 1
        max_dist = min(max_dist, dist)
        fill_sites!(stack, max_dist, sites)
        jdac!(grid, get_sites(stack), jdac_aux2c!, p, depth - 1, rect, stack)
        pop!(stack)
    end
    grid
end

@inbounds function jdac_aux3a!(grid, sites, p, depth, rect, stack)
    (t, l), (N, M) = rect
    if all(>(0), (N, M)) && any(==(0), @view grid[t:t+N-1, l:l+M-1])
        center = (t+(N-1)/2, l+(M-1)/2) 
        min_site = sites[findmin(site -> distance(center, site[2], p), sites)[2]][2]
        corners = (t, l), (t + N - 1, l), (t + N - 1, l + M - 1), (t, l + M - 1)
        stack_sites = filter(site -> any(corner -> distance(corner, site[2], p) <= distance(corner, min_site, p), corners), sites)
        jdac!(grid, stack_sites, jdac_aux3a!, p, depth - 1, rect, stack)
    end
    grid
end

@inbounds function jdac_aux3b!(grid, sites, p, depth, rect, stack)
    (t, l), (N, M) = rect
    if all(>(0), (N, M)) && any(==(0), @view grid[t:t+N-1, l:l+M-1])
        center = (t+(N-1)/2, l+(M-1)/2) 
        min_site = sites[findmin(site -> distance(center, site[2], p), sites)[2]][2]
        corners = (t, l), (t + N - 1, l), (t + N - 1, l + M - 1), (t, l + M - 1)
        dists = map(corner -> distance(corner, min_site, p), corners)
        stack_sites = filter(site -> any(corner -> distance(corner[2], site[2], p) <= dists[corner[1]], enumerate(corners)), sites)
        jdac!(grid, stack_sites, jdac_aux3b!, p, depth - 1, rect, stack)
    end
    grid
end

@inbounds function jdac_aux3c!(grid, sites, p, depth, rect, stack)
    (t, l), (N, M) = rect
    if all(>(0), (N, M)) && any(==(0), @view grid[t:t+N-1, l:l+M-1])
        center = (t+(N-1)/2, l+(M-1)/2) 
        push_empty!(stack)
        min_site = sites[findmin(site -> distance(center, site[2], p), sites)[2]][2]
        corners = (t, l), (t + N - 1, l), (t + N - 1, l + M - 1), (t, l + M - 1)
        fill_dists!(stack, min_site, enumerate(corners), p)
        fill_sites!(stack, corners, sites, p)
        jdac!(grid, get_sites(stack), jdac_aux3c!, p, depth - 1, rect, stack)
        pop!(stack)
    end
    grid
end

@inbounds function jdac_aux4a!(grid, sites, p, depth, rect, stack)
    (t, l), (N, M) = rect
    if all(>(0), (N, M)) && any(==(0), @view grid[t:t+N-1, l:l+M-1])
        mean_site = reduce(.+, site[2] for site in sites) ./ length(sites)
        anchor = closest(rect, mean_site)
        min_site = sites[findmin(site -> distance(anchor, site[2], p), sites)[2]][2]
        corners = (t, l), (t + N - 1, l), (t + N - 1, l + M - 1), (t, l + M - 1)
        stack_sites = filter(site -> any(corner -> distance(corner, site[2], p) <= distance(corner, min_site, p), corners), sites)
        jdac!(grid, stack_sites, jdac_aux4a!, p, depth - 1, rect, stack)
    end
    grid
end

@inbounds function jdac_aux4c!(grid, sites, p, depth, rect, stack)
    (t, l), (N, M) = rect
    if all(>(0), (N, M)) && any(==(0), @view grid[t:t+N-1, l:l+M-1])
        mean_site = get_mean_site(stack)
        push_empty!(stack)
        anchor = closest(rect, mean_site)
        min_site = sites[findmin(site -> distance(anchor, site[2], p), sites)[2]][2]
        corners = (t, l), (t + N - 1, l), (t + N - 1, l + M - 1), (t, l + M - 1)
        fill_dists!(stack, min_site, enumerate(corners), p)
        fill_sites!(stack, corners, sites, p)
        jdac!(grid, get_sites(stack), jdac_aux4c!, p, depth - 1, rect, stack)
        pop!(stack)
    end
    grid
end

@inbounds function jdac_aux5a!(grid, sites, p, depth, rect, stack)
    (t, l), (N, M) = rect
    if all(>(0), (N, M)) && any(==(0), @view grid[t:t+N-1, l:l+M-1])
        sites_geometric_mean = reduce(.+, site[2] for site in sites) ./ length(sites)
        inner_anchor = closest_anchor_to_rectangle(rect, sites_geometric_mean)
        min_site = sites[findmin(site -> distance(inner_anchor, site[2]), sites)[2]][2]
        corners = (t, l), (t + N - 1, l), (t, l + M - 1), (t + N - 1, l + M - 1)
        stack_sites = filter(site -> any(corner -> distance(corner, site[2], p) <= distance(corner, min_site, p), corners), sites)
        jdac!(grid, stack_sites, jdac_aux5a!, p, depth - 1, rect, stack)
    end
    grid
end

@inbounds function jdac_aux5c!(grid, sites, p, depth, rect, stack)
    (t, l), (N, M) = rect
    if all(>(0), (N, M)) && any(==(0), @view grid[t:t+N-1, l:l+M-1])
        sites_geometric_mean = get_mean_site(stack)
        push_empty!(stack)
        inner_anchor = closest_anchor_to_rectangle(rect, sites_geometric_mean)
        min_site = sites[findmin(site -> distance(inner_anchor, site[2], p), sites)[2]][2]
        corners = (t, l), (t + N - 1, l), (t + N - 1, l + M - 1), (t, l + M - 1)
        fill_dists!(stack, min_site, enumerate(corners), p)
        fill_sites!(stack, corners, sites, p)
        jdac!(grid, get_sites(stack), jdac_aux4c!, p, depth - 1, rect, stack)
        pop!(stack)
    end
    grid
end

@inbounds function jdac!(grid, sites, aux!, p::Real = 2, depth::Int = 1, rect = ((1, 1), size(grid)), stack = SiteStack{Float64,eltype(sites)}(sites))
    (t, l), (N, M) = rect
    if (N == 1 && M == 1) || length(sites) == 1
        center = (t+(N-1)/2, l+(M-1)/2) 
        _, min_index = findmin(site -> distance(center, site[2], p), sites)
        (@view grid[t:t+N-1, l:l+M-1]) .= convert(eltype(grid), sites[min_index][1])
    else
        Nd = N Ć· 2 + N % 2
        Md = M Ć· 2 + M % 2
        tl_rect = (t, l), (Nd, Md)
        bl_rect = (t + Nd, l), (N - Nd, Md)
        tr_rect = (t, l + Md), (Nd, M - Md)
        br_rect = (t + Nd, l + Md), (N - Nd, M - Md)
        sub_rects = (tl_rect, bl_rect, tr_rect, br_rect)
        if depth > 0
            Threads.@threads for i in 1:4
                aux!(grid, sites, p, depth, sub_rects[i], SiteStack{Float64,eltype(sites)}(sites))
            end
        else
            for i in 1:4
                aux!(grid, sites, p, depth, sub_rects[i], stack)
            end
        end
    end
    grid
end

function jdacx!(grid, sites, aux!, p::Real = 2, depth::Int = 1)
    for (color, point) in sites
        closest_point = round_tuple(point)
        if 1 <= closest_point[1] <= size(grid, 1) && 1 <= closest_point[2] <= size(grid, 2)
            grid[closest_point...] = convert(eltype(grid), color)
        end
    end
    jdac!(grid, sites, aux!, p, depth)
end
