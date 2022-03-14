export jdac!, jdacx!, jdac_aux0!, jdac_aux1a!, jdac_aux1b!, jdac_aux1c!, jdac_aux2a!, jdac_aux2b!, jdac_aux2c!, jdac_aux3a!, jdac_aux3b!, jdac_aux3c!


"""
    exact_site_filter(sites, corners, p)
Returns {s ∈ 𝒮 : ∄ s' ∈ 𝒮 s.t. ∀ C ∈ {C₁, C₂, C₃, C₄}, d(s', C) < d(s, C)}
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
        corners = (t, l), (t+N-1, l), (t, l+M-1), (t+N-1, l+M-1)
        stack_sites = exact_site_filter(sites, corners, p)
        jdac!(grid, stack_sites, jdac_aux1a!, p, depth - 1, rect, stack)
    end
    grid
end

@inbounds function jdac_aux1a!(grid, sites, p, depth, rect, stack)
    (t, l), (N, M) = rect
    if all(>(0), (N, M)) && any(==(0), @view grid[t:t+N-1, l:l+M-1])
        center = (t+N/2, l+M/2) 
        max_dist = findmin(site -> distance(center, site[2], p), sites)[1] + _norm((N, M), p) + 1
        stack_sites = filter(site -> distance(center, site[2], p) <= max_dist, sites)
        jdac!(grid, stack_sites, jdac_aux1a!, p, depth - 1, rect, stack)
    end
    grid
end

@inbounds function jdac_aux1b!(grid, sites, p, depth, rect, stack)
    (t, l), (N, M) = rect
    if all(>(0), (N, M)) && any(==(0), @view grid[t:t+N-1, l:l+M-1])
        center = (t+N/2, l+M/2) 
        dists = map(site -> distance(center, site[2], p), sites)
        max_dist = findmin(dists)[1] + _norm((N, M), p) + 1
        stack_sites = [site for (dist, site) in zip(dists, sites) if dist <= max_dist]
        jdac!(grid, stack_sites, jdac_aux1b!, p, depth - 1, rect, stack)
    end
    grid
end

@inbounds function jdac_aux1c!(grid, sites, p, depth, rect, stack)
    (t, l), (N, M) = rect
    if all(>(0), (N, M)) && any(==(0), @view grid[t:t+N-1, l:l+M-1])
        center = (t+N/2, l+M/2) 
        push_empty!(stack)
        fill_dists!(stack, center, sites, p)
        max_dist = findmin(get_dists(stack))[1] + _norm((N, M), p) + 1
        fill_sites!(stack, max_dist, sites)
        jdac!(grid, get_sites(stack), jdac_aux1c!, p, depth - 1, rect, stack)
        pop!(stack)
    end
    grid
end

@inbounds function jdac_aux2a!(grid, sites, p, depth, rect, stack)
    (t, l), (N, M) = rect
    if all(>(0), (N, M)) && any(==(0), @view grid[t:t+N-1, l:l+M-1])
        center = (t+N/2, l+M/2) 
        dist, index = findmin(site -> distance(center, site[2], p), sites)
        max_dist = dist + _norm((N, M), p) + 1
        min_site = sites[index][2]
        corners = (t, l), (t+N-1, l), (t, l+M-1), (t+N-1, l+M-1)
        dist = maximum(corner -> distance(corner, min_site, p), corners) + _norm(center, p) + 1
        max_dist = min(max_dist, dist)
        stack_sites = filter(site -> distance(center, site[2], p) <= max_dist, sites)
        jdac!(grid, stack_sites, jdac_aux2a!, p, depth - 1, rect, stack)
    end
    grid
end

@inbounds function jdac_aux2b!(grid, sites, p, depth, rect, stack)
    (t, l), (N, M) = rect
    if all(>(0), (N, M)) && any(==(0), @view grid[t:t+N-1, l:l+M-1])
        center = (t+N/2, l+M/2) 
        dists = map(site -> distance(center, site[2], p), sites)
        dist, index = findmin(dists)
        max_dist = dist + _norm((N, M), p) + 1
        min_site = sites[index][2]
        corners = (t, l), (t+N-1, l), (t, l+M-1), (t+N-1, l+M-1)
        dist = maximum(corner -> distance(corner, min_site, p), corners) + _norm(center, p) + 1
        max_dist = min(max_dist, dist)
        stack_sites = [site for (dist, site) in zip(dists, sites) if dist <= max_dist]
        jdac!(grid, stack_sites, jdac_aux2b!, p, depth - 1, rect, stack)
    end
    grid
end

@inbounds function jdac_aux2c!(grid, sites, p, depth, rect, stack)
    (t, l), (N, M) = rect
    if all(>(0), (N, M)) && any(==(0), @view grid[t:t+N-1, l:l+M-1])
        center = (t+N/2, l+M/2) 
        push_empty!(stack)
        fill_dists!(stack, center, sites, p)
        dist, index = findmin(get_dists(stack))
        max_dist = dist + _norm((N, M), p) + 1
        min_site = sites[index][2]
        corners = (t, l), (t+N-1, l), (t, l+M-1), (t+N-1, l+M-1)
        dist = maximum(corner -> distance(corner, min_site, p), corners) + _norm(center, p) + 1
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
        center = (t+N/2, l+M/2) 
        min_site = sites[findmin(site -> distance(center, site[2], p), sites)[2]][2]
        corners = (t, l), (t+N-1, l), (t, l+M-1), (t+N-1, l+M-1)
        stack_sites = filter(site -> any(corner -> distance(corner, site[2], p) <= distance(corner, min_site, p), corners), sites)
        jdac!(grid, stack_sites, jdac_aux3a!, p, depth - 1, rect, stack)
    end
    grid
end

@inbounds function jdac_aux3b!(grid, sites, p, depth, rect, stack)
    (t, l), (N, M) = rect
    if all(>(0), (N, M)) && any(==(0), @view grid[t:t+N-1, l:l+M-1])
        center = (t+N/2, l+M/2) 
        min_site = sites[findmin(site -> distance(center, site[2], p), sites)[2]][2]
        corners = (t, l), (t+N-1, l), (t, l+M-1), (t+N-1, l+M-1)
        dists = map(corner -> distance(corner, min_site, p), corners)
        stack_sites = filter(site -> any(corner -> distance(corner[2], site[2], p) <= dists[corner[1]], enumerate(corners)), sites)
        jdac!(grid, stack_sites, jdac_aux3b!, p, depth - 1, rect, stack)
    end
    grid
end

@inbounds function jdac_aux3c!(grid, sites, p, depth, rect, stack)
    (t, l), (N, M) = rect
    if all(>(0), (N, M)) && any(==(0), @view grid[t:t+N-1, l:l+M-1])
        center = (t+N/2, l+M/2) 
        push_empty!(stack)
        min_site = sites[findmin(site -> distance(center, site[2], p), sites)[2]][2]
        corners = (t, l), (t+N-1, l), (t, l+M-1), (t+N-1, l+M-1)
        fill_dists!(stack, min_site, enumerate(corners), p)
        fill_sites!(stack, min_site, corners, sites, p)
        jdac!(grid, get_sites(stack), jdac_aux3c!, p, depth - 1, rect, stack)
        pop!(stack)
    end
    grid
end

@inbounds function jdac!(grid, sites, aux!, p::Real = 2, depth::Int = 1, rect = ((1,1), size(grid)), stack = SiteStack{Float64,eltype(sites)}())
    (t, l), (N, M) = rect
    if (N == 1 && M == 1) || length(sites) == 1
        _, min_index = findmin(site -> distance((t+N-1, l+M-1), site[2], p), sites)
        (@view grid[t:t+N-1, l:l+M-1]) .= convert(eltype(grid), sites[min_index][1])
    else
        Nd = N ÷ 2 + N % 2
        Md = M ÷ 2 + M % 2
        tl_rect = (t, l), (Nd, Md)
        bl_rect = (t+Nd, l), (N-Nd, Md)
        tr_rect = (t, l+Md), (Nd, M-Md)
        br_rect = (t+Nd, l+Md), (N-Nd, M-Md)
        sub_rects = (tl_rect, bl_rect, tr_rect, br_rect)
        if depth > 0
            Threads.@threads for i in 1:4
                aux!(grid, sites, p, depth, sub_rects[i], SiteStack{Float64,eltype(sites)}())
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
        nearest_point = nearest(point)
        if 1 <= nearest_point[1] <= size(grid, 1) && 1 <= nearest_point[2] <= size(grid, 2)
            grid[nearest_point...] = convert(eltype(grid), color)
        end
    end
    jdac!(grid, sites, aux!, p, depth)
end