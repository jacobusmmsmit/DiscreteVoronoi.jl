export jdac!, jdacx!, jdac_aux0!, jdac_aux1a!, jdac_aux1b!, jdac_aux1c!, jdac_aux2a!, jdac_aux2b!, jdac_aux2c!, jdac_aux3!


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

function jdac_aux0!(grid, sites, p, depth, stack)
    if all(.>(0), size(grid)) && any(==(0), grid)
        N, M = size(grid)
        corners = ((1, 1), (N, 1), (1, M), (N, M))
        stack_sites = exact_site_filter(sites, corners, p)
        jdac!(grid, stack_sites, jdac_aux1a!, p, depth - 1, stack)
    end
    grid
end

function jdac_aux1a!(grid, sites, p, depth, stack)
    if all(.>(0), size(grid)) && any(==(0), grid)
        center = size(grid) ./ 2
        max_dist = findmin(site -> distance(center, site[2], p), sites)[1] + norm(size(grid), p) + 1
        stack_sites = filter(site -> distance(center, site[2], p) <= max_dist, sites)
        jdac!(grid, stack_sites, jdac_aux1a!, p, depth - 1, stack)
    end
    grid
end

function jdac_aux1b!(grid, sites, p, depth, stack)
    if all(.>(0), size(grid)) && any(==(0), grid)
        center = size(grid) ./ 2
        dists = map(site -> distance(center, site[2], p), sites)
        max_dist = findmin(dists)[1] + norm(size(grid), p) + 1
        stack_sites = [site for (dist, site) in zip(dists, sites) if dist <= max_dist]
        jdac!(grid, stack_sites, jdac_aux1b!, p, depth - 1, stack)
    end
    grid
end

function jdac_aux1c!(grid, sites, p, depth, stack)
    if all(.>(0), size(grid)) && any(==(0), grid)
        center = size(grid) ./ 2
        push_empty!(stack)
        fill_dists!(stack, center, sites, p)
        max_dist = findmin(get_dists(stack))[1] + norm(size(grid), p) + 1
        fill_sites!(stack, max_dist, sites)
        jdac!(grid, get_sites(stack), jdac_aux1c!, p, depth - 1, stack)
        pop!(stack)
    end
    grid
end

function jdac_aux2a!(grid, sites, p, depth, stack)
    if all(.>(0), size(grid)) && any(==(0), grid)
        center = size(grid) ./ 2
        dist, index = findmin(site -> distance(center, site[2], p), sites)
        max_dist = dist + norm(size(grid), p) + 1
        min_site = sites[index][2]
        N, M = size(grid)
        corners = ((1, 1), (N, 1), (1, M), (N, M))
        dist = maximum(corner -> distance(corner, min_site, p), corners) + norm(center, p) + 1
        max_dist = min(max_dist, dist)
        stack_sites = filter(site -> distance(center, site[2], p) <= max_dist, sites)
        jdac!(grid, stack_sites, jdac_aux2a!, p, depth - 1, stack)
    end
    grid
end

function jdac_aux2b!(grid, sites, p, depth, stack)
    if all(.>(0), size(grid)) && any(==(0), grid)
        center = size(grid) ./ 2
        dists = map(site -> distance(center, site[2], p), sites)
        dist, index = findmin(dists)
        max_dist = dist + norm(size(grid), p) + 1
        min_site = sites[index][2]
        N, M = size(grid)
        corners = ((1, 1), (N, 1), (1, M), (N, M))
        dist = maximum(corner -> distance(corner, min_site, p), corners) + norm(center, p) + 1
        max_dist = min(max_dist, dist)
        stack_sites = [site for (dist, site) in zip(dists, sites) if dist <= max_dist]
        jdac!(grid, stack_sites, jdac_aux2b!, p, depth - 1, stack)
    end
    grid
end

function jdac_aux2c!(grid, sites, p, depth, stack)
    if all(.>(0), size(grid)) && any(==(0), grid)
        push_empty!(stack)
        center = size(grid) ./ 2
        fill_dists!(stack, center, sites, p)
        dist, index = findmin(get_dists(stack))
        max_dist = dist + norm(size(grid), p) + 1
        min_site = sites[index][2]
        N, M = size(grid)
        corners = ((1, 1), (N, 1), (1, M), (N, M))
        dist = maximum(corner -> distance(corner, min_site, p), corners) + norm(center, p) + 1
        max_dist = min(max_dist, dist)
        fill_sites!(stack, max_dist, sites)
        jdac!(grid, get_sites(stack), jdac_aux2c!, p, depth - 1, stack)
        pop!(stack)
    end
    grid
end

function jdac_aux3!(grid, sites, p, depth, stack)
    if all(.>(0), size(grid)) && any(==(0), grid)
        center = size(grid) ./ 2
        min_site = sites[findmin(site -> distance(center, site[2], p), sites)[2]][2]
        N, M = size(grid)
        corners = ((1, 1), (N, 1), (1, M), (N, M))
        stack_sites = filter(site -> mapreduce(corner -> distance(corner, site[2], p) <= distance(corner, min_site, p), (x, y) -> x || y, corners), sites)
        jdac!(grid, stack_sites, jdac_aux3!, p, depth - 1, stack)
    end
    grid
end

function jdac!(grid, sites, aux!, p::Real = 2, depth::Int = 1, stack = SiteStack{Float64,eltype(sites)}())
    N, M = size(grid)
    if (N == 1 && M == 1) || length(sites) == 1
        min_dist, min_index = findmin(site -> distance(size(grid), site[2], p), sites)
        grid .= convert(eltype(grid), sites[min_index][1])
    else
        Nd = N ÷ 2 + N % 2
        Md = M ÷ 2 + M % 2
        offsets = ((0, 0), (Nd, 0), (0, Md), (Nd, Md))
        tl_grid = @view grid[1:Nd, 1:Md]
        bl_grid = @view grid[Nd+1:N, 1:Md]
        tr_grid = @view grid[1:Nd, Md+1:M]
        br_grid = @view grid[Nd+1:N, Md+1:M]
        sub_grids = (tl_grid, bl_grid, tr_grid, br_grid)
        if depth > 0
            Threads.@threads for i in 1:4
                thread_sites = map(site -> (site[1], site[2] .- offsets[i]), sites)
                aux!(sub_grids[i], thread_sites, p, depth, SiteStack{Float64,eltype(sites)}())
            end
        else
            for i in 1:4
                replace!(site -> (site[1], site[2] .- offsets[i]), sites)
                aux!(sub_grids[i], sites, p, depth, stack)
                replace!(site -> (site[1], site[2] .+ offsets[i]), sites)
            end
        end
    end
    grid
end

function jdacx!(grid, sites, aux!, p::Real = 2, depth::Int = 1, stack = SiteStack{Float64,eltype(sites)}())
    for (color, point) in sites
        nearest_point = nearest(point)
        if 1 <= nearest_point[1] <= size(grid, 1) && 1 <= nearest_point[2] <= size(grid, 2)
            grid[nearest_point...] = convert(eltype(grid), color)
        end
    end
    jdac!(grid, sites, aux!, p, depth, stack)
end