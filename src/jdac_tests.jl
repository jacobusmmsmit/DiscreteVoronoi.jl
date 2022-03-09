export jdac!, jdacx!, jdac_aux1!, jdac_aux2!, jdac_aux3!, jdac_aux4!, jdac_aux5!, jdac_aux6!

function jdac_aux1!(grid, sites, depth, stack)
    if all(.>(0), size(grid)) && any(==(0), grid)
        center = size(grid) .÷ 2 .+ size(grid) .% 2
        max_dist = findmin(site -> distance(center, site[2]), sites)[1] + distance((0, 0), size(grid)) + 1
        stack_sites = filter(site -> distance(center, site[2]) <= max_dist, sites)
        jdac!(grid, stack_sites, jdac_aux1!, depth - 1, stack)
    end
    grid
end

function jdac_aux2!(grid, sites, depth, stack)
    if all(.>(0), size(grid)) && any(==(0), grid)
        center = size(grid) .÷ 2 .+ size(grid) .% 2
        dists = map(site -> distance(center, site[2]), sites)
        max_dist = findmin(dists)[1] + distance((0, 0), size(grid)) + 1
        stack_sites = [site for (dist, site) in zip(dists, sites) if dist <= max_dist]
        jdac!(grid, stack_sites, jdac_aux2!, depth - 1, stack)
    end
    grid
end

function jdac_aux3!(grid, sites, depth, stack)
    if all(.>(0), size(grid)) && any(==(0), grid)
        center = size(grid) .÷ 2 .+ size(grid) .% 2
        push_empty!(stack)
        fill_dists!(center, sites, stack)
        max_dist = findmin(get_dists(stack))[1] + distance((0, 0), size(grid)) + 1
        fill_sites!(max_dist, sites, stack)
        jdac!(grid, get_sites(stack), jdac_aux3!, depth - 1, stack)
        pop!(stack)
    end
    grid
end

function jdac_aux4!(grid, sites, depth, stack)
    if all(.>(0), size(grid)) && any(==(0), grid)
        N, M = size(grid)
        center = size(grid) .÷ 2 .+ size(grid) .% 2
        min_site = sites[findmin(site -> distance(center, site[2]), sites)[2]][2]
        corners = ((1, 1), (N, 1), (1, M), (N, M))
        max_dist = maximum(corner -> distance(corner, min_site), corners) + distance((0, 0), center) + 1
        stack_sites = filter(site -> distance(center, site[2]) <= max_dist, sites)
        jdac!(grid, stack_sites, jdac_aux4!, depth - 1, stack)
    end
    grid
end

function jdac_aux5!(grid, sites, depth, stack)
    if all(.>(0), size(grid)) && any(==(0), grid)
        N, M = size(grid)
        center = size(grid) .÷ 2 .+ size(grid) .% 2
        dists = map(site -> distance(center, site[2]), sites)
        min_site = sites[findmin(dists)[2]][2]
        corners = ((1, 1), (N, 1), (1, M), (N, M))
        max_dist = maximum(corner -> distance(corner, min_site), corners) + distance((0, 0), center) + 1
        stack_sites = [site for (dist, site) in zip(dists, sites) if dist <= max_dist]
        jdac!(grid, stack_sites, jdac_aux5!, depth - 1, stack)
    end
    grid
end

function jdac_aux6!(grid, sites, depth, stack)
    if all(.>(0), size(grid)) && any(==(0), grid)
        N, M = size(grid)
        center = size(grid) .÷ 2 .+ size(grid) .% 2
        push_empty!(stack)
        fill_dists!(center, sites, stack)
        min_site = sites[findmin(get_dists(stack))[2]]
        corners = ((1, 1), (N, 1), (1, M), (N, M))
        max_dist = maximum(corner -> distance(corner, min_site[2]), corners) + distance((0, 0), center) + 1
        fill_sites!(max_dist, sites, stack)
        jdac!(grid, get_sites(stack), jdac_aux6!, depth - 1, stack)
        pop!(stack)
    end
    grid
end

function jdac!(grid, sites, aux!, depth=1, stack=SiteStack())
    N, M = size(grid)
    corners = ((1, 1), (N, 1), (1, M), (N, M))
    nearest_colors = map(corner -> findmin(site -> distance(corner, site[2]), sites)[2], corners)
    color = sites[nearest_colors[1]][1]
    if all(nearest_color -> sites[nearest_color][1] == color, nearest_colors[2:end])
        grid .= convert(eltype(grid), color)
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
                aux!(sub_grids[i], thread_sites, depth, SiteStack())
            end
        else
            for i in 1:4
                replace!(site -> (site[1], site[2] .- offsets[i]), sites)
                aux!(sub_grids[i], sites, depth, stack)
                replace!(site -> (site[1], site[2] .+ offsets[i]), sites)
            end
        end
    end
    grid
end

function jdacx!(grid, sites, aux!, depth=1, stack=SiteStack())
    for (color, (x, y)) in sites
        grid[x, y] = convert(eltype(grid), color)
    end
    jdac!(grid, sites, aux!, depth, stack)
end