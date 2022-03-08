export jdac!, jdacx!, jdac_aux1!, jdac_aux2!, jdac_aux3!, jdac_aux4!, jdac_aux5!, jdac_aux6!

function jdac_aux1!(grid, seeds, depth, stack)
    if all(.>(0), size(grid)) && any(==(0), grid)
        center = size(grid) .÷ 2 .+ size(grid) .% 2
        max_dist = findmin(seed -> distance(center, seed[2]), seeds)[1] + distance((0, 0), size(grid)) + 1
        stack_seeds = filter(seed -> distance(center, seed[2]) <= max_dist, seeds)
        jdac!(grid, stack_seeds, jdac_aux1!, depth - 1, stack)
    end
    grid
end

function jdac_aux2!(grid, seeds, depth, stack)
    if all(.>(0), size(grid)) && any(==(0), grid)
        center = size(grid) .÷ 2 .+ size(grid) .% 2
        dists = map(seed -> distance(center, seed[2]), seeds)
        max_dist = findmin(dists)[1] + distance((0, 0), size(grid)) + 1
        stack_seeds = [seed for (dist, seed) in zip(dists, seeds) if dist <= max_dist]
        jdac!(grid, stack_seeds, jdac_aux2!, depth - 1, stack)
    end
    grid
end

function jdac_aux3!(grid, seeds, depth, stack)
    if all(.>(0), size(grid)) && any(==(0), grid)
        center = size(grid) .÷ 2 .+ size(grid) .% 2
        push_empty!(stack)
        fill_dists!(center, seeds, stack)
        max_dist = findmin(get_dists(stack))[1] + distance((0, 0), size(grid)) + 1
        fill_seeds!(max_dist, seeds, stack)
        jdac!(grid, get_seeds(stack), jdac_aux3!, depth - 1, stack)
        pop!(stack)
    end
    grid
end

function jdac_aux4!(grid, seeds, depth, stack)
    if all(.>(0), size(grid)) && any(==(0), grid)
        N, M = size(grid)
        center = size(grid) .÷ 2 .+ size(grid) .% 2
        min_seed = seeds[findmin(seed -> distance(center, seed[2]), seeds)[2]][2]
        corners = ((1, 1), (N, 1), (1, M), (N, M))
        max_dist = maximum(corner -> distance(corner, min_seed), corners) + distance((0, 0), center) + 1
        stack_seeds = filter(seed -> distance(center, seed[2]) <= max_dist, seeds)
        jdac!(grid, stack_seeds, jdac_aux4!, depth - 1, stack)
    end
    grid
end

function jdac_aux5!(grid, seeds, depth, stack)
    if all(.>(0), size(grid)) && any(==(0), grid)
        N, M = size(grid)
        center = size(grid) .÷ 2 .+ size(grid) .% 2
        dists = map(seed -> distance(center, seed[2]), seeds)
        min_seed = seeds[findmin(dists)[2]][2]
        corners = ((1, 1), (N, 1), (1, M), (N, M))
        max_dist = maximum(corner -> distance(corner, min_seed), corners) + distance((0, 0), center) + 1
        stack_seeds = [seed for (dist, seed) in zip(dists, seeds) if dist <= max_dist]
        jdac!(grid, stack_seeds, jdac_aux5!, depth - 1, stack)
    end
    grid
end

function jdac_aux6!(grid, seeds, depth, stack)
    if all(.>(0), size(grid)) && any(==(0), grid)
        N, M = size(grid)
        center = size(grid) .÷ 2 .+ size(grid) .% 2
        push_empty!(stack)
        fill_dists!(center, seeds, stack)
        min_seed = seeds[findmin(get_dists(stack))[2]]
        corners = ((1, 1), (N, 1), (1, M), (N, M))
        max_dist = maximum(corner -> distance(corner, min_seed[2]), corners) + distance((0, 0), center) + 1
        fill_seeds!(max_dist, seeds, stack)
        jdac!(grid, get_seeds(stack), jdac_aux6!, depth - 1, stack)
        pop!(stack)
    end
    grid
end

function jdac!(grid, seeds, aux!, depth=1, stack=SiteStack())
    N, M = size(grid)
    corners = ((1, 1), (N, 1), (1, M), (N, M))
    nearest_colors = map(corner -> findmin(seed -> distance(corner, seed[2]), seeds)[2], corners)
    color = seeds[nearest_colors[1]][1]
    if all(nearest_color -> seeds[nearest_color][1] == color, nearest_colors[2:end])
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
                thread_seeds = map(seed -> (seed[1], seed[2] .- offsets[i]), seeds)
                aux!(sub_grids[i], thread_seeds, depth, SiteStack())
            end
        else
            for i in 1:4
                replace!(seed -> (seed[1], seed[2] .- offsets[i]), seeds)
                aux!(sub_grids[i], seeds, depth, stack)
                replace!(seed -> (seed[1], seed[2] .+ offsets[i]), seeds)
            end
        end
    end
    grid
end

function jdacx!(grid, seeds, aux!, depth=1, stack=SiteStack())
    for (color, (x, y)) in seeds
        grid[x, y] = convert(eltype(grid), color)
    end
    jdac!(grid, seeds, aux!, depth, stack)
end