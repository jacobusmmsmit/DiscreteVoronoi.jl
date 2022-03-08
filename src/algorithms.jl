export naive_voronoi, jfa!, dac_aux!, dac!, dacx!

naive_voronoi(grid, seeds) = map(cell -> findmin(seed -> distance(cell, seed), seeds)[2], grid)

function jfa!(grid, seeds)
    for (color, (x, y)) in enumerate(seeds)
        grid[x, y] = convert(eltype(grid), color)
    end
    k = max(size(grid)...)
    while k > 1
        k = k ÷ 2 + k % 2
        @inbounds for y in 1:size(grid, 2), x in 1:size(grid, 1)
            for j in (-k, 0, k), i in (-k, 0, k)
                checkbounds(Bool, grid, x + i, y + j) || continue
                colorq = grid[x + i, y + j]
                colorq !== 0 || continue
                colorp = grid[x, y]
                if colorp === 0
                    grid[x, y] = colorq
                elseif distance(seeds[colorp], (x, y)) > distance(seeds[colorq], (x, y))
                    grid[x, y] = colorq
                end
            end
        end
    end
    return grid
end

function dac_aux!(grid, seeds, depth)
    if all(.>(0), size(grid)) && any(==(0), grid)
        dac!(grid, seeds, depth - 1)
    end
    return grid
end

function dac!(grid, seeds, depth=1)
    N, M = size(grid)
    corners = ((1, 1), (N, 1), (1, M), (N, M))
    nearest_colors = map(corner -> findmin(seed -> distance(corner, seed), seeds)[2], corners)
    color = nearest_colors[1]
    if all(nearest_color -> nearest_color == color, nearest_colors[2:end])
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
                thread_seeds = map(seed -> seed .- offsets[i], seeds)
                dac_aux!(sub_grids[i], thread_seeds, depth)
            end
        else
            for i in 1:4
                replace!(seed -> seed .- offsets[i], seeds)
                dac_aux!(sub_grids[i], seeds, depth)
                replace!(seed -> seed .+ offsets[i], seeds)
            end
        end
    end
    return grid
end

function dacx!(grid, seeds, depth=1)
    for (color, (x, y)) in enumerate(seeds)
        grid[x, y] = convert(eltype(grid), color)
    end
    dac!(grid, seeds, depth)
end