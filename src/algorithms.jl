export naive_voronoi, jfa!, dac_aux!, dac!, dacx!

naive_voronoi(grid, sites) = map(cell -> findmin(site -> distance(cell, site), sites)[2], grid)

function jfa!(grid, sites)
    for (color, (x, y)) in enumerate(sites)
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
                elseif distance(sites[colorp], (x, y)) > distance(sites[colorq], (x, y))
                    grid[x, y] = colorq
                end
            end
        end
    end
    return grid
end

function dac_aux!(grid, sites, depth)
    if all(.>(0), size(grid)) && any(==(0), grid)
        dac!(grid, sites, depth - 1)
    end
    return grid
end

function dac!(grid, sites, depth=1)
    N, M = size(grid)
    corners = ((1, 1), (N, 1), (1, M), (N, M))
    nearest_colors = map(corner -> findmin(site -> distance(corner, site), sites)[2], corners)
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
                thread_sites = map(site -> site .- offsets[i], sites)
                dac_aux!(sub_grids[i], thread_sites, depth)
            end
        else
            for i in 1:4
                replace!(site -> site .- offsets[i], sites)
                dac_aux!(sub_grids[i], sites, depth)
                replace!(site -> site .+ offsets[i], sites)
            end
        end
    end
    return grid
end

function dacx!(grid, sites, depth=1)
    for (color, (x, y)) in enumerate(sites)
        grid[x, y] = convert(eltype(grid), color)
    end
    dac!(grid, sites, depth)
end