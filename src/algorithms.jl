using OffsetArrays

export naive_voronoi, jfa!, dac_aux!, dac!, dacx!

naive_voronoi(grid, sites, p::Real=2) = map(cell -> findmin(site -> distance(cell, site, p), sites)[2], grid)

function jfa!(grid, sites, p::Real=2)
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
                elseif distance(sites[colorp], (x, y), p) > distance(sites[colorq], (x, y), p)
                    grid[x, y] = colorq
                end
            end
        end
    end
    return grid
end

function dac_aux!(grid, sites, p, depth)
    if all(.>(0), size(grid)) && any(==(0), grid)
        dac!(grid, sites, p, depth - 1)
    end
    return grid
end

function dac!(grid, sites, p::Real=2, depth::Int=1)
    N, M = size(grid)
    if (N == 1 && M == 1) || length(sites) == 1
        min_dist, min_index = findmin(site -> distance(size(grid), site, p), sites)
        grid .= convert(eltype(grid), min_index)
    else
        center = size(grid) ./ 2
        min_dist, min_index = findmin(site -> distance(center, site, p), sites)
        dist, _ = findmin(site -> distance(center, site, p), @view sites[1:end .!= min_index])
        if dist > min_dist + norm(size(grid), p) + 1
            grid .= convert(eltype(grid), min_index)
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
                    dac_aux!(sub_grids[i], thread_sites, p, depth)
                end
            else
                for i in 1:4
                    replace!(site -> site .- offsets[i], sites)
                    dac_aux!(sub_grids[i], sites, p, depth)
                    replace!(site -> site .+ offsets[i], sites)
                end
            end
        end
        return grid
    end
end

function dacx!(grid, sites, p::Real=2, depth::Int=1)
    for (color, point) in enumerate(sites)
        nearest_point = nearest(point)
        if 1 <= nearest_point[1] <= size(grid, 1) && 1 <= nearest_point[2] <= size(grid, 2)
            grid[nearest_point...] = convert(eltype(grid), color)
        end
    end
    dac!(grid, sites, p, depth)
end