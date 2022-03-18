using Polyester

export naive_voronoi, jfa!, jfa_par!, dac_aux!, dac!, dacx!

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

function jfa_par!(grid, sites, p::Real=2)
    for (color, (x, y)) in enumerate(sites)
        grid[x, y] = convert(eltype(grid), color)
    end
    gridp = copy(grid)
    k = max(size(grid)...)
    while k > 1
        k = k ÷ 2 + k % 2
        @inbounds @batch for y in 1:size(grid, 2), x in 1:size(grid, 1)
            for j in (-k, 0, k), i in (-k, 0, k)
                checkbounds(Bool, grid, x + i, y + j) || continue
                colorq = grid[x + i, y + j]
                colorq !== 0 || continue
                colorp = grid[x, y]
                if colorp === 0
                    gridp[x, y] = colorq
                elseif distance(sites[colorp], (x, y), p) > distance(sites[colorq], (x, y), p)
                    gridp[x, y] = colorq
                end
            end
        end
        grid .= gridp
    end
    return grid
end

@inbounds function dac_aux!(grid, sites, p, depth, rect)
    (t, l), (N, M) = rect
    if all(>(0), (N, M)) && any(==(0), @view grid[t:t+N-1, l:l+M-1])
        dac!(grid, sites, p, depth - 1, rect)
    end
    return grid
end

@inbounds function dac!(grid, sites, p::Real = 2, depth::Int = 1, rect = ((1,1), size(grid)))
    (t, l), (N, M) = rect
    center = (t+(N-1)/2, l+(M-1)/2) 
    if (N == 1 && M == 1) || length(sites) == 1
        _, min_index = findmin(site -> distance(center, site, p), sites)
        (@view grid[t:t+N-1, l:l+M-1]) .= convert(eltype(grid), min_index)
    else
        min_dist, min_index = findmin(site -> distance(center, site, p), sites)
        dist, _ = findmin(site -> distance(center, site, p), @view sites[1:end .!= min_index])
        if dist > min_dist + distance((1, 1), (N, M), p) + 1
            (@view grid[t:t+N-1, l:l+M-1]) .= convert(eltype(grid), min_index)
        else
            Nd = N ÷ 2 + N % 2
            Md = M ÷ 2 + M % 2
            tl_rect = (t, l), (Nd, Md)
            bl_rect = (t+Nd, l), (N-Nd, Md)
            tr_rect = (t, l+Md), (Nd, M-Md)
            br_rect = (t+Nd, l+Md), (N-Nd, M-Md)
            sub_rects = (tl_rect, bl_rect, tr_rect, br_rect)
            #= if depth > 0
                Threads.@threads for i in 1:4
                    dac_aux!(grid, sites, p, depth, sub_rects[i])
                end
            else
                for i in 1:4
                    dac_aux!(grid, sites, p, depth, sub_rects[i])
                end
            end =#
            @batch for i in 1:4
                dac_aux!(grid, sites, p, depth, sub_rects[i])
            end
        end
        return grid
    end
end

function dacx!(grid, sites, p::Real=2, depth::Int=1)
    for (color, point) in enumerate(sites)
        closest_point = round_tuple(point)
        if 1 <= closest_point[1] <= size(grid, 1) && 1 <= closest_point[2] <= size(grid, 2)
            grid[closest_point...] = convert(eltype(grid), color)
        end
    end
    dac!(grid, sites, p, depth)
end
