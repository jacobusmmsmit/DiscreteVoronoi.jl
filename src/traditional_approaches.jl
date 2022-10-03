"""
    naive_voronoi!(grid, sites::T, p=2) where {T<:Vector{SVector{2,Int}}}

Construct in-place a Voronoi diagram in the most basic way possible: check every cell and every combination.
"""
function naive_voronoi!(grid, sites::T, p=2) where {T<:Vector{SVector{2,Int}}}
    for I in CartesianIndices(grid)
        grid[I] = find_closest_site(Tuple(I), sites, p)
    end
    return nothing
end

"""
    jfa_voronoi!(grid, sites::T, p=2) where {T<:Vector{SVector{2,Int}}}

Construct in-place a Voronoi diagram using the [jump flooding algorithm](https://en.wikipedia.org/wiki/Jump_flooding_algorithm).
"""
function jfa_voronoi!(grid, sites::T, p=2) where {T<:Vector{SVector{2,Int}}}
    for site in sites
        # Splatting (grid[site...] = site) causes allocations?
        x, y = site
        grid[x, y] = site
    end
    k = max(size(grid)...)
    while k > 1
        k = k ÷ 2 + k % 2
        # @inbounds for y in 1:size(grid, 2), x in 1:size(grid, 1)
        @inbounds for I in CartesianIndices(grid)
            x, y = Tuple(I)
            for j in (-k, 0, k), i in (-k, 0, k)
                checkbounds(Bool, grid, x + i, y + j) || continue
                i == j == 0 && continue
                siteq = grid[x+i, y+j]
                siteq !== SVector(0, 0) || continue
                sitep = grid[x, y]
                if sitep == SVector(0, 0)
                    grid[x, y] = siteq
                elseif distance(sitep, (x, y), p) > distance(siteq, (x, y), p)
                    grid[x, y] = siteq
                end
            end
        end
    end
    return nothing
end

"""
    dac_voronoi!(grid, sites::T, p=2) where {T<:Vector{SVector{2,Int}}}

Construct in-place a Voronoi diagram using the [divide-and-conquer algorithm](https://doi.org/10.1109%2Feit48999.2020.9208270).
"""
function dac_voronoi!(grid, sites::T, p=2) where {T<:Vector{SVector{2,Int}}}
    TL = (1, 1)
    BR = size(grid)
    _dac_voronoi!(grid, TL, BR, sites, p)
    return nothing
end

function _dac_voronoi!(grid, TL, BR, sites, p)
    # First, if the grid is a single cell then we are done
    if all((BR .- TL) .== 0)
        grid[TL...] = find_closest_site(TL, sites, p)
    else
        # Otherwise we check if all corners have the same closest site
        corners = get_corners(TL, BR)
        closest_corners = (find_closest_site(corner, sites, p) for corner in corners)
        if allequal(closest_corners)
            grid[TL[1]:BR[1], TL[2]:BR[2]] .= Ref(first(closest_corners))
        else
            # And if not we divide the grid into quadrants and "conquer" each one
            for (quadrant_TL, quadrant_BR) in get_quadrants(TL, BR)
                _dac_voronoi!(grid, quadrant_TL, quadrant_BR, sites, p)
            end
        end
    end
    return nothing
end