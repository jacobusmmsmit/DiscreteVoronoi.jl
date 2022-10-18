# This is the core functionality for the methods of producing discrete Voronoi
# diagrams. The JFA is most commonly used on a GPU where every cell is run in
# parallel, but it is known to produce (small) errors. The most efficient
# algorithm for single-core performance is redac_voronoi!, but it cannot be
# multithreaded in its current state as it heavily mutates its arguments.

"""
    naive_voronoi!(grid::Matrix{T}, sites::Vector{T}; distance=euclidean) where {T<:SVector{2,Int}}

Construct in-place a Voronoi diagram in the most basic way possible: check every cell and every combination.
"""
function naive_voronoi!(grid::Matrix{T}, sites::Vector{T}; distance=euclidean) where {T<:SVector{2,Int}}
    for I in CartesianIndices(grid)
        @inbounds grid[I] = find_closest_site(Tuple(I), sites; distance=distance)
    end
    return nothing
end

"""
    jfa_voronoi!(grid::Matrix{T}, sites::Vector{T}; distance=euclidean) where {T<:SVector{2,Int}}

Construct in-place a Voronoi diagram using the [jump flooding algorithm](https://en.wikipedia.org/wiki/Jump_flooding_algorithm).
The algorithm assumes that a blank cell in the grid has value `SVector(0, 0)` and that sites are inside the grid.
"""
function jfa_voronoi!(grid::Matrix{T}, sites::Vector{T}; distance=euclidean) where {T<:SVector{2,Int}}
    for site in sites
        # Splatting (grid[site...] = site) causes allocations?
        x, y = site
        grid[x, y] = site
    end
    k = max(size(grid)...)
    while k > 1
        k = k ÷ 2 + k % 2
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
                elseif distance(sitep, (x, y)) > distance(siteq, (x, y))
                    grid[x, y] = siteq
                end
            end
        end
    end
    return nothing
end

"""
    dac_voronoi!(grid::Matrix{T}, sites::Vector{T}; distance=euclidean) where {T<:SVector{2,Int}}

Construct in-place a Voronoi diagram using the [divide-and-conquer algorithm](https://doi.org/10.1109%2Feit48999.2020.9208270).
"""
function dac_voronoi!(grid::Matrix{T}, sites::Vector{T}; distance=euclidean) where {T<:SVector{2,Int}}
    TL = (1, 1)
    BR = size(grid)
    _dac_voronoi!(grid, TL, BR, sites, distance=distance)
    return nothing
end

@inbounds function _dac_voronoi!(grid, TL, BR, sites; distance)
    any(TL .> BR) && return nothing
    if all(BR .== TL)
        grid[TL...] = find_closest_site(TL, sites, distance=distance)
    elseif length(sites) == 1 # Same if there is a single site
        view(grid, TL[1]:BR[1], TL[2]:BR[2]) .= Ref(first(sites))
    else
        # Otherwise we check if all corners have the same closest site
        corners = get_corners(TL, BR)
        closest_corners = (find_closest_site(corner, sites, distance=distance) for corner in corners)
        if allequal(closest_corners)
            grid[TL[1]:BR[1], TL[2]:BR[2]] .= Ref(first(closest_corners))
        else
            # And if not we divide the grid into quadrants and "conquer" each one
            for (quadrant_TL, quadrant_BR) in get_quadrants(TL, BR)
                _dac_voronoi!(grid, quadrant_TL, quadrant_BR, sites; distance=distance)
            end
        end
    end
    return nothing
end

"""
    redac_voronoi!(grid::Matrix{T}, sites::Vector{T}; distance=euclidean, auxiliary=exact_aux) where {T<:SVector{2,Int}}

Performs a divide-and-conquer method similar to `dac_voronoi!` but has an additional site-elimination
step which aims to reduce the work of subsequent steps.
"""
function redac_voronoi!(grid::Matrix{T}, sites::Vector{T}; distance=euclidean, auxiliary=exact_aux) where {T<:SVector{2,Int}}
    TL = 1, 1
    BR = size(grid)
    _redac_voronoi!(grid, TL, BR, EarlyStopper(sites); distance=distance, auxiliary=auxiliary)
    return nothing
end

@inbounds function _redac_voronoi!(grid, TL, BR, sites::ES; distance, auxiliary) where {ES<:EarlyStopper}
    any(TL .> BR) && return nothing
    # Then, if the grid is a single cell then we are done
    if all(BR .== TL)
        grid[TL...] = find_closest_site(TL, sites, distance=distance)
    elseif length(sites) == 1 # Same if there is a single site
        view(grid, TL[1]:BR[1], TL[2]:BR[2]) .= Ref(first(sites))
    else
        # Otherwise we check if all corners have the same closest site
        corners = get_corners(TL, BR)
        closest_corners = (find_closest_site(corner, sites, distance=distance) for corner in corners)
        if allequal(closest_corners)
            view(grid, TL[1]:BR[1], TL[2]:BR[2]) .= Ref(first(closest_corners))
        else
            # And if not we eliminate faraway seeds from subsequent steps
            # `auxiliary` sorts sites by whether the predicate is true and stores how many are true.
            local_sites = auxiliary(sites, TL, BR; distance=distance)
            # then divide the grid into quadrants and "conquer" each one
            for (quadrant_TL, quadrant_BR) in get_quadrants(TL, BR)
                _redac_voronoi!(grid, quadrant_TL, quadrant_BR, local_sites; distance=distance, auxiliary=auxiliary)
            end
        end
    end
    return nothing
end