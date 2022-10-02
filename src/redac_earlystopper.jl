using LinearAlgebra: norm
using StaticArrays
using BenchmarkTools

include("helper_functions.jl")
include("elimination_methods.jl")
include("traditional_approaches.jl")

# A key part of the implementation is the monadic (?) EarlyStopper type which
# wraps information about how far an iterator should go in an object.
include("EarlyStopper.jl")

# Reduce, Divide, and Conquer
function early_stop_sort!(early_stopper::ES, predicate) where {ES<:EarlyStopper}
    N = length(early_stopper)
    mask = SVector{N,Bool}([predicate(el) for el in early_stopper]) # Apply predicate
    n_trues = count(mask) # Get the index of the last "true" after sorting
    perm = sortperm(mask, rev=true) # Return permutation that sorts the mask
    permute!(view(early_stopper.obj, 1:N), perm) # Apply permutation to sites
    return EarlyStopper(early_stopper.obj, n_trues)
end

early_stop_sort!(non_early_stopper, predicate) = early_stop_sort!(EarlyStopper(non_early_stopper), predicate)

function site_sort!(sites::ES, keep_predicate, TL, BR) where {ES<:EarlyStopper}
    predicate(site) = keep_predicate(site, sites, TL, BR)
    return early_stop_sort!(sites, predicate)
end

site_sort!(sites, keep_predicate) = site_sort!(EarlyStopper(sites), keep_predicate)

function redac_voronoi!(grid, sites; p=2, predicate)
    TL = 1, 1
    BR = size(grid)
    _redac_voronoi!(grid, TL, BR, EarlyStopper(sites), p, predicate)
    return nothing
end

function _redac_voronoi!(grid, TL, BR, sites::ES, p, predicate) where {ES<:EarlyStopper}
    # Again, first, if the grid is a single cell, we are done
    if all((BR .- TL) .== 0)
        grid[TL...] = find_closest_site(TL, sites, p)
    else
        # Otherwise we check if all corners have the same closest site
        corners = (TL, (TL[1], BR[2]), (BR[1], TL[2]), BR)
        closest_corners = (find_closest_site(corner, sites, p) for corner in corners)
        if allequal(closest_corners)
            view(grid, TL[1]:BR[1], TL[2]:BR[2]) .= Ref(first(closest_corners))
        else
            # And if not we eliminate faraway seeds from subsequent steps
            # `site_sort` sorts sites by whether the predicate is true and stores how many are true.
            local_sites = site_sort!(sites, predicate, TL, BR)
            # then divide the grid into quadrants and "conquer" each one
            for (quadrant_TL, quadrant_BR) in get_quadrants(TL, BR)
                _redac_voronoi!(grid, quadrant_TL, quadrant_BR, local_sites, p, predicate)
            end
        end
    end
    return nothing
end

### Example ###
n = 50
l = 3
TL = (1, 1)
BR = (n, n)
locs = sort([SVector{2}(rand(1:n, 2)) for _ in 1:l])
# sites = [Site(loc, lab) for (lab, loc) in enumerate(locs)]
grid = zeros(Int, (n, n))
loc_grid = zeros(SVector{2,Int}, (n, n))

begin
    dac_time = Float64[]
    redac_time = Float64[]
    redac_anchor_time = Float64[]
    naive_time = Float64[]
    ls = 1:100
    for l in ls
        locs = sort([SVector{2}(rand(1:n, 2)) for _ in 1:l])
        for (alg, timer) in zip((dac_voronoi!, (x, y) -> redac_voronoi!(x, y, predicate=exact_condition), (x, y) -> redac_voronoi!(x, y, predicate=centre_anchor_condition), naive_voronoi!), (dac_time, redac_time, redac_anchor_time, naive_time))
            push!(timer, @elapsed(alg(loc_grid, locs)))
        end
    end
end
using Plots
plot(xlabel="Number of sites", ylabel="Time (s)", legendposition=:topleft)
plot!(ls, dac_time, label="DAC")
plot!(ls[2:end], redac_time[2:end], label="REDAC + Exact condition")
plot!(ls[2:end], redac_anchor_time[2:end], label="REDAC + Centre-Anchor condition")
plot!(ls, naive_time, label="Naive")

@belapsed redac_voronoi!($loc_grid, $locs, predicate=exact_condition)
@belapsed redac_voronoi!($loc_grid, $locs, predicate=centre_anchor_condition)
@belapsed dac_voronoi!($loc_grid, $locs)
@belapsed naive_voronoi!($loc_grid, $locs)
@belapsed jfa_voronoi!($loc_grid, $locs)