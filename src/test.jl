using StaticArrays
using BenchmarkTools

include("helper_functions.jl")

struct Site
    loc::SVector{2, Int}
    label::Int
end

isless(s::Site, t::Site) = isless(s.loc, t.loc)

function site_in_rect(s::Site, TL, BR)
    return all(TL .<= s.loc .<= BR)
end

function site_in_rect(v::V, TL, BR) where {V<:AbstractVector}
    return all(TL .<= v .<= BR)
end

function in_TL_quadrant(s, TL, BR)
    TL_quadrant = (TL, (TL .+ BR) .÷ 2)
    return site_in_rect(s, TL_quadrant...)
end

# function main()
n = 500
l = 50
TL = (1, 1)
BR = (n, n)
locs = sort([SVector{2}(rand(1:n, 2)) for _ in 1:l])
sites = [Site(loc, lab) for (lab, loc) in enumerate(locs)]
grid = zeros(Int, (n, n))
@benchmark predicate_sort!($sites, x -> in_TL_quadrant(x, $TL, $BR))
@benchmark predicate_sort!($locs, x -> in_TL_quadrant(x, $TL, $BR))

mycount = predicate_sort!(sites, x -> in_TL_quadrant(x, TL, BR))
mycount = predicate_sort!(locs, x -> in_TL_quadrant(x, TL, BR))

locs
sites

for (i, loc) in enumerate(locs)
    i > mycount && break
    println(loc)
end

function site_sort!(indices, sites, predicate; stop=length(sites))
    v = @views sites[indices] # Get the subset of v that we're working with
    mask = SVector{length(v), Bool}([predicate(el) for el in v]) # Apply predicate
    count = sum(mask) # Get the index of the last "true" after sorting
    perm = sortperm(mask, rev=true) # Return permutation that sorts the mask
    permute!(v, perm) # Apply permutation to our subset
    return count
end

mypred(x) = in_TL_quadrant(x, TL, BR)

@btime site_sort!(1:length($locs), $locs, mypred)
# end

main()