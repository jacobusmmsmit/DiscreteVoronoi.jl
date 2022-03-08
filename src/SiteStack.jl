#######################################################################################
# %% Fundamental type definitions
#######################################################################################
export SiteStack

mutable struct SiteStack{T}
    depth::Int
    dists::Vector{Vector{T}}
    seeds::Vector{Vector{Tuple{Int, Tuple{Int, Int}}}}
    SiteStack() = new{Float64}(0, Vector{Vector{Float64}}(), Vector{Vector{Tuple{Int, Tuple{Int, Int}}}}())
end

SiteStack(::T) where T = new(0, Vector{Vector{T}}(), Vector{Vector{Tuple{Int, Tuple{Int, Int}}}}())

function push_empty!(stack::SiteStack{T}) where T
    stack.depth += 1
    if stack.depth > length(stack.dists)
        push!(stack.dists, Vector{T}())
    end
    if stack.depth > length(stack.seeds)
        push!(stack.seeds, Vector{Tuple{Int, Tuple{Int, Int}}}())
    end
    stack.depth
end

function append_dist!(stack::SiteStack, len, dist)
    if len <= length(stack.dists[stack.depth])
        @inbounds stack.dists[stack.depth][len] = dist
    else
        @inbounds push!(stack.dists[stack.depth], dist)
    end
end
function resize_dists!(stack::SiteStack, len)
    @inbounds resize!(stack.dists[stack.depth], len)
end

function get_dists(stack::SiteStack)
    @inbounds stack.dists[stack.depth]
end

function append_seed!(stack::SiteStack, len, seed)
    if len <= length(stack.seeds[stack.depth])
        @inbounds stack.seeds[stack.depth][len] = seed
    else
        @inbounds push!(stack.seeds[stack.depth], seed)
    end
end

function resize_seeds!(stack::SiteStack, len)
    @inbounds resize!(stack.seeds[stack.depth], len)
end

function get_seeds(stack::SiteStack)
    @inbounds stack.seeds[stack.depth]
end

function pop!(stack::SiteStack)
    stack.depth -= 1
end

function fill_dists!(center, seeds, stack::SiteStack)
    #= push_empty!(stack)
    len = 0
    for seed in seeds
        len += 1
        append_dist!(stack, len, distance(center, seed[2]))
    end
    resize_dists!(stack, len) =#
    resize_dists!(stack, length(seeds))
    copyto!(get_dists(stack), distance(center, seed[2]) for seed in seeds)
    # get_dists(stack)
end

function fill_seeds!(max_dist, seeds, stack::SiteStack)
    #= len = 0
    for (dist, seed) in zip(get_dists(stack), seeds)
        if dist <= max_dist
            len += 1
            append_seed!(stack, len, seed)
        end
    end
    resize_seeds!(stack, len) =#
    len = sum(1 for dist in get_dists(stack) if dist < max_dist)
    resize_seeds!(stack, len)
    copyto!(get_seeds(stack),
        seed for (dist, seed) in zip(get_dists(stack), seeds) if dist < max_dist)
    # get_seeds(stack)
end