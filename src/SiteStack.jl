#######################################################################################
# %% Fundamental type definitions
#######################################################################################
export SiteStack

mutable struct SiteStack{DT,ST<:Tuple{T1,Tuple{T2,T2}} where {T1,T2}}
    depth::Int
    dists::Vector{Vector{DT}}
    sites::Vector{Vector{ST}}
    SiteStack{DT,ST}() where {DT,ST} = new{DT,ST}(0, Vector{Vector{DT}}(), Vector{Vector{ST}}())
end

SiteStack(::DT, ::ST) where {DT,ST} = new(0, Vector{Vector{DT}}(), Vector{Vector{ST}}())

disttype(::SiteStack{DT,ST}) where {DT,ST} = DT
sitetype(::SiteStack{DT,ST}) where {DT,ST} = ST

function push_empty!(stack::SiteStack{DT,ST}) where {DT,ST}
    stack.depth += 1
    if stack.depth > length(stack.dists)
        push!(stack.dists, Vector{DT}())
    end
    if stack.depth > length(stack.sites)
        push!(stack.sites, Vector{ST}())
    end
    stack.depth
end

@inbounds function append_dist!(stack::SiteStack, len, dist)
    if len <= length(stack.dists[stack.depth])
        stack.dists[stack.depth][len] = dist
    else
        push!(stack.dists[stack.depth], dist)
    end
end
@inbounds function resize_dists!(stack::SiteStack, len)
    resize!(stack.dists[stack.depth], len)
end

@inbounds function get_dists(stack::SiteStack)
    stack.dists[stack.depth]
end

@inbounds function append_site!(stack::SiteStack, len, site)
    if len <= length(stack.sites[stack.depth])
        stack.sites[stack.depth][len] = site
    else
        push!(stack.sites[stack.depth], site)
    end
end

@inbounds function resize_sites!(stack::SiteStack, len)
    resize!(stack.sites[stack.depth], len)
end

@inbounds function get_sites(stack::SiteStack)
    stack.sites[stack.depth]
end

function pop!(stack::SiteStack)
    stack.depth -= 1
end

function fill_dists!(stack::SiteStack, reference, sites, p::Real = 2)
    #= push_empty!(stack)
    len = 0
    for site in sites
        len += 1
        append_dist!(stack, len, distance(reference, site[2], p))
    end
    resize_dists!(stack, len) =#
    resize_dists!(stack, length(sites))
    copyto!(get_dists(stack), distance(reference, site[2], p) for site in sites)
    # get_dists(stack)
end

function fill_sites!(stack::SiteStack, max_dist, sites)
    #= len = 0
    for (dist, site) in zip(get_dists(stack), sites)
        if dist <= max_dist
            len += 1
            append_site!(stack, len, site)
        end
    end
    resize_sites!(stack, len) =#
    len = sum(1 for dist in get_dists(stack) if dist < max_dist)
    resize_sites!(stack, len)
    copyto!(get_sites(stack),
        site for (dist, site) in zip(get_dists(stack), sites) if dist < max_dist)
    # get_sites(stack)
end

function fill_sites!(stack::SiteStack, min_site, corners, sites, p::Real = 2)
    dists = get_dists(stack)
    len = 0
    for site in sites
        if any(corner -> distance(corner[2], site[2], p) <= dists[corner[1]], enumerate(corners))
            len += 1
            append_site!(stack, len, site)
        end
    end
    resize_sites!(stack, len) 
    #= len = sum(1 for site in sites if any(corner -> distance(corner[2], site[2], p) <= dists[corner[1]], enumerate(corners)))
    resize_sites!(stack, len)
    copyto!(get_sites(stack),
        site for site in sites if any(corner -> distance(corner[2], site[2], p) <= dists[corner[1]], enumerate(corners))) =#
    # get_sites(stack)
end
