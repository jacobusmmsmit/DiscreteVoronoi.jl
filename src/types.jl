mutable struct SiteStack{ST<:Tuple{T1,Tuple{T2,T2}} where {T1,T2}}
    depth::Int
    sites::Vector{Vector{ST}}
    SiteStack{ST}() where {ST} = new{ST}(0, Vector{Vector{ST}}())
end

SiteStack(::ST) where {ST} = new(0, Vector{Vector{ST}}())

sitetype(::SiteStack{ST}) where {ST} = ST

function push_empty!(stack::SiteStack{ST}) where {ST}
    stack.depth += 1
    if stack.depth > length(stack.sites)
        push!(stack.sites, Vector{ST}())
    end
    stack.depth
end

@inbounds function append_site!(stack::SiteStack, index, site)
    if index <= length(stack.sites[stack.depth])
        stack.sites[stack.depth][index] = site
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

function set_sites!(stack::SiteStack, sites)
    len = 0
    for (index, site) in enumerate(sites)
        append_site!(stack, index, site)
        len += 1
    end
    resize_sites!(stack, len)
    #= len = length(sites)
    resize_sites!(stack, len)
    copyto!(get_sites(stack), sites) =#
    get_sites(stack)
end

function pop!(stack::SiteStack)
    stack.depth -= 1
end
