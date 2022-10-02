"""
    EarlyStopper{T}
A simple container that contains information about how far to iterate through the object.
 
"""
struct EarlyStopper{T}
    obj::T
    stop_at::Int
end

EarlyStopper(obj::T) where {T} = EarlyStopper(obj, length(eachindex(obj)))

Base.iterate(ES::EarlyStopper{T}, state=1) where {T} = state > ES.stop_at ? nothing : (ES.obj[state], state + 1)
Base.length(ES::EarlyStopper{T}) where {T} = ES.stop_at
Base.eltype(::Type{EarlyStopper{T}}) where {T} = eltype(T)
Base.iterate(rES::Iterators.Reverse{EarlyStopper{T}}, state=rES.itr.stop_at) where {T} = state < 1 ? nothing : (rES.itr.obj[state], state - 1)

function Base.getindex(ES::EarlyStopper{T}, i::Int) where {T}
    1 <= i <= ES.stop_at || throw(BoundsError(ES, i))
    return ES.obj[i]
end
Base.getindex(ES::EarlyStopper{T}, i::Number) where {T} = ES[convert(Int, i)]
Base.firstindex(ES::EarlyStopper{T}) where {T} = 1
Base.lastindex(ES::EarlyStopper{T}) where {T} = length(ES)