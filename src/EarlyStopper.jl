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
    if !(1 <= i <= ES.stop_at)
        throw(BoundsError(ES, i))
    end
    return ES.obj[i]
end
Base.getindex(ES::EarlyStopper{T}, i::Number) where {T} = ES[convert(Int, i)]

function Base.setindex!(ES::EarlyStopper{T}, val, i::Int) where {T}
    if !(1 <= i <= ES.stop_at)
        throw(BoundsError(ES, i))
    end
    (ES.obj[i] = val)
end

Base.firstindex(ES::EarlyStopper{T}) where {T} = 1
Base.lastindex(ES::EarlyStopper{T}) where {T} = length(ES)

"""
    early_stop_sort!(early_stopper::ES, predicate::F)::ES where {ES<:EarlyStopper, F<:Function}

    Sorts `early_stopper.obj` based on `predicate` and returns an `EarlyStopper` of the same type
    that contains updated `stop_at` information. The sorting is unstable. 
"""
function early_stop_sort!(early_stopper::ES, predicate::F)::ES where {ES<:EarlyStopper, F<:Function}
    n_trues = early_stopper.stop_at
    for (i, element) in enumerate(early_stopper)
        # n_trues is updated during iteration so we need to manually break.
        if i > n_trues
            break
        end
        # n_trues is the index of the last predicate-truthy element in the vector, swapping
        # causes the current value to be the first predicate-falsey element.
        if !predicate(element)
            counter = 0
            while !predicate(early_stopper[n_trues - counter])
                counter += 1
            end
            swap!(early_stopper, i, n_trues - counter)
            n_trues -= counter
        end
    end
    return EarlyStopper(early_stopper.obj, n_trues)
end

early_stop_sort!(non_early_stopper, predicate) = early_stop_sort!(EarlyStopper(non_early_stopper), predicate)