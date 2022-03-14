function _norm(xy, p::Real=2)::Float64
    if p == 1
        @fastmath abs(xy[1]) + abs(xy[2])
    elseif p == 2
        @fastmath sqrt(xy[1]^2 + xy[2]^2)
    elseif p == Inf
        @fastmath max(abs(xy[1]), abs(xy[2]))
    else
        @assert false
    end
end

function distance(x, y, p::Real=2)
    _norm((x[1] - y[1], x[2] - y[2]), p)
end

function nearest(x)
    (round(Int, x[1]), round(Int, x[2]))
end