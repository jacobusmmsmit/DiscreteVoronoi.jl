using LinearAlgebra

function distance(x, y, p::Real=2)
    norm((x[1] - y[1], x[2] - y[2]), p)
end

function nearest(x)
    (round(Int, x[1]), round(Int, x[2]))
end