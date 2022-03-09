function distance(p, q)
    @fastmath sqrt((p[1] - q[1])^2 + (p[2] - q[2])^2)
end

function nearest(p)
    (round(Int, p[1]), round(Int, p[2]))
end