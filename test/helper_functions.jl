function rand_points(::Type{Int}, N, M, K)
    Coord.(shuffle!([Iterators.product(1:N, 1:M)...])[1:min(M * N, K)])
end

function rand_sites(::Type{Int}, N, M, K)
    [(color, point) for (color, point) in enumerate(rand_points(Int, N, M, K))]
end

