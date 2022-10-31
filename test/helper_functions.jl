function random_coordinates(N, K)
    Coord.(shuffle!([Iterators.product(1:N, 1:N)...])[1:min(N * N, K)])
end

function random_coordinates(N, M, K)
    Coord.(shuffle!([Iterators.product(1:N, 1:M)...])[1:min(N * M, K)])
end