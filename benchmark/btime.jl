using DiscreteVoronoi
using BenchmarkTools
using Random
using StaticArrays

using DiscreteVoronoi: exact_aux, centre_anchor_aux

const Coord = SVector{2,Int}
Base.zero(::Type{Coord}) = Coord(0, 0)

function random_coordinates(N, M, K)
    Coord.(shuffle!([Iterators.product(1:N, 1:M)...])[1:min(M * N, K)])
end

for n in [100, 1000]
    for s in [isqrt(n), n, n * isqrt(n)] # , n * n]
        @show n, s

        println("jfa_voronoi!")
        @btime jfa_voronoi!(grid, sites) setup = (
            Random.seed!(42);
            grid = zeros(Coord, ($n, $n));
            sites = random_coordinates(size(grid)..., $s)) evals=1

        if s <= n
            println("dac_voronoi!")
            @btime dac_voronoi!(grid, sites) setup=(
                Random.seed!(42);
                grid = zeros(Coord, ($n, $n));
                sites = random_coordinates(size(grid)..., $s)) evals=1
        end
    
        println("redac_voronoi!")
        @show exact_aux
        @btime redac_voronoi!(grid, sites, auxiliary=exact_aux) setup=(
            Random.seed!(42);
            grid = zeros(Coord, ($n, $n));
            sites = random_coordinates(size(grid)..., $s)) evals=1

        println("redac_voronoi!")
        @show centre_anchor_aux
        @btime redac_voronoi!(grid, sites, auxiliary=centre_anchor_aux) setup=(
            Random.seed!(42);
            grid = zeros(Coord, ($n, $n));
            sites = random_coordinates(size(grid)..., $s)) evals=1
    end
end
