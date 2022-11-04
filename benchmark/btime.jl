using DiscreteVoronoi
using BenchmarkTools
using Random
using StaticArrays

using DiscreteVoronoi: Coord, exact_aux, centre_anchor_aux, exact_aux_es, centre_anchor_aux_es

function random_coordinates(N, M, K)
    Coord.(shuffle!([Iterators.product(1:N, 1:M)...])[1:min(N * M, K)])
end

for n in [400]
    for s in [isqrt(n), n, n * isqrt(n), n * n]
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

        @show centre_anchor_aux
        @btime redac_voronoi!(grid, sites, auxiliary=centre_anchor_aux) setup=(
            Random.seed!(42);
            grid = zeros(Coord, ($n, $n));
            sites = random_coordinates(size(grid)..., $s)) evals=1

        println("redac_voronoi_es!")
        @show exact_aux_es
        @btime redac_voronoi_es!(grid, sites, auxiliary=exact_aux_es) setup=(
            Random.seed!(42);
            grid = zeros(Coord, ($n, $n));
            sites = random_coordinates(size(grid)..., $s)) evals=1

        @show centre_anchor_aux_es
        @btime redac_voronoi_es!(grid, sites, auxiliary=centre_anchor_aux_es) setup=(
            Random.seed!(42);
            grid = zeros(Coord, ($n, $n));
            sites = random_coordinates(size(grid)..., $s)) evals=1
    end
end
