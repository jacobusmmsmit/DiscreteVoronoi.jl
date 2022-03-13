include("../src/DiscreteVoronoi.jl")
using .DiscreteVoronoi
using Random
using BenchmarkTools

function get_sites(::Type{Int}, N, M, K)
    idx = collect(Iterators.product(1:N, 1:M))
    shuffle!(idx)
    idx[1:K]
end


for n in [10, 100, 1000]
    for s in [isqrt(n), n, n * isqrt(n), n * n]
        @show n, s
        println("jfa")
        @btime jfa!(grid, sites) setup = (
            Random.seed!(42);
            grid = zeros(Int, $n, $n);
            sites = get_sites(Int, size(grid)..., $s))

        #= println("dac!")
        @btime dac!(grid, sites) setup=(
            Random.seed!(42);
            grid = zeros(Int32, $n, $n);
            sites = get_sites(size(grid)..., $s)) =#

        println("jdac")
        # Not testing jdac_aux0! due to poor scaling in s (number of sites)
        for aux! in [jdac_aux1a!, jdac_aux1b!, jdac_aux1c!, jdac_aux2a!, jdac_aux2b!, jdac_aux2c!, jdac_aux3!]
            @show aux!
            @btime jdacx!(grid, sites, $aux!) setup = (
                Random.seed!(42);
                grid = zeros(Int, $n, $n);
                sites = collect(enumerate(get_sites(Int, size(grid)..., $s))))
        end
    end
end
