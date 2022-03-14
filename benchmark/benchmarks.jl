include("../src/DiscreteVoronoi.jl")
using .DiscreteVoronoi
using Random
using BenchmarkTools

function rand_sites(::Type{Int}, N, M, K)
    idx = collect(Iterators.product(1:N, 1:M))
    shuffle!(idx)
    idx[1:K]
end


const SUITE = BenchmarkGroup()

for n in [10, 100, 1000]
    SUITE[string("grid ", n, "x", n)] = BenchmarkGroup()
    for s in [isqrt(n), n, n * isqrt(n), n * n]
        SUITE[string("grid ", n, "x", n)][string(s, " sites")] = BenchmarkGroup()

        SUITE[string("grid ", n, "x", n)][string(s, " sites")]["jfa!"] = @benchmarkable jfa!(grid, sites) setup=(
            Random.seed!(42);
            grid = zeros(Int, $n, $n);
            sites = rand_sites(Int, size(grid)..., $s))

        for aux! in [jdac_aux1a!, jdac_aux1b!, jdac_aux1c!, jdac_aux2a!, jdac_aux2b!, jdac_aux2c!, jdac_aux3a!, jdac_aux3b!, jdac_aux3c!]
            SUITE[string("grid ", n, "x", n)][string(s, " sites")][aux!] = @benchmarkable jdacx!(grid, sites, $aux!) setup=(
                Random.seed!(42);
                grid = zeros(Int, $n, $n);
                sites = collect(enumerate(rand_sites(Int, size(grid)..., $s))))
        end
    end
end