include("../src/DiscreteVoronoi.jl")
using .DiscreteVoronoi
using Random
using BenchmarkTools

function get_sites(N, M, K)
    idx = [(n, m) for n in 1:N, m in 1:M]
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
            sites = get_sites(size(grid)..., $s))

        for aux! in [jdac_aux1!, jdac_aux2!, jdac_aux3!, jdac_aux4!, jdac_aux5!, jdac_aux6!]
            SUITE[string("grid ", n, "x", n)][string(s, " sites")][aux!] = @benchmarkable jdacx!(grid, sites, $aux!) setup=(
                Random.seed!(42);
                grid = zeros(Int, $n, $n);
                sites = collect(enumerate(get_sites(size(grid)..., $s))))
        end
    end
end 

