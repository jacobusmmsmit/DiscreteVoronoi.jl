include("../src/DiscreteVoronoi.jl")
using .DiscreteVoronoi
using Random
using BenchmarkTools

function rand_sites(::Type{Int}, N, M, K)
    idx = collect(Iterators.product(1:N, 1:M))
    shuffle!(idx)
    idx[1:K]
end


Random.seed!(42)
grid = zeros(Int, rand(1:100), rand(1:100))
sites = rand_sites(Int, size(grid)..., 10)
jfa!(grid, sites)
@show grid

Random.seed!(42)
grid = zeros(Int, rand(1:100), rand(1:100))
sites = rand_sites(Int, size(grid)..., 10)
dac!(grid, sites)
@show grid

Random.seed!(42)
grid = zeros(Int, rand(1:100), rand(1:100))
sites = collect(enumerate(rand_sites(Int, size(grid)..., 10)))
jdac!(grid, sites, jdac_aux1a!)
@show grid