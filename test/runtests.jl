include("../src/DiscreteVoronoi.jl")
using .DiscreteVoronoi
using Test
using Random

function get_seeds(N, M, K)
    idx = [(n, m) for n in 1:N, m in 1:M]
    shuffle!(idx)
    idx[1:K]
end

@testset "dac matching results" begin
    Random.seed!(42)
    for i in 1:100
        N, M = rand(1:1000, 2)
        seeds = collect(enumerate(get_seeds(N, M, rand(1:100))))

        grid1 = naive_voronoi(CartesianIndices((1:N, 1:M)), map(seed -> seed[2], seeds))
        grid2 = zeros(Int, N, M)
        dac!(grid2, map(seed -> seed[2], seeds))
        @test grid2 == grid1
    end
end

@testset "dacx matching results" begin
    Random.seed!(42)
    for i in 1:100
        N, M = rand(1:1000, 2)
        seeds = collect(enumerate(get_seeds(N, M, rand(1:100))))

        grid1 = naive_voronoi(CartesianIndices((1:N, 1:M)), map(seed -> seed[2], seeds))
        grid2 = zeros(Int, N, M)
        dacx!(grid2, map(seed -> seed[2], seeds))
        @test grid2 == grid1
    end
end

@testset "jdac matching results" begin
    Random.seed!(42)
    for i in 1:100
        N, M = rand(1:1000, 2)
        seeds = collect(enumerate(get_seeds(N, M, rand(1:100))))

        grid1 = naive_voronoi(CartesianIndices((1:N, 1:M)), map(seed -> seed[2], seeds))
        for aux! in [jdac_aux1!, jdac_aux2!, jdac_aux3!, jdac_aux4!, jdac_aux5!, jdac_aux6!]
            grid2 = zeros(Int, N, M)
            jdac!(grid2, seeds, aux!)
            @test grid2 == grid1
        end
    end
end

@testset "jdacx matching results" begin
    Random.seed!(42)
    for i in 1:100
        N, M = rand(1:1000, 2)
        seeds = collect(enumerate(get_seeds(N, M, rand(1:100))))

        grid1 = naive_voronoi(CartesianIndices((1:N, 1:M)), map(seed -> seed[2], seeds))
        for aux! in [jdac_aux1!, jdac_aux2!, jdac_aux3!, jdac_aux4!, jdac_aux5!, jdac_aux6!]
            grid2 = zeros(Int, N, M)
            jdacx!(grid2, seeds, aux!)
            @test grid2 == grid1
        end
    end
end

