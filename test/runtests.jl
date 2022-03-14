include("../src/DiscreteVoronoi.jl")
using .DiscreteVoronoi
using Test
using Random
using Distributions

const REPETITIONS = 10

function rand_sites(::Type{Int}, N, M, K)
    idx = collect(Iterators.product(1:N, 1:M))
    shuffle!(idx)
    idx[1:K]
end

function rand_sites(::Type{T}, N, M, K) where {T<:AbstractFloat}
    [(rand(Uniform(0, N)), rand(Uniform(0, M))) for k in 1:K]
end

@testset "jfa! working for Int sites" begin
    Random.seed!(42)
    for p in [1, 2, Inf]
        for i in 1:REPETITIONS
            N, M = rand(1:1000, 2)
            sites = collect(enumerate(rand_sites(Int, N, M, rand(1:100))))

            grid1 = naive_voronoi(CartesianIndices((1:N, 1:M)), map(site -> site[2], sites), p)
            grid2 = zeros(Int, N, M)
            @test isa(jfa!(grid2, map(site -> site[2], sites), p), Matrix)
        end
    end
end

@testset "dac! matching results for Int sites" begin
    Random.seed!(42)
    for p in [1, 2, Inf]
        for i in 1:REPETITIONS
            N, M = rand(1:1000, 2)
            sites = collect(enumerate(rand_sites(Int, N, M, rand(1:100))))

            grid1 = naive_voronoi(CartesianIndices((1:N, 1:M)), map(site -> site[2], sites), p)
            grid2 = zeros(Int, N, M)
            dac!(grid2, map(site -> site[2], sites), p)
            @test grid2 == grid1
        end
    end
end

@testset "dac! matching results for Float64 sites" begin
    Random.seed!(42)
    for p in [1, 2, Inf]
        for i in 1:REPETITIONS
            N, M = rand(1:1000, 2)
            sites = collect(enumerate(rand_sites(Float64, N, M, rand(1:100))))

            grid1 = naive_voronoi(CartesianIndices((1:N, 1:M)), map(site -> site[2], sites), p)
            grid2 = zeros(Int, N, M)
            dac!(grid2, map(site -> site[2], sites), p)
            @test grid2 == grid1
        end
    end
end

@testset "dacx! matching results for Int sites" begin
    Random.seed!(42)
    for p in [1, 2, Inf]
        for i in 1:REPETITIONS
            N, M = rand(1:1000, 2)
            sites = collect(enumerate(rand_sites(Int, N, M, rand(1:100))))

            grid1 = naive_voronoi(CartesianIndices((1:N, 1:M)), map(site -> site[2], sites), p)
            grid2 = zeros(Int, N, M)
            dacx!(grid2, map(site -> site[2], sites), p)
            @test grid2 == grid1
        end
    end
end

@testset "dacx! working for Float64 sites" begin
    Random.seed!(42)
    for p in [1, 2, Inf]
        for i in 1:REPETITIONS
            N, M = rand(1:1000, 2)
            sites = collect(enumerate(rand_sites(Float64, N, M, rand(1:100))))

            grid1 = naive_voronoi(CartesianIndices((1:N, 1:M)), map(site -> site[2], sites), p)
            grid2 = zeros(Int, N, M)
            @test isa(dacx!(grid2, map(site -> site[2], sites), p), Matrix)
        end
    end
end 

@testset "jdac! matching results for Int sites" begin
    Random.seed!(42)
    for p in [1, 2, Inf]
        for i in 1:REPETITIONS
            N, M = rand(1:1000, 2)
            sites = collect(enumerate(rand_sites(Int, N, M, rand(1:100))))

            grid1 = naive_voronoi(CartesianIndices((1:N, 1:M)), map(site -> site[2], sites), p)
            for aux! in [jdac_aux0!, jdac_aux1a!, jdac_aux1b!, jdac_aux1c!, jdac_aux2a!, jdac_aux2b!, jdac_aux2c!, jdac_aux3a!, jdac_aux3b!, jdac_aux3c!]
                grid2 = zeros(Int, N, M)
                jdac!(grid2, sites, aux!, p)
                @test grid2 == grid1
            end
        end
    end
end

@testset "jdac! matching results for Float64 sites" begin
    Random.seed!(42)
    for p in [1, 2, Inf]
        for i in 1:REPETITIONS
            N, M = rand(1:1000, 2)
            sites = collect(enumerate(rand_sites(Float64, N, M, rand(1:100))))

            grid1 = naive_voronoi(CartesianIndices((1:N, 1:M)), map(site -> site[2], sites), p)
            for aux! in [jdac_aux0!, jdac_aux1a!, jdac_aux1b!, jdac_aux1c!, jdac_aux2a!, jdac_aux2b!, jdac_aux2c!, jdac_aux3a!, jdac_aux3b!, jdac_aux3c!]
                grid2 = zeros(Int, N, M)
                jdac!(grid2, sites, aux!, p)
                @test grid2 == grid1
            end
        end
    end
end

@testset "jdacx! matching results for Int sites" begin
    Random.seed!(42)
    for p in [1, 2, Inf]
        for i in 1:REPETITIONS
            N, M = rand(1:1000, 2)
            sites = collect(enumerate(rand_sites(Int, N, M, rand(1:100))))

            grid1 = naive_voronoi(CartesianIndices((1:N, 1:M)), map(site -> site[2], sites), p)
            for aux! in [jdac_aux0!, jdac_aux1a!, jdac_aux1b!, jdac_aux1c!, jdac_aux2a!, jdac_aux2b!, jdac_aux2c!, jdac_aux3a!, jdac_aux3b!, jdac_aux3c!]
                grid2 = zeros(Int, N, M)
                jdacx!(grid2, sites, aux!, p)
                @test grid2 == grid1
            end
        end
    end
end

@testset "jdacx! working for Float64 sites" begin
    Random.seed!(42)
    for p in [1, 2, Inf]
        for i in 1:REPETITIONS
            N, M = rand(1:1000, 2)
            sites = collect(enumerate(rand_sites(Float64, N, M, rand(1:100))))

            grid1 = naive_voronoi(CartesianIndices((1:N, 1:M)), map(site -> site[2], sites), p)
            for aux! in [jdac_aux0!, jdac_aux1a!, jdac_aux1b!, jdac_aux1c!, jdac_aux2a!, jdac_aux2b!, jdac_aux2c!, jdac_aux3a!, jdac_aux3b!, jdac_aux3c!]
                grid2 = zeros(Int, N, M)
                @test isa(jdacx!(grid2, sites, aux!), Matrix)
            end
        end
    end
end 
