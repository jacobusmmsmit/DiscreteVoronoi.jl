include("../src/DiscreteVoronoi.jl")
using .DiscreteVoronoi
using Test
using Random
using Distributions

function get_sites(::Type{Int}, N, M, K)
    idx = [(n, m) for n in 1:N, m in 1:M]
    shuffle!(idx)
    idx[1:K]
end

function get_sites(::Type{T}, N, M, K) where T<:AbstractFloat
    [(rand(Uniform(0, N)), rand(Uniform(0, M))) for k in 1:K]
end

@testset "jfa! working for Int sites" begin
    Random.seed!(42)
    for i in 1:100
        N, M = rand(1:1000, 2)
        sites = collect(enumerate(get_sites(Int, N, M, rand(1:100))))

        grid1 = naive_voronoi(CartesianIndices((1:N, 1:M)), map(site -> site[2], sites))
        grid2 = zeros(Int, N, M)
        @test isa(jfa!(grid2, map(site -> site[2], sites)), Matrix)
    end
end

@testset "dac! matching results for Int sites" begin
    Random.seed!(42)
    for i in 1:100
        N, M = rand(1:1000, 2)
        sites = collect(enumerate(get_sites(Int, N, M, rand(1:100))))

        grid1 = naive_voronoi(CartesianIndices((1:N, 1:M)), map(site -> site[2], sites))
        grid2 = zeros(Int, N, M)
        dac!(grid2, map(site -> site[2], sites))
        @test grid2 == grid1
    end
end

@testset "dac! matching results for Float64 sites" begin
    Random.seed!(42)
    for i in 1:100
        N, M = rand(1:1000, 2)
        sites = collect(enumerate(get_sites(Float64, N, M, rand(1:100))))

        grid1 = naive_voronoi(CartesianIndices((1:N, 1:M)), map(site -> site[2], sites))
        grid2 = zeros(Int, N, M)
        dac!(grid2, map(site -> site[2], sites))
        @test grid2 == grid1
    end
end

@testset "dacx! matching results for Int sites" begin
    Random.seed!(42)
    for i in 1:100
        N, M = rand(1:1000, 2)
        sites = collect(enumerate(get_sites(Int, N, M, rand(1:100))))

        grid1 = naive_voronoi(CartesianIndices((1:N, 1:M)), map(site -> site[2], sites))
        grid2 = zeros(Int, N, M)
        dacx!(grid2, map(site -> site[2], sites))
        @test grid2 == grid1
    end
end 

@testset "dacx! working for Float64 sites" begin
    Random.seed!(42)
    for i in 1:100
        N, M = rand(1:1000, 2)
        sites = collect(enumerate(get_sites(Float64, N, M, rand(1:100))))

        grid1 = naive_voronoi(CartesianIndices((1:N, 1:M)), map(site -> site[2], sites))
        grid2 = zeros(Int, N, M)
        @test isa(dacx!(grid2, map(site -> site[2], sites)), Matrix)
    end
end

@testset "jdac! matching results for Int sites" begin
    Random.seed!(42)
    for i in 1:100
        N, M = rand(1:1000, 2)
        sites = collect(enumerate(get_sites(Int, N, M, rand(1:100))))

        grid1 = naive_voronoi(CartesianIndices((1:N, 1:M)), map(site -> site[2], sites))
        for aux! in [jdac_aux1!, jdac_aux2!, jdac_aux3!, jdac_aux4!, jdac_aux5!, jdac_aux6!]
            grid2 = zeros(Int, N, M)
            jdac!(grid2, sites, aux!)
            @test grid2 == grid1
        end
    end
end

@testset "jdac! matching results for Float64 sites" begin
    Random.seed!(42)
    for i in 1:100
        N, M = rand(1:1000, 2)
        sites = collect(enumerate(get_sites(Float64, N, M, rand(1:100))))

        grid1 = naive_voronoi(CartesianIndices((1:N, 1:M)), map(site -> site[2], sites))
        for aux! in [jdac_aux1!, jdac_aux2!, jdac_aux3!, jdac_aux4!, jdac_aux5!, jdac_aux6!]
            grid2 = zeros(Int, N, M)
            jdac!(grid2, sites, aux!)
            @test grid2 == grid1
        end
    end
end

@testset "jdacx! matching results for Int sites" begin
    Random.seed!(42)
    for i in 1:100
        N, M = rand(1:1000, 2)
        sites = collect(enumerate(get_sites(Int, N, M, rand(1:100))))

        grid1 = naive_voronoi(CartesianIndices((1:N, 1:M)), map(site -> site[2], sites))
        for aux! in [jdac_aux1!, jdac_aux2!, jdac_aux3!, jdac_aux4!, jdac_aux5!, jdac_aux6!]
            grid2 = zeros(Int, N, M)
            jdacx!(grid2, sites, aux!)
            @test grid2 == grid1
        end
    end
end

@testset "jdacx working for Float64 sites" begin
    Random.seed!(42)
    for i in 1:100
        N, M = rand(1:1000, 2)
        sites = collect(enumerate(get_sites(Float64, N, M, rand(1:100))))

        grid1 = naive_voronoi(CartesianIndices((1:N, 1:M)), map(site -> site[2], sites))
        for aux! in [jdac_aux1!, jdac_aux2!, jdac_aux3!, jdac_aux4!, jdac_aux5!, jdac_aux6!]
            grid2 = zeros(Int, N, M)
            @test isa(jdacx!(grid2, sites, aux!), Matrix)
        end
    end
end
