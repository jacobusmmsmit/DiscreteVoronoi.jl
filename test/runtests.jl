using DiscreteVoronoi
using Distances
using Test
using Random

const REPETITIONS = 100

@testset "jfa_voronoi! working for Int sites" begin
    Random.seed!(42)
    for distance in [cityblock, euclidean, chebyshev]
        for i in 1:REPETITIONS
            N, M = rand(1:100, 2)
            points = rand_points(Int, N, M, rand(1:100))

            # grid1 = naive_voronoi(Tuple.(CartesianIndices((1:N, 1:M))), points, distance)
            grid2 = zeros(Int, N, M)
            @test !any(==(0), jfa_voronoi!(grid2, points, distance)) 
        end
    end
end

@testset "jfa_voronoi_parallel! working for Int sites" begin
    Random.seed!(42)
    for distance in [cityblock, euclidean, chebyshev]
        for i in 1:REPETITIONS
            N, M = rand(1:100, 2)
            points = rand_points(Int, N, M, rand(1:100))

            # grid1 = naive_voronoi(Tuple.(CartesianIndices((1:N, 1:M))), points, distance)
            grid2 = zeros(Int, N, M)
            @test !any(==(0), jfa_voronoi_parallel!(grid2, points, distance)) 
        end
    end
end

@testset "dac_voronoi! matching results for Int sites" begin
    Random.seed!(42)
    for distance in [cityblock, euclidean, chebyshev]
        for i in 1:REPETITIONS
            N, M = rand(1:100, 2)
            sites = rand_sites(Int, N, M, rand(1:100))

            grid1 = naive_voronoi(Tuple.(CartesianIndices((1:N, 1:M))), map(site -> site[2], sites), distance)
            for site_find in [original_site_find, center_site_find]
                grid2 = preset_voronoi!(zeros(Int, N, M), sites)
                @test dac_voronoi!(grid2, sites, site_find, distance) == grid1
            end
        end
    end
end

@testset "dac_voronoi! matching results for Float64 sites" begin
    Random.seed!(42)
    for distance in [cityblock, euclidean, chebyshev]
        for i in 1:REPETITIONS
            N, M = rand(1:100, 2)
            sites = rand_sites(Float64, N, M, rand(1:100))

            grid1 = naive_voronoi(Tuple.(CartesianIndices((1:N, 1:M))), map(site -> site[2], sites), distance)
            for site_find in [original_site_find, center_site_find]
                grid2 = zeros(Int, N, M)
                @test dac_voronoi!(grid2, sites, site_find, distance) == grid1
            end
        end
    end
end

@testset "dac_voronoi! working for Float64 sites" begin
    Random.seed!(42)
    for distance in [cityblock, euclidean, chebyshev]
        for i in 1:REPETITIONS
            N, M = rand(1:100, 2)
            sites = rand_sites(Float64, N, M, rand(1:100))

            for site_find in [original_site_find, center_site_find]
                grid = preset_voronoi_rounded!(zeros(Int, N, M), sites)
                @test !any(==(0), dac_voronoi!(grid, sites, site_find, distance))
            end
         end
    end
end

@testset "redac_voronoi! matching results for Int sites" begin
    Random.seed!(42)
    for distance in [cityblock, euclidean, chebyshev]
        for i in 1:REPETITIONS
            N, M = rand(1:100, 2)
            sites = rand_sites(Int, N, M, rand(1:100))

            grid1 = naive_voronoi(Tuple.(CartesianIndices((1:N, 1:M))), map(site -> site[2], sites), distance)
            for site_filter in [center_site_filter, anchor_site_filter]
                grid2 = preset_voronoi!(zeros(Int, N, M), sites)
                @test redac_voronoi!(grid2, sites, site_filter, distance) == grid1
            end
        end
    end
end

@testset "redac_voronoi! matching results for Float64 sites" begin
    Random.seed!(42)
    for distance in [cityblock, euclidean, chebyshev]
        for i in 1:REPETITIONS
            N, M = rand(1:100, 2)
            sites = rand_sites(Float64, N, M, rand(1:100))

            grid1 = naive_voronoi(Tuple.(CartesianIndices((1:N, 1:M))), map(site -> site[2], sites), distance)
            for site_filter in [center_site_filter, anchor_site_filter]
                grid2 = zeros(Int, N, M)
                @test redac_voronoi!(grid2, sites, site_filter, distance) == grid1
            end
        end
    end
end

@testset "redac_voronoi! working for Float64 sites" begin
    Random.seed!(42)
    for distance in [cityblock, euclidean, chebyshev]
        for i in 1:REPETITIONS
            N, M = rand(1:100, 2)
            sites = rand_sites(Float64, N, M, rand(1:100))

            for site_filter in [center_site_filter, anchor_site_filter]
                grid = preset_voronoi_rounded!(zeros(Int, N, M), sites)
                @test !any(==(0), redac_voronoi!(grid, sites, site_filter, distance))
            end
        end
    end
end

@testset "redac_voronoi_optimized! matching results for Int sites" begin
    Random.seed!(42)
    for distance in [cityblock, euclidean, chebyshev]
        for i in 1:REPETITIONS
            N, M = rand(1:100, 2)
            sites = rand_sites(Int, N, M, rand(1:100))

            grid1 = naive_voronoi(Tuple.(CartesianIndices((1:N, 1:M))), map(site -> site[2], sites), distance)
            for site_filter! in [center_site_filter!, anchor_site_filter!]
                grid2 = preset_voronoi!(zeros(Int, N, M), sites)
                @test redac_voronoi_optimized!(grid2, sites, site_filter!, distance) == grid1
            end
        end
    end
end

@testset "redac_voronoi_optimized! matching results for Float64 sites" begin
    Random.seed!(42)
    for distance in [cityblock, euclidean, chebyshev]
        for i in 1:REPETITIONS
            N, M = rand(1:100, 2)
            sites = rand_sites(Float64, N, M, rand(1:100))

            grid1 = naive_voronoi(Tuple.(CartesianIndices((1:N, 1:M))), map(site -> site[2], sites), distance)
            for site_filter! in [center_site_filter!, anchor_site_filter!]
                grid2 = zeros(Int, N, M)
                @test redac_voronoi_optimized!(grid2, sites, site_filter!, distance) == grid1
            end
        end
    end
end

@testset "redac_voronoi_optimized! working for Float64 sites" begin
    Random.seed!(42)
    for distance in [cityblock, euclidean, chebyshev]
        for i in 1:REPETITIONS
            N, M = rand(1:100, 2)
            sites = rand_sites(Float64, N, M, rand(1:100))

            for site_filter! in [center_site_filter!, anchor_site_filter!]
                grid = preset_voronoi_rounded!(zeros(Int, N, M), sites)
                @test !any(==(0), redac_voronoi_optimized!(grid, sites, site_filter!, distance))
            end
        end
    end
end
