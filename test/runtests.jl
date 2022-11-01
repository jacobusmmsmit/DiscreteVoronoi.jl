using DiscreteVoronoi
using Distances
using Test
using Random

const REPETITIONS = 100

@testset "jfa_voronoi! working for Int sites" begin
    Random.seed!(42)
    for distance in [cityblock, euclidean, chebyshev]
        @testset verbose=true "$distance" begin
            for i in 1:REPETITIONS
                N, M = rand(1:100, 2)
                points = rand_points(Int, N, M, rand(1:100))

                # grid1 = zeros(Int, N, M)
                # naive_voronoi!(grid1, points, distance)
                grid2 = zeros(Int, N, M)
                @test !any(==(0), jfa_voronoi!(grid2, points, distance))
            end
        end
    end
end

@testset "jfa_voronoi_parallel! working for Int sites" begin
    Random.seed!(42)
    for distance in [cityblock, euclidean, chebyshev]
        @testset verbose=true "$distance" begin
            for i in 1:REPETITIONS
                N, M = rand(1:100, 2)
                points = rand_points(Int, N, M, rand(1:100))

                # grid1 = zeros(Int, N, M)
                # naive_voronoi!(grid1, points, distance)
                grid2 = zeros(Int, N, M)
                @test !any(==(0), jfa_voronoi_parallel!(grid2, points, distance))
            end
        end
    end
end

@testset "dac_voronoi! matching results for Int sites" begin
    Random.seed!(42)
    for site_find in [original_site_find, center_site_find]
        @testset verbose=true "$site_find" begin
            for distance in [cityblock, euclidean, chebyshev]
                @testset verbose=true "$distance" begin
                    for i in 1:REPETITIONS
                        N, M = rand(1:100, 2)
                        sites = rand_sites(Int, N, M, rand(1:100))

                        grid1 = zeros(Int, N, M)
                        naive_voronoi!(grid1, map(site -> site[2], sites), distance)
                        grid2 = preset_voronoi!(zeros(Int, N, M), sites)
                        @test dac_voronoi!(grid2, sites, site_find, distance) == grid1
                    end
                end
            end
        end
    end
end

@testset "dac_voronoi! matching results for Float64 sites" begin
    Random.seed!(42)
    for site_find in [original_site_find, center_site_find]
        @testset verbose=true "$site_find" begin
            for distance in [cityblock, euclidean, chebyshev]
                @testset verbose=true "$distance" begin
                    for i in 1:REPETITIONS
                        N, M = rand(1:100, 2)
                        sites = rand_sites(Float64, N, M, rand(1:100))

                        grid1 = zeros(Int, N, M)
                        naive_voronoi!(grid1, map(site -> site[2], sites), distance)
                        grid2 = zeros(Int, N, M)
                        @test dac_voronoi!(grid2, sites, site_find, distance) == grid1
                    end
                end
            end
        end
    end
end

@testset "dac_voronoi! working for Float64 sites" begin
    Random.seed!(42)
    for site_find in [original_site_find, center_site_find]
        @testset verbose=true "$site_find" begin
            for distance in [cityblock, euclidean, chebyshev]
                @testset verbose=true "$distance" begin
                    for i in 1:REPETITIONS
                        N, M = rand(1:100, 2)
                        sites = rand_sites(Float64, N, M, rand(1:100))

                        grid = preset_voronoi_rounded!(zeros(Int, N, M), sites)
                        @test !any(==(0), dac_voronoi!(grid, sites, site_find, distance))
                    end
                end
            end
        end
    end
end

@testset "redac_voronoi! matching results for Int sites" begin
    Random.seed!(42)
    for method in [:filter, :partition]
        @testset verbose=true "$method" begin
            for site_filter in [naive_site_filter, center_site_filter, anchor_site_filter, corner_site_filter]
                @testset verbose=true "$site_filter" begin
                    for distance in [cityblock, euclidean, chebyshev]
                        @testset verbose=true "$distance" begin
                            for i in 1:REPETITIONS
                                N, M = rand(1:100, 2)
                                sites = rand_sites(Int, N, M, rand(1:100))

                                grid1 = zeros(Int, N, M)
                                naive_voronoi!(grid1, map(site -> site[2], sites), distance)
                                grid2 = preset_voronoi!(zeros(Int, N, M), sites)
                                if method === :filter
                                    @test redac_voronoi!(Val(method), grid2, sites, site_filter, distance) == grid1
                                elseif method === :partition
                                    @test equal_grid_sites(redac_voronoi!(Val(method), grid2, sites, site_filter, distance), grid1, sites, distance)
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

@testset "redac_voronoi! matching results for Float64 sites" begin
    Random.seed!(42)
    for site_filter in [naive_site_filter, center_site_filter, anchor_site_filter, corner_site_filter]
        @testset verbose=true "$site_filter" begin
            for distance in [cityblock, euclidean, chebyshev]
                @testset verbose=true "$distance" begin
                    for i in 1:REPETITIONS
                        N, M = rand(1:100, 2)
                        sites = rand_sites(Float64, N, M, rand(1:100))

                        grid1 = zeros(Int, N, M)
                        naive_voronoi!(grid1, map(site -> site[2], sites), distance)
                        grid2 = zeros(Int, N, M)
                        @test redac_voronoi!(Val(:filter), grid2, sites, site_filter, distance) == grid1
                    end
                end
            end
        end
    end
end

@testset "redac_voronoi! working for Float64 sites" begin
    Random.seed!(42)
    for site_filter in [naive_site_filter, center_site_filter, anchor_site_filter, corner_site_filter]
        @testset verbose=true "$site_filter" begin
            for distance in [cityblock, euclidean, chebyshev]
                @testset verbose=true "$distance" begin
                    for i in 1:REPETITIONS
                        N, M = rand(1:100, 2)
                        sites = rand_sites(Float64, N, M, rand(1:100))

                        grid = preset_voronoi_rounded!(zeros(Int, N, M), sites)
                        @test !any(==(0), redac_voronoi!(Val(:filter), grid, sites, site_filter, distance))
                    end
                end
            end
        end
    end
end
