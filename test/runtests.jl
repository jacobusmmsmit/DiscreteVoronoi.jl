using DiscreteVoronoi
using StaticArrays
using Test
using BenchmarkTools
using Random
using Distances

include("helper_functions.jl")

using DiscreteVoronoi: get_corners, voronoi_equality, exact_condition
using DiscreteVoronoi: exact_aux, centre_anchor_aux, exact_aux_es, centre_anchor_aux_es
using DiscreteVoronoi: EarlyStopper, early_stop_sort!

const Coord = SVector{2,Int}

@testset "Elimination methods" begin
    # Setup: a grid with sites on and directly around each corner.
    # Expected result: after exact elimination, we should see only the sites inside the grid.
    function sites_on_and_around(coord)
        a, b = coord
        dirs = SVector((-1, 0), (0, 0), (1, 0), (0, -1), (0, 1))
        return [Coord(a + i, b + j) for (i, j) in dirs]
    end
    TL, BR = (1, 1), (5, 5)
    sites = mapreduce(sites_on_and_around, vcat, get_corners(TL, BR))
    exact_newsites = sites[map(s -> exact_condition(s, sites, TL, BR), sites)]
    @test all(s -> all(TL .<= s .<= BR), exact_newsites)  # All `exact_newsites` are inside rectangle TL-BR
end

@testset "EarlyStopper.jl" begin
    v = [1, 2, 3, 4, 5]
    es = EarlyStopper(v, 4) # Specifically define to not satisfy the predicate below
    @testset "Basic Functionality" begin
        @test sum(i for i in es) == 10
        @test length(es) == 4
        @test eltype(es) == eltype(v)
        @test all(es[i] == v[i] for (i, _) in enumerate(es))
    end

    @testset "Sorting" begin
        # Basic test
        predicate(x) = 2 <= x <= 3
        new_es = early_stop_sort!(es, predicate)
        @test length(new_es) == count(predicate.(v))

        # Allocation test
        w = shuffle(1:100)
        es2 = EarlyStopper(w)
        predicate2(x) = x <= 50
        @test @ballocated(
            early_stop_sort!($es2, $predicate2), seconds = 0.2, samples = 100
        ) == 0

        # Representative workload correctness test
        n = 500
        l = 30
        TL, BR = (1, 1), (n, n)
        sites = random_coordinates(n, l)
        # locs = sort([Coord(rand(1:n, 2)) for _ in 1:l]) # Possibly overlapping sites.
        grid = zeros(Coord, (n, n))
        es3 = EarlyStopper(sites)
        predicate3(x) = exact_condition(x, sites, TL, BR)
        new_es3 = early_stop_sort!(es3, predicate3)
        @test new_es3.obj == sort(sites; by=predicate3)
    end
end

@testset "DiscreteVoronoi.jl" begin
    n = 50
    l = 30
    TL, BR = (1, 1), (n, n)

    #TODO: Test with overlapping sites.
    #TODO: Concrete tests to ensure `naive_voronoi` is correct.
    sites = random_coordinates(n, l)
    grid = zeros(Coord, (n, n))

    @testset "DAC correctness" begin
        for metric in (euclidean, cityblock, chebyshev)
            naive_grid = copy(grid)
            naive_voronoi!(naive_grid, sites) # Baseline for correctness
            dac_grid = copy(grid)
            dac_voronoi!(dac_grid, sites)
            @test voronoi_equality(naive_grid, dac_grid; distance=metric)
        end
    end

    @testset "REDAC correctness" begin
        for metric in (euclidean, cityblock, chebyshev)
            naive_grid = copy(grid)
            naive_voronoi!(naive_grid, sites; distance=metric) # Baseline for correctness
            for aux in (exact_aux, centre_anchor_aux)
                redac_grid = copy(grid)
                redac_voronoi!(redac_grid, sites; distance=metric, auxiliary=aux)
                @test voronoi_equality(naive_grid, redac_grid; distance=metric)
            end
        end
    end

    @testset "Auxiliary function allocations" begin
        @test @ballocated(
            exact_aux($grid, $sites, $TL, $(BR .÷ 2)), seconds = 1.0
        ) == 0
        @test @ballocated(
            centre_anchor_aux($grid, $sites, $TL, $(BR .÷ 2)), seconds = 1.0
        ) == 0
    end

    @testset "Main function allocations" begin
        @test @ballocated(naive_voronoi!($grid, $sites), seconds = 1.0) == 0
        @test @ballocated(dac_voronoi!($grid, $sites), seconds = 1.0) == 0
        @test @ballocated(jfa_voronoi!($grid, $sites), seconds = 1.0) == 0
        @test @ballocated(
            redac_voronoi!($grid, $sites, auxiliary=exact_aux), seconds = 1.0
        ) == 0
        @test @ballocated(
            redac_voronoi!($grid, $sites, auxiliary=centre_anchor_aux), seconds = 1.0
        ) == 0
    end
end

const REPETITIONS = 100

@testset verbose=true "jfa_voronoi! working for Int sites" begin
    Random.seed!(42)
    for distance in [cityblock, euclidean, chebyshev]
        test_passes = true
        @testset "$distance" begin 
            for i in 1:REPETITIONS
                N, M = rand(1:100, 2)
                points = random_coordinates(N, M, rand(1:100))

                # grid1 = zeros(Coord, (N, M))
                # naive_voronoi!(grid1, points, distance=distance)
                grid2 = zeros(Coord, (N, M))
                jfa_voronoi!(grid2, points, distance=distance)
                test_passes = test_passes && !any(==(0), grid2) 
                test_passes || break
            end
            test_passes ? (@test test_passes) : (@test_broken test_passes)
        end
    end
end

@testset verbose=true "dac_voronoi! matching results for Int sites" begin
    Random.seed!(42)
    for distance in [cityblock, euclidean, chebyshev]
        @testset "$distance" begin
            test_passes = true
            for i in 1:REPETITIONS
                N, M = rand(1:100, 2)
                points = random_coordinates(N, M, rand(1:100))

                grid1 = zeros(Coord, (N, M))
                naive_voronoi!(grid1, points, distance=distance)
                grid2 = zeros(Coord, (N, M))
                dac_voronoi!(grid2, points, distance=distance)
                test_passes = test_passes && grid2 == grid1
                test_passes || break
            end
            test_passes ? (@test test_passes) : (@test_broken test_passes)
        end
    end
end

@testset verbose=true "redac_voronoi! matching results for Int sites" begin
    Random.seed!(42)
    for distance in [cityblock, euclidean, chebyshev]
        @testset verbose=true "$distance" begin
            for auxiliary in [exact_aux, centre_anchor_aux]
                @testset "$auxiliary" begin
                    test_passes = true
                    for i in 1:REPETITIONS
                        N, M = rand(1:100, 2)
                        points = random_coordinates(N, M, rand(1:100))

                        grid1 = zeros(Coord, (N, M))
                        naive_voronoi!(grid1, points, distance=distance)
                        grid2 = zeros(Coord, (N, M))
                        redac_voronoi!(grid2, points, distance=distance, auxiliary=auxiliary)
                        test_passes = test_passes && voronoi_equality(grid2, grid1; distance=distance)
                        test_passes || break
                    end
                    test_passes ? (@test test_passes) : (@test_broken test_passes)
                end
            end
        end
    end
end

@testset verbose=true "redac_voronoi_es! matching results for Int sites" begin
    Random.seed!(42)
    for distance in [cityblock, euclidean, chebyshev]
        @testset verbose=true "$distance" begin
            for auxiliary in [exact_aux_es, centre_anchor_aux_es]
                @testset "$auxiliary" begin
                    test_passes = true
                    for i in 1:REPETITIONS
                        N, M = rand(1:100, 2)
                        points = random_coordinates(N, M, rand(1:100))

                        grid1 = zeros(Coord, (N, M))
                        naive_voronoi!(grid1, points, distance=distance)
                        grid2 = zeros(Coord, (N, M))
                        redac_voronoi_es!(grid2, points, distance=distance, auxiliary=auxiliary)
                        test_passes = test_passes && voronoi_equality(grid2, grid1; distance=distance)
                        test_passes || break
                    end
                    test_passes ? (@test test_passes) : (@test_broken test_passes)
                end
            end
        end
    end
end
