using DiscreteVoronoi
using StaticArrays
using Test
using BenchmarkTools
using Random

using DiscreteVoronoi: get_corners, label, voronoi_equality
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
        @test @ballocated(early_stop_sort!($es2, $predicate2), seconds = 0.2, samples = 100) == 0

        # Representative workload correctness test
        n = 500
        l = 30
        TL, BR = (1, 1), (n, n)
        locs = sort([Coord(rand(1:n, 2)) for _ in 1:l])
        grid = zeros(Coord, (n, n))
        es3 = EarlyStopper(locs)
        predicate3(x) = exact_condition(x, locs, TL, BR)
        new_es3 = early_stop_sort!(es3, predicate3)
        @test new_es3.obj == sort(locs, by=predicate3)
    end
end

@testset "DiscreteVoronoi.jl" begin
    n = 11
    l = 3
    TL, BR = (1, 1), (n, n)
    locs = sort([Coord(rand(1:n, 2)) for _ in 1:l])
    grid = zeros(Coord, (n, n))

    ### Test for correctness
    naive_grid = deepcopy(grid)
    naive_voronoi!(naive_grid, locs) # Baseline for correctness
    dac_grid = deepcopy(grid)
    dac_voronoi!(dac_grid, locs)
    redac_grid = deepcopy(grid)
    redac_voronoi!(redac_grid, locs, predicate=exact_condition)

    @testset "Correctness" begin
        @test voronoi_equality(dac_grid, naive_grid)
        @test voronoi_equality(redac_grid, naive_grid)
    end

    @testset "Allocations" begin
        @test @ballocated(redac_voronoi!($grid, $locs, predicate=exact_condition), seconds = 1.0) == 0
        @test @ballocated(redac_voronoi!($grid, $locs, predicate=centre_anchor_condition), seconds = 1.0) == 0
        @test @ballocated(dac_voronoi!($grid, $locs), seconds = 1.0) == 0
        @test @ballocated(naive_voronoi!($grid, $locs), seconds = 1.0) == 0
        @test @ballocated(jfa_voronoi!($grid, $locs), seconds = 1.0) == 0
    end
end
