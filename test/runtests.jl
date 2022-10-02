using DiscreteVoronoi
using StaticArrays
using Test
using BenchmarkTools
using Random

@testset "EarlyStopper.jl" begin
    iv = [1, 2, 3, 4, 5]
    es = EarlyStopper(v, 4) # Specifically define to not satisfy the predicate below
    @testset "Basic Functionality" begin
        @test sum(i for i in es) == 10
        @test length(es) == 4
        @test eltype(es) == eltype(v)
        @test all(es[i] == v[i] for (i, _) in enumerate(es)) 
    end
    
    @testset "Sorting" begin
        predicate(x) = 2 <= x <= 3
        new_es = early_stop_sort!(es, predicate)
        w = shuffle(1:100)
        es2 = EarlyStopper(w)
        predicate2(x) = x <= 50
        @test length(new_es) == count(predicate.(v))
        @test @ballocated(early_stop_sort!($es2, $predicate2), seconds=0.2, samples=100) == 0 
    end
end

@testset "DiscreteVoronoi.jl" begin
    n = 50
    l = 3
    TL = (1, 1)
    BR = (n, n)
    locs = sort([SVector{2}(rand(1:n, 2)) for _ in 1:l])
    grid = zeros(SVector{2,Int}, (n, n))

    ### Test for correctness
    naive_grid = deepcopy(grid); naive_voronoi!(naive_grid, locs) # Baseline for correctness
    dac_grid = deepcopy(grid); dac_voronoi!(dac_grid, locs)
    redac_grid = deepcopy(grid); redac_voronoi!(redac_grid, locs, predicate=exact_condition)

    @testset "Correctness" begin
        #TODO: These will sometimes fail because both answers are technically correct
        # I should implement a function which compares distances only where the two differ.
        # i.e. @test voronoi_equality(grids..., sites)
        @test dac_grid == naive_grid
        @test redac_grid == naive_grid
    end

    @testset "Allocations" begin
        @test (@ballocated redac_voronoi!($grid, $locs, predicate=exact_condition)) == 0
        @test (@ballocated redac_voronoi!($grid, $locs, predicate=centre_anchor_condition)) == 0
        @test (@ballocated dac_voronoi!($grid, $locs)) == 0
        @test (@ballocated naive_voronoi!($grid, $locs)) == 0
        @test (@ballocated jfa_voronoi!($grid, $locs)) == 0
    end
end
