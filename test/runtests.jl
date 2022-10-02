using DiscreteVoronoi
using StaticArrays
using Test

@testset "DiscreteVoronoi.jl" begin
    ### Example ###
    n = 50
    l = 3
    TL = (1, 1)
    BR = (n, n)
    locs = sort([SVector{2}(rand(1:n, 2)) for _ in 1:l])
    grid = zeros(SVector{2,Int}, (n, n))

    ### Test for correctness
    naive_grid = deepcopy(grid); naive_voronoi!(naive_grid, locs) # Baseline for correctness
    dac_grid = deepcopy(grid); dac_voronoi!(dac_grid, locs)
    redac_grid = deepcopy(grid); redac_voronoi!(redac_grid, locs, exact_condition)

    @test dac_grid == naive_grid
    @test redac_grid == naive_grid

    ### Test for allocations
    @test @ballocated redac_voronoi!($grid, $locs, predicate=exact_condition) == 0
    @test @ballocated redac_voronoi!($grid, $locs, predicate=centre_anchor_condition) == 0
    @test @ballocated dac_voronoi!($grid, $locs) == 0
    @test @ballocated naive_voronoi!($grid, $locs) == 0
    @test @ballocated jfa_voronoi!($grid, $locs) == 0
end
