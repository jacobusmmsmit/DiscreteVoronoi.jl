using DiscreteVoronoi
using Distances
using Random
using Plots

Random.seed!(42)
grid = zeros(Int, rand(1:100), rand(1:100))
sites = rand_points(Int, size(grid)..., 10)
jfa_voronoi!(grid, sites)
heatmap(grid)

Random.seed!(42)
grid = zeros(Int, rand(1:100), rand(1:100))
sites = rand_sites(Int, size(grid)..., 10)
dac_voronoi!(grid, sites, original_site_find)
heatmap(grid)

Random.seed!(42)
grid = zeros(Int, rand(1:100), rand(1:100))
sites = rand_sites(Int, size(grid)..., 10)
redac_voronoi!(grid, sites, corner_site_filter)
heatmap(grid)

