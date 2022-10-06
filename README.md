# DiscreteVoronoi

A package for computing discrete approximations of Voronoi diagrams. All Voronoi diagram calculating functions are in-place.

Currently, to use this package you need to declare your sites as a `Vector{SVector{2,Int}}` and your grid as a `Matrix{SVector{2,Int}}` and hence there is a strong dependency on `StaticArrays`.

```julia
using DiscreteVoronoi
using StaticArrays

# Testing
n = 50
p = 10
grid = zeros(SVector{2,Int}, n, n)
sites = [SVector{2,Int}(rand(1:n, 2)) for _ in 1:p]
redac_voronoi!(grid, sites)
```

There are currently three ways of computing discrete Voronoi diagrams exported, they are all completely allocation free and called in the same way:
* `*_voronoi!(grid::Matrix{SVector{2,Int}}, sites::Vector{SVector{2,Int}}` is the generic function call for all of the following functions,
* `naive_voronoi!` simply compares all cells to all sites and chooses the closest
* `jfa_voronoi!` uses the jump flooding algorithm, which is explained in more detail [here](https://en.wikipedia.org/wiki/Jump_flooding_algorithm)
* `dac_voronoi!` employs a divide-and-conquer method first detailed [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7840081/).
* `redac_voronoi!` (reduce-divide-and-conquer) additionally runs a site elimination algorithm before each recursive step. This elimination algorithm aims to reduce the amount of unnecessary work performed by the algorithm in subsequent recursions by removing seeds that are far away from the corners.

Which algorithm you should use for the fastest execution time depends somewhat on the task at hand. The best thing to do is to try all above for your usecase and decide from benchmarks. As a rule of thumb, the larger the grid is the better the divide-and-conquer methods will be in comparison. Particularly if the number of sites scales with the size of the grid (`n^2`) faster than `log(n)` (natural log), then `redac_voronoi!` will be much faster than `dac_voronoi!`. 

Additionally, the package exports some helper functions for analysing Voronoi diagrams and writing your own algorithms:
* `find_closest_site` finds the closest site to a specified cell in the Lp sense.
* `get_corners` and `get_quadrants` take the top-left (TL) and bottom-right (BR) corners of a rectangle and return the TL and BR corners, and non-overlapping quadrants (calculated by integer division) respectively.
* `label_voronoi_grid` takes a grid of `SVector{2, Int}` and labels each unique value with an integer in a new grid of the same size so it can be visualised.
* `voronoi_equality` can be used to test equality of resulting Voronoi diagrams taking into account that some sites may be the same distance from certain cells and so there are multiple valid/correct diagrams that could be produced.

## Work in progress:
 * Support for arbitrary distance functions.
 * Implementing a hybrid version of `redac_voronoi!` and `dac_voronoi!` that switches to `naive_voronoi!` once a certain size is reached.
 * Currently, I have not implemented multithreaded (or GPU in the case of JFA) versions of these algorithms on the main branch, but the legacy branch contains versions of the algorithms that do have this capability.

## Contributions:
* [@jacobusmmsmit](github.com/jacobusmmsmit) - Author and maintainer
* [@goerch](github.com/goerch) - Author
* [@marcelroed](github.com/marcelroed) - Algorithmic improvements `early_stop_sort!`
