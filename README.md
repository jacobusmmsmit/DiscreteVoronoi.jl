# DiscreteVoronoi

A package for computing discrete approximations of Voronoi diagrams. All Voronoi diagram calculating functions are in-place.

```julia
using DiscreteVoronoi
using StaticVectors

# Testing
n = 50
p = 10
grid = zeros(SVector{2,Int}, n, n)
sites = [SVector{2,Int}(rand(1:n, 2)) for _ in 1:p]
naive_voronoi!(grid, sites)
```

There are currently three ways of computing discrete Voronoi diagrams exported, the first three of which are completely allocation free:
* A naive method which simply compares all cells to all sites and chooses the closest
* The jump flooding algorithm, which is explained in more detail [here](https://en.wikipedia.org/wiki/Jump_flooding_algorithm)
* A divide-and-conquer method first detailed [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7840081/).
* A reduce-divide-and-conquer method additionally runs a site elimination algorithm before each recursive step. This elimination algorithm aims to reduce the amount of unnecessary work performed by the algorithm in subsequent recursions by removing seeds that are far away from the corners.
