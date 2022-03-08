# DiscreteVoronoi.jl

This is a package for the efficient calculation of discrete voronoi diagrams.

It implements the jump-flood algorithm (on the CPU, maybe we can do GPU in `CUDA.jl`),
the divide-and-conquer algorith (link to 2020 paper), and a "seed justified" divide-
and-conquer algorithm.

```
using DiscreteVoronoi.jl

grid = zeros(10, 10)
nsites = 4
sites = [rand(1:10, 1:10) for _ in 1:nsites]
jfa!(grid, sites)
dacx!(grid, sites)
jdacx!(grid, sites)
```
