# DiscreteVoronoi.jl

This is a package for the efficient calculation of discrete voronoi diagrams.

It implements the jump-flood algorithm (on the CPU, maybe we can do GPU in `CUDA.jl`),
the divide-and-conquer algorithm ([Smith et al. 2020](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7840081/)), and a "site justified" divide-
and-conquer algorithm.

To use this package, install it with this command in the REPL:
```
using Pkg; Pkg.add(url="https://github.com/jacobusmmsmit/DiscreteVoronoi.jl")
```
and we can use the functionality as such:
```
using DiscreteVoronoi

grid = zeros(10, 10)
nsites = 4
sites = [rand(1:10, 1:10) for _ in 1:nsites]
jfa!(grid, sites)
dacx!(grid, sites)
jdacx!(grid, sites)
```
