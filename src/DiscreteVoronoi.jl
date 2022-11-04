module DiscreteVoronoi

#TODO: Implement hybrid algorithms for certain gridsize and number of sites.
#TODO: Add nice UnicodePlot recipe :)

export find_closest_site, find_closest_site!, get_corners, get_quadrants, label_voronoi_grid, voronoi_equality # Helper functions
export naive_voronoi!, jfa_voronoi!, dac_voronoi!, redac_voronoi!, redac_voronoi_es! # Core functionality

using LinearAlgebra: norm
using Random: shuffle
using StaticArrays
using Distances: euclidean

const Coord = SVector{2,Int}

# Completely non-exported files with core definitions
include("EarlyStopper.jl") # Because every package needs at least one struct...
include("elimination_methods.jl")

include("helper_functions.jl")
include("core_algorithms.jl")

end # module
