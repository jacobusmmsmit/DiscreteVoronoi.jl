module DiscreteVoronoi

#TODO: Fix bug with redac_voronoi! not returning the right answer
#TODO: Implement hybrid algorithms for certain gridsize and number of sites

export find_closest_site # Helper functions
export naive_voronoi!, jfa_voronoi!, dac_voronoi! # Traditional approaches

# Novel work (Joerg Walter and Jacobus Smit):
export exact_condition, centre_anchor_condition, exact_elimination, centre_anchor_elimination # Elimination methods
export redac_voronoi! # Reduce-Divide-and-Conquer

using LinearAlgebra: norm
using Random: shuffle
using StaticArrays

include("helper_functions.jl")
include("traditional_approaches.jl")
include("elimination_methods.jl")
include("EarlyStopper.jl")
include("redac_earlystopper.jl")

end # module
