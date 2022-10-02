module DiscreteVoronoi

#TODO: Sort sites for better cache locality?
#TODO: ] remove BenchmarkTools PProf Plots
#TODO: Implement hybrid DAC+Naive for certain size

export find_closest_site # Helper functions
export naive_voronoi!, jfa_voronoi!, dac_voronoi! # Traditional approaches

# Novel work (Joerg Walter and Jacobus Smit):
export exact_condition, centre_anchor_condition, exact_elimination, centre_anchor_elimination # Elimination methods
export redac_voronoi! # Reduce-Divide-and-Conquer

using LinearAlgebra: norm
using StaticArrays

include("helper_functions.jl")
include("traditional_approaches.jl")
include("elimination_methods.jl")
include("EarlyStopper.jl")
include("redac_earlystopper.jl")

end # module
