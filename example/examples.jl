include("../src/DiscreteVoronoi.jl")
using .DiscreteVoronoi
using Random
using BenchmarkTools

function rand_sites(::Type{Int}, N, M, K)
    idx = collect(Iterators.product(1:N, 1:M))
    shuffle!(idx)
    idx[1:K]
end


Random.seed!(42)
grid = zeros(Int, rand(1:100), rand(1:100))
sites = rand_sites(Int, size(grid)..., 10)
jfa!(grid, sites)
@show grid

Random.seed!(42)
grid = zeros(Int, rand(1:100), rand(1:100))
sites = rand_sites(Int, size(grid)..., 10)
dac!(grid, sites)
@show grid

Random.seed!(42)
grid = zeros(Int, rand(1:100), rand(1:100))
sites = collect(enumerate(rand_sites(Int, size(grid)..., 10)))
jdac!(grid, sites, jdac_aux1a!)
@show grid

#= L0(A, B) = sum(A .!= B)

begin
    # time_naive = Float64[]
    time_jfa = Float64[]
    # time_dac = Float64[]
    time_jdac = Float64[]

    nruns = 3
    nsites_list = 100:100:1000
    for nsites in nsites_list
        N, M = 1000, 1000
        @show N, M, nsites
        for times in 1:nruns
            local sites = [(rand(1:N), rand(1:M)) for i in 1:nsites]

            #= grid = zeros(Int, N, M)
            push!(time_naive, @belapsed naive_voronoi(CartesianIndices((1:$N, 1:$M)), $sites))
            # grid1 = @time naive_voronoi(CartesianIndices((1:N, 1:M)), map(site -> site[2], sites)) =#

            grid2 = zeros(Int, N, M)
            push!(time_jfa, @belapsed jfa!($grid2, $sites))

            #= grid3 = zeros(Int, N, M)
            push!(time_dac, @belapsed dacx!($grid3, $sites, 0)) =#

            grid4 = zeros(Int, N, M)
            sites2 = collect(enumerate(sites))
            push!(time_jdac, @belapsed jdac!($grid4, $sites2, jdac_aux2c!))
            # println("N. errors in JFA: ", L0(grid3, grid2))
            # println("N. errors in jdac: ", L0(grid3, grid4))
        end
    end
    # time_naive = reshape(time_naive, nruns, length(nsites_list))
    time_jfa = reshape(time_jfa, nruns, length(nsites_list))
    # time_dac = reshape(time_dac, nruns, length(nsites_list))
    time_jdac = reshape(time_jdac, nruns, length(nsites_list))
end

using Plots

begin
    # median_naive = map(median, eachcol(time_naive))
    median_jfa = map(median, eachcol(time_jfa))
    # median_dac = map(median, eachcol(time_dac))
    median_jdac = map(median, eachcol(time_jdac))

    begin
        times_plot = plot(xlabel = "Number of Sites", ylabel = "Median Runtime", title = "Runtime varying number of sites", legend = :topleft)
        # plot!(nsites_list, median_naive, label = "Naive")
        plot!(nsites_list, median_jfa, label = "Jump-flood")
        # plot!(nsites_list, median_dac, label = "DAC")
        plot!(nsites_list, median_jdac, label = "JDAC")

        times_plot_fast = plot(xlabel = "Number of Sites", ylabel = "Median Runtime", title = "Runtime varying number of sites", legend = :left)
        # plot!(nsites_list, median_naive, label = "Naive")
        plot!(nsites_list, median_jfa, label = "Jump-flood")
        # plot!(nsites_list, median_dac, label = "DAC")
        plot!(nsites_list, median_jdac, label = "JDAC")
    end

    # plot(times_plot, times_plot_fast, size = (1000, 500), margin = 3mm)
    plot(times_plot, times_plot_fast, size = (1000, 500))
end =#

