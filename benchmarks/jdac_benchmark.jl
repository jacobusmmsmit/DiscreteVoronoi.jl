include("../src/DiscreteVoronoi.jl")
using .DiscreteVoronoi
using Random
using BenchmarkTools


function get_sites(N, M, K)
    idx = [(n, m) for n in 1:N, m in 1:M]
    shuffle!(idx)
    idx[1:K]
end

Random.seed!(42)
grid = zeros(Int, rand(1:100), rand(1:100))
sites = get_sites(size(grid)..., 10)
jfa!(grid, sites)
@show grid

Random.seed!(42)
grid = zeros(Int, rand(1:100), rand(1:100))
sites = get_sites(size(grid)..., 10)
dac!(grid, sites, 0)
@show grid

Random.seed!(42)
grid = zeros(Int, rand(1:100), rand(1:100))
sites = collect(enumerate(get_sites(size(grid)..., 10)))
jdac!(grid, sites, jdac_aux1!, 0)
@show grid

println("jfa!")
@btime jfa!(grid, sites) setup=(
    Random.seed!(42);
    grid = zeros(Int, 100, 40);
    sites = get_sites(size(grid)..., 30))

println("dac!")
@btime dac!(grid, sites) setup=(
    Random.seed!(42);
    grid = zeros(Int, 100, 40);
    sites = get_sites(size(grid)..., 30))

println("jdac!")
for aux! in [jdac_aux1!, jdac_aux2!, jdac_aux3!, jdac_aux4!, jdac_aux5!, jdac_aux6!]
    @show aux!
    @btime jdac!(grid, sites, $aux!) setup=(
        Random.seed!(42);
        grid = zeros(Int, 100, 40);
        sites = collect(enumerate(get_sites(size(grid)..., 30))))
end =#

for n in 10:10:50
    for s in [isqrt(n), n, n * isqrt(n), n * n]
        @show n, s
        println("jfa")
        @btime jfa!(grid, sites) setup=(
            Random.seed!(42);
            grid = zeros(Int, $n, $n);
            sites = get_sites(size(grid)..., $s))

        #= println("dac!")
        @btime dac!(grid, sites) setup=(
            Random.seed!(42);
            grid = zeros(Int32, $n, $n);
            sites = get_sites(size(grid)..., $s)) =#

        println("jdac")
        for aux! in [jdac_aux1!, jdac_aux2!, jdac_aux3!, jdac_aux4!, jdac_aux5!, jdac_aux6!]
            @show aux!
            @btime jdacx!(grid, sites, $aux!) setup=(
                Random.seed!(42);
                grid = zeros(Int, $n, $n);
                sites = collect(enumerate(get_sites(size(grid)..., $s))))
        end
    end
end 

#= using Profile
Profile.clear()
grid = zeros(Int, 1000, 1000)
sites = collect(enumerate(get_sites(size(grid)..., 1000 * isqrt(1000))))
jdac!(grid, sites, jdac_aux6!, 0)
grid = zeros(Int, 1000, 1000)
sites = collect(enumerate(get_sites(size(grid)..., 1000 * isqrt(1000))))
@profile jdac!(grid, sites, jdac_aux6!, 0);
# Profile.print()
Juno.profiler() =#

# L0(A, B) = sum(A .!= B)

#= begin
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
            push!(time_dac, @belapsed dac!($grid3, $sites, 0)) =#

            grid4 = zeros(Int, N, M)
            sites2 = collect(enumerate(sites))
            push!(time_jdac, @belapsed jdac!($grid4, $sites2, 0))
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

