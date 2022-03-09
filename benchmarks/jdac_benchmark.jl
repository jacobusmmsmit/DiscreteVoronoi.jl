include("../src/DiscreteVoronoi.jl")
using .DiscreteVoronoi
using Random
using BenchmarkTools


function get_seeds(N, M, K)
    idx = [(n, m) for n in 1:N, m in 1:M]
    shuffle!(idx)
    idx[1:K]
end

Random.seed!(42)
grid = zeros(Int, rand(1:100), rand(1:100))
seeds = get_seeds(size(grid)..., 10)
jfa!(grid, seeds)
@show grid

Random.seed!(42)
grid = zeros(Int, rand(1:100), rand(1:100))
seeds = get_seeds(size(grid)..., 10)
dac!(grid, seeds, 0)
@show grid

Random.seed!(42)
grid = zeros(Int, rand(1:100), rand(1:100))
seeds = collect(enumerate(get_seeds(size(grid)..., 10)))
jdac!(grid, seeds, jdac_aux1!, 0)
@show grid

println("jfa!")
@btime jfa!(grid, seeds) setup=(
    Random.seed!(42);
    grid = zeros(Int, 100, 40);
    seeds = get_seeds(size(grid)..., 30))

println("dac!")
@btime dac!(grid, seeds) setup=(
    Random.seed!(42);
    grid = zeros(Int, 100, 40);
    seeds = get_seeds(size(grid)..., 30))

println("jdac!")
for aux! in [jdac_aux1!, jdac_aux2!, jdac_aux3!, jdac_aux4!, jdac_aux5!, jdac_aux6!]
    @show aux!
    @btime jdac!(grid, seeds, $aux!) setup=(
        Random.seed!(42);
        grid = zeros(Int, 100, 40);
        seeds = collect(enumerate(get_seeds(size(grid)..., 30))))
end =#

for n in 10:10:50
    for s in [isqrt(n), n, n * isqrt(n), n * n]
        @show n, s
        println("jfa")
        @btime jfa!(grid, seeds) setup=(
            Random.seed!(42);
            grid = zeros(Int, $n, $n);
            seeds = get_seeds(size(grid)..., $s))

        #= println("dac!")
        @btime dac!(grid, seeds) setup=(
            Random.seed!(42);
            grid = zeros(Int32, $n, $n);
            seeds = get_seeds(size(grid)..., $s)) =#

        println("jdac")
        for aux! in [jdac_aux1!, jdac_aux2!, jdac_aux3!, jdac_aux4!, jdac_aux5!, jdac_aux6!]
            @show aux!
            @btime jdacx!(grid, seeds, $aux!) setup=(
                Random.seed!(42);
                grid = zeros(Int, $n, $n);
                seeds = collect(enumerate(get_seeds(size(grid)..., $s))))
        end
    end
end 

#= using Profile
Profile.clear()
grid = zeros(Int, 1000, 1000)
seeds = collect(enumerate(get_seeds(size(grid)..., 1000 * isqrt(1000))))
jdac!(grid, seeds, jdac_aux6!, 0)
grid = zeros(Int, 1000, 1000)
seeds = collect(enumerate(get_seeds(size(grid)..., 1000 * isqrt(1000))))
@profile jdac!(grid, seeds, jdac_aux6!, 0);
# Profile.print()
Juno.profiler() =#

# L0(A, B) = sum(A .!= B)

#= begin
    # time_naive = Float64[]
    time_jfa = Float64[]
    # time_dac = Float64[]
    time_jdac = Float64[]

    nruns = 3
    nseeds_list = 100:100:1000
    for nseeds in nseeds_list
        N, M = 1000, 1000
        @show N, M, nseeds
        for times in 1:nruns
            local seeds = [(rand(1:N), rand(1:M)) for i in 1:nseeds]

            #= grid = zeros(Int, N, M)
            push!(time_naive, @belapsed naive_voronoi(CartesianIndices((1:$N, 1:$M)), $seeds))
            # grid1 = @time naive_voronoi(CartesianIndices((1:N, 1:M)), map(seed -> seed[2], seeds)) =#

            grid2 = zeros(Int, N, M)
            push!(time_jfa, @belapsed jfa!($grid2, $seeds))

            #= grid3 = zeros(Int, N, M)
            push!(time_dac, @belapsed dac!($grid3, $seeds, 0)) =#

            grid4 = zeros(Int, N, M)
            seeds2 = collect(enumerate(seeds))
            push!(time_jdac, @belapsed jdac!($grid4, $seeds2, 0))
            # println("N. errors in JFA: ", L0(grid3, grid2))
            # println("N. errors in jdac: ", L0(grid3, grid4))
        end
    end
    # time_naive = reshape(time_naive, nruns, length(nseeds_list))
    time_jfa = reshape(time_jfa, nruns, length(nseeds_list))
    # time_dac = reshape(time_dac, nruns, length(nseeds_list))
    time_jdac = reshape(time_jdac, nruns, length(nseeds_list))
end

using Plots

begin
    # median_naive = map(median, eachcol(time_naive))
    median_jfa = map(median, eachcol(time_jfa))
    # median_dac = map(median, eachcol(time_dac))
    median_jdac = map(median, eachcol(time_jdac))

    begin
        times_plot = plot(xlabel = "Number of Sites", ylabel = "Median Runtime", title = "Runtime varying number of sites", legend = :topleft)
        # plot!(nseeds_list, median_naive, label = "Naive")
        plot!(nseeds_list, median_jfa, label = "Jump-flood")
        # plot!(nseeds_list, median_dac, label = "DAC")
        plot!(nseeds_list, median_jdac, label = "JDAC")

        times_plot_fast = plot(xlabel = "Number of Sites", ylabel = "Median Runtime", title = "Runtime varying number of sites", legend = :left)
        # plot!(nseeds_list, median_naive, label = "Naive")
        plot!(nseeds_list, median_jfa, label = "Jump-flood")
        # plot!(nseeds_list, median_dac, label = "DAC")
        plot!(nseeds_list, median_jdac, label = "JDAC")
    end

    # plot(times_plot, times_plot_fast, size = (1000, 500), margin = 3mm)
    plot(times_plot, times_plot_fast, size = (1000, 500))
end =#

