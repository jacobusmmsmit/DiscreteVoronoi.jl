using Random, BenchmarkTools, TimerOutputs

function distance(p, q)
    @fastmath sqrt((p[1] - q[1])^2 + (p[2] - q[2])^2)
end

function naive(grid, seeds)
    map(cell -> findmin(seed -> distance(cell, seed), seeds)[2], grid)
end

function jfa!(grid, seeds)
    for (color, (x, y)) in enumerate(seeds)
        grid[x, y] = convert(eltype(grid), color)
    end
    k = max(size(grid)...)
    while k > 1
        k = k ÷ 2 + k % 2
        @inbounds for y in 1:size(grid, 2), x in 1:size(grid, 1)
            for j in (-k, 0, k), i in (-k, 0, k)
                checkbounds(Bool, grid, x + i, y + j) || continue
                colorq = grid[x + i, y + j]
                colorq !== 0 || continue
                colorp = grid[x, y]
                if colorp === 0
                    grid[x, y] = colorq
                elseif distance(seeds[colorp], (x, y)) > distance(seeds[colorq], (x, y))
                    grid[x, y] = colorq
                end
            end
        end
    end
    grid
end

function dac_aux!(grid, seeds, depth)
    if all(.>(0), size(grid)) && any(==(0), grid)
        dac!(grid, seeds, depth - 1)
    end
    grid
end
function dac!(grid, seeds, depth=1)
    N, M = size(grid)
    corners = ((1, 1), (N, 1), (1, M), (N, M))
    nearest_colors = map(corner -> findmin(seed -> distance(corner, seed), seeds)[2], corners)
    color = nearest_colors[1]
    if all(nearest_color -> nearest_color == color, nearest_colors[2:end])
        grid .= convert(eltype(grid), color)
    else
        Nd = N ÷ 2 + N % 2
        Md = M ÷ 2 + M % 2
        offsets = ((0, 0), (Nd, 0), (0, Md), (Nd, Md))
        tl_grid = @view grid[1:Nd, 1:Md]
        bl_grid = @view grid[Nd+1:N, 1:Md]
        tr_grid = @view grid[1:Nd, Md+1:M]
        br_grid = @view grid[Nd+1:N, Md+1:M]
        sub_grids = (tl_grid, bl_grid, tr_grid, br_grid)
        if depth > 0
            Threads.@threads for i in 1:4
                thread_seeds = map(seed -> seed .- offsets[i], seeds)
                dac_aux!(sub_grids[i], thread_seeds, depth)
            end
        else
            for i in 1:4
                replace!(seed -> seed .- offsets[i], seeds)
                dac_aux!(sub_grids[i], seeds, depth)
                replace!(seed -> seed .+ offsets[i], seeds)
            end
        end
    end
    grid
end
function dacx!(grid, seeds, depth=1)
    for (color, (x, y)) in enumerate(seeds)
        grid[x, y] = convert(eltype(grid), color)
    end
    dac!(grid, seeds, depth)
end

mutable struct Stack
    depth::Int
    dists::Vector{Vector{Float64}}
    seeds::Vector{Vector{Tuple{Int, Tuple{Int, Int}}}}
    Stack() = new(0, Vector{Vector{Float64}}(), Vector{Vector{Tuple{Int, Tuple{Int, Int}}}}())
end

function _push_empty!(stack)
    stack.depth += 1
    if stack.depth > length(stack.dists)
        push!(stack.dists, Vector{Float64}())
    end
    if stack.depth > length(stack.seeds)
        push!(stack.seeds, Vector{Tuple{Int, Tuple{Int, Int}}}())
    end
    stack.depth
end
function _append_dist!(stack, len, dist)
    if len <= length(stack.dists[stack.depth])
        @inbounds stack.dists[stack.depth][len] = dist
    else
        @inbounds push!(stack.dists[stack.depth], dist)
    end
end
function _resize_dists!(stack, len)
    @inbounds resize!(stack.dists[stack.depth], len)
end
function _get_dists(stack)
    @inbounds stack.dists[stack.depth]
end
function _append_seed!(stack, len, seed)
    if len <= length(stack.seeds[stack.depth])
        @inbounds stack.seeds[stack.depth][len] = seed
    else
        @inbounds push!(stack.seeds[stack.depth], seed)
    end
end
function _resize_seeds!(stack, len)
    @inbounds resize!(stack.seeds[stack.depth], len)
end
function _get_seeds(stack)
    @inbounds stack.seeds[stack.depth]
end
function _pop!(stack)
    stack.depth -= 1
end

function _fill_dists!(center, seeds, stack)
    #= _push_empty!(stack)
    len = 0
    for seed in seeds
        len += 1
        _append_dist!(stack, len, distance(center, seed[2]))
    end
    _resize_dists!(stack, len) =#
    _resize_dists!(stack, length(seeds))
    copyto!(_get_dists(stack), distance(center, seed[2]) for seed in seeds)
    # _get_dists(stack)
end

function _fill_seeds!(max_dist, seeds, stack)
    #= len = 0
    for (dist, seed) in zip(_get_dists(stack), seeds)
        if dist <= max_dist
            len += 1
            _append_seed!(stack, len, seed)
        end
    end
    _resize_seeds!(stack, len) =#
    len = sum(1 for dist in _get_dists(stack) if dist < max_dist)
    _resize_seeds!(stack, len)
    copyto!(_get_seeds(stack),
        seed for (dist, seed) in zip(_get_dists(stack), seeds) if dist < max_dist)
    # _get_seeds(stack)
end

function jdac_aux1!(grid, seeds, depth, stack)
    if all(.>(0), size(grid)) && any(==(0), grid)
        center = size(grid) .÷ 2 .+ size(grid) .% 2
        max_dist = findmin(seed -> distance(center, seed[2]), seeds)[1] + distance((0, 0), size(grid)) + 1
        stack_seeds = filter(seed -> distance(center, seed[2]) <= max_dist, seeds)
        jdac!(grid, stack_seeds, jdac_aux1!, depth - 1, stack)
    end
    grid
end
function jdac_aux2!(grid, seeds, depth, stack)
    if all(.>(0), size(grid)) && any(==(0), grid)
        center = size(grid) .÷ 2 .+ size(grid) .% 2
        dists = map(seed -> distance(center, seed[2]), seeds)
        max_dist = findmin(dists)[1] + distance((0, 0), size(grid)) + 1
        stack_seeds = [seed for (dist, seed) in zip(dists, seeds) if dist <= max_dist]
        jdac!(grid, stack_seeds, jdac_aux2!, depth - 1, stack)
    end
    grid
end
function jdac_aux3!(grid, seeds, depth, stack)
    if all(.>(0), size(grid)) && any(==(0), grid)
        center = size(grid) .÷ 2 .+ size(grid) .% 2
        _push_empty!(stack)
        _fill_dists!(center, seeds, stack)
        max_dist = findmin(_get_dists(stack))[1] + distance((0, 0), size(grid)) + 1
        _fill_seeds!(max_dist, seeds, stack)
        jdac!(grid, _get_seeds(stack), jdac_aux3!, depth - 1, stack)
        _pop!(stack)
    end
    grid
end
function jdac_aux4!(grid, seeds, depth, stack)
    if all(.>(0), size(grid)) && any(==(0), grid)
        N, M = size(grid)
        center = size(grid) .÷ 2 .+ size(grid) .% 2
        min_seed = seeds[findmin(seed -> distance(center, seed[2]), seeds)[2]][2]
        corners = ((1, 1), (N, 1), (1, M), (N, M))
        max_dist = maximum(corner -> distance(corner, min_seed), corners) + distance((0, 0), center) + 1
        stack_seeds = filter(seed -> distance(center, seed[2]) <= max_dist, seeds)
        jdac!(grid, stack_seeds, jdac_aux4!, depth - 1, stack)
    end
    grid
end
function jdac_aux5!(grid, seeds, depth, stack)
    if all(.>(0), size(grid)) && any(==(0), grid)
        N, M = size(grid)
        center = size(grid) .÷ 2 .+ size(grid) .% 2
        dists = map(seed -> distance(center, seed[2]), seeds)
        min_seed = seeds[findmin(dists)[2]][2]
        corners = ((1, 1), (N, 1), (1, M), (N, M))
        max_dist = maximum(corner -> distance(corner, min_seed), corners) + distance((0, 0), center) + 1
        stack_seeds = [seed for (dist, seed) in zip(dists, seeds) if dist <= max_dist]
        jdac!(grid, stack_seeds, jdac_aux5!, depth - 1, stack)
    end
    grid
end
function jdac_aux6!(grid, seeds, depth, stack)
    if all(.>(0), size(grid)) && any(==(0), grid)
        N, M = size(grid)
        center = size(grid) .÷ 2 .+ size(grid) .% 2
        _push_empty!(stack)
        _fill_dists!(center, seeds, stack)
        min_seed = seeds[findmin(_get_dists(stack))[2]]
        corners = ((1, 1), (N, 1), (1, M), (N, M))
        max_dist = maximum(corner -> distance(corner, min_seed[2]), corners) + distance((0, 0), center) + 1
        _fill_seeds!(max_dist, seeds, stack)
        jdac!(grid, _get_seeds(stack), jdac_aux6!, depth - 1, stack)
        _pop!(stack)
    end
    grid
end
function jdac!(grid, seeds, aux!, depth=1, stack=Stack())
    N, M = size(grid)
    corners = ((1, 1), (N, 1), (1, M), (N, M))
    nearest_colors = map(corner -> findmin(seed -> distance(corner, seed[2]), seeds)[2], corners)
    color = seeds[nearest_colors[1]][1]
    if all(nearest_color -> seeds[nearest_color][1] == color, nearest_colors[2:end])
        grid .= convert(eltype(grid), color)
    else
        Nd = N ÷ 2 + N % 2
        Md = M ÷ 2 + M % 2
        offsets = ((0, 0), (Nd, 0), (0, Md), (Nd, Md))
        tl_grid = @view grid[1:Nd, 1:Md]
        bl_grid = @view grid[Nd+1:N, 1:Md]
        tr_grid = @view grid[1:Nd, Md+1:M]
        br_grid = @view grid[Nd+1:N, Md+1:M]
        sub_grids = (tl_grid, bl_grid, tr_grid, br_grid)
        if depth > 0
            Threads.@threads for i in 1:4
                thread_seeds = map(seed -> (seed[1], seed[2] .- offsets[i]), seeds)
                aux!(sub_grids[i], thread_seeds, depth, Stack())
            end
        else
            for i in 1:4
                replace!(seed -> (seed[1], seed[2] .- offsets[i]), seeds)
                aux!(sub_grids[i], seeds, depth, stack)
                replace!(seed -> (seed[1], seed[2] .+ offsets[i]), seeds)
            end
        end
    end
    grid
end
function jdacx!(grid, seeds, aux!, depth=1, stack=Stack())
    for (color, (x, y)) in seeds
        grid[x, y] = convert(eltype(grid), color)
    end
    jdac!(grid, seeds, aux!, depth, stack)
end

function get_seeds(N, M, K)
    idx = [(n, m) for n in 1:N, m in 1:M]
    shuffle!(idx)
    idx[1:K]
end

#= Random.seed!(42)
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
@show grid =#

#= Random.seed!(42)
for i in 1:1000
    N, M = rand(1:1000, 2)
    local seeds = collect(enumerate(get_seeds(N, M, rand(1:100))))
    @show N, M, length(seeds)

    println("naive")
    grid1 = @time naive(CartesianIndices((1:N, 1:M)), map(seed -> seed[2], seeds))

    println("jfa!")
    grid2 = zeros(Int, N, M)
    @time jfa!(grid2, map(seed -> seed[2], seeds))
    # @assert grid2 == grid1

    println("dac!")
    grid3 = zeros(Int, N, M)
    @time dac!(grid3, map(seed -> seed[2], seeds))
    @assert grid3 == grid1

    println("jdac!")
    for aux! in [jdac_aux1!, jdac_aux2!, jdac_aux3!, jdac_aux4!, jdac_aux5!, jdac_aux6!]
        @show aux!
        grid4 = zeros(Int, N, M)
        @time jdac!(grid4, seeds, aux!, 0)
        @assert grid4 == grid3
    end
end =#

#= println("jfa!")
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

for n in 1000:1000:5000
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
            push!(time_naive, @belapsed naive(CartesianIndices((1:$N, 1:$M)), $seeds))
            # grid1 = @time naive(CartesianIndices((1:N, 1:M)), map(seed -> seed[2], seeds)) =#

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
