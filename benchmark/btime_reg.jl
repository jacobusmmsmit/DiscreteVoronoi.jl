using DiscreteVoronoi
using Distances
using BenchmarkTools

function reg_points(::Type{Int}, N, M, K)
    [Iterators.product(2:3:N, 2:3:M)...]
end

function reg_sites(::Type{Int}, N, M, K)
    [(color, point) for (color, point) in enumerate(reg_points(Int, N, M, K))]
end

for s in [1, 2, 4, 8, 16, 32, 64]
    n = 3 * s
    @show n, s * s

    println("jfa_voronoi!")
    @btime jfa_voronoi!(grid, sites) setup=(
        sites = reg_points(Int, $n, $n, $s);
        grid = zeros(Int, $n, $n)) evals=1 samples=1
            

    println("dac_voronoi!")
    for site_find in [original_site_find, center_site_find]
        @show site_find
        @btime dac_voronoi!(grid, sites, site_find, euclidean, 0) setup=(
            sites = reg_sites(Int, $n, $n, $s);
            grid = preset_voronoi!(zeros(Int, $n, $n), sites);
            site_find = $site_find) evals=1 samples=1
    end

    println("redac_voronoi!")
    for site_find in [no_site_find, original_site_find, center_site_find]
        @show site_find
        for site_filter in [original_site_filter, center_site_filter, anchor_site_filter, corner_site_filter]
            @show site_filter
            @btime redac_voronoi!(grid, sites, site_find, site_filter, euclidean, 0) setup=(
                sites = reg_sites(Int, $n, $n, $s);
                grid = preset_voronoi!(zeros(Int, $n, $n), sites);
                site_find = $site_find;
                site_filter = $site_filter) evals=1 samples=1
        end
    end
end

