using DiscreteVoronoi
using Distances
using BenchmarkTools
using Random

for n in [100, 1000]
    for s in [isqrt(n), n, n * isqrt(n), n * n]
        @show n, s

        println("jfa_voronoi!")
        @btime jfa_voronoi!(grid, sites) setup=(
            Random.seed!(42);
            sites = rand_points(Int, $n, $n, $s);
            grid = zeros(Int, $n, $n)) evals=1
             

        println("jfa_voronoi_parallel!")
        @btime jfa_voronoi_parallel!(grid, sites) setup=(
            Random.seed!(42);
            sites = rand_points(Int, $n, $n, $s);
            grid = zeros(Int, $n, $n)) evals=1
    
        if s <= n
            println("dac_voronoi!")
            for site_find in [original_site_find, center_site_find]
                @show site_find
                @btime dac_voronoi!(grid, sites, site_find) setup=(
                    Random.seed!(42);
                    sites = rand_sites(Int, $n, $n, $s);
                    grid = preset_voronoi!(zeros(Int, $n, $n), sites);
                    site_find = $site_find) evals=1
            end
        end

        println("redac_voronoi!")
        for method in [:filter, :partition]
            @show method
            for site_find in [no_site_find, original_site_find, center_site_find]
                @show site_find
                for site_filter in [naive_site_filter, center_site_filter, anchor_site_filter, corner_site_filter]
                    @show site_filter
                    if method === :filter
                        @btime redac_voronoi!(Val(:filter), grid, sites, site_find, site_filter, euclidean) setup=(
                            Random.seed!(42);
                            sites = rand_sites(Int, $n, $n, $s);
                            grid = preset_voronoi!(zeros(Int, $n, $n), sites);
                            site_find = $site_find;
                            site_filter = $site_filter) evals=1
                    elseif method === :partition
                        @btime redac_voronoi!(Val(:partition), grid, sites, site_find, site_filter, euclidean) setup=(
                            Random.seed!(42);
                            sites = rand_sites(Int, $n, $n, $s);
                            grid = preset_voronoi!(zeros(Int, $n, $n), sites);
                            site_find = $site_find;
                            site_filter = $site_filter) evals=1
                    end
                end
            end
        end
    end
end

