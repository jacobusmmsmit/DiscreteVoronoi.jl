using DiscreteVoronoi
using Random
using BenchmarkTools

const SUITE = BenchmarkGroup()

for n in [10, 100, 1000]
    SUITE[string("grid ", n, "x", n)] = BenchmarkGroup()
    for s in [isqrt(n), n, n * isqrt(n), n * n]
        SUITE[string("grid ", n, "x", n)][string(s, " sites")] = BenchmarkGroup()

        SUITE[string("grid ", n, "x", n)][string(s, " sites")]["jfa_voronoi!"] = @benchmarkable jfa_voronoi!(grid, points) setup=(
            Random.seed!(42);            
            points = rand_points(Int, $n, $n, $s);
            grid = zeros(Int, $n, $n)) evals=1

        for site_filter in [center_site_filter, anchor_site_filter, corner_site_filter]
            SUITE[string("grid ", n, "x", n)][string(s, " sites")][site_find][site_filter] = @benchmarkable redac_voronoi!(Val(:filter), grid, sites, site_filter) setup=(
                Random.seed!(42);
                sites = rand_sites(Int, $n, $n, $s);
                grid = preset_voronoi!(zeros(Int, $n, $n), sites);
                site_filter = $site_filter) evals=1
        end
    end
end