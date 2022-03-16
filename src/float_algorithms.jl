using Distributions
using BenchmarkTools
using VoronoiCells
using GeometryBasics
using Plots
using Random

function distance(p, q)
    @fastmath sqrt((p[1] - q[1])^2 + (p[2] - q[2])^2)
end

allequal(itr) = isempty(itr) ? true : all(isequal(first(itr)), itr)

struct Rect{T}
    BL::Tuple{T,T}
    TR::Tuple{T,T}
end
Rect(a::T, b::T) where {T} = Rect((zero(T), zero(T)), (a, b))
get_corners(r::Rect) = (r.BL, (r.BL[1], r.TR[2]), r.TR, (r.TR[1], r.BL[2]))
function quadrants(r::Rect)
    BL, TR = r.BL, r.TR
    C = (BL[1] + TR[1]) / 2, (BL[2] + TR[2]) / 2
    BC = (BL[1] + TR[1]) / 2, BL[2]
    TC = (BL[1] + TR[1]) / 2, TR[2]
    CL = BL[1], (BL[2] + TR[2]) / 2
    CR = TR[1], (BL[2] + TR[2]) / 2
    return (Rect(BL, C), Rect(BC, CR), Rect(CL, TC), Rect(C, TR))
end
area(r::Rect) = reduce(*, (r.TR .- r.BL))
Base.show(io::IO, r::Rect) = print(io, "BL = $(r.BL), TR = $(r.TR)")

function fDAC(rectangle::Rect{T}, sites, site_area = zeros(T, length(sites)); nrecursions::Int = 0, total_recursions::Int = 3) where {T}
    corners = get_corners(rectangle)
    closest_sites = map(corner -> findmin(site -> distance(corner, site), sites)[2], corners)
    if allequal(closest_sites)
        site_area[first(closest_sites)] += area(rectangle)
    elseif nrecursions < total_recursions
        for quadrant in quadrants(rectangle)
            _fDAC(quadrant, sites, site_area; nrecursions = nrecursions + 1, total_recursions = total_recursions)
        end
    end
    # After all recursions are finished, distribute the extra area evenly to all sites.
    # This is better than rescaling as otherwise the more densely packed sites would be
    # strongly underrepresented
    site_area .+= (area(rectangle) - sum(site_area)) / length(sites)
end

function _fDAC(rectangle, sites, site_area; nrecursions, total_recursions)
    corners = get_corners(rectangle)
    closest_sites = map(corner -> findmin(site -> distance(corner, site), sites)[2], corners)
    if allequal(closest_sites)
        site_area[first(closest_sites)] += area(rectangle)
    elseif nrecursions < total_recursions
        for quadrant in quadrants(rectangle)
            _fDAC(quadrant, sites, site_area; nrecursions = nrecursions + 1, total_recursions = total_recursions)
        end
    end
    site_area
end

function main()
    N1, M1 = 0.0, 0.0
    N2, M2 = 10.0, 100.0
    r = Main.Rect((N1, M1), (N2, M2))
    sites = [(rand(Uniform(N1, N2)), rand(Uniform(M1, M2))) for i in 1:100]
    site_area = zeros(length(sites))
    # @time fDAC(r, sites, site_area, total_recursions = 10)
    # @btime fDAC($r, $sites, $site_area, total_recursions = 10) 

    rect = Rectangle(Point2(N1, M1), Point2(N2, M2))
    points = [Point2(site...) for site in sites]
    # @btime fDAC($r, $sites, $site_area, total_recursions = 10)
    # @btime voronoiarea(voronoicells($points, $rect))
    n = 10
    diff = zeros(n)
    for (i, tot) in enumerate(1:n)
        fd, va = fDAC(r, sites, site_area, total_recursions = tot), voronoiarea(voronoicells(points, rect))
        diff[i] = sum(abs.(fd .- va))
    end
    plot(diff)
end

main()
