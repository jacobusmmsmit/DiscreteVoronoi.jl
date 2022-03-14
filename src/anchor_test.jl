## In this file I am testing whether the choice of the centre of the rectangle
## can be improved upon by choosing a different "anchor" from which to calculate
## the closest site. The anchor will hopefully choose a more important closest site
## than the centre in some circumstances such as when sites are heavily skewed in
## one direction.

## Note: the anchor *must* (I think, proof by contradiction exists in my head)
## be inside/on the edge of the rectangle otherwise the anchor may choose a site
## that results in too large circles.

## Geometric median is probably better at eliminating the problem with multimodality
## but harder to calculate

using BenchmarkTools
using Statistics
using Plots
using Random

function distance(p, q)
    @fastmath sqrt((p[1] - q[1])^2 + (p[2] - q[2])^2)
end


"""
Calculate the intersection of two lines defined by 4 points,
2 per line
TODO: remake for line segments using wikipedia article
https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection#Given_two_points_on_each_line_segment
"""
function intersection_two_lines((x1, y1), (x2, y2), (x3, y3), (x4, y4))
    # Numerators
    x_num = (x1 * y2 - x2 * y1) * (x3 - x4) - (x3 * y4 - x4 * y3) * (x1 - x2)
    y_num = (x1 * y2 - x2 * y1) * (y3 - y4) - (x3 * y4 - x4 * y3) * (y1 - y2)
    # Denominator
    denom = (x1 - x2) * (y3 - y4) - (x3 - x4) * (y1 - y2)
    return (x_num, y_num) ./ denom
end


"""
    intersection_two_line_segments((x1, y1), (x2, y2), (x3, y3), (x4, y4))
Returns the intersection point of the lines defined by the points and a tuple
of Boolean values that say if the intersection is within the segment or outside.
"""
function intersection_two_line_segments((x1, y1), (x2, y2), (x3, y3), (x4, y4))
    t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / ((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4))
    u = ((x1 - x3) * (y1 - y2) - (y1 - y3) * (x1 - x2)) / ((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4))
    # 0.0 ≤ t ≤ 1.0 ⟹ intersection within first line segment
    # 0.0 ≤ u ≤ 1.0 ⟹ intersection second first line segment
    return (x1 + t * (x2 - x1), y1 + t * (y2 - y1)), (t, u)
end

function closest_anchor_in_rectangle(rect, anchor)
    (T, L), (N, M) = rect
    centre = (T + N / 2, L + M / 2)
    corners = (T, L), (T + N, L), (T + N, L + M), (T, L + M)
    j = 4
    for i in 1:4
        intersection_point, (t, u) = intersection_two_line_segments(
            corners[i], corners[j],
            centre, anchor
        )
        # if the centre to anchor (C2A) direction collides with rectangle edge segment
        if 0.0 ≤ t ≤ 1.0 && u ≥ 0
            # if C2A intersects rectangle edge segment
            if u ≤ 1.0
                inner_anchor = intersection_point
            else # i.e. the anchor is already inside
                inner_anchor = anchor
            end
            return inner_anchor
        end
        j = i
    end
end


# A = [(i, i + 2) for i in 1:5]
# reduce(.+, A)
# @time reduce(.+, A) ./ length(A)
# @btime reduce(.+, $A) ./ length($A)


begin
    Random.seed!(5)
    scene = plot(xlims = (0, 1), ylims = (0, 1))
    sites = [rand(2) for i in 1:10]
    N, M = 0.5, 0.5
    corners = (0.0, 0.0), (N, 0.0), (N, M), (0.0, M)
    rectangle_shape = Shape([0.0, N, N, 0.0], [0.0, 0.0, M, M])
    plot!(rectangle_shape, alpha = 0.5, label = "Rectangle")
    scatter!(first.(sites), last.(sites), label = "Sites", shape = :star)
    anchor = reduce(.+, sites) ./ length(sites)
    scatter!([centre[1]], [centre[2]], label = "Centre")
    plot!(scene, size = (500, 400), aspect_ratio = 1.0, legend = :outerright)
end

cl_ia = sites[findmin(site -> distance(inner_anchor, site), sites)[2]]
cl_ce = sites[findmin(site -> distance(centre, site), sites)[2]]

radii_cl_ia = map(corner -> distance(cl_ia, corner), corners)
radii_cl_ce = map(corner -> distance(cl_ce, corner), corners)
function circleShape(centre, r)
    x, y = centre
    θ = range(start = 0, stop = 2π, length = 500)
    x .+ r * sin.(θ), y .+ r * cos.(θ)
end

# scene_ce = scene
# for i in 1:4    
#     circ = circleShape(corners[i], radii_cl_ce[i])
#     plot!(scene_ce, circ, fillalpha = 0.2)
# end
# scene_ce
scene_ia = deepcopy(scene)
begin
    for i in 1:4
        plot!(scene_ia, circleShape(corners[i], radii_cl_ia[i]), seriestype = [:shape,], lw = 0.5, c = :yellow, legend = false, fillalpha = 0.1)
    end
    scatter!(scene_ia, [anchor[1]], [anchor[2]], label = "Anchor")
    plot!(scene_ia, [centre[1], anchor[1]], [centre[2], anchor[2]], label = "", lw = 2.0, lc = :black, ls = :dash)
    rect = ((0, 0), (N, M))
    inner_anchor = closest_anchor_in_rectangle(rect, anchor)
    plot!(scene_ia, [cl_ia[1], inner_anchor[1]], [cl_ia[2], inner_anchor[2]], label = "", lw = 2.0, lc = :black)
    scatter!(scene_ia, [inner_anchor[1]], [inner_anchor[2]], label = "Inner Anchor")
end

scene_ce = deepcopy(scene)
begin
    for i in 1:4
        plot!(scene_ce, circleShape(corners[i], radii_cl_ce[i]), seriestype = [:shape,], lw = 0.5, c = :red, legend = false, fillalpha = 0.1)
    end
    plot!(scene_ce, [cl_ce[1], centre[1]], [cl_ce[2], centre[2]], label = "", lw = 2.0, lc = :black)
end





plot(scene_ia, scene_ce, size = (800, 500), xlims = (-0.5, 1), ylims = (-0.5, 1))