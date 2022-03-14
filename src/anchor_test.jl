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
using LinearAlgebra: ⋅
using Distributions

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

function closest_anchor_from_centre(rect, anchor)
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


"""
    closest_point(point, segment_start, segment_end)
Calculate the closest point from a point to a line segment defined by two points
Arguments   
===
+ `point`: a euclidean co-ordinate,
+ `segment_start`: the start of the line segment
+ `segment_end`: the end of the line segment
Result
===
An elementwise convex combination of `segment_start` and `segment_end` closest to `point`
"""
function closest_point(point, segment_start, segment_end)
    segment = segment_end .- segment_start
    t0 = clamp(((point .- segment_start) ⋅ segment) / (segment ⋅ segment), 0, 1)
    segment_start .+ t0 .* segment
end
begin
    # i = 8
    i += 1
    @show i
    Random.seed!(i)

    scene = plot(xlims = (0, 1), ylims = (0, 1))
    sites = [(rand(Normal(0.5, 1)), rand(Exponential(1))) for i in 1:100]
    N, M = 0.5, 0.5
    corners = (0.0, 0.0), (N, 0.0), (N, M), (0.0, M)
    centre = (corners[1] .+ corners[3]) ./ 2
    rectangle_shape = Shape([0.0, N, N, 0.0], [0.0, 0.0, M, M])
    plot!(rectangle_shape, alpha = 0.5, label = "Rectangle")
    scatter!(first.(sites), last.(sites), label = "Sites", shape = :star, mlw = 0.5)
    anchor = reduce(.+, sites) ./ length(sites)
    scatter!([centre[1]], [centre[2]], label = "Centre")
    plot!(scene, size = (500, 400), aspect_ratio = 1.0, legend = :outerright)
end

function closest_anchor_to_rectangle(rect, anchor)
    (T, L), (N, M) = rect
    corners = (T, L), (T + N, L), (T + N, L + M), (T, L + M)
    # If inside, we are done
    if (T ≤ anchor[1] ≤ T + N) && (L ≤ anchor[2] ≤ L + M)
        min_point = anchor
    else
        min_dist = Inf
        min_point = (-Inf, -Inf) # break noisily if it doesn't work 
        j = 4
        for i in 1:4
            curr_point = closest_point(anchor, corners[i], corners[j])
            curr_dist = distance(curr_point, anchor)
            if curr_dist < min_dist
                min_point = curr_point
                min_dist = curr_dist
            end
            j = i
        end
    end
    return min_point
end

function circleShape(centre, r)
    x, y = centre
    θ = range(start = 0, stop = 2π, length = 500)
    x .+ r * sin.(θ), y .+ r * cos.(θ)
end

scene_ia = deepcopy(scene)
begin
    rect = ((0, 0), (N, M))
    inner_anchor = closest_anchor_to_rectangle(rect, anchor)
    cl_ia = sites[findmin(site -> distance(inner_anchor, site), sites)[2]]
    radii_cl_ia = map(corner -> distance(cl_ia, corner), corners)
    for i in 1:3
        plot!(scene_ia, circleShape(corners[i], radii_cl_ia[i]), seriestype = [:shape,], lw = 0.5, c = :green, fillalpha = 0.1, label = "")
    end
    plot!(scene_ia, circleShape(corners[4], radii_cl_ia[4]), seriestype = [:shape,], lw = 0.5, c = :green, fillalpha = 0.1, label = "Inclusion Area")
    scatter!(scene_ia, [anchor[1]], [anchor[2]], label = "Anchor")
    plot!(scene_ia, [inner_anchor[1], anchor[1]], [inner_anchor[2], anchor[2]], label = "", lw = 2.0, lc = :black, ls = :dash)

    plot!(scene_ia, [cl_ia[1], inner_anchor[1]], [cl_ia[2], inner_anchor[2]], label = "", lw = 2.0, lc = :black)
    scatter!(scene_ia, [inner_anchor[1]], [inner_anchor[2]], label = "Inner Anchor")
end

scene_fc = deepcopy(scene)
begin
    rect = ((0, 0), (N, M))
    inner_anchor_fc = closest_anchor_from_centre(rect, anchor)
    cl_fc = sites[findmin(site -> distance(inner_anchor_fc, site), sites)[2]]
    radii_cl_fc = map(corner -> distance(cl_fc, corner), corners)
    for i in 1:3
        plot!(scene_fc, circleShape(corners[i], radii_cl_fc[i]), seriestype = [:shape,], lw = 0.5, c = :yellow, fillalpha = 0.1, label = "")
    end
    plot!(scene_fc, circleShape(corners[4], radii_cl_fc[4]), seriestype = [:shape,], lw = 0.5, c = :yellow, fillalpha = 0.1, label = "Inclusion Area")
    scatter!(scene_fc, [anchor[1]], [anchor[2]], label = "Anchor")
    plot!(scene_fc, [centre[1], anchor[1]], [centre[2], anchor[2]], label = "", lw = 2.0, lc = :black, ls = :dash)

    plot!(scene_fc, [cl_fc[1], inner_anchor_fc[1]], [cl_fc[2], inner_anchor_fc[2]], label = "", lw = 2.0, lc = :black)
    scatter!(scene_fc, [inner_anchor_fc[1]], [inner_anchor_fc[2]], label = "Inner Anchor")
end

scene_ce = deepcopy(scene)
begin
    cl_ce = sites[findmin(site -> distance(centre, site), sites)[2]]
    radii_cl_ce = map(corner -> distance(cl_ce, corner), corners)
    for i in 1:3
        plot!(scene_ce, circleShape(corners[i], radii_cl_ce[i]), seriestype = [:shape,], lw = 0.5, c = :red, fillalpha = 0.1, label = "")
    end
    plot!(scene_ce, circleShape(corners[4], radii_cl_ce[4]), seriestype = [:shape,], lw = 0.5, c = :red, fillalpha = 0.1, label = "Inclusion Area")
    plot!(scene_ce, [cl_ce[1], centre[1]], [cl_ce[2], centre[2]], label = "", lw = 2.0, lc = :black)
end

scene_fc
plot(scene_ia, scene_ce, scene_fc, layout = (3, 1), size = (800, 1000), xlims = (-0.5, 1), ylims = (-0.5, 1))