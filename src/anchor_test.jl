## In this file I am testing whether the choice of the centre of the rectangle
## can be improved upon by choosing a different "anchor" from which to calculate
## the closest site. The anchor will hopefully choose a more important closest site
## than the centre in some circumstances such as when sites are heavily skewed in
## one direction.

## Note: the anchor *must* (I think, proof by contradiction exists in my head)
## be inside/on the edge of the rectangle otherwise the anchor may choose a site
## that results in too large circles.

using BenchmarkTools
using Statistics
using Plots

function distance(p, q)
    @fastmath sqrt((p[1] - q[1])^2 + (p[2] - q[2])^2)
end


"""
Calculate the intersection of two lines defined by 4 points,
2 per line
TODO: remake for line segments using wikipedia article
"""
function intersection_two_lines((x1, y1), (x2, y2), (x3, y3), (x4, y4))
    # Numerators
    pxn = (x1 * y2 - x2 * y1) * (x3 - x4) - (x3 * y4 - x4 * y3) * (x1 - x2)
    pyn = (x1 * y2 - x2 * y1) * (y3 - y4) - (x3 * y4 - x4 * y3) * (y1 - y2)
    # Denominator
    pd = (x1 - x2) * (y3 - y4) - (x3 - x4) * (y1 - y2)
    return pxn / pd, pyn / pd
end

function closest_anchor_in_rectangle(rectangle, anchor)
    # calculate centre of rectangle
    # calculate intersection of rectangle and line segment from centre to anchor
end


A = [(i, i + 2) for i in 1:100]
reduce(.+, A)
@time reduce(.+, A) ./ length(A)
@btime reduce(.+, $A) ./ length($A)

scene = plot(xlims = (0, 1), ylims = (0, 1))
sites = [rand(2) for i in 1:10]
N, M = 0.5, 0.5
corners = (0.0, 0.0), (N, 0.0), (N, M), (0.0, M)
anchor = reduce(.+, sites) ./ length(sites)
centre = N / 2, M / 2
sites[findmin(site -> distance(anchor, site), sites)[2]]
sites[findmin(site -> distance(centre, site), sites)[2]]

dists_a = map(site -> distance(centre, site), sites)
areas_a = map(r -> π * r^2, dists)

