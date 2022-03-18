using LinearAlgebra # Not sure where to put this

function _norm(xy, p::Real = 2)::Float64
    if p == 1
        @fastmath abs(xy[1]) + abs(xy[2])
    elseif p == 2
        @fastmath sqrt(xy[1]^2 + xy[2]^2)
    elseif p == Inf
        @fastmath max(abs(xy[1]), abs(xy[2]))
    else
        @assert false
    end
end

function _dot(x, y)
    sum(x .* y)
end

function distance(x, y, p::Real = 2)
    _norm((x[1] - y[1], x[2] - y[2]), p)
end

function segment(s1, s2, t)
    s1 .+ t .* (s2 .- s1)
end
function closest(s1, s2, x)
    s = s2 .- s1
    segment(s1, s2, clamp(_dot(x .- s1, s) / _dot(s, s), 0.0, 1.0))
end
function closest(rect, mean)::Tuple{Float64, Float64}
    (t, l), (N, M) = rect
    if N == 1 && M == 1 
        (t, l)
    elseif t <= mean[1] <= t + N - 1 && l <= mean[2] <= l + M - 1
        mean
    else
        corners = (t, l), (t + N - 1, l), (t + N - 1, l + M - 1), (t, l + M - 1)
        segments = (1, 2), (2, 3), (3, 4), (4, 1)
        # _, index = findmin(segment -> distance(closest(corners[segment[1]], corners[segment[2]], mean), mean, 2), segments)
        # closest(corners[segments[index][1]], corners[segments[index][2]], mean)
        anchor = (NaN, NaN)
        anchor_dist = Inf
        for segment in segments
            point = closest(corners[segment[1]], corners[segment[2]], mean)
            point_dist = distance(point, mean, 2)
            if point_dist < anchor_dist
                anchor = point
                anchor_dist = point_dist
            end
        end
        anchor
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

function closest_anchor_to_rectangle(rect, anchor)::Tuple{Float64, Float64}
    (T, L), (N, M) = rect
    corners = (T, L), (T + N - 1, L), (T + N - 1, L + M - 1), (T, L + M - 1)
    # If rect is empty, we are done
    if N == 1 && M == 1 
        min_point = T, L
    # If inside, we are done
    elseif (T ≤ anchor[1] ≤ T + N - 1) && (L ≤ anchor[2] ≤ L + M - 1)
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

function round_tuple(x)
    (round(Int, x[1]), round(Int, x[2]))
end
