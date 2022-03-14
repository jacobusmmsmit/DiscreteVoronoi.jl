using LinearAlgebra: ⋅ # Not sure where to put this

function _norm(xy, p::Real=2)::Float64
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


function distance(x, y, p::Real=2)
    _norm((x[1] - y[1], x[2] - y[2]), p)
end

function nearest(x)
    (round(Int, x[1]), round(Int, x[2]))
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