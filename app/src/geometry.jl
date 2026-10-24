# =============================================================================
# Geometry for the 2-D / 3-D plots of the web interface.
#
# Everything is computed exactly in Julia (vertex enumeration of small
# polytopes) and sent to the browser as plain coordinates; Plotly only draws.
# =============================================================================

const GEOM_TOL = 1e-7

"""
    inequality_system(model) -> (G, h)

The LP relaxation of `model` as `G x <= h` (rows, then variable bounds).
"""
function inequality_system(model::Model)
    data = lp_matrix_data(model)
    A = Matrix(data.A)
    n = size(A, 2)
    rows, rhs = Vector{Vector{Float64}}(), Float64[]
    add!(a, b) = (push!(rows, a); push!(rhs, b))
    for i in 1:size(A, 1)
        isfinite(data.b_upper[i]) && add!(A[i, :], data.b_upper[i])
        isfinite(data.b_lower[i]) && add!(-A[i, :], -data.b_lower[i])
    end
    for j in 1:n
        e = zeros(n); e[j] = 1.0
        isfinite(data.x_upper[j]) && add!(e, data.x_upper[j])
        isfinite(data.x_lower[j]) && add!(-e, -data.x_lower[j])
    end
    return permutedims(reduce(hcat, rows)), rhs
end

"""
    box_system(lo, hi) -> (G, h)

The box `lo <= x <= hi` as `G x <= h`.
"""
function box_system(lo, hi)
    n = length(lo)
    return vcat(Matrix(1.0I, n, n), -Matrix(1.0I, n, n)), vcat(hi, -lo)
end

"""
    plot_box(relaxed_model, ω) -> (lo, hi)

Axis ranges for the plot: from the variable bounds (or 0) to the largest value
each variable reaches over the relaxation, padded by half a unit. Unbounded
directions are clipped a few units past the LP optimum `ω`.
"""
function plot_box(relaxed_model::Model, ω)
    data = lp_matrix_data(relaxed_model)
    n = length(ω)
    lo, hi = zeros(n), zeros(n)
    for j in 1:n
        lo[j] = isfinite(data.x_lower[j]) ? data.x_lower[j] : min(0.0, ω[j] - 3)
        m = copy(relaxed_model)
        @objective(m, Max, all_variables(m)[j])
        optimize_model!(m)
        hi[j] = termination_status(m) == MOI.OPTIMAL ? objective_value(m) : ω[j] + 3
        hi[j] = max(hi[j], lo[j] + 1)
    end
    return floor.(lo) .- 0.5, ceil.(hi) .+ 0.5
end

"""All `d`-element subsets of `1:m`."""
function combinations(m::Int, d::Int)
    d == 0 && return [Int[]]
    return [vcat(c, j) for c in combinations(m, d - 1) for j in (isempty(c) ? 1 : last(c) + 1):m]
end

"""
    vertices(G, h; eq = nothing) -> Vector{Vector{Float64}}

Vertices of `{x : G x <= h}` (intersected with the hyperplane `eq = (a, b)`,
i.e. `a'x == b`, when given), by enumerating every square subsystem.
"""
function vertices(G::AbstractMatrix, h::AbstractVector; eq = nothing)
    n = size(G, 2)
    k = eq === nothing ? n : n - 1
    found = Vector{Vector{Float64}}()
    for idx in combinations(size(G, 1), k)
        M, r = G[idx, :], h[idx]
        if eq !== nothing
            M, r = vcat(M, eq[1]'), vcat(r, eq[2])
        end
        abs(det(M)) < 1e-10 && continue
        v = M \ r
        all(G * v .<= h .+ GEOM_TOL) || continue
        any(w -> norm(w - v) < 1e-6, found) || push!(found, v)
    end
    return found
end

"""Order the vertices of a convex polygon (in 2-D, or on a plane in 3-D) counter-clockwise."""
function order_polygon(pts::Vector{Vector{Float64}}, normal = nothing)
    length(pts) < 3 && return pts
    c = sum(pts) / length(pts)
    if length(c) == 2
        u, w = [1.0, 0.0], [0.0, 1.0]
    else                                    # orthonormal basis of the plane
        nrm = normalize(normal)
        u = normalize(abs(nrm[1]) < 0.9 ? cross(nrm, [1.0, 0, 0]) : cross(nrm, [0, 1.0, 0]))
        w = cross(nrm, u)
    end
    return sort(pts; by = p -> atan(dot(p - c, w), dot(p - c, u)))
end

"""2-D convex hull (Andrew's monotone chain), counter-clockwise."""
function convex_hull_2d(pts::Vector{Vector{Float64}})
    P = sort(unique(pts))
    length(P) <= 2 && return P
    turn(o, a, b) = (a[1] - o[1]) * (b[2] - o[2]) - (a[2] - o[2]) * (b[1] - o[1])
    half(seq) = foldl(seq; init = Vector{Vector{Float64}}()) do H, p
        while length(H) >= 2 && turn(H[end-1], H[end], p) <= 0
            pop!(H)
        end
        push!(H, p)
    end
    lower, upper = half(P), half(reverse(P))
    return vcat(lower[1:end-1], upper[1:end-1])
end

"""
    integer_points(G, h, lo, hi; limit = 20_000) -> (feasible, infeasible)

Integer points of the plot box, split by feasibility for `G x <= h`.
Returns empty lists when the box holds more than `limit` points.
"""
function integer_points(G, h, lo, hi; limit::Int = 20_000)
    ranges = [ceil(Int, lo[j]):floor(Int, hi[j]) for j in eachindex(lo)]
    prod(length, ranges) > limit && return Vector{Vector{Float64}}(), Vector{Vector{Float64}}()
    feas, infeas = Vector{Vector{Float64}}(), Vector{Vector{Float64}}()
    for p in Iterators.product(ranges...)
        x = collect(Float64, p)
        push!(all(G * x .<= h .+ GEOM_TOL) ? feas : infeas, x)
    end
    return feas, infeas
end

"""
    region_geometry(relaxed_model, ω, cuts; all_integer = true) -> Dict

Everything the browser needs to draw the feasible region of a 2- or
3-variable problem: plot box, LP polytope, integer points and integer hull
(only when every variable is integer), and for each cut `(α, β)` its
line/plane and the polytope that remains after it.
"""
function region_geometry(relaxed_model::Model, ω, cuts; all_integer::Bool = true)
    G, h = inequality_system(relaxed_model)
    n = size(G, 2)
    lo, hi = plot_box(relaxed_model, ω)
    Gb, hb = box_system(lo, hi)
    Gc, hc = vcat(G, Gb), vcat(h, hb)                 # relaxation clipped to the box

    lp_region = vertices(Gc, hc)
    n == 2 && (lp_region = order_polygon(lp_region))

    feas, infeas = all_integer ? integer_points(G, h, lo, hi) :
                                 (Vector{Vector{Float64}}(), Vector{Vector{Float64}}())
    hull = n == 2 ? convex_hull_2d(feas) : feas

    cut_geoms = map(cuts) do cut
        α, β = cut.α, cut.β
        plane = vertices(Gb, hb; eq = (α, β))             # cut line / plane inside the box
        plane = order_polygon(plane, n == 3 ? α : nothing)
        after = vertices(vcat(Gc, -α'), vcat(hc, -β))     # region left after the cut
        n == 2 && (after = order_polygon(after))
        Dict("plane" => plane, "region" => after)
    end

    return Dict("dim" => n, "lo" => lo, "hi" => hi, "lp_region" => lp_region,
                "integer_feasible" => feas, "integer_infeasible" => infeas,
                "integer_hull" => hull, "cuts" => cut_geoms)
end

"""
    iteration_geometry(relaxed_model, cuts, n_rounds; all_integer = true, max_cuts = 150) -> Dict or nothing

Geometry for replaying an iterative run in a 2- or 3-variable plot: the base
picture of `region_geometry` (box, LP polytope, integer points, hull) plus
each cut's line/plane (with its round and family) and the region left after
each round `0:n_rounds` (cumulative). `cuts` are `(α, β, round, method)`.
Returns `nothing` when there are more than `max_cuts` cuts (too slow to draw).
"""
function iteration_geometry(relaxed_model::Model, cuts, n_rounds::Int; all_integer::Bool = true, max_cuts::Int = 150)
    length(cuts) > max_cuts && return nothing
    G, h = inequality_system(relaxed_model)
    n = size(G, 2)
    ω = value.(all_variables(relaxed_model))
    lo, hi = plot_box(relaxed_model, ω)
    Gb, hb = box_system(lo, hi)
    Gc, hc = vcat(G, Gb), vcat(h, hb)
    order(pts) = n == 2 ? order_polygon(pts) : pts

    feas, infeas = all_integer ? integer_points(G, h, lo, hi) : (Vector{Vector{Float64}}(), Vector{Vector{Float64}}())
    planes = map(cuts) do c
        Dict("plane" => order_polygon(vertices(Gb, hb; eq = (c.α, c.β)), n == 3 ? c.α : nothing),
             "round" => c.round, "method" => String(c.method))
    end
    regions = map(0:n_rounds) do r
        active = [c for c in cuts if c.round <= r]
        Ga = isempty(active) ? Gc : vcat(Gc, reduce(vcat, [-c.α' for c in active]))
        ha = isempty(active) ? hc : vcat(hc, [-c.β for c in active])
        order(vertices(Ga, ha))
    end
    return Dict("dim" => n, "lo" => lo, "hi" => hi, "lp_region" => regions[1],
                "integer_feasible" => feas, "integer_infeasible" => infeas,
                "integer_hull" => n == 2 ? convex_hull_2d(feas) : feas,
                "cuts" => planes, "regions" => regions)
end
