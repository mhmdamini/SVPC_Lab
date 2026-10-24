# Checks on the toy example: the generated VPC must cut off the LP optimum and
# must not cut off any integer-feasible point.

using CuttingPlane
using JuMP
using Test

include(joinpath(@__DIR__, "..", "examples", "toy_model.jl"))

@testset "VPC on toy model" begin
    model = toy_model()
    optimize_model!(model)
    z_int = objective_value(model)

    relaxed = relax_model(model)
    optimize_model!(relaxed)
    z_lp = objective_value(relaxed)
    ω, _ = get_variable_values(relaxed)
    @test z_lp < z_int

    integer_vars = integer_variable_indices(model)
    k = first(fractional_integer_variables(relaxed, integer_vars))

    std = create_standard_form_model(relaxed)
    subproblems, _ = create_disjunctive_subproblems(std, k, integer_vars)
    results = [process_subproblem(sp) for sp in subproblems]
    _, α, β = solve_prlp(ω, [r.point for r in results], [r.rays for r in results]; normalization = :auto)

    # The cut separates the LP optimum ...
    @test α' * ω < β - 1e-6

    # ... and is satisfied by every feasible integer point (brute force).
    feasible(x) = -2x[1] - 2x[2] <= -0.5 && 2x[1] - x[2] <= 1.5 && x[1] <= 1 && x[2] <= 1
    for x1 in 0:1, x2 in 0:1
        feasible((x1, x2)) && @test α' * [x1, x2] >= β - 1e-6
    end

    # Adding the cut never weakens the bound and never passes the MIP optimum.
    cut_model = add_cut(relaxed, α, β; verbose = false)
    optimize_model!(cut_model)
    z_cut = objective_value(cut_model)
    @test z_lp - 1e-6 <= z_cut <= z_int + 1e-6
    @test 0 <= gap_closed(z_cut, z_int, z_lp) <= 100 + 1e-6
end

@testset "All cut families are valid on the toy model" begin
    model = toy_model()
    relaxed = relax_model(model)
    optimize_model!(relaxed)
    ω, _ = get_variable_values(relaxed)
    integer_vars = integer_variable_indices(model)
    feasible(x) = -2x[1] - 2x[2] <= -0.5 && 2x[1] - x[2] <= 1.5 && x[1] <= 1 && x[2] <= 1

    for method in keys(CUT_METHODS)
        cut = generate_cut(relaxed, 1, integer_vars, method)
        @test cut.α' * ω < cut.β - 1e-6                       # separates the LP optimum
        for x1 in 0:1, x2 in 0:1                              # keeps every integer point
            feasible((x1, x2)) && @test cut.α' * [x1, x2] >= cut.β - 1e-6
        end
    end
end

# 3-D knapsack-style example from the web interface:
#   max 5x + 4y + 3z  s.t.  2x + 3y + z <= 5.5,  4x + y + 2z <= 11.3,
#                           3x + 4y + 2z <= 10.7,  x, y, z ∈ ℤ₊
function knapsack3d_model()
    model = Model()
    @variable(model, x[1:3] >= 0, Int)
    @objective(model, Max, 5x[1] + 4x[2] + 3x[3])
    @constraint(model, 2x[1] + 3x[2] + x[3] <= 5.5)
    @constraint(model, 4x[1] + x[2] + 2x[3] <= 11.3)
    @constraint(model, 3x[1] + 4x[2] + 2x[3] <= 10.7)
    return model
end
knapsack3d_feasible(x) = 2x[1] + 3x[2] + x[3] <= 5.5 && 4x[1] + x[2] + 2x[3] <= 11.3 &&
                         3x[1] + 4x[2] + 2x[3] <= 10.7

@testset "Monoidal strengthening is valid ($label)" for (label, model, feasible, box) in (
        ("toy model", toy_model(),
         x -> -2x[1] - 2x[2] <= -0.5 && 2x[1] - x[2] <= 1.5 && x[1] <= 1 && x[2] <= 1, 0:1),
        ("3-D knapsack", knapsack3d_model(), knapsack3d_feasible, 0:6))
    relaxed = relax_model(model)
    optimize_model!(relaxed)
    ω, _ = get_variable_values(relaxed)
    integer_vars = integer_variable_indices(model)
    n = length(ω)
    points = [collect(p) for p in Iterators.product(fill(box, n)...) if feasible(p)]
    @test !isempty(points)

    for k in fractional_integer_variables(relaxed, integer_vars)
        vpc = generate_cut(relaxed, k, integer_vars, :vpc)
        mono = generate_cut(relaxed, k, integer_vars, :vpc_mono)
        @test mono.β ≈ vpc.β
        @test all(mono.α .<= vpc.α .+ 1e-9)                 # never weaker than the VPC
        @test mono.α' * ω < mono.β - 1e-6                   # still separates ω
        @test all(p -> mono.α' * p >= mono.β - 1e-6, points) # keeps every integer point
    end
end

@testset "Monoidal strengthening on the 3-D example, split on x" begin
    model = knapsack3d_model()
    relaxed = relax_model(model)
    optimize_model!(relaxed)
    integer_vars = integer_variable_indices(model)
    vpc = generate_cut(relaxed, 1, integer_vars, :vpc)
    mono = generate_cut(relaxed, 1, integer_vars, :vpc_mono)
    @test mono.α[2] < vpc.α[2] - 1e-3                       # the y coefficient is strengthened
    points = [collect(p) for p in Iterators.product(0:6, 0:6, 0:6) if knapsack3d_feasible(p)]
    @test length(points) == 16
    @test count(p -> mono.α' * p < mono.β - 1e-6, points) == 0
end

@testset "Intersection cut is valid and dominated by GMI ($label)" for (label, model, feasible, box) in (
        ("toy model", toy_model(),
         x -> -2x[1] - 2x[2] <= -0.5 && 2x[1] - x[2] <= 1.5 && x[1] <= 1 && x[2] <= 1, 0:1),
        ("3-D knapsack", knapsack3d_model(), knapsack3d_feasible, 0:6))
    relaxed = relax_model(model)
    optimize_model!(relaxed)
    z_lp = objective_value(relaxed)
    ω, _ = get_variable_values(relaxed)
    integer_vars = integer_variable_indices(model)
    points = [collect(p) for p in Iterators.product(fill(box, length(ω))...) if feasible(p)]
    improvement(cut) = (m = add_cut(relaxed, cut.α, cut.β; verbose = false);
                        optimize_model!(m); abs(objective_value(m) - z_lp))

    for k in fractional_integer_variables(relaxed, integer_vars)
        sic = generate_cut(relaxed, k, integer_vars, :sic)
        gmi = generate_cut(relaxed, k, integer_vars, :gmi)
        @test sic.α' * ω < sic.β - 1e-6                       # separates ω
        @test all(p -> sic.α' * p >= sic.β - 1e-6, points)    # keeps every integer point
        @test improvement(gmi) >= improvement(sic) - 1e-6     # GMI is at least as strong
    end
end

# ---- Simplex: the textbook two-phase simplex must agree with HiGHS ----------

function lp_from(build)
    m = Model()
    build(m)
    return m
end

const SIMPLEX_CASES = [
    "toy (min, Phase I)" => lp_from(m -> begin
        @variable(m, 0 <= x[1:2] <= 1)
        @objective(m, Min, x[1] + 2x[2])
        @constraint(m, -2x[1] - 2x[2] <= -0.5); @constraint(m, 2x[1] - x[2] <= 1.5)
    end),
    "3-D knapsack" => relax_model(knapsack3d_model()),
    "≥ and = rows (Phase I)" => lp_from(m -> begin
        @variable(m, x[1:3] >= 0)
        @objective(m, Min, 2x[1] + 3x[2] + x[3])
        @constraint(m, x[1] + x[2] + x[3] >= 4); @constraint(m, x[1] - x[2] == 1); @constraint(m, x[2] + 2x[3] >= 3)
    end),
    # Beale (1955): cycles under the largest-coefficient rule, not under Bland's rule.
    "Beale's cycling example" => lp_from(m -> begin
        @variable(m, x[1:4] >= 0)
        @objective(m, Max, 0.75x[1] - 150x[2] + 0.02x[3] - 6x[4])
        @constraint(m, 0.25x[1] - 60x[2] - 0.04x[3] + 9x[4] <= 0)
        @constraint(m, 0.5x[1] - 90x[2] - 0.02x[3] + 3x[4] <= 0)
        @constraint(m, x[3] <= 1)
    end),
]

@testset "Two-phase simplex agrees with HiGHS ($label)" for (label, lp) in SIMPLEX_CASES
    std = create_standard_form_model(lp)
    optimize_model!(std)
    run = tableau_simplex(std)
    @test run.status == :optimal
    @test run.objective ≈ objective_value(std) atol = 1e-7
    data = lp_matrix_data(std)
    @test all(Matrix(data.A) * run.x .<= data.b_upper .+ 1e-7) && all(run.x .>= -1e-9)   # feasible

    normal = create_normal_form_model(std)
    optimize_model!(normal)
    for c in verify_optimal_tableau(normal)
        @test c.ok
    end
end

@testset "Two-phase simplex detects infeasible and unbounded LPs" begin
    infeasible = lp_from(m -> (@variable(m, x >= 0); @objective(m, Max, x); @constraint(m, x <= 1); @constraint(m, x >= 2)))
    unbounded = lp_from(m -> (@variable(m, x[1:2] >= 0); @objective(m, Max, x[1]); @constraint(m, x[2] - x[1] <= 1)))
    @test tableau_simplex(create_standard_form_model(infeasible)).status == :infeasible
    @test tableau_simplex(create_standard_form_model(unbounded)).status == :unbounded
end

# ---- Gomory fractional cut ---------------------------------------------------

@testset "Gomory fractional cut is valid ($label)" for (label, model, feasible, box) in (
        ("toy model", toy_model(),
         x -> -2x[1] - 2x[2] <= -0.5 && 2x[1] - x[2] <= 1.5 && x[1] <= 1 && x[2] <= 1, 0:1),
        ("3-D knapsack", knapsack3d_model(), knapsack3d_feasible, 0:6))
    relaxed = relax_model(model)
    optimize_model!(relaxed)
    ω, _ = get_variable_values(relaxed)
    integer_vars = integer_variable_indices(model)
    points = [collect(p) for p in Iterators.product(fill(box, length(ω))...) if feasible(p)]
    for k in fractional_integer_variables(relaxed, integer_vars)
        gfc = generate_cut(relaxed, k, integer_vars, :gfc)
        @test gfc.α' * ω < gfc.β - 1e-6                       # separates ω
        @test all(p -> gfc.α' * p >= gfc.β - 1e-6, points)    # keeps every integer point
    end
end

@testset "Gomory fractional cut refuses mixed-integer problems" begin
    m = Model()
    @variable(m, x >= 0, Int); @variable(m, y >= 0)
    @objective(m, Max, x + y)
    @constraint(m, 2x + y <= 3.5); @constraint(m, y <= 0.7)
    relaxed = relax_model(m); optimize_model!(relaxed)
    @test_throws ErrorException generate_cut(relaxed, 1, [1], :gfc)
end

# ---- Iterative methods: cutting-plane loop, branch and bound, branch and cut ----

function production_model()
    m = Model(); @variable(m, x >= 0, Int); @variable(m, y >= 0, Int)
    @objective(m, Max, 5x + 4y)
    @constraint(m, 6x + 4y <= 24.5); @constraint(m, x + 2y <= 6.5); @constraint(m, -x + y <= 1)
    return m
end
function knapsack_model()                      # 0/1 variables
    w = [12, 7, 11, 8, 9, 6, 5, 14]; v = [24, 13, 23, 15, 16, 10, 9, 25]
    m = Model(); @variable(m, x[1:8], Bin); @objective(m, Max, v' * x); @constraint(m, w' * x <= 35)
    return m
end
function mixed_model()                         # y is continuous
    m = Model(); @variable(m, x >= 0, Int); @variable(m, y >= 0); @variable(m, z >= 0, Int)
    @objective(m, Min, -3x - 2y - 4z)
    @constraint(m, x + y + 2z <= 6.6); @constraint(m, 2x - y + z <= 4.3); @constraint(m, x + 3z <= 7.5)
    return m
end
const ITERATIVE_CASES = ["toy" => toy_model, "3-D knapsack" => knapsack3d_model, "production" => production_model,
                         "0/1 knapsack" => knapsack_model, "mixed integer" => mixed_model]
mip_optimum(make) = (m = make(); optimize_model!(m); objective_value(m))

@testset "Branch and bound / cut find the MIP optimum ($label)" for (label, make) in ITERATIVE_CASES
    z = mip_optimum(make)
    for kw in ((;), (node_selection = :depth_first,), (node_selection = :breadth_first, branching = :first_fractional),
               (cut_methods = [:gmi], root_cut_rounds = 3),
               (cut_methods = [:gmi, :vpc], root_cut_rounds = 2, node_cut_rounds = 1))
        r = branch_and_bound(make(); kw...)
        @test r.status == :optimal
        @test r.objective ≈ z atol = 1e-6
    end
end

@testset "Cutting-plane loop never passes the MIP optimum ($label)" for (label, make) in ITERATIVE_CASES
    z = mip_optimum(make)
    minimize = objective_sense(make()) == MIN_SENSE
    pure = length(integer_variable_indices(make())) == length(all_variables(make()))
    for methods in ([:sic], [:gmi], [:vpc], [:vpc_mono], [:gmi, :sic, :vpc], pure ? [:gfc] : Symbol[])
        isempty(methods) && continue
        r = cutting_plane_loop(make(); methods, max_rounds = 15, split_rule = :all)
        @test minimize ? r.objective <= z + 1e-6 : r.objective >= z - 1e-6     # every cut was valid
        r.status == :optimal && @test r.objective ≈ z atol = 1e-6
    end
end

@testset "Gomory's pure cutting-plane method solves pure-integer examples" begin
    for make in (toy_model, knapsack3d_model, production_model)
        r = cutting_plane_loop(make(); methods = [:gfc], max_rounds = 30, split_rule = :all)
        @test r.status == :optimal
        @test r.objective ≈ mip_optimum(make) atol = 1e-6
    end
end
