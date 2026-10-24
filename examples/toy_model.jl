# Small two-variable integer program used as the default example and in the tests.
#
#   min   x₁ + 2x₂
#   s.t. -2x₁ - 2x₂ <= -0.5
#         2x₁ -  x₂ <=  1.5
#          x₁       <=  1
#               x₂  <=  1
#          x₁, x₂ ∈ ℤ₊
#
# LP optimum: (0.25, 0) with value 0.25.  Integer optimum: (0, 1) with value 2.

using JuMP

function toy_model()
    model = Model()
    @variable(model, x[1:2] >= 0, Int)
    @objective(model, Min, x[1] + 2x[2])
    @constraint(model, C1, -2x[1] - 2x[2] <= -0.5)
    @constraint(model, C2, 2x[1] - x[2] <= 1.5)
    @constraint(model, C3, x[1] <= 1)
    @constraint(model, C4, x[2] <= 1)
    return model
end
