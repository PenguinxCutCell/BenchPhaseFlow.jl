using Test
import Pkg

# Activate project environment for tests
Pkg.activate(joinpath(@__DIR__, ".."))

using BenchPhaseFlow

@testset "Basic sanity tests" begin
    # test count_inside_cells with a simple named tuple having cell_types
    capacity = (cell_types = [1, 0, 1, 2])
    @test count_inside_cells(capacity) == 2

    # test compute_pairwise_orders for a simple halving-error sequence
    h = [0.5, 0.25, 0.125]
    errs = [1.0, 0.5, 0.25]
    orders = compute_pairwise_orders(h, errs)
    @test isapprox(orders[2], 1.0; atol=1e-12)
    @test isapprox(orders[3], 1.0; atol=1e-12)
end
