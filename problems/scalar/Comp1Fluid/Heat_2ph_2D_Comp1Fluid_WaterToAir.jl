include(joinpath(@__DIR__, "Heat_2ph_2D_Comp1Fluid.jl"))

"""
Water-to-air Comp1Fluid benchmark (He = 1/30, Dg/Dl = 0.1).
"""
function run_water_to_air(; csv_path=nothing, nx_list=nothing, norm::Real=2, relative::Bool=false)
    return run_case(:water_to_air; csv_path=csv_path, nx_list=nx_list, norm=norm, relative=relative)
end

    results = run_water_to_air()

    @testset "Comp1Fluid water_to_air" begin
        orders = results.data.orders
        @test !isnan(orders.all)
        @test length(results.data.h_vals) == length(results.data.err_vals)
        @test results.data.h_vals[1] > results.data.h_vals[end]
        @test isfile(results.csv_path)
    end
