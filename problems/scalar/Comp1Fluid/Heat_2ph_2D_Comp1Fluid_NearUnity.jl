include(joinpath(@__DIR__, "Heat_2ph_2D_Comp1Fluid.jl"))

"""
Near-unity Henry jump Comp1Fluid benchmark (He = 2, Dg/Dl = 10).
"""
function run_near_unity(; csv_path=nothing, nx_list=nothing, norm::Real=2, relative::Bool=false)
    return run_case(:near_unity; csv_path=csv_path, nx_list=nx_list, norm=norm, relative=relative)
end

    results = run_near_unity()

    @testset "Comp1Fluid near_unity" begin
        orders = results.data.orders
        @test !isnan(orders.all)
        @test length(results.data.h_vals) == length(results.data.err_vals)
        @test results.data.h_vals[1] > results.data.h_vals[end]
        @test isfile(results.csv_path)
    end

