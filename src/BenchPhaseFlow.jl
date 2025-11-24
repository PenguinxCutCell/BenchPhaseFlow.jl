module BenchPhaseFlow

using LinearAlgebra
using SparseArrays

# Try to load commonly used packages; warn if missing so the environment can be prepared
for pkg in (:LsqFit, :DataFrames, :IterativeSolvers, :CSV, :Test, :Roots, :SpecialFunctions, :Penguin)
    try
        @eval using $(pkg)
    catch err
        @warn "Package $(pkg) not available in the active environment. Add it with Pkg.add or Pkg.develop for local packages." error=err
    end
end

# Include repository utilities if present
_utils_path = joinpath(@__DIR__, "..", "utils", "convergence.jl")
if isfile(_utils_path)
    include(_utils_path)
else
    @warn "`utils/convergence.jl` not found. Some helper functions will be unavailable until that file exists." path=_utils_path
end

export count_inside_cells, compute_orders, compute_pairwise_orders, make_convergence_dataframe, make_diphasic_convergence_dataframe

end # module
