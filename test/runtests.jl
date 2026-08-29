using Test
using TWA

@testset "TWA.jl" begin
    include("test_geometry.jl")
    include("test_interactions.jl")
    include("test_spin_model.jl")
    include("test_pair_couplings.jl")
    include("test_fields.jl")
    include("test_hamiltonians.jl")
end