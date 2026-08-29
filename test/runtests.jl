using Test
using TWA

@testset "TWA.jl" begin
    include("test_geometry.jl")
    include("test_interactions.jl")
    include("test_spin_model.jl")
    include("test_pair_couplings.jl")
    include("test_fields.jl")
    include("test_hamiltonians.jl")

    include("test_states.jl")
    include("test_dtwa_sampling.jl")
    include("test_methods.jl")
    include("test_dtwa_ensemble.jl")
    include("test_dtwa_dynamics.jl")
    include("test_dtwa_compiled.jl")
    include("test_simulation.jl")
    include("test_observables.jl")
    include("test_collective_observables.jl")
    include("test_correlation_profiles.jl")
    include("test_structure_factors.jl")

    include("test_pauli_basis.jl")
    include("test_pauli_algebra.jl")
    include("test_clustering.jl")
    include("test_ctwa_hamiltonian.jl")
    include("test_ctwa_dynamics.jl")
    include("test_ctwa_compiled_algebra.jl")
    include("test_ctwa_sampling.jl")
    include("test_ctwa_compiled.jl")
    include("test_ctwa_simulation.jl")
    include("test_ctwa_observables.jl")
    include("test_ctwa_collective_observables.jl")
    include("test_ctwa_spatial_observables.jl")
end