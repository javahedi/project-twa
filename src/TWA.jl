module TWA

import Random

export Chain
export OpenBoundary, PeriodicBoundary
export bonds, distance

export XXZ, NearestNeighbor, PowerLaw

export Field, field_components

export Hamiltonian

export SpinModel
export geometry, hamiltonian, terms

export coupling_strength, coupling_components


export AbstractState
export Polarized, Up, Down, DomainWall
export state_direction

export sample_dtwa

export AbstractApproximation
export DTWA

export sample_initial
export dtwa_rhs!

export CompiledDTWA
export compile

export AbstractSimulationResult, DTWAResult
export simulate
export times, ntrajectories, trajectory


export AbstractObservable
export LocalMagnetization, Magnetization
export Correlation, ConnectedCorrelation
export expectation
export SpinSecondMoment, SpinVariance
export CorrelationProfile
export correlation_profile
export StructureFactor
export structure_factor



export PauliStringBasis
export basis_size
export pauli_code, pauli_index, pauli_digits
export pauli_symbol, pauli_symbols


export pauli_product_digit
export pauli_product
export pauli_commutator
export pauli_anticommutes


export Clustering
export cluster_count, cluster_sites
export site_cluster, site_position, site_cluster_position
export same_cluster
export local_pauli_digits, local_pauli_index

export CTWALocalTerm
export CTWAInterClusterTerm
export CompiledCTWAHamiltonian
export compile_ctwa_hamiltonian
export local_terms, intercluster_terms

export ctwa_rhs!

export CTWACommutatorEntry
export CTWACommutatorAction
export CompiledCTWAAlgebra
export compile_ctwa_algebra
export commutator_action


export CTWASamplingPlan
export compile_ctwa_sampling
export sample_ctwa
export sample_ctwa!

export CTWA
export CompiledCTWA
export clustering
export basis

export CTWAResult

include("models/geometry.jl")
include("models/hamiltonians.jl")
include("models/interactions.jl")
include("models/fields.jl")
include("models/spin_model.jl")
include("models/pair_couplings.jl")


include("states/states.jl")

include("methods/methods.jl")
include("methods/dtwa.jl")

include("sampling/dtwa.jl")
include("sampling/dtwa_ensemble.jl")

include("dynamics/dtwa.jl")
include("dynamics/dtwa_compiled.jl")

include("simulation/results.jl")
include("simulation/dtwa.jl")



include("algebra/pauli_basis.jl")
include("algebra/pauli_algebra.jl")

include("clustering/ctwa.jl")

include("compilation/ctwa_hamiltonian.jl")
include("compilation/ctwa_algebra.jl")

include("sampling/ctwa.jl")
include("methods/ctwa.jl")
include("dynamics/ctwa.jl")

include("compilation/ctwa.jl")
include("simulation/ctwa.jl")

include("observables/observables.jl")
include("observables/collective.jl")
include("observables/correlation_profiles.jl")
include("observables/structure_factors.jl")




end
