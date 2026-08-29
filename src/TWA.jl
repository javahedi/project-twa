module TWA

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



end
