"""
    AbstractSimulationResult

Abstract supertype for results returned by `simulate`.
"""
abstract type AbstractSimulationResult end


"""
    DTWAResult

Result of a closed-system DTWA simulation.

Fields:

- `t`: saved simulation times.
- `trajectories`: array with shape `(nsites, 3, ntimes, ntrajectories)`.

The Cartesian component dimension is ordered `(x, y, z)`.

A single trajectory at a saved time can be accessed without allocation as

    @view result.trajectories[:, :, time_index, trajectory_index]
"""
struct DTWAResult{
    T<:Real,
    A<:AbstractArray{<:Real,4},
} <: AbstractSimulationResult
    t::Vector{T}
    trajectories::A
end


"""
    times(result)

Return the saved simulation times.
"""
times(result::DTWAResult) = result.t


"""
    ntrajectories(result)

Return the number of DTWA trajectories stored in a result.
"""
ntrajectories(
    result::DTWAResult,
) = size(result.trajectories, 4)


"""
    trajectory(result, trajectory_index)

Return a non-allocating view of one complete trajectory.

The returned view has shape `(nsites, 3, ntimes)`.
"""
function trajectory(
    result::DTWAResult,
    trajectory_index::Integer,
)
    1 <= trajectory_index <= ntrajectories(result) ||
        throw(BoundsError(
            1:ntrajectories(result),
            trajectory_index,
        ))

    return @view result.trajectories[
        :,
        :,
        :,
        trajectory_index,
    ]
end
