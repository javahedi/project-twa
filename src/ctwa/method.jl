"""
    CTWA(; cluster_size=2, trajectories=1000)

Cluster truncated-Wigner approximation configuration.

`cluster_size` controls the number of physical spin-1/2 sites grouped into
each CTWA cluster. A cluster of size `n` has `4^n - 1` traceless Pauli-string
phase-space coordinates.

`trajectories` is the number of independently sampled classical trajectories.

Solver choices, RNGs, precision, and environments are intentionally not stored
in this type.
"""
struct CTWA <: AbstractApproximation
    cluster_size::Int
    trajectories::Int

    function CTWA(;
        cluster_size::Integer=2,
        trajectories::Integer=1000,
    )
        cluster_size > 0 ||
            throw(ArgumentError(
                "cluster_size must be positive",
            ))

        trajectories > 0 ||
            throw(ArgumentError(
                "trajectories must be positive",
            ))

        new(
            Int(cluster_size),
            Int(trajectories),
        )
    end
end
