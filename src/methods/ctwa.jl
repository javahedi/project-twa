abstract type AbstractCTWASampling end

struct GaussianSampling <: AbstractCTWASampling end
struct DiscreteSampling <: AbstractCTWASampling end


"""
    CTWA(; cluster_size=2, trajectories=1000,
           sampling=DiscreteSampling())

Cluster truncated-Wigner approximation configuration.

`cluster_size` controls the number of physical spin-1/2 sites grouped into
each CTWA cluster. A cluster of size `n` has `4^n - 1` traceless Pauli-string
phase-space coordinates.

`trajectories` is the number of independently sampled classical trajectories.

`sampling` selects the initial cluster phase-space sampling prescription.
"""
struct CTWA{S<:AbstractCTWASampling} <: AbstractApproximation
    cluster_size::Int
    trajectories::Int
    sampling::S

    function CTWA(;
        cluster_size::Integer=2,
        trajectories::Integer=1000,
        sampling::AbstractCTWASampling=DiscreteSampling(),
    )
        cluster_size > 0 ||
            throw(ArgumentError(
                "cluster_size must be positive",
            ))

        trajectories > 0 ||
            throw(ArgumentError(
                "trajectories must be positive",
            ))

        new{typeof(sampling)}(
            Int(cluster_size),
            Int(trajectories),
            sampling,
        )
    end
end