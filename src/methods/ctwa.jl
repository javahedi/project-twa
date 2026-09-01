"""
    AbstractCTWASampling

Abstract supertype for CTWA initial phase-space sampling strategies.

Sampling objects control how the initial cluster phase-space coordinates
are drawn, independently of the cluster size and number of trajectories.

Concrete sampling strategies currently provided by `TWA.jl` are
[`GaussianSampling`](@ref) and [`DiscreteSampling`](@ref).
"""
abstract type AbstractCTWASampling end


"""
    GaussianSampling()

Gaussian initial sampling for the cluster truncated-Wigner approximation.

This strategy samples the initial cluster phase-space distribution from
the Gaussian approximation used in conventional cluster TWA (gcTWA).

Use it with [`CTWA`](@ref), for example:

```julia
CTWA(
    cluster_size=2,
    trajectories=1000,
    sampling=GaussianSampling(),
)
```
"""
struct GaussianSampling <: AbstractCTWASampling end


"""
    DiscreteSampling()

Discrete initial sampling for the cluster truncated-Wigner approximation.

The discrete cluster sampler first draws the physical single-spin phase-space
variables using the DTWA prescription and then constructs every cluster
Pauli-string coordinate as the corresponding product of the sampled local
spin components.

This is the default sampling strategy used by [`CTWA`](@ref).

For example:

```julia
CTWA(
    cluster_size=2,
    trajectories=1000,
    sampling=DiscreteSampling(),
)
```
"""
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
The default is [`DiscreteSampling`](@ref). Use [`GaussianSampling`](@ref)
for Gaussian cluster sampling.
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
