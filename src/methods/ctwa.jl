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

    CTWA(; clustering, trajectories=1000,
           sampling=DiscreteSampling())

Cluster truncated-Wigner approximation configuration.

`cluster_size` selects the number of physical spin-1/2 sites per cluster using
the default clustering associated with the model geometry.

Alternatively, `clustering` may be supplied explicitly to specify the physical
site membership of every cluster. This is useful for spatially compact
clusters in multidimensional geometries.

A cluster of size `n` has `4^n - 1` traceless Pauli-string phase-space
coordinates.

`sampling` selects the initial cluster phase-space sampling prescription.
The default is [`DiscreteSampling`](@ref).
"""
struct CTWA{
    S<:AbstractCTWASampling,
    C,
} <: AbstractApproximation
    cluster_size::Int
    trajectories::Int
    sampling::S
    clustering::C

    function CTWA(;
        cluster_size::Union{Nothing,Integer}=nothing,
        clustering::Union{Nothing,Clustering}=nothing,
        trajectories::Integer=1000,
        sampling::AbstractCTWASampling=DiscreteSampling(),
    )
        trajectories > 0 ||
            throw(ArgumentError(
                "trajectories must be positive",
            ))

        if clustering === nothing
            k =
                cluster_size === nothing ?
                2 :
                Int(cluster_size)

            k > 0 ||
                throw(ArgumentError(
                    "cluster_size must be positive",
                ))

            return new{
                typeof(sampling),
                Nothing,
            }(
                k,
                Int(trajectories),
                sampling,
                nothing,
            )
        end

        if cluster_size !== nothing
            Int(cluster_size) == clustering.cluster_size ||
                throw(ArgumentError(
                    "cluster_size=$(Int(cluster_size)) does not match " *
                    "explicit clustering size $(clustering.cluster_size)",
                ))
        end

        return new{
            typeof(sampling),
            typeof(clustering),
        }(
            clustering.cluster_size,
            Int(trajectories),
            sampling,
            clustering,
        )
    end
end