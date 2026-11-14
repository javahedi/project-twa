"""
    LocalDissipation(; decay=0.0, dephasing=0.0)

Uniform local spontaneous emission and pure dephasing for spin-½ systems.
The corresponding Lindblad operators at each site are
`sqrt(decay) * σ⁻` and `sqrt(dephasing / 2) * σᶻ`. Each active
channel at each site has an independent Wiener increment.
"""
struct LocalDissipation{T<:AbstractFloat}
    decay::T
    dephasing::T

    function LocalDissipation(; decay::Real=0.0, dephasing::Real=0.0)
        all(rate -> isfinite(rate) && rate >= 0, (decay, dephasing)) ||
            throw(ArgumentError("dissipation rates must be finite and nonnegative"))
        γ, γφ = promote(float(decay), float(dephasing))
        return new{typeof(γ)}(γ, γφ)
    end
end

Base.iszero(rates::LocalDissipation) =
    iszero(rates.decay) && iszero(rates.dephasing)

# Include the step size in the drift and sample one increment per channel.
# All strings involving this site share those same increments.
@inline function _local_step(rng, h, rates)
    γ, γφ = rates.decay, rates.dephasing
    decay_noise = iszero(γ) ? zero(h) : sqrt(γ * h) * randn(rng, typeof(h))
    phase_noise = iszero(γφ) ? zero(h) : sqrt(2 * γφ * h) * randn(rng, typeof(h))
    return (γ / 2 + γφ) * h, decay_noise + phase_noise, decay_noise - γ * h
end

function _dissipative_step!(next, u, compiled::CompiledDTWA, rates, h, rng)
    @inbounds for site in axes(u, 1)
        transverse, rotation, longitudinal = _local_step(rng, h, rates)
        x, y, z = u[site, 1], u[site, 2], u[site, 3]
        next[site, 1] -= transverse * x + rotation * y
        next[site, 2] += rotation * x - transverse * y
        next[site, 3] += longitudinal * (z + one(z))
    end
    return nothing
end

function _dissipative_step!(next, u, compiled::CompiledCTWA, rates, h, rng)
    partition = clustering(compiled)
    dimension = basis_size(basis(compiled))
    @inbounds for site in 1:partition.nsites
        cluster, position = site_cluster_position(partition, site)
        transverse, rotation, longitudinal = _local_step(rng, h, rates)
        weight = 4^(partition.cluster_size - position)
        # Iterate strings with I at this site, then substitute X, Y and Z.
        for block in 0:(4 * weight):dimension
            for offset in 0:(weight - 1)
                identity = block + offset
                ix, iy, iz = identity + weight, identity + 2 * weight, identity + 3 * weight
                x, y, z = u[ix, cluster], u[iy, cluster], u[iz, cluster]
                reduced = iszero(identity) ? one(z) : u[identity, cluster]
                next[ix, cluster] -= transverse * x + rotation * y
                next[iy, cluster] += rotation * x - transverse * y
                next[iz, cluster] += longitudinal * (z + reduced)
            end
        end
    end
    return nothing
end
