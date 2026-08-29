"""
    sample_dtwa(state, nsites; rng=nothing, T=Float64)

Generate one discrete truncated-Wigner (DTWA) phase-space sample for a
physical product state.

The returned array has shape `(nsites, 3)`, with columns ordered as
`(x, y, z)`.

For a spin-1/2 Cartesian product state, the component parallel to the
physical polarization is fixed to its eigenvalue `±1`, while the two
transverse Pauli components are sampled independently from `±1` with
equal probability.

This function generates only the initial phase-space point. It does not
perform time evolution and does not create an ensemble of trajectories.
"""
function sample_dtwa(
    state::AbstractState,
    nsites::Integer;
    rng=nothing,
    T::Type{<:AbstractFloat}=Float64,
)
    nsites > 0 ||
        throw(ArgumentError("nsites must be positive"))

    sample = Matrix{T}(undef, Int(nsites), 3)

    @inbounds for site in 1:nsites
        direction = state_direction(state, site, nsites)
        _sample_dtwa_site!(
            sample,
            site,
            direction;
            rng=rng,
        )
    end

    return sample
end


"""
    sample_dtwa(state, chain; kwargs...)

Convenience overload using the number of sites stored in a `Chain`.
"""
function sample_dtwa(
    state::AbstractState,
    chain::Chain;
    kwargs...,
)
    return sample_dtwa(
        state,
        chain.nsites;
        kwargs...,
    )
end


@inline function _sample_dtwa_site!(
    sample::AbstractMatrix{T},
    site::Integer,
    direction::NTuple{3,<:Integer};
    rng=nothing,
) where {T<:AbstractFloat}

    @inbounds for component in 1:3
        value = direction[component]

        sample[site, component] =
            value == 0 ?
            T(_random_dtwa_sign(rng)) :
            T(value)
    end

    return nothing
end


@inline function _random_dtwa_sign(rng)
    draw = rng === nothing ?
           rand(Bool) :
           rand(rng, Bool)

    return draw ? 1 : -1
end
