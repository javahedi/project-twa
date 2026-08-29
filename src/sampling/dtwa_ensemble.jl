"""
    sample_initial(model, state, method::DTWA; rng=nothing, T=Float64)

Generate the full ensemble of DTWA initial phase-space trajectories.

The returned array has shape

    (nsites, 3, trajectories)

where the second dimension stores Cartesian components `(x, y, z)` and
the third dimension indexes trajectories.

Each trajectory is therefore available as

    @view samples[:, :, trajectory]

without allocating a separate matrix.

The physical initial state is specified independently of DTWA. This
method only applies the DTWA phase-space sampling rule to that state.

# Arguments

- `model`: spin model whose geometry determines the number of sites.
- `state`: physical initial state.
- `method`: `DTWA` approximation containing the trajectory count.
- `rng`: optional random-number generator.
- `T`: floating-point element type of the returned array.

# Example

    model = SpinModel(Chain(50), XXZ())
    state = DomainWall()
    method = DTWA(trajectories=1000)

    samples = sample_initial(
        model,
        state,
        method;
        rng=Xoshiro(1234),
    )
"""
function sample_initial(
    model::SpinModel,
    state::AbstractState,
    method::DTWA;
    rng=nothing,
    T::Type{<:AbstractFloat}=Float64,
)
    nsites = model.geometry.nsites
    ntrajectories = method.trajectories

    samples = Array{T,3}(
        undef,
        nsites,
        3,
        ntrajectories,
    )

    @inbounds for trajectory in 1:ntrajectories
        for site in 1:nsites
            direction = state_direction(
                state,
                site,
                nsites,
            )

            for component in 1:3
                value = direction[component]

                samples[
                    site,
                    component,
                    trajectory,
                ] = value == 0 ?
                    T(_random_dtwa_sign(rng)) :
                    T(value)
            end
        end
    end

    return samples
end
