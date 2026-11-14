using OrdinaryDiffEq


"""
    simulate(model, state, method::DTWA; kwargs...)

Run a discrete truncated-Wigner ensemble simulation.

This first simulation backend performs four steps:

1. compile the high-level spin model into DTWA dynamics data;
2. sample the requested initial DTWA ensemble;
3. evolve each trajectory with `OrdinaryDiffEq`;
4. collect all trajectories on a common output-time grid.

# Required keyword

- `tspan`: two-element tuple `(t_initial, t_final)`.

# Optional keywords

- `saveat`: output times or output spacing accepted by SciML. If omitted,
  101 equally spaced output times spanning `tspan` are used.
- `rng`: optional random-number generator for initial-state sampling.
- `T`: floating-point type used for sampled phase-space coordinates.
- `solver`: OrdinaryDiffEq algorithm. Defaults to `Tsit5()`.
- additional keyword arguments are forwarded to `solve`.

# Example

    result = simulate(
        model,
        DomainWall(),
        DTWA(trajectories=1000);
        tspan=(0.0, 12.0),
    )

Use `dissipation=LocalDissipation(...)` for open-system dynamics. Nonzero
rates select Itô Euler–Maruyama with maximum step `dt` (default `1e-3`),
shortened to land on each saved time. The supplied RNG controls both
sampling and noise. Solver overrides and additional solver keywords apply
to coherent dynamics only.
"""
function simulate(
    model::SpinModel,
    state::AbstractState,
    method::DTWA;
    tspan,
    saveat=nothing,
    rng=nothing,
    T::Type{<:AbstractFloat}=Float64,
    solver=nothing,
    dissipation::LocalDissipation=LocalDissipation(),
    dt=nothing,
    kwargs...,
)
    _validate_open_options(dissipation, solver, kwargs)
    open_system = !iszero(dissipation)
    algorithm = isnothing(solver) ? Tsit5() : solver
    solver_options = isnothing(dt) ? kwargs : (; kwargs..., dt)
    open_dt = isnothing(dt) ? 1e-3 : dt

    t0, tf = _validate_tspan(tspan)

    output_times = _simulation_saveat(
        t0,
        tf,
        saveat,
        T,
    )

    compiled = compile(
        model,
        method,
    )

    noise_rng = open_system ? Random.Xoshiro(rand(
        isnothing(rng) ? Random.default_rng() : rng, UInt64,
    )) : nothing

    initial_ensemble = sample_initial(
        model,
        state,
        method;
        rng=rng,
        T=T,
    )

    N = nsites(model.geometry)
    ntimes = length(output_times)
    ntraj = method.trajectories

    values = Array{T,4}(
        undef,
        N,
        3,
        ntimes,
        ntraj,
    )

    for k in 1:ntraj
        u0 = copy(
            @view initial_ensemble[:, :, k]
        )

        if open_system
            _open_trajectory!(
                view(values, :, :, :, k), u0, dtwa_rhs!,
                compiled, dissipation, T(t0), output_times, open_dt, noise_rng,
            )
            continue
        end

        problem = ODEProblem(
            dtwa_rhs!,
            u0,
            (T(t0), T(tf)),
            compiled,
        )

        solution = solve(
            problem,
            algorithm;
            saveat=output_times,
            solver_options...,
        )

        length(solution.u) == ntimes ||
            throw(ErrorException(
                "solver returned $(length(solution.u)) saved states; " *
                "expected $ntimes",
            ))

        @inbounds for time_index in 1:ntimes
            values[
                :,
                :,
                time_index,
                k,
            ] .= solution.u[time_index]
        end
    end

    return DTWAResult(
        collect(T, output_times),
        values,
    )
end


"""
    _validate_tspan(tspan)

Validate and return the two endpoints of a simulation time interval.
"""
function _validate_tspan(tspan)
    length(tspan) == 2 ||
        throw(ArgumentError(
            "tspan must contain exactly two values",
        ))

    t0, tf = tspan

    tf > t0 ||
        throw(ArgumentError(
            "tspan must satisfy t_final > t_initial; got $tspan",
        ))

    return t0, tf
end


"""
    _simulation_saveat(t0, tf, saveat, T)

Construct the common output grid used by all ensemble trajectories.
"""
function _simulation_saveat(
    t0,
    tf,
    ::Nothing,
    T::Type{<:AbstractFloat},
)
    return collect(
        range(
            T(t0),
            T(tf);
            length=101,
        ),
    )
end


function _simulation_saveat(
    t0,
    tf,
    saveat,
    T::Type{<:AbstractFloat},
)
    if saveat isa Real
        saveat > 0 ||
            throw(ArgumentError(
                "saveat spacing must be positive",
            ))

        times = collect(
            T(t0):T(saveat):T(tf),
        )

        if isempty(times) || times[end] != T(tf)
            push!(times, T(tf))
        end

        return times
    end

    times = T.(collect(saveat))

    isempty(times) &&
        throw(ArgumentError(
            "saveat must not be empty",
        ))

    first(times) >= T(t0) ||
        throw(ArgumentError(
            "saveat contains times before tspan",
        ))

    last(times) <= T(tf) ||
        throw(ArgumentError(
            "saveat contains times after tspan",
        ))

    issorted(times) ||
        throw(ArgumentError(
            "saveat times must be sorted",
        ))

    return times
end
