"""
    CTWAResult

CTWA simulation result.

`trajectories` has layout

    (generator, cluster, time, trajectory)

where `generator` is interpreted using `basis` and the physical site-to-cluster
mapping is interpreted using `clustering`.

Unlike a DTWA result, CTWA generator indices are not physically meaningful
without this metadata, so both are retained explicitly.
"""
struct CTWAResult{
    T<:Real,
    A<:AbstractArray{<:Real,4},
    B<:PauliStringBasis,
    C<:Clustering,
} <: AbstractSimulationResult
    t::Vector{T}
    trajectories::A
    basis::B
    clustering::C
end


"""
    simulate(model, state, method::CTWA; kwargs...)

Run a CTWA ensemble simulation.

Keyword arguments:

* `tspan`: required two-element time interval;
* `saveat`: saved times. If omitted, 101 uniformly spaced times are used;
* `rng`: optional RNG used for initial-state sampling;
* `T`: floating-point state type, default `Float64`;
* `solver`: OrdinaryDiffEq solver, default `Tsit5()`;
* additional keyword arguments are forwarded to `solve`.

The compiled CTWA Hamiltonian, sparse commutator algebra, and sampling plan are
constructed once and reused by all trajectories.

Use `dissipation=LocalDissipation(...)` for open-system dynamics. Nonzero
rates use RK4 for Hamiltonian evolution and an Itô Euler–Maruyama
dissipative update with maximum step `dt` (default `1e-3`),
shortened to land on each saved time. The supplied RNG controls both
sampling and noise. Solver overrides and additional solver keywords apply
to coherent dynamics only.
"""
function simulate(
    model::SpinModel,
    state::AbstractState,
    method::CTWA;
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

    t0, tf =
        _ctwa_validate_tspan(
            tspan,
        )

    saved_times =
    _simulation_saveat(
        t0,
        tf,
        saveat,
        T,
    )

    compiled =
        compile(
            model,
            state,
            method;
            T=T,
        )

    dimension =
        basis_size(
            basis(compiled),
        )

    nclusters =
        cluster_count(
            clustering(compiled),
        )

    nsave =
        length(
            saved_times,
        )

    ntrajectories =
        method.trajectories

    trajectories =
        Array{T,4}(
            undef,
            dimension,
            nclusters,
            nsave,
            ntrajectories,
        )

    random =
        isnothing(rng) ?
        Random.default_rng() :
        rng

    noise_rng = open_system ? Random.Xoshiro(rand(random, UInt64)) : nothing

    x0 =
        Matrix{T}(
            undef,
            dimension,
            nclusters,
        )

    for trajectory in 1:ntrajectories
        sample_ctwa!(
            x0,
            compiled;
            rng=random,
        )

        if open_system
            _open_trajectory!(
                view(trajectories, :, :, :, trajectory), copy(x0), ctwa_rhs!,
                compiled, dissipation, T(t0), saved_times, open_dt, noise_rng,
            )
            continue
        end

        problem =
            ODEProblem(
                ctwa_rhs!,
                copy(x0),
                (
                    T(t0),
                    T(tf),
                ),
                compiled,
            )

        solution =
            solve(
                problem,
                algorithm;
                saveat=saved_times,
                solver_options...,
            )

        length(solution.u) == nsave ||
            throw(ArgumentError(
                "solver returned $(length(solution.u)) saved states, " *
                "but $nsave were requested",
            ))

        @inbounds for time_index in 1:nsave
            trajectories[
                :,
                :,
                time_index,
                trajectory,
            ] .=
                solution.u[
                    time_index
                ]
        end
    end

    return CTWAResult(
        collect(
            T,
            saved_times,
        ),
        trajectories,
        basis(compiled),
        clustering(compiled),
    )
end


function _ctwa_validate_tspan(
    tspan,
)
    length(tspan) == 2 ||
        throw(ArgumentError(
            "tspan must contain exactly two times",
        ))

    t0, tf =
        tspan

    tf > t0 ||
        throw(ArgumentError(
            "tspan must satisfy tf > t0",
        ))

    return t0, tf
end
