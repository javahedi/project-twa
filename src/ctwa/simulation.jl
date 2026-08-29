"""
    CTWAResult

Closed-system CTWA simulation result.

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

Run a closed-system CTWA ensemble simulation.

Keyword arguments:

* `tspan`: required two-element time interval;
* `saveat`: saved times. If omitted, 101 uniformly spaced times are used;
* `rng`: optional RNG used for initial-state sampling;
* `T`: floating-point state type, default `Float64`;
* `solver`: OrdinaryDiffEq solver, default `Tsit5()`;
* additional keyword arguments are forwarded to `solve`.

The compiled CTWA Hamiltonian, sparse commutator algebra, and sampling plan are
constructed once and reused by all trajectories.
"""
function simulate(
    model::SpinModel,
    state::AbstractState,
    method::CTWA;
    tspan,
    saveat=nothing,
    rng=nothing,
    T::Type{<:AbstractFloat}=Float64,
    solver=Tsit5(),
    kwargs...,
)
    t0, tf =
        _ctwa_validate_tspan(
            tspan,
        )

    saved_times =
        isnothing(saveat) ?
        collect(
            range(
                T(t0),
                T(tf);
                length=101,
            ),
        ) :
        collect(
            T,
            saveat,
        )

    isempty(saved_times) &&
        throw(ArgumentError(
            "saveat must contain at least one time",
        ))

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
                solver;
                saveat=saved_times,
                kwargs...,
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
        T.(
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
