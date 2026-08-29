"""
    AbstractObservable

Abstract supertype for physical observables evaluated from simulation results.
"""
abstract type AbstractObservable end


"""
    LocalMagnetization(axis, site)

Single-site Pauli expectation value

    ⟨σᵢᵅ⟩

for Cartesian axis `:x`, `:y`, or `:z`.

Evaluating this observable on a `DTWAResult` returns one value per saved time.
"""
struct LocalMagnetization <: AbstractObservable
    axis::Symbol
    site::Int

    function LocalMagnetization(axis::Symbol, site::Integer)
        _check_observable_axis(axis)

        site > 0 ||
            throw(ArgumentError(
                "site must be positive; got $site",
            ))

        new(axis, Int(site))
    end
end


"""
    Magnetization(axis)

Collective magnetization density

    ⟨Sᵅ⟩ / N
    = (1/N) Σᵢ ⟨σᵢᵅ⟩.

Evaluating this observable on a `DTWAResult` returns one value per saved time.
"""
struct Magnetization <: AbstractObservable
    axis::Symbol

    function Magnetization(axis::Symbol)
        _check_observable_axis(axis)
        new(axis)
    end
end


"""
    Correlation(axis_i, site_i, axis_j, site_j)

Two-point Pauli correlation

    ⟨σᵢᵅ σⱼᵝ⟩.

For distinct sites, DTWA evaluates this as the ensemble average of the
trajectory-wise product of the two phase-space coordinates.

For the same site, Pauli operator identities are used instead of naively
multiplying phase-space coordinates:

    σᵅ σᵅ = I.

The current observable layer therefore supports same-site correlations only
when `axis_i == axis_j`.
"""
struct Correlation <: AbstractObservable
    axis_i::Symbol
    site_i::Int
    axis_j::Symbol
    site_j::Int

    function Correlation(
        axis_i::Symbol,
        site_i::Integer,
        axis_j::Symbol,
        site_j::Integer,
    )
        _check_observable_axis(axis_i)
        _check_observable_axis(axis_j)

        site_i > 0 ||
            throw(ArgumentError(
                "site_i must be positive; got $site_i",
            ))

        site_j > 0 ||
            throw(ArgumentError(
                "site_j must be positive; got $site_j",
            ))

        if site_i == site_j && axis_i != axis_j
            throw(ArgumentError(
                "same-site correlations with different axes are not " *
                "supported yet because the ordered Pauli product can be complex",
            ))
        end

        new(axis_i, Int(site_i), axis_j, Int(site_j))
    end
end


"""
    ConnectedCorrelation(axis_i, site_i, axis_j, site_j)

Connected two-point correlation

    Cᵅᵝᵢⱼ
    = ⟨σᵢᵅ σⱼᵝ⟩
      - ⟨σᵢᵅ⟩⟨σⱼᵝ⟩.
"""
struct ConnectedCorrelation <: AbstractObservable
    correlation::Correlation
end


function ConnectedCorrelation(
    axis_i::Symbol,
    site_i::Integer,
    axis_j::Symbol,
    site_j::Integer,
)
    return ConnectedCorrelation(
        Correlation(
            axis_i,
            site_i,
            axis_j,
            site_j,
        ),
    )
end


"""
    expectation(result, observable)

Evaluate a physical observable on a simulation result.

The returned vector is aligned with `times(result)`.
"""
function expectation(
    result::DTWAResult,
    observable::LocalMagnetization,
)
    axis = _axis_index(observable.axis)
    site = observable.site

    _check_observable_site(result, site)

    ntimes = length(result.t)
    ntraj = ntrajectories(result)

    Tout = promote_type(eltype(result.trajectories), Float64)
    values = Vector{Tout}(undef, ntimes)

    @inbounds for time_index in 1:ntimes
        total = zero(Tout)

        for trajectory_index in 1:ntraj
            total += result.trajectories[
                site,
                axis,
                time_index,
                trajectory_index,
            ]
        end

        values[time_index] = total / ntraj
    end

    return values
end


function expectation(
    result::DTWAResult,
    observable::Magnetization,
)
    axis = _axis_index(observable.axis)

    nsites = size(result.trajectories, 1)
    ntimes = length(result.t)
    ntraj = ntrajectories(result)

    Tout = promote_type(eltype(result.trajectories), Float64)
    values = Vector{Tout}(undef, ntimes)

    normalization = nsites * ntraj

    @inbounds for time_index in 1:ntimes
        total = zero(Tout)

        for trajectory_index in 1:ntraj
            for site in 1:nsites
                total += result.trajectories[
                    site,
                    axis,
                    time_index,
                    trajectory_index,
                ]
            end
        end

        values[time_index] = total / normalization
    end

    return values
end


function expectation(
    result::DTWAResult,
    observable::Correlation,
)
    _check_observable_site(result, observable.site_i)
    _check_observable_site(result, observable.site_j)

    ntimes = length(result.t)
    Tout = promote_type(eltype(result.trajectories), Float64)

    if observable.site_i == observable.site_j
        return ones(Tout, ntimes)
    end

    axis_i = _axis_index(observable.axis_i)
    axis_j = _axis_index(observable.axis_j)
    ntraj = ntrajectories(result)

    values = Vector{Tout}(undef, ntimes)

    @inbounds for time_index in 1:ntimes
        total = zero(Tout)

        for trajectory_index in 1:ntraj
            total += (
                result.trajectories[
                    observable.site_i,
                    axis_i,
                    time_index,
                    trajectory_index,
                ] *
                result.trajectories[
                    observable.site_j,
                    axis_j,
                    time_index,
                    trajectory_index,
                ]
            )
        end

        values[time_index] = total / ntraj
    end

    return values
end


function expectation(
    result::DTWAResult,
    observable::ConnectedCorrelation,
)
    correlation = observable.correlation

    pair = expectation(result, correlation)

    left = expectation(
        result,
        LocalMagnetization(
            correlation.axis_i,
            correlation.site_i,
        ),
    )

    right = expectation(
        result,
        LocalMagnetization(
            correlation.axis_j,
            correlation.site_j,
        ),
    )

    return pair .- left .* right
end


@inline function _check_observable_axis(axis::Symbol)
    axis in (:x, :y, :z) ||
        throw(ArgumentError(
            "axis must be :x, :y, or :z; got $axis",
        ))

    return nothing
end


@inline function _axis_index(axis::Symbol)
    axis === :x && return 1
    axis === :y && return 2
    axis === :z && return 3

    throw(ArgumentError(
        "axis must be :x, :y, or :z; got $axis",
    ))
end


@inline function _check_observable_site(
    result::DTWAResult,
    site::Integer,
)
    nsites = size(result.trajectories, 1)

    1 <= site <= nsites ||
        throw(BoundsError(
            1:nsites,
            site,
        ))

    return nothing
end
