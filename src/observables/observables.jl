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
"""
Primitive observables for `CTWAResult`.

The CTWA result stores trajectories as

    x[generator, cluster, time, trajectory].

Physical spin operators must therefore be mapped through the retained
`Clustering` and `PauliStringBasis`.

For two operators on different sites:

* same cluster:
      σ_i^a σ_j^b
  is represented by one cluster Pauli-string generator;

* different clusters:
      σ_i^a σ_j^b
  is represented by the product of two one-site cluster coordinates.

This distinction is essential to CTWA.
"""


"""
    expectation(result::CTWAResult, observable::LocalMagnetization)

Return the trajectory-averaged local Pauli expectation as a function of time.
"""
function expectation(
    result::CTWAResult,
    observable::LocalMagnetization,
)
    axis =
        _ctwa_axis_digit(
            observable.axis,
        )

    site =
        observable.site

    generator, cluster =
        _ctwa_local_generator(
            result,
            axis,
            site,
        )

    ntimes =
        length(
            result.t,
        )

    ntrajectories =
        size(
            result.trajectories,
            4,
        )

    values =
        zeros(
            eltype(result.trajectories),
            ntimes,
        )

    @inbounds for trajectory in 1:ntrajectories
        for time_index in 1:ntimes
            values[time_index] +=
                result.trajectories[
                    generator,
                    cluster,
                    time_index,
                    trajectory,
                ]
        end
    end

    values ./=
        ntrajectories

    return values
end


"""
    expectation(result::CTWAResult, observable::Magnetization)

Return the magnetization density

    <S_a>/N = (1/N) sum_i <σ_i^a>

as a function of time.
"""
function expectation(
    result::CTWAResult,
    observable::Magnetization,
)
    axis =
        _ctwa_axis_digit(
            observable.axis,
        )

    nsites =
        result.clustering.nsites

    ntimes =
        length(
            result.t,
        )

    ntrajectories =
        size(
            result.trajectories,
            4,
        )

    values =
        zeros(
            eltype(result.trajectories),
            ntimes,
        )

    @inbounds for site in 1:nsites
        generator, cluster =
            _ctwa_local_generator(
                result,
                axis,
                site,
            )

        for trajectory in 1:ntrajectories
            for time_index in 1:ntimes
                values[time_index] +=
                    result.trajectories[
                        generator,
                        cluster,
                        time_index,
                        trajectory,
                    ]
            end
        end
    end

    values ./=
        nsites *
        ntrajectories

    return values
end


"""
    expectation(result::CTWAResult, observable::Correlation)

Return the equal-time two-site Pauli correlation

    <σ_i^a σ_j^b>

as a function of time.

For the same physical site:

* same axis uses the exact identity σ_a² = I;
* different axes are rejected because the ordered Pauli product is generally
  complex, while the current observable API is real-valued.
"""
function expectation(
    result::CTWAResult,
    observable::Correlation,
)
    axis_i =
        _ctwa_axis_digit(
            observable.axis_i,
        )

    axis_j =
        _ctwa_axis_digit(
            observable.axis_j,
        )

    site_i =
        observable.site_i

    site_j =
        observable.site_j

    _ctwa_check_site(
        result,
        site_i,
    )

    _ctwa_check_site(
        result,
        site_j,
    )

    ntimes =
        length(
            result.t,
        )

    T =
        eltype(
            result.trajectories,
        )

    if site_i == site_j
        axis_i == axis_j ||
            throw(ArgumentError(
                "same-site correlations with different Pauli axes are " *
                "complex ordered products and are not supported",
            ))

        return ones(
            T,
            ntimes,
        )
    end

    cluster_i =
        site_cluster(
            result.clustering,
            site_i,
        )

    cluster_j =
        site_cluster(
            result.clustering,
            site_j,
        )

    if cluster_i == cluster_j
        return _ctwa_same_cluster_correlation(
            result,
            axis_i,
            site_i,
            axis_j,
            site_j,
            cluster_i,
        )
    end

    return _ctwa_intercluster_correlation(
        result,
        axis_i,
        site_i,
        axis_j,
        site_j,
        cluster_i,
        cluster_j,
    )
end


"""
    expectation(result::CTWAResult, observable::ConnectedCorrelation)

Return

    <σ_i^a σ_j^b>
    - <σ_i^a><σ_j^b>

as a function of time.
"""
function expectation(
    result::CTWAResult,
    observable::ConnectedCorrelation,
)
    correlation =
        observable.correlation

    left =
        LocalMagnetization(
            correlation.axis_i,
            correlation.site_i,
        )

    right =
        LocalMagnetization(
            correlation.axis_j,
            correlation.site_j,
        )

    return expectation(
        result,
        correlation,
    ) .-
           expectation(
        result,
        left,
    ) .*
           expectation(
        result,
        right,
    )
end


@inline function _ctwa_axis_digit(
    axis::Symbol,
)
    axis === :x &&
        return 1

    axis === :y &&
        return 2

    axis === :z &&
        return 3

    throw(ArgumentError(
        "Pauli axis must be :x, :y, or :z; got $axis",
    ))
end


@inline function _ctwa_check_site(
    result::CTWAResult,
    site::Integer,
)
    1 <= site <= result.clustering.nsites ||
        throw(BoundsError(
            1:result.clustering.nsites,
            site,
        ))

    return nothing
end


@inline function _ctwa_local_generator(
    result::CTWAResult,
    axis::Int,
    site::Integer,
)
    _ctwa_check_site(
        result,
        site,
    )

    cluster, position =
        site_cluster_position(
            result.clustering,
            site,
        )

    generator =
        _ctwa_single_site_generator(
            result.basis,
            position,
            axis,
        )

    return generator, cluster
end


@inline function _ctwa_single_site_generator(
    basis::PauliStringBasis,
    position::Int,
    axis::Int,
)
    n =
        basis.cluster_size

    1 <= position <= n ||
        throw(BoundsError(
            1:n,
            position,
        ))

    code =
        axis *
        4^(
            n -
            position
        )

    return code
end


@inline function _ctwa_two_site_generator(
    basis::PauliStringBasis,
    position_i::Int,
    axis_i::Int,
    position_j::Int,
    axis_j::Int,
)
    position_i != position_j ||
        throw(ArgumentError(
            "two-site cluster generator requires distinct positions",
        ))

    n =
        basis.cluster_size

    code =
        axis_i *
        4^(
            n -
            position_i
        ) +
        axis_j *
        4^(
            n -
            position_j
        )

    return code
end


function _ctwa_same_cluster_correlation(
    result::CTWAResult,
    axis_i::Int,
    site_i::Int,
    axis_j::Int,
    site_j::Int,
    cluster::Int,
)
    position_i =
        site_position(
            result.clustering,
            site_i,
        )

    position_j =
        site_position(
            result.clustering,
            site_j,
        )

    generator =
        _ctwa_two_site_generator(
            result.basis,
            position_i,
            axis_i,
            position_j,
            axis_j,
        )

    ntimes =
        length(
            result.t,
        )

    ntrajectories =
        size(
            result.trajectories,
            4,
        )

    values =
        zeros(
            eltype(result.trajectories),
            ntimes,
        )

    @inbounds for trajectory in 1:ntrajectories
        for time_index in 1:ntimes
            values[time_index] +=
                result.trajectories[
                    generator,
                    cluster,
                    time_index,
                    trajectory,
                ]
        end
    end

    values ./=
        ntrajectories

    return values
end


function _ctwa_intercluster_correlation(
    result::CTWAResult,
    axis_i::Int,
    site_i::Int,
    axis_j::Int,
    site_j::Int,
    cluster_i::Int,
    cluster_j::Int,
)
    position_i =
        site_position(
            result.clustering,
            site_i,
        )

    position_j =
        site_position(
            result.clustering,
            site_j,
        )

    generator_i =
        _ctwa_single_site_generator(
            result.basis,
            position_i,
            axis_i,
        )

    generator_j =
        _ctwa_single_site_generator(
            result.basis,
            position_j,
            axis_j,
        )

    ntimes =
        length(
            result.t,
        )

    ntrajectories =
        size(
            result.trajectories,
            4,
        )

    values =
        zeros(
            eltype(result.trajectories),
            ntimes,
        )

    @inbounds for trajectory in 1:ntrajectories
        for time_index in 1:ntimes
            values[time_index] +=
                result.trajectories[
                    generator_i,
                    cluster_i,
                    time_index,
                    trajectory,
                ] *
                result.trajectories[
                    generator_j,
                    cluster_j,
                    time_index,
                    trajectory,
                ]
        end
    end

    values ./=
        ntrajectories

    return values
end
