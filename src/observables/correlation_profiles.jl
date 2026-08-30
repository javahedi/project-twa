"""
    CorrelationProfile

Time-dependent, distance-resolved correlation data.

Fields:

- `distances`: geometric separations included in the profile.
- `values`: matrix with shape `(ntimes, ndistances)`.
- `axis`: Cartesian spin axis.
- `connected`: whether connected correlations were computed.

The row index of `values` follows `times(result)`, while the column index
corresponds to `distances`.
"""
struct CorrelationProfile{T<:Real,A<:AbstractMatrix{T}}
    distances::Vector{Int}
    values::A
    axis::Symbol
    connected::Bool
end


"""
    correlation_profile(result, geometry, axis; connected=false)

Compute the spatially averaged same-axis correlation profile

    Cᵅ(r, t)
    = (1 / N_r) Σ_{d(i,j)=r} ⟨σᵢᵅ σⱼᵅ⟩

for all geometric distances supported by `geometry`.

If `connected=true`, compute instead

    Cᵅ_conn(r, t)
    = (1 / N_r) Σ_{d(i,j)=r}
      [⟨σᵢᵅ σⱼᵅ⟩
       - ⟨σᵢᵅ⟩⟨σⱼᵅ⟩].

The profile includes `r = 0`. For equal Cartesian components the on-site
Pauli identity gives

    ⟨(σᵢᵅ)²⟩ = 1.

For an open chain the returned distances are `0:(N-1)`. For a periodic
chain they are `0:fld(N, 2)`.

The geometry is passed explicitly because `DTWAResult` currently stores
trajectory data and times, but not the originating model metadata.
"""
function correlation_profile(
    result::DTWAResult,
    geometry::Chain,
    axis::Symbol;
    connected::Bool=false,
)
    _check_observable_axis(axis)

    nsites_result = size(result.trajectories, 1)

    geometry.nsites == nsites_result ||
        throw(DimensionMismatch(
            "geometry has $(geometry.nsites) sites, but result has " *
            "$nsites_result sites",
        ))

    distances = collect(
        0:_maximum_profile_distance(geometry),
    )

    Tout = promote_type(
        eltype(result.trajectories),
        Float64,
    )

    values = Matrix{Tout}(
        undef,
        length(result.t),
        length(distances),
    )

    for (column, r) in pairs(distances)
        _distance_correlation!(
            @view(values[:, column]),
            result,
            geometry,
            axis,
            r,
            connected,
        )
    end

    return CorrelationProfile(
        distances,
        values,
        axis,
        connected,
    )
end


@inline function _maximum_profile_distance(
    chain::Chain{<:OpenBoundary},
)
    return chain.nsites - 1
end


@inline function _maximum_profile_distance(
    chain::Chain{<:PeriodicBoundary},
)
    return fld(chain.nsites, 2)
end


function _distance_correlation!(
    destination::AbstractVector,
    result::DTWAResult,
    geometry::Chain,
    axis::Symbol,
    r::Integer,
    connected::Bool,
)
    axis_index = _axis_index(axis)
    ntimes = length(result.t)
    ntraj = ntrajectories(result)
    nsites = geometry.nsites

    Tout = eltype(destination)

    if r == 0
        @inbounds for time_index in 1:ntimes
            if connected
                total = zero(Tout)

                for site in 1:nsites
                    mean = zero(Tout)

                    for trajectory_index in 1:ntraj
                        mean += result.trajectories[
                            site,
                            axis_index,
                            time_index,
                            trajectory_index,
                        ]
                    end

                    mean /= ntraj
                    total += one(Tout) - mean * mean
                end

                destination[time_index] =
                    total / nsites
            else
                destination[time_index] =
                    one(Tout)
            end
        end

        return nothing
    end

    site_pairs = _site_pairs_at_distance(
        geometry,
        r,
    )

    npairs = length(site_pairs)

    npairs > 0 ||
        throw(ArgumentError(
            "no site pairs exist at distance $r",
        ))

    @inbounds for time_index in 1:ntimes
        total = zero(Tout)

        for (i, j) in site_pairs
            pair_mean = zero(Tout)

            for trajectory_index in 1:ntraj
                pair_mean += (
                    result.trajectories[
                        i,
                        axis_index,
                        time_index,
                        trajectory_index,
                    ] *
                    result.trajectories[
                        j,
                        axis_index,
                        time_index,
                        trajectory_index,
                    ]
                )
            end

            pair_mean /= ntraj

            if connected
                mean_i = zero(Tout)
                mean_j = zero(Tout)

                for trajectory_index in 1:ntraj
                    mean_i += result.trajectories[
                        i,
                        axis_index,
                        time_index,
                        trajectory_index,
                    ]

                    mean_j += result.trajectories[
                        j,
                        axis_index,
                        time_index,
                        trajectory_index,
                    ]
                end

                mean_i /= ntraj
                mean_j /= ntraj

                pair_mean -= mean_i * mean_j
            end

            total += pair_mean
        end

        destination[time_index] =
            total / npairs
    end

    return nothing
end


"""
    _site_pairs_at_distance(chain, r)

Return each unordered physical site pair once at geometric distance `r`.
"""
function _site_pairs_at_distance(
    chain::Chain,
    r::Integer,
)
    r >= 0 ||
        throw(ArgumentError(
            "distance must be nonnegative; got $r",
        ))

    pairs = Tuple{Int,Int}[]

    if r == 0
        return pairs
    end

    nsites = chain.nsites

    for i in 1:(nsites - 1)
        for j in (i + 1):nsites
            distance(chain, i, j) == r ||
                continue

            push!(pairs, (i, j))
        end
    end

    return pairs
end






"""
Correlation profiles and static structure factors for `CTWAResult`.

These methods reuse the primitive CTWA observable layer rather than
reimplementing cluster/operator mapping.

The public conventions match the DTWA observable API:

    correlation_profile(result, geometry, axis; connected=false)

and

    structure_factor(
        result,
        geometry,
        axis;
        connected=false,
        momenta=nothing,
    )

For correlation profiles, every unordered physical pair is counted once.
For periodic chains the distance range is 0:floor(N/2); for open chains it
is 0:N-1.

The structure factor is computed directly from physical site separations,

    S^a(q) = (1/N) sum_{i,j} cos(q (j-i)) <sigma_i^a sigma_j^a>,

rather than Fourier transforming a distance-averaged profile.
"""


function correlation_profile(
    result::CTWAResult,
    geometry::Chain,
    axis::Symbol;
    connected::Bool=false,
)
    _ctwa_validate_geometry(
        result,
        geometry,
    )

    distances =
        _ctwa_profile_distances(
            geometry,
        )

    ntimes =
        length(
            result.t,
        )

    T =
        eltype(
            result.trajectories,
        )

    values =
        zeros(
            T,
            ntimes,
            length(distances),
        )

    counts =
        zeros(
            Int,
            length(distances),
        )

    distance_to_column =
        Dict(
            distance => column
            for (column, distance) in enumerate(distances)
        )

    nsites =
        geometry.nsites

    # r = 0 contains the exact same-site identity contribution.
    zero_column =
        distance_to_column[0]

    @inbounds for site in 1:nsites
        values[:, zero_column] .+=
            expectation(
                result,
                connected ?
                    ConnectedCorrelation(
                        axis,
                        site,
                        axis,
                        site,
                    ) :
                    Correlation(
                        axis,
                        site,
                        axis,
                        site,
                    ),
            )

        counts[zero_column] +=
            1
    end

    for site_i in 1:(nsites - 1)
        for site_j in (site_i + 1):nsites
            distance =
                TWA.distance(
                    geometry,
                    site_i,
                    site_j,
                )

            column =
                distance_to_column[
                    distance
                ]

            observable =
                connected ?
                ConnectedCorrelation(
                    axis,
                    site_i,
                    axis,
                    site_j,
                ) :
                Correlation(
                    axis,
                    site_i,
                    axis,
                    site_j,
                )

            values[:, column] .+=
                expectation(
                    result,
                    observable,
                )

            counts[column] +=
                1
        end
    end

    for column in eachindex(distances)
        count =
            counts[column]

        count > 0 ||
            continue

        values[:, column] ./=
            T(count)
    end

    return CorrelationProfile(
        collect(distances),
        values,
        axis,
        connected,
    )
end






function _ctwa_profile_distances(
    geometry::Chain{<:OpenBoundary},
)
    return 0:(
        geometry.nsites -
        1
    )
end


function _ctwa_profile_distances(
    geometry::Chain{<:PeriodicBoundary},
)
    return 0:fld(
        geometry.nsites,
        2,
    )
end


@inline function _ctwa_validate_geometry(
    result::CTWAResult,
    geometry::Chain,
)
    geometry.nsites ==
    result.clustering.nsites ||
        throw(ArgumentError(
            "geometry has $(geometry.nsites) sites, but CTWAResult has " *
            "$(result.clustering.nsites)",
        ))

    return nothing
end


