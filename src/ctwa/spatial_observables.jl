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


function structure_factor(
    result::CTWAResult,
    geometry::Chain,
    axis::Symbol;
    connected::Bool=false,
    momenta=nothing,
)
    _ctwa_validate_geometry(
        result,
        geometry,
    )

    nsites =
        geometry.nsites

    qs =
        momenta === nothing ?
        _ctwa_default_momenta(
            nsites,
            eltype(result.trajectories),
        ) :
        collect(momenta)

    ntimes =
        length(
            result.t,
        )

    T =
        promote_type(
            eltype(result.trajectories),
            eltype(qs),
        )

    values =
        zeros(
            T,
            ntimes,
            length(qs),
        )

    # Diagonal i=j terms.
    for site in 1:nsites
        diagonal =
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

        for q_index in eachindex(qs)
            values[:, q_index] .+=
                diagonal
        end
    end

    # Off-diagonal terms.  Because same-axis correlations on different
    # physical sites commute, the (i,j) and (j,i) terms combine into
    # 2*cos(q*(j-i)).
    for site_i in 1:(nsites - 1)
        for site_j in (site_i + 1):nsites
            correlation =
                expectation(
                    result,
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
                        ),
                )

            separation =
                site_j -
                site_i

            for (q_index, q) in pairs(qs)
                weight =
                    T(2) *
                    cos(
                        T(q) *
                        T(separation),
                    )

                values[:, q_index] .+=
                    weight .*
                    correlation
            end
        end
    end

    values ./=
        T(nsites)

    return StructureFactor(
        qs,
        values,
        axis,
        connected,
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


function _ctwa_default_momenta(
    nsites::Int,
    ::Type{T},
) where {T<:Real}
    two_pi =
        T(2) *
        T(pi)

    return T[
        two_pi *
        T(m) /
        T(nsites)
        for m in 0:(nsites - 1)
    ]
end
