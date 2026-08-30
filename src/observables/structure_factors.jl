"""
    StructureFactor

Time-dependent static structure-factor data.

Fields:

- `momenta`: momentum values `q`.
- `values`: matrix with shape `(ntimes, nq)`.
- `axis`: Cartesian spin axis.
- `connected`: whether connected correlations were used.

The row index of `values` follows `times(result)`, while the column index
corresponds to `momenta`.
"""
struct StructureFactor{T<:Real,A<:AbstractMatrix{T}}
    momenta::Vector{T}
    values::A
    axis::Symbol
    connected::Bool
end


"""
    structure_factor(result, geometry, axis; connected=false, momenta=nothing)

Compute the same-axis static structure factor

    Sᵅ(q, t)
    = (1/N) Σᵢⱼ exp[i q (j-i)]
      ⟨σᵢᵅ σⱼᵅ⟩.

For same-axis Hermitian correlations the pair contribution is real, so the
implementation uses the equivalent cosine form

    Sᵅ(q, t)
    = (1/N) [
        Σᵢ ⟨(σᵢᵅ)²⟩
        + 2 Σᵢ<ⱼ cos(q(j-i))
          ⟨σᵢᵅ σⱼᵅ⟩
      ].

The on-site Pauli identity gives

    ⟨(σᵢᵅ)²⟩ = 1.

If `connected=true`, each two-point function is replaced by the connected
correlator

    ⟨σᵢᵅ σⱼᵅ⟩
    - ⟨σᵢᵅ⟩⟨σⱼᵅ⟩,

including the on-site contribution

    1 - ⟨σᵢᵅ⟩².

By default the momentum grid is

    qₘ = 2πm/N,    m = 0, ..., N-1.

This is the natural reciprocal grid for a periodic chain and also provides a
convenient discrete Fourier grid for an open chain. Custom momentum values can
be supplied with the `momenta` keyword.

The normalization is `1/N`, so a perfectly polarized state has

    Sᵅ(q=0) = N

along its polarization axis.
"""
function structure_factor(
    result::DTWAResult,
    geometry::Chain,
    axis::Symbol;
    connected::Bool=false,
    momenta=nothing,
)
    _check_observable_axis(axis)

    nsites_result = size(result.trajectories, 1)

    geometry.nsites == nsites_result ||
        throw(DimensionMismatch(
            "geometry has $(geometry.nsites) sites, but result has " *
            "$nsites_result sites",
        ))

    qs = _structure_factor_momenta(
        geometry,
        momenta,
        eltype(result.trajectories),
    )

    Tout = promote_type(
        eltype(result.trajectories),
        eltype(qs),
        Float64,
    )

    values = Matrix{Tout}(
        undef,
        length(result.t),
        length(qs),
    )

    _structure_factor_values!(
        values,
        result,
        axis,
        qs,
        connected,
    )

    return StructureFactor(
        collect(Tout, qs),
        values,
        axis,
        connected,
    )
end


function _structure_factor_momenta(
    geometry::Chain,
    ::Nothing,
    ::Type{T},
) where {T<:Real}
    Tout = promote_type(T, Float64)
    nsites = geometry.nsites

    return Tout[
        2 * pi * m / nsites
        for m in 0:(nsites - 1)
    ]
end


function _structure_factor_momenta(
    ::Chain,
    momenta,
    ::Type{T},
) where {T<:Real}
    qs = collect(momenta)

    isempty(qs) &&
        throw(ArgumentError(
            "momenta must not be empty",
        ))

    all(q -> q isa Real, qs) ||
        throw(ArgumentError(
            "all momentum values must be real",
        ))

    Tout = promote_type(
        T,
        mapreduce(typeof, promote_type, qs),
        Float64,
    )

    return Tout.(qs)
end


function _structure_factor_values!(
    values::AbstractMatrix,
    result::DTWAResult,
    axis::Symbol,
    momenta::AbstractVector,
    connected::Bool,
)
    axis_index = _axis_index(axis)

    nsites = size(result.trajectories, 1)
    ntimes = length(result.t)
    ntraj = ntrajectories(result)

    Tout = eltype(values)

    @inbounds for time_index in 1:ntimes
        site_means = connected ?
            Vector{Tout}(undef, nsites) :
            nothing

        if connected
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

                site_means[site] =
                    mean / ntraj
            end
        end

        for momentum_index in eachindex(momenta)
            q = momenta[momentum_index]

            total = zero(Tout)

            # On-site contribution.
            if connected
                for site in 1:nsites
                    mean = site_means[site]
                    total += one(Tout) - mean * mean
                end
            else
                total = Tout(nsites)
            end

            # Distinct-site contribution. We use the actual signed lattice
            # displacement j-i before taking the cosine, rather than a
            # distance-averaged correlation profile.
            for i in 1:(nsites - 1)
                for j in (i + 1):nsites
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
                        pair_mean -=
                            site_means[i] *
                            site_means[j]
                    end

                    total += (
                        2 *
                        cos(q * (j - i)) *
                        pair_mean
                    )
                end
            end

            values[
                time_index,
                momentum_index,
            ] = total / nsites
        end
    end

    return nothing
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
