"""
    CTWADiscreteSamplingPlan

Discrete cluster-TWA sampling data.

Each physical spin is sampled from the single-spin discrete phase space.
Every cluster Pauli-string coordinate is then obtained as the product of the
corresponding sampled local Pauli components.

For the canonical cluster basis,

    0 = I
    1 = X
    2 = Y
    3 = Z,

site 1 is the most-significant base-4 digit.
"""
struct CTWADiscreteSamplingPlan{T<:AbstractFloat}
    basis::PauliStringBasis
    clustering::Clustering
    axes::Matrix{Int}
    signs::Matrix{Int}
end


function compile_ctwa_sampling(
    state::AbstractState,
    clustering::Clustering,
    ::DiscreteSampling;
    T::Type{<:AbstractFloat}=Float64,
)
    basis =
        PauliStringBasis(
            clustering.cluster_size,
        )

    nclusters =
        cluster_count(
            clustering,
        )

    axes =
        Matrix{Int}(
            undef,
            clustering.cluster_size,
            nclusters,
        )

    signs =
        Matrix{Int}(
            undef,
            clustering.cluster_size,
            nclusters,
        )

    for cluster in 1:nclusters
        sites =
            cluster_sites(
                clustering,
                cluster,
            )

        for (
            position,
            site,
        ) in enumerate(sites)
            direction =
                state_direction(
                    state,
                    site,
                    clustering.nsites,
                )

            axis, sign =
                _state_axis_sign(
                    direction,
                )

            axes[
                position,
                cluster,
            ] =
                axis

            signs[
                position,
                cluster,
            ] =
                sign
        end
    end

    return CTWADiscreteSamplingPlan{T}(
        basis,
        clustering,
        axes,
        signs,
    )
end


function sample_ctwa(
    plan::CTWADiscreteSamplingPlan{T};
    rng=nothing,
) where {T<:AbstractFloat}
    x =
        Matrix{T}(
            undef,
            basis_size(plan.basis),
            cluster_count(plan.clustering),
        )

    sample_ctwa!(
        x,
        plan;
        rng=rng,
    )

    return x
end


function sample_ctwa!(
    x::AbstractMatrix{T},
    plan::CTWADiscreteSamplingPlan{T};
    rng=nothing,
) where {T<:AbstractFloat}
    expected =
        (
            basis_size(plan.basis),
            cluster_count(plan.clustering),
        )

    size(x) == expected ||
        throw(DimensionMismatch(
            "CTWA sample must have shape $expected; got $(size(x))",
        ))

    random =
        isnothing(rng) ?
        Random.default_rng() :
        rng

    n =
        plan.clustering.cluster_size

    local_components =
        Matrix{T}(
            undef,
            3,
            n,
        )

    for cluster in 1:cluster_count(plan.clustering)

        for position in 1:n
            axis =
                plan.axes[
                    position,
                    cluster,
                ]

            sign =
                plan.signs[
                    position,
                    cluster,
                ]

            _sample_discrete_local_spin!(
                local_components,
                position,
                axis,
                sign,
                random,
            )
        end

        @inbounds for generator in 1:basis_size(plan.basis)
            code =
                generator

            value =
                one(T)

            # Site 1 is the most-significant base-4 digit.
            for position in n:-1:1
                digit =
                    code % 4

                code ÷= 4

                if digit != 0
                    value *=
                        local_components[
                            digit,
                            position,
                        ]
                end
            end

            x[
                generator,
                cluster,
            ] =
                value
        end
    end

    return x
end


@inline function _sample_discrete_local_spin!(
    components::AbstractMatrix{T},
    position::Int,
    axis::Int,
    sign::Int,
    rng,
) where {T<:AbstractFloat}

    # Two independent ±1 variables for the two transverse directions.
    r1 =
        rand(rng, Bool) ?
        one(T) :
        -one(T)

    r2 =
        rand(rng, Bool) ?
        one(T) :
        -one(T)

    if axis == 1
        components[1, position] = T(sign)
        components[2, position] = r1
        components[3, position] = r2

    elseif axis == 2
        components[1, position] = r1
        components[2, position] = T(sign)
        components[3, position] = r2

    elseif axis == 3
        components[1, position] = r1
        components[2, position] = r2
        components[3, position] = T(sign)

    else
        throw(ArgumentError(
            "Pauli axis digit must be 1, 2, or 3; got $axis",
        ))
    end

    return nothing
end