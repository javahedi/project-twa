"""
Matrix-free Gaussian sampling for closed-system CTWA.

For an initial cluster product state |ψ>, CTWA samples phase-space variables
with the exact first moments and symmetrized second moments of the cluster
Pauli-string generators:

    E[x_μ] = <P_μ>

and

    Cov(x_μ, x_ν)
        = 1/2 <{P_μ, P_ν}> - <P_μ><P_ν>.

For product eigenstates of local Pauli operators, a Pauli string acting on
|ψ> maps to exactly one product-basis flip pattern. This lets us factor the
covariance combinatorially instead of constructing a dense
(4^n - 1) × (4^n - 1) matrix.

Each nonzero flip mask gets two independent standard Gaussian variables.
Thus a cluster of size n requires only

    2 * (2^n - 1)

Gaussian random variables, while still producing all 4^n - 1 CTWA
coordinates.
"""


"""
    CTWASamplingPlan

Precompiled matrix-free Gaussian sampling data for one physical product state
and one equal-size `Clustering`.

Arrays have shape

    (basis_dimension, nclusters).

For generator μ in cluster c:

* `means[μ,c]` is the exact quantum expectation <P_μ>;
* `flip_masks[μ,c] == 0` means the generator is deterministic;
* otherwise its fluctuation is
      amplitude_re * ξ_mask + amplitude_im * η_mask.
"""
struct CTWASamplingPlan{T<:Real}
    basis::PauliStringBasis
    clustering::Clustering
    means::Matrix{T}
    flip_masks::Matrix{Int}
    amplitude_re::Matrix{T}
    amplitude_im::Matrix{T}
end


"""
    compile_ctwa_sampling(state, clustering; T=Float64)

Compile a matrix-free CTWA Gaussian sampler for a product state described by
`state_direction`.

The state must be a product of local Pauli eigenstates, which is the state
family currently represented by the package's state layer.
"""
function compile_ctwa_sampling(
    state::AbstractState,
    clustering::Clustering;
    T::Type{<:AbstractFloat}=Float64,
)
    basis =
        PauliStringBasis(
            clustering.cluster_size,
        )

    dimension =
        basis_size(
            basis,
        )

    nclusters =
        cluster_count(
            clustering,
        )

    means =
        zeros(
            T,
            dimension,
            nclusters,
        )

    flip_masks =
        zeros(
            Int,
            dimension,
            nclusters,
        )

    amplitude_re =
        zeros(
            T,
            dimension,
            nclusters,
        )

    amplitude_im =
        zeros(
            T,
            dimension,
            nclusters,
        )

    for cluster in 1:nclusters
        first_site =
            first(
                cluster_sites(
                    clustering,
                    cluster,
                ),
            )

        axes =
            Vector{Int}(
                undef,
                clustering.cluster_size,
            )

        signs =
            Vector{Int}(
                undef,
                clustering.cluster_size,
            )

        for position in 1:clustering.cluster_size
            site =
                first_site +
                position -
                1

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

            axes[position] =
                axis

            signs[position] =
                sign
        end

        for generator in 1:dimension
            code = generator

            mask = 0
            amplitude =
                Complex{T}(
                    one(T),
                    zero(T),
                )

            # Decode from the least-significant base-4 digit while mapping
            # back to physical position. Site 1 is the most-significant digit.
            for position in clustering.cluster_size:-1:1
                digit =
                    code % 4

                code ÷= 4

                flips, local_amplitude =
                    _local_pauli_on_eigenstate(
                        T,
                        axes[position],
                        signs[position],
                        digit,
                    )

                amplitude *=
                    local_amplitude

                if flips
                    mask |=
                        1 << (
                            position - 1
                        )
                end
            end

            if iszero(mask)
                # A no-flip Pauli string leaves the product state invariant
                # up to its real eigenvalue ±1.
                means[
                    generator,
                    cluster,
                ] =
                    real(amplitude)
            else
                flip_masks[
                    generator,
                    cluster,
                ] =
                    mask

                amplitude_re[
                    generator,
                    cluster,
                ] =
                    real(amplitude)

                amplitude_im[
                    generator,
                    cluster,
                ] =
                    imag(amplitude)
            end
        end
    end

    return CTWASamplingPlan{T}(
        basis,
        clustering,
        means,
        flip_masks,
        amplitude_re,
        amplitude_im,
    )
end


"""
    sample_ctwa(plan; rng=nothing)

Draw one CTWA phase-space sample with shape

    (basis_dimension, nclusters).

The output element type is the numeric type used to compile `plan`.
"""
function sample_ctwa(
    plan::CTWASamplingPlan{T};
    rng=nothing,
) where {T<:AbstractFloat}
    random =
        isnothing(rng) ?
        Random.default_rng() :
        rng

    x =
        similar(
            plan.means,
        )

    sample_ctwa!(
        x,
        plan;
        rng=random,
    )

    return x
end


"""
    sample_ctwa!(x, plan; rng=nothing)

Fill a preallocated CTWA phase-space state in place.
"""
function sample_ctwa!(
    x::AbstractMatrix{T},
    plan::CTWASamplingPlan{T};
    rng=nothing,
) where {T<:AbstractFloat}
    expected =
        size(
            plan.means,
        )

    size(x) == expected ||
        throw(DimensionMismatch(
            "CTWA sample must have shape $expected; got $(size(x))",
        ))

    random =
        isnothing(rng) ?
        Random.default_rng() :
        rng

    nmasks =
        (1 << plan.clustering.cluster_size) -
        1

    ξ =
        Vector{T}(
            undef,
            nmasks,
        )

    η =
        Vector{T}(
            undef,
            nmasks,
        )

    for cluster in 1:cluster_count(plan.clustering)
        Random.randn!(
            random,
            ξ,
        )

        Random.randn!(
            random,
            η,
        )

        @inbounds for generator in 1:basis_size(plan.basis)
            mask =
                plan.flip_masks[
                    generator,
                    cluster,
                ]

            if iszero(mask)
                x[
                    generator,
                    cluster,
                ] =
                    plan.means[
                        generator,
                        cluster,
                    ]
            else
                x[
                    generator,
                    cluster,
                ] =
                    plan.means[
                        generator,
                        cluster,
                    ] +
                    plan.amplitude_re[
                        generator,
                        cluster,
                    ] *
                    ξ[mask] +
                    plan.amplitude_im[
                        generator,
                        cluster,
                    ] *
                    η[mask]
            end
        end
    end

    return x
end


"""
Convert a cardinal state direction `(sx,sy,sz)` into Pauli axis digit and
eigenvalue sign.
"""
@inline function _state_axis_sign(
    direction,
)
    sx, sy, sz =
        direction

    if sx == 1 && sy == 0 && sz == 0
        return 1, 1
    elseif sx == -1 && sy == 0 && sz == 0
        return 1, -1
    elseif sx == 0 && sy == 1 && sz == 0
        return 2, 1
    elseif sx == 0 && sy == -1 && sz == 0
        return 2, -1
    elseif sx == 0 && sy == 0 && sz == 1
        return 3, 1
    elseif sx == 0 && sy == 0 && sz == -1
        return 3, -1
    end

    throw(ArgumentError(
        "CTWA product-state sampling requires a local Pauli eigenstate " *
        "direction; got $direction",
    ))
end


"""
Return `(flips, amplitude)` for σ_digit acting on a local Pauli eigenstate.

The local orthogonal basis vector is chosen as the opposite eigenstate of the
same Pauli axis. Its phase is arbitrary; this convention is fixed so that
generators sharing the same flip mask have the correct relative phase.

Digits:
    0 = I
    1 = X
    2 = Y
    3 = Z
"""
@inline function _local_pauli_on_eigenstate(
    ::Type{T},
    axis::Int,
    sign::Int,
    digit::Int,
) where {T<:AbstractFloat}
    0 <= digit <= 3 ||
        throw(ArgumentError(
            "Pauli digit must be in 0:3; got $digit",
        ))

    sign == 1 || sign == -1 ||
        throw(ArgumentError(
            "Pauli eigenvalue sign must be ±1; got $sign",
        ))

    one_complex =
        Complex{T}(
            one(T),
            zero(T),
        )

    if digit == 0
        return false, one_complex
    end

    if digit == axis
        return false, Complex{T}(
            T(sign),
            zero(T),
        )
    end

    # Off-axis matrix elements in the local basis
    # {|axis,sign>, |axis,-sign>}.
    if axis == 1
        # X eigenstate:
        #   <−s_x|Y|s_x> = -i s
        #   <−s_x|Z|s_x> = 1
        if digit == 2
            return true, Complex{T}(
                zero(T),
                T(-sign),
            )
        else
            return true, one_complex
        end
    elseif axis == 2
        # Y eigenstate:
        #   <−s_y|X|s_y> = +i s
        #   <−s_y|Z|s_y> = 1
        if digit == 1
            return true, Complex{T}(
                zero(T),
                T(sign),
            )
        else
            return true, one_complex
        end
    elseif axis == 3
        # Z eigenstate:
        #   <−s_z|X|s_z> = 1
        #   <−s_z|Y|s_z> = +i s
        if digit == 1
            return true, one_complex
        else
            return true, Complex{T}(
                zero(T),
                T(sign),
            )
        end
    end

    throw(ArgumentError(
        "Pauli axis must be in 1:3; got $axis",
    ))
end
