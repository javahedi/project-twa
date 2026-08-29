"""
Sparse CTWA representation of a Hamiltonian after a physical spin model has
been partitioned into equal-size clusters.

The compiler implements the central CTWA mapping

    H_W =
        sum_c,μ B[c,μ] x[c,μ]
        +
        sum_{c<c',μ,ν} J[c,c',μ,ν] x[c,μ] x[c',ν].

Physical terms are classified as follows:

* on-site fields become linear cluster terms;
* pair interactions whose sites lie in the same cluster also become linear
  cluster terms, but on multi-site Pauli-string generators;
* pair interactions crossing a cluster boundary become sparse bilinear terms
  connecting two cluster generators.

No dense cluster coupling matrices are constructed.
"""


"""
    CTWALocalTerm(cluster, generator, coefficient)

One nonzero linear term

    coefficient * x_cluster^generator

in the CTWA Weyl Hamiltonian.
"""
struct CTWALocalTerm{T<:Real}
    cluster::Int
    generator::Int
    coefficient::T
end


"""
    CTWAInterClusterTerm(left_cluster, left_generator,
                         right_cluster, right_generator, coefficient)

One nonzero bilinear inter-cluster term

    coefficient *
    x_left_cluster^left_generator *
    x_right_cluster^right_generator.

The compiler always stores `left_cluster < right_cluster`.
"""
struct CTWAInterClusterTerm{T<:Real}
    left_cluster::Int
    left_generator::Int
    right_cluster::Int
    right_generator::Int
    coefficient::T
end


"""
    CompiledCTWAHamiltonian

Sparse cluster representation of a physical Hamiltonian.

`local_terms` contains fields and intra-cluster interactions.
`intercluster_terms` contains only interactions crossing cluster boundaries.
"""
struct CompiledCTWAHamiltonian{T<:Real}
    basis::PauliStringBasis
    clustering::Clustering
    local_terms::Vector{CTWALocalTerm{T}}
    intercluster_terms::Vector{CTWAInterClusterTerm{T}}
end


"""
    local_terms(compiled)

Return the sparse linear cluster terms of a compiled CTWA Hamiltonian.
"""
local_terms(
    compiled::CompiledCTWAHamiltonian,
) = compiled.local_terms


"""
    intercluster_terms(compiled)

Return the sparse bilinear cluster terms of a compiled CTWA Hamiltonian.
"""
intercluster_terms(
    compiled::CompiledCTWAHamiltonian,
) = compiled.intercluster_terms


"""
    compile_ctwa_hamiltonian(model, clustering; T=Float64)

Compile a normalized physical `SpinModel` into the sparse CTWA cluster
Hamiltonian associated with `clustering`.

The coefficient type can be selected explicitly with `T`. This keeps the
compiler compatible with later Float32/Float64 simulation choices while
avoiding abstractly typed sparse term containers.
"""
function compile_ctwa_hamiltonian(
    model::SpinModel,
    clustering::Clustering;
    T::Type{<:Real}=Float64,
)
    model.geometry.nsites == clustering.nsites ||
        throw(DimensionMismatch(
            "model has $(model.geometry.nsites) sites but clustering has " *
            "$(clustering.nsites)",
        ))

    basis = PauliStringBasis(
        clustering.cluster_size,
    )

    local_accumulator =
        Dict{Tuple{Int,Int},T}()

    inter_accumulator =
        Dict{NTuple{4,Int},T}()

    _compile_ctwa_fields!(
        local_accumulator,
        model,
        basis,
        clustering,
        T,
    )

    _compile_ctwa_pairs!(
        local_accumulator,
        inter_accumulator,
        model,
        basis,
        clustering,
        T,
    )

    local_sparse =
        CTWALocalTerm{T}[]

    for ((cluster, generator), coefficient) in local_accumulator
        iszero(coefficient) && continue

        push!(
            local_sparse,
            CTWALocalTerm{T}(
                cluster,
                generator,
                coefficient,
            ),
        )
    end

    sort!(
        local_sparse;
        by = term -> (
            term.cluster,
            term.generator,
        ),
    )

    inter_sparse =
        CTWAInterClusterTerm{T}[]

    for (
        (
            left_cluster,
            left_generator,
            right_cluster,
            right_generator,
        ),
        coefficient,
    ) in inter_accumulator
        iszero(coefficient) && continue

        push!(
            inter_sparse,
            CTWAInterClusterTerm{T}(
                left_cluster,
                left_generator,
                right_cluster,
                right_generator,
                coefficient,
            ),
        )
    end

    sort!(
        inter_sparse;
        by = term -> (
            term.left_cluster,
            term.right_cluster,
            term.left_generator,
            term.right_generator,
        ),
    )

    return CompiledCTWAHamiltonian{T}(
        basis,
        clustering,
        local_sparse,
        inter_sparse,
    )
end


function _compile_ctwa_fields!(
    accumulator::Dict{Tuple{Int,Int},T},
    model::SpinModel,
    basis::PauliStringBasis,
    clustering::Clustering,
    ::Type{T},
) where {T<:Real}
    for term in terms(model)
        term isa AbstractField || continue

        hx, hy, hz = field_components(term)

        _compile_uniform_field_component!(
            accumulator,
            basis,
            clustering,
            1,
            T(hx),
        )

        _compile_uniform_field_component!(
            accumulator,
            basis,
            clustering,
            2,
            T(hy),
        )

        _compile_uniform_field_component!(
            accumulator,
            basis,
            clustering,
            3,
            T(hz),
        )
    end

    return nothing
end


function _compile_uniform_field_component!(
    accumulator::Dict{Tuple{Int,Int},T},
    basis::PauliStringBasis,
    clustering::Clustering,
    digit::Int,
    coefficient::T,
) where {T<:Real}
    iszero(coefficient) &&
        return nothing

    for site in 1:clustering.nsites
        cluster, position =
            site_cluster_position(
                clustering,
                site,
            )

        generator =
            _single_site_pauli_code(
                basis,
                position,
                digit,
            )

        key = (
            cluster,
            generator,
        )

        accumulator[key] =
            get(
                accumulator,
                key,
                zero(T),
            ) + coefficient
    end

    return nothing
end


function _compile_ctwa_pairs!(
    local_accumulator::Dict{Tuple{Int,Int},T},
    inter_accumulator::Dict{NTuple{4,Int},T},
    model::SpinModel,
    basis::PauliStringBasis,
    clustering::Clustering,
    ::Type{T},
) where {T<:Real}
    nsites = clustering.nsites

    for i in 1:(nsites - 1)
        cluster_i, position_i =
            site_cluster_position(
                clustering,
                i,
            )

        for j in (i + 1):nsites
            Jx, Jy, Jz =
                coupling_components(
                    model,
                    i,
                    j,
                )

            _compile_ctwa_pair_component!(
                local_accumulator,
                inter_accumulator,
                basis,
                clustering,
                cluster_i,
                position_i,
                j,
                1,
                T(Jx),
            )

            _compile_ctwa_pair_component!(
                local_accumulator,
                inter_accumulator,
                basis,
                clustering,
                cluster_i,
                position_i,
                j,
                2,
                T(Jy),
            )

            _compile_ctwa_pair_component!(
                local_accumulator,
                inter_accumulator,
                basis,
                clustering,
                cluster_i,
                position_i,
                j,
                3,
                T(Jz),
            )
        end
    end

    return nothing
end


function _compile_ctwa_pair_component!(
    local_accumulator::Dict{Tuple{Int,Int},T},
    inter_accumulator::Dict{NTuple{4,Int},T},
    basis::PauliStringBasis,
    clustering::Clustering,
    cluster_i::Int,
    position_i::Int,
    j::Int,
    digit::Int,
    coefficient::T,
) where {T<:Real}
    iszero(coefficient) &&
        return nothing

    cluster_j, position_j =
        site_cluster_position(
            clustering,
            j,
        )

    if cluster_i == cluster_j
        generator =
            _two_site_pauli_code(
                basis,
                position_i,
                digit,
                position_j,
                digit,
            )

        key = (
            cluster_i,
            generator,
        )

        local_accumulator[key] =
            get(
                local_accumulator,
                key,
                zero(T),
            ) + coefficient

        return nothing
    end

    generator_i =
        _single_site_pauli_code(
            basis,
            position_i,
            digit,
        )

    generator_j =
        _single_site_pauli_code(
            basis,
            position_j,
            digit,
        )

    # Physical site ordering normally implies cluster_i < cluster_j for
    # contiguous clusters, but canonicalize explicitly so this invariant
    # remains true if the compiler is generalized later.
    if cluster_i < cluster_j
        key = (
            cluster_i,
            generator_i,
            cluster_j,
            generator_j,
        )
    else
        key = (
            cluster_j,
            generator_j,
            cluster_i,
            generator_i,
        )
    end

    inter_accumulator[key] =
        get(
            inter_accumulator,
            key,
            zero(T),
        ) + coefficient

    return nothing
end


"""
Construct the canonical basis code for a one-site Pauli operator without
creating an intermediate tuple.

Site 1 is the most-significant base-4 digit.
"""
@inline function _single_site_pauli_code(
    basis::PauliStringBasis,
    position::Int,
    digit::Int,
)
    1 <= position <= basis.cluster_size ||
        throw(BoundsError(
            1:basis.cluster_size,
            position,
        ))

    1 <= digit <= 3 ||
        throw(ArgumentError(
            "non-identity Pauli digit must be in 1:3; got $digit",
        ))

    power =
        basis.cluster_size -
        position

    return digit * 4^power
end


"""
Construct the canonical basis code for a two-site Pauli string without
creating an intermediate tuple.
"""
@inline function _two_site_pauli_code(
    basis::PauliStringBasis,
    position_a::Int,
    digit_a::Int,
    position_b::Int,
    digit_b::Int,
)
    position_a != position_b ||
        throw(ArgumentError(
            "two-site Pauli operator requires distinct local positions",
        ))

    return (
        _single_site_pauli_code(
            basis,
            position_a,
            digit_a,
        ) +
        _single_site_pauli_code(
            basis,
            position_b,
            digit_b,
        )
    )
end
