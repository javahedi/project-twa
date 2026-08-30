"""
CTWA deterministic equations of motion using sparse precompiled Pauli algebra.

State layout:

    x[generator, cluster]

with shape

    (4^cluster_size - 1, nclusters).
"""


"""
    ctwa_rhs!(dx, x, compiled, algebra, t)

Evaluate deterministic CTWA dynamics in place using a precompiled sparse
commutator cache.
"""
function ctwa_rhs!(
    dx::AbstractMatrix,
    x::AbstractMatrix,
    compiled::CompiledCTWAHamiltonian,
    algebra::CompiledCTWAAlgebra,
    t,
)
    _check_ctwa_state_shape(
        dx,
        compiled,
    )

    _check_ctwa_state_shape(
        x,
        compiled,
    )

    size(dx) == size(x) ||
        throw(DimensionMismatch(
            "dx and x must have identical shapes; got " *
            "$(size(dx)) and $(size(x))",
        ))

    fill!(
        dx,
        zero(eltype(dx)),
    )

    @inbounds for term in compiled.local_terms
        action =
            commutator_action(
                algebra,
                term.generator,
            )

        _apply_ctwa_action!(
            dx,
            x,
            term.cluster,
            action,
            term.coefficient,
        )
    end

    @inbounds for term in compiled.intercluster_terms
        left_field =
            term.coefficient *
            x[
                term.right_generator,
                term.right_cluster,
            ]

        right_field =
            term.coefficient *
            x[
                term.left_generator,
                term.left_cluster,
            ]

        left_action =
            commutator_action(
                algebra,
                term.left_generator,
            )

        right_action =
            commutator_action(
                algebra,
                term.right_generator,
            )

        _apply_ctwa_action!(
            dx,
            x,
            term.left_cluster,
            left_action,
            left_field,
        )

        _apply_ctwa_action!(
            dx,
            x,
            term.right_cluster,
            right_action,
            right_field,
        )
    end

    return nothing
end


"""
    ctwa_rhs!(dx, x, compiled, t)

Reference-compatible CTWA RHS.

This fallback computes Pauli commutators directly and is retained for testing
and correctness comparisons. Production trajectory evolution should compile
the algebra once and call the five-argument method.
"""
function ctwa_rhs!(
    dx::AbstractMatrix,
    x::AbstractMatrix,
    compiled::CompiledCTWAHamiltonian,
    t,
)
    _check_ctwa_state_shape(
        dx,
        compiled,
    )

    _check_ctwa_state_shape(
        x,
        compiled,
    )

    size(dx) == size(x) ||
        throw(DimensionMismatch(
            "dx and x must have identical shapes; got " *
            "$(size(dx)) and $(size(x))",
        ))

    fill!(
        dx,
        zero(eltype(dx)),
    )

    @inbounds for term in compiled.local_terms
        _apply_ctwa_generator_reference!(
            dx,
            x,
            compiled.basis,
            term.cluster,
            term.generator,
            term.coefficient,
        )
    end

    @inbounds for term in compiled.intercluster_terms
        left_field =
            term.coefficient *
            x[
                term.right_generator,
                term.right_cluster,
            ]

        right_field =
            term.coefficient *
            x[
                term.left_generator,
                term.left_cluster,
            ]

        _apply_ctwa_generator_reference!(
            dx,
            x,
            compiled.basis,
            term.left_cluster,
            term.left_generator,
            left_field,
        )

        _apply_ctwa_generator_reference!(
            dx,
            x,
            compiled.basis,
            term.right_cluster,
            term.right_generator,
            right_field,
        )
    end

    return nothing
end


@inline function _apply_ctwa_action!(
    dx::AbstractMatrix,
    x::AbstractMatrix,
    cluster::Int,
    action::CTWACommutatorAction,
    field_coefficient,
)
    iszero(field_coefficient) &&
        return nothing

    @inbounds for entry in action.entries
        dx[
            entry.state_generator,
            cluster,
        ] +=
            2 *
            entry.sign *
            field_coefficient *
            x[
                entry.result_generator,
                cluster,
            ]
    end

    return nothing
end


@inline function _apply_ctwa_generator_reference!(
    dx::AbstractMatrix,
    x::AbstractMatrix,
    basis::PauliStringBasis,
    cluster::Int,
    hamiltonian_generator::Int,
    field_coefficient,
)
    iszero(field_coefficient) &&
        return nothing

    dimension = basis_size(
        basis,
    )

    @inbounds for state_generator in 1:dimension
        coefficient, result_generator =
            pauli_commutator(
                basis,
                state_generator,
                hamiltonian_generator,
            )

        iszero(coefficient) &&
            continue

        structure_sign =
            coefficient == 2im ?
            1 :
            -1

        dx[
            state_generator,
            cluster,
        ] +=
            2 *
            structure_sign *
            field_coefficient *
            x[
                result_generator,
                cluster,
            ]
    end

    return nothing
end


@inline function _check_ctwa_state_shape(
    state::AbstractMatrix,
    compiled::CompiledCTWAHamiltonian,
)
    expected =
        (
            basis_size(compiled.basis),
            cluster_count(compiled.clustering),
        )

    size(state) == expected ||
        throw(DimensionMismatch(
            "CTWA state must have shape $expected; got $(size(state))",
        ))

    return nothing
end
