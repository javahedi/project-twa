"""
Sparse precompiled Pauli-commutator actions for CTWA dynamics.

For a Hamiltonian generator ν we need repeated evaluations of

    [P_μ, P_ν] = 2 i f_{μνρ} P_ρ.

For Pauli strings, each μ either commutes with ν or maps to exactly one ρ with
structure sign ±1. We therefore cache only the nonzero action entries.
"""


"""
    CTWACommutatorEntry(state_generator, result_generator, sign)

One nonzero structure-constant action

    [P_state_generator, P_hamiltonian_generator]
        = 2im * sign * P_result_generator,

with `sign ∈ {-1, +1}`.
"""
struct CTWACommutatorEntry
    state_generator::Int
    result_generator::Int
    sign::Int8
end


"""
    CTWACommutatorAction(generator, entries)

Sparse commutator action associated with one Hamiltonian generator.
"""
struct CTWACommutatorAction
    generator::Int
    entries::Vector{CTWACommutatorEntry}
end


"""
    CompiledCTWAAlgebra

Collection of sparse commutator actions used by a compiled CTWA Hamiltonian.

`actions[g]` is either `nothing` if generator `g` never appears in the
Hamiltonian, or a `CTWACommutatorAction`.
"""
struct CompiledCTWAAlgebra
    actions::Vector{Union{Nothing,CTWACommutatorAction}}
end


"""
    compile_ctwa_algebra(compiled)

Precompute sparse commutator actions for every generator that actually appears
in `compiled`.

This is setup work. The resulting cache is intended to remove repeated Pauli
algebra from the hot CTWA RHS loop.
"""
function compile_ctwa_algebra(
    compiled::CompiledCTWAHamiltonian,
)
    basis = compiled.basis
    dimension = basis_size(basis)

    used = falses(
        dimension,
    )

    @inbounds for term in compiled.local_terms
        used[
            term.generator,
        ] = true
    end

    @inbounds for term in compiled.intercluster_terms
        used[
            term.left_generator,
        ] = true

        used[
            term.right_generator,
        ] = true
    end

    actions =
        Vector{
            Union{
                Nothing,
                CTWACommutatorAction,
            }
        }(
            undef,
            dimension,
        )

    fill!(
        actions,
        nothing,
    )

    for generator in 1:dimension
        used[generator] ||
            continue

        entries =
            CTWACommutatorEntry[]

        # For any non-identity Pauli string, exactly half of the full Pauli
        # operator space anticommutes with it. Identity never contributes, so
        # this is a useful capacity hint and avoids repeated vector growth.
        sizehint!(
            entries,
            (dimension + 1) ÷ 2,
        )

        @inbounds for state_generator in 1:dimension
            coefficient, result_generator =
                pauli_commutator(
                    basis,
                    state_generator,
                    generator,
                )

            iszero(coefficient) &&
                continue

            sign =
                coefficient == 2im ?
                Int8(1) :
                Int8(-1)

            push!(
                entries,
                CTWACommutatorEntry(
                    state_generator,
                    result_generator,
                    sign,
                ),
            )
        end

        actions[generator] =
            CTWACommutatorAction(
                generator,
                entries,
            )
    end

    return CompiledCTWAAlgebra(
        actions,
    )
end


"""
    commutator_action(algebra, generator)

Return the cached sparse action for `generator`.

Throws if the generator was not present in the compiled Hamiltonian.
"""
@inline function commutator_action(
    algebra::CompiledCTWAAlgebra,
    generator::Integer,
)
    1 <= generator <= length(algebra.actions) ||
        throw(BoundsError(
            algebra.actions,
            generator,
        ))

    action =
        @inbounds algebra.actions[
            Int(generator)
        ]

    action === nothing &&
        throw(ArgumentError(
            "no cached commutator action for generator $generator",
        ))

    return action::CTWACommutatorAction
end
