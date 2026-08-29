"""
    pauli_product_digit(a, b)

Multiply two single-site Pauli operators encoded by

    0 => I
    1 => X
    2 => Y
    3 => Z

and return `(phase, digit)` such that

    σₐ σ_b = phase * σ_digit.

The phase is one of `1`, `-1`, `im`, `-im`.
"""
@inline function pauli_product_digit(
    a::Integer,
    b::Integer,
)
    _check_pauli_digit(a)
    _check_pauli_digit(b)

    ai = Int(a)
    bi = Int(b)

    ai == 0 && return (1, bi)
    bi == 0 && return (1, ai)
    ai == bi && return (1, 0)

    if ai == 1
        bi == 2 && return (im, 3)   # X Y = +i Z
        return (-im, 2)             # X Z = -i Y
    elseif ai == 2
        bi == 3 && return (im, 1)   # Y Z = +i X
        return (-im, 3)             # Y X = -i Z
    else
        bi == 1 && return (im, 2)   # Z X = +i Y
        return (-im, 1)             # Z Y = -i X
    end
end


"""
    _pauli_product_digit_unchecked(a, b)

Internal allocation-free single-site Pauli multiplication.

Inputs must already be valid Pauli digits in `0:3`.
"""
@inline function _pauli_product_digit_unchecked(
    a::Int,
    b::Int,
)
    a == 0 && return (1, b)
    b == 0 && return (1, a)
    a == b && return (1, 0)

    if a == 1
        b == 2 && return (im, 3)
        return (-im, 2)
    elseif a == 2
        b == 3 && return (im, 1)
        return (-im, 3)
    else
        b == 1 && return (im, 2)
        return (-im, 1)
    end
end


"""
    pauli_product(basis, left, right)

Multiply two traceless Pauli-string basis elements.

Returns `(phase, code)` such that

    P_left P_right = phase * P_code.

`code` is the raw canonical Pauli code and may be zero when the product is the
all-identity string.

The implementation works directly on the base-4 codes. It does not construct
Pauli-digit tuples or dense matrices.
"""
function pauli_product(
    basis::PauliStringBasis,
    left::Integer,
    right::Integer,
)
    _check_basis_index(basis, left)
    _check_basis_index(basis, right)

    left_code = Int(left)
    right_code = Int(right)

    # We process from the least-significant base-4 digit (last cluster site)
    # toward the most-significant digit. `place` restores the canonical code.
    place = 1
    result_code = 0
    phase = 1 + 0im

    @inbounds for _ in 1:basis.cluster_size
        left_digit = left_code % 4
        right_digit = right_code % 4

        local_phase, result_digit =
            _pauli_product_digit_unchecked(
                left_digit,
                right_digit,
            )

        phase *= local_phase
        result_code += result_digit * place

        left_code ÷= 4
        right_code ÷= 4
        place *= 4
    end

    return (phase, result_code)
end


"""
    pauli_commutator(basis, left, right)

Return the sparse commutator of two traceless Pauli strings.

The result is represented as

    (coefficient, code)

with

    [P_left, P_right] = coefficient * P_code.

For Pauli strings:

- commuting strings return `(0im, 0)`;
- anticommuting strings return `(±2im, code)`.

Only one pass through the encoded strings is required.
"""
function pauli_commutator(
    basis::PauliStringBasis,
    left::Integer,
    right::Integer,
)
    _check_basis_index(basis, left)
    _check_basis_index(basis, right)

    left_code = Int(left)
    right_code = Int(right)

    place = 1
    result_code = 0
    phase = 1 + 0im

    # Two Pauli strings anticommute iff the number of sites on which both
    # operators are non-identity and different is odd.
    anticommute_parity = false

    @inbounds for _ in 1:basis.cluster_size
        left_digit = left_code % 4
        right_digit = right_code % 4

        local_phase, result_digit =
            _pauli_product_digit_unchecked(
                left_digit,
                right_digit,
            )

        phase *= local_phase
        result_code += result_digit * place

        if left_digit != 0 &&
           right_digit != 0 &&
           left_digit != right_digit
            anticommute_parity = !anticommute_parity
        end

        left_code ÷= 4
        right_code ÷= 4
        place *= 4
    end

    anticommute_parity ||
        return (0im, 0)

    return (2 * phase, result_code)
end


"""
    pauli_anticommutes(basis, left, right)

Return `true` if two Pauli strings anticommute and `false` if they commute.

This checks only the local anticommutation parity and does not construct the
product string.
"""
function pauli_anticommutes(
    basis::PauliStringBasis,
    left::Integer,
    right::Integer,
)
    _check_basis_index(basis, left)
    _check_basis_index(basis, right)

    left_code = Int(left)
    right_code = Int(right)
    parity = false

    @inbounds for _ in 1:basis.cluster_size
        left_digit = left_code % 4
        right_digit = right_code % 4

        if left_digit != 0 &&
           right_digit != 0 &&
           left_digit != right_digit
            parity = !parity
        end

        left_code ÷= 4
        right_code ÷= 4
    end

    return parity
end
