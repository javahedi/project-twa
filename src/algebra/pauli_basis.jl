"""
    PauliStringBasis(cluster_size)

Canonical traceless Pauli-string basis for a cluster of spin-1/2 sites.

Each local operator is encoded by a base-4 digit

    0 => I
    1 => X
    2 => Y
    3 => Z

For a cluster of `n` sites, a Pauli string

    (a₁, a₂, ..., aₙ)

is assigned the integer code

    code = a₁ 4^(n-1) + a₂ 4^(n-2) + ... + aₙ.

Site 1 is therefore the most-significant base-4 digit.

The all-identity string has code `0` and is excluded from the traceless
dynamical basis. Consequently the basis contains

    4^n - 1

operators, indexed by exactly the same integers as their nonzero codes:

    basis index = Pauli-string code,    1 <= index <= 4^n - 1.

For example, for a two-site cluster the basis begins

    1  => (0, 1)  = IX
    2  => (0, 2)  = IY
    3  => (0, 3)  = IZ
    4  => (1, 0)  = XI
    5  => (1, 1)  = XX

This direct arithmetic mapping avoids dictionary lookup and provides one
canonical ordering shared by CTWA Hamiltonian compilation, algebra, sampling,
dynamics, and observables.
"""
struct PauliStringBasis
    cluster_size::Int
    dimension::Int

    function PauliStringBasis(cluster_size::Integer)
        cluster_size > 0 ||
            throw(ArgumentError(
                "cluster_size must be positive; got $cluster_size",
            ))

        n = Int(cluster_size)

        # checked multiplication avoids silently overflowing Int for an
        # impossible cluster size.
        operator_space_dimension = Base.checked_pow(4, n)
        dimension = operator_space_dimension - 1

        new(n, dimension)
    end
end


"""
    basis_size(basis)

Number of traceless Pauli strings in `basis`, equal to `4^n - 1`.
"""
basis_size(basis::PauliStringBasis) = basis.dimension


Base.length(basis::PauliStringBasis) = basis.dimension


"""
    pauli_code(basis, digits)

Return the base-4 integer code of a Pauli string.

`digits` must contain exactly `basis.cluster_size` entries, each in `0:3`.

The identity string is allowed here and has code `0`. Use `pauli_index` when
the string must belong to the traceless CTWA basis.
"""
function pauli_code(
    basis::PauliStringBasis,
    digits,
)
    length(digits) == basis.cluster_size ||
        throw(DimensionMismatch(
            "expected $(basis.cluster_size) Pauli digits, got " *
            "$(length(digits))",
        ))

    code = 0

    for digit in digits
        _check_pauli_digit(digit)
        code = 4 * code + Int(digit)
    end

    return code
end


"""
    pauli_index(basis, digits)

Return the canonical basis index of a non-identity Pauli string.

Because the traceless basis excludes only the all-identity string, the basis
index is identical to the nonzero base-4 code.
"""
function pauli_index(
    basis::PauliStringBasis,
    digits,
)
    code = pauli_code(basis, digits)

    code != 0 ||
        throw(ArgumentError(
            "the all-identity Pauli string is not part of the traceless basis",
        ))

    return code
end


"""
    pauli_digits(basis, index)

Decode a basis index into its canonical tuple of Pauli digits.

The returned tuple has length `basis.cluster_size`, with site 1 first.
"""
function pauli_digits(
    basis::PauliStringBasis,
    index::Integer,
)
    _check_basis_index(basis, index)

    value = Int(index)
    n = basis.cluster_size

    return ntuple(n) do site
        power = Base.checked_pow(4, n - site)
        (value ÷ power) % 4
    end
end


"""
    pauli_symbol(digit)

Convert an integer Pauli digit to its symbolic label.
"""
function pauli_symbol(digit::Integer)
    digit == 0 && return :I
    digit == 1 && return :X
    digit == 2 && return :Y
    digit == 3 && return :Z

    throw(ArgumentError(
        "Pauli digit must be one of 0, 1, 2, 3; got $digit",
    ))
end


"""
    pauli_symbols(basis, index)

Decode a basis index into symbolic Pauli labels such as `(:I, :X)`.
"""
function pauli_symbols(
    basis::PauliStringBasis,
    index::Integer,
)
    digits = pauli_digits(basis, index)
    return map(pauli_symbol, digits)
end


@inline function _check_pauli_digit(digit)
    digit isa Integer ||
        throw(ArgumentError(
            "Pauli digits must be integers in 0:3; got $digit",
        ))

    0 <= digit <= 3 ||
        throw(ArgumentError(
            "Pauli digit must be in 0:3; got $digit",
        ))

    return nothing
end


@inline function _check_basis_index(
    basis::PauliStringBasis,
    index::Integer,
)
    1 <= index <= basis.dimension ||
        throw(BoundsError(
            1:basis.dimension,
            index,
        ))

    return nothing
end
