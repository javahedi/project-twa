"""
    AbstractField

Abstract supertype for one-body fields in a spin Hamiltonian.

Fields are kept separate from pair interactions because they represent
one-site Hamiltonian terms. This distinction will later allow DTWA and
CTWA backends to compile one-body and two-body contributions efficiently.
"""
abstract type AbstractField <: AbstractHamiltonianTerm end


"""
    Field(axis, strength)

Uniform local field acting along one Cartesian spin axis.

# Arguments

- `axis`: one of `:x`, `:y`, or `:z`.
- `strength`: field strength.

# Examples

    Field(:z, 0.3)
    Field(:x, 1.0)

The field represents a one-body Hamiltonian contribution of the form

    H_field = h * sum_i S_i^axis

with the precise sign convention fixed consistently by the package's
Hamiltonian and equations-of-motion conventions.
"""
struct Field{T<:Real} <: AbstractField
    axis::Symbol
    strength::T

    function Field(
        axis::Symbol,
        strength::T,
    ) where {T<:Real}
        _check_axis(axis)
        new{T}(axis, strength)
    end
end


"""
    _check_axis(axis)

Internal validation helper for Cartesian spin-component symbols.
"""
@inline function _check_axis(axis::Symbol)
    axis in (:x, :y, :z) ||
        throw(ArgumentError(
            "axis must be one of :x, :y, or :z; got :$axis",
        ))

    return nothing
end


"""
    field_components(field::Field)

Return the Cartesian field tuple `(hx, hy, hz)`.

The numeric type of the field strength is preserved.
"""
@inline function field_components(field::Field{T}) where {T}
    z = zero(T)
    h = field.strength

    if field.axis === :x
        return (h, z, z)
    elseif field.axis === :y
        return (z, h, z)
    else
        return (z, z, h)
    end
end
