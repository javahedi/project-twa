"""
    AbstractInteraction

Abstract supertype for two-body spin interactions.

Bare pair interactions are Hamiltonian terms. When passed to `SpinModel`,
they are normalized to nearest-neighbor coupling profiles by default.
"""
abstract type AbstractInteraction <: AbstractHamiltonianTerm end


"""
    AbstractCouplingProfile

Abstract supertype for spatially resolved pair-interaction terms.

A coupling profile is itself a Hamiltonian term because it represents a
complete two-body contribution: both operator structure and spatial
dependence.
"""
abstract type AbstractCouplingProfile <: AbstractHamiltonianTerm end


"""
    XXZ(; J=1.0, Δ=1.0)

XXZ pair interaction with transverse coupling `J` and anisotropy `Δ`.
"""
struct XXZ{T<:Real} <: AbstractInteraction
    J::T
    Δ::T
end


function XXZ(; J::Real=1.0, Δ::Real=1.0)
    Jp, Δp = promote(J, Δ)
    return XXZ(Jp, Δp)
end


"""
    NearestNeighbor(interaction)

Nearest-neighbor spatial coupling profile for a pair interaction.
"""
struct NearestNeighbor{I<:AbstractInteraction} <: AbstractCouplingProfile
    interaction::I
end


"""
    PowerLaw(interaction; α)

Long-range coupling profile with spatial weight `1 / r^α`.

No Kac normalization is applied implicitly.
"""
struct PowerLaw{
    I<:AbstractInteraction,
    T<:Real,
} <: AbstractCouplingProfile
    interaction::I
    α::T

    function PowerLaw(
        interaction::I;
        α::T,
    ) where {
        I<:AbstractInteraction,
        T<:Real,
    }
        α >= zero(T) ||
            throw(ArgumentError("α must be non-negative"))

        new{I,T}(interaction, α)
    end
end





"""
    Ising(axis; J=1.0)

Ising pair interaction along one Cartesian spin axis.

The interaction contributes

    H = J * σᵢᵃ σⱼᵃ

for `axis = :x`, `:y`, or `:z`, with the interacting pairs determined
separately by the coupling profile.

# Examples

    Ising(:x; J=1.0)
    PowerLaw(Ising(:x; J=1.0); α=3.0)
"""
struct Ising{T<:Real} <: AbstractInteraction
    axis::Symbol
    J::T

    function Ising(
        axis::Symbol;
        J::T=1.0,
    ) where {T<:Real}
        axis in (:x, :y, :z) ||
            throw(ArgumentError(
                "axis must be one of :x, :y, or :z; got :$axis",
            ))

        new{T}(axis, J)
    end
end