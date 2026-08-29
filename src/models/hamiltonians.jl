"""
    AbstractHamiltonianTerm

Common interface for terms that can appear in a spin Hamiltonian.

Concrete terms currently include pair interactions such as `XXZ` and
one-body terms such as `Field`.
"""
abstract type AbstractHamiltonianTerm end


"""
    Hamiltonian(terms...)

A lightweight immutable collection of Hamiltonian terms.

Users normally construct this implicitly with `+`, for example

    XXZ(J=1.0, Δ=0.5) + Field(:z, 0.3)

The container stores the physics description only. Backend-specific
representations are constructed later by the simulation compiler.
"""
struct Hamiltonian{T<:Tuple}
    terms::T
end

Hamiltonian(terms::AbstractHamiltonianTerm...) = Hamiltonian(terms)


# Pair interactions and fields participate in Hamiltonian composition.
Base.:+(
    a::AbstractHamiltonianTerm,
    b::AbstractHamiltonianTerm,
) = Hamiltonian((a, b))

Base.:+(
    h::Hamiltonian,
    term::AbstractHamiltonianTerm,
) = Hamiltonian((h.terms..., term))

Base.:+(
    term::AbstractHamiltonianTerm,
    h::Hamiltonian,
) = Hamiltonian((term, h.terms...))

Base.:+(
    a::Hamiltonian,
    b::Hamiltonian,
) = Hamiltonian((a.terms..., b.terms...))


"""
    Base.length(h::Hamiltonian)

Return the number of terms in the Hamiltonian.
"""
Base.length(h::Hamiltonian) = length(h.terms)


"""
    Base.iterate(h::Hamiltonian[, state])

Iterate over Hamiltonian terms without allocating an intermediate array.
"""
Base.iterate(h::Hamiltonian, state...) = iterate(h.terms, state...)


Base.IteratorSize(::Type{<:Hamiltonian}) = Base.HasLength()
Base.eltype(::Type{Hamiltonian{T}}) where {T} = eltype(T)
