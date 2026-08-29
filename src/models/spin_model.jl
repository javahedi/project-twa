"""
    SpinModel(geometry, hamiltonian)

Physical spin model defined by a geometry and a normalized Hamiltonian.

Normalization rules:

- `XXZ(...)` becomes `NearestNeighbor(XXZ(...))`;
- explicit coupling profiles such as `NearestNeighbor(...)` and
  `PowerLaw(...)` are preserved;
- one-body terms such as `Field(...)` are preserved.

The model contains only the physical description. Numerical
representations for DTWA and CTWA are compiled later.
"""
struct SpinModel{
    G<:AbstractGeometry,
    H<:Hamiltonian,
}
    geometry::G
    hamiltonian::H

    # Internal raw constructor.
    #
    # The extra Val argument deliberately gives this constructor a
    # different signature from the public constructors below. This avoids
    # overwriting Julia's ordinary two-argument constructor and prevents
    # recursive constructor calls during Hamiltonian normalization.
    function SpinModel(
        geometry::G,
        hamiltonian::H,
        ::Val{:raw},
    ) where {
        G<:AbstractGeometry,
        H<:Hamiltonian,
    }
        new{G,H}(geometry, hamiltonian)
    end
end


"""
    SpinModel(geometry, term::AbstractHamiltonianTerm)

Construct a spin model from a single Hamiltonian term.

Bare two-body interactions are interpreted as nearest-neighbor terms.
"""
function SpinModel(
    geometry::AbstractGeometry,
    term::AbstractHamiltonianTerm,
)
    normalized = _normalize_term(term)
    hamiltonian = Hamiltonian((normalized,))

    return SpinModel(
        geometry,
        hamiltonian,
        Val(:raw),
    )
end


"""
    SpinModel(geometry, hamiltonian::Hamiltonian)

Construct a spin model from a composed Hamiltonian.

Each term is normalized independently before storage.
"""
function SpinModel(
    geometry::AbstractGeometry,
    hamiltonian::Hamiltonian,
)
    normalized_terms = map(
        _normalize_term,
        hamiltonian.terms,
    )

    normalized_hamiltonian = Hamiltonian(
        normalized_terms,
    )

    return SpinModel(
        geometry,
        normalized_hamiltonian,
        Val(:raw),
    )
end


# ---------------------------------------------------------------------------
# Hamiltonian-term normalization
# ---------------------------------------------------------------------------

"""
    _normalize_term(interaction::AbstractInteraction)

A bare pair interaction defaults to nearest-neighbor coupling.
"""
@inline _normalize_term(
    interaction::AbstractInteraction,
) = NearestNeighbor(interaction)


"""
    _normalize_term(profile::AbstractCouplingProfile)

Explicit spatial coupling profiles are already complete pair terms and
are therefore preserved.
"""
@inline _normalize_term(
    profile::AbstractCouplingProfile,
) = profile


"""
    _normalize_term(field::AbstractField)

One-body fields are preserved unchanged.
"""
@inline _normalize_term(
    field::AbstractField,
) = field


# ---------------------------------------------------------------------------
# Public accessors
# ---------------------------------------------------------------------------

"""
    geometry(model::SpinModel)

Return the geometry associated with the model.
"""
geometry(model::SpinModel) = model.geometry


"""
    hamiltonian(model::SpinModel)

Return the normalized Hamiltonian associated with the model.
"""
hamiltonian(model::SpinModel) = model.hamiltonian


"""
    terms(model::SpinModel)

Return the normalized Hamiltonian terms as a tuple.
"""
terms(model::SpinModel) = model.hamiltonian.terms
