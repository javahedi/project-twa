"""
    coupling_components(model::SpinModel, i::Integer, j::Integer)

Return the total Cartesian pair-coupling components `(Jx, Jy, Jz)`
between sites `i` and `j`.

All pair terms are summed. One-body fields contribute zero.

The numeric type is inherited from the Hamiltonian terms rather than
being forced to `Float64`.
"""
function coupling_components(
    model::SpinModel,
    i::Integer,
    j::Integer,
)
    _check_site(model.geometry, i)
    _check_site(model.geometry, j)

    return _sum_pair_components(
        model.geometry,
        terms(model),
        i,
        j,
    )
end


"""
    coupling_strength(model::SpinModel, i::Integer, j::Integer)

Return the sum of spatial pair-profile weights between two sites.

This is mainly a model-inspection helper. Dynamics backends will later
use compiled pair representations instead.

The numeric type is inherited from the pair-coupling terms.
"""
function coupling_strength(
    model::SpinModel,
    i::Integer,
    j::Integer,
)
    _check_site(model.geometry, i)
    _check_site(model.geometry, j)

    return _sum_pair_strengths(
        model.geometry,
        terms(model),
        i,
        j,
    )
end


# ---------------------------------------------------------------------------
# Contribution of an individual Hamiltonian term
# ---------------------------------------------------------------------------

# Fields are one-body terms and therefore make no pair contribution.
# Integer zero is intentional: Julia promotes it to the numeric type of
# the actual pair term when the tuple reduction below performs addition.
@inline _pair_strength(
    ::AbstractGeometry,
    ::AbstractField,
    ::Integer,
    ::Integer,
) = 0

@inline _pair_strength(
    geometry::AbstractGeometry,
    profile::AbstractCouplingProfile,
    i::Integer,
    j::Integer,
) = _coupling_strength(geometry, profile, i, j)


@inline _pair_components(
    ::AbstractGeometry,
    ::AbstractField,
    ::Integer,
    ::Integer,
) = (0, 0, 0)

@inline _pair_components(
    geometry::AbstractGeometry,
    profile::AbstractCouplingProfile,
    i::Integer,
    j::Integer,
) = _coupling_components(geometry, profile, i, j)


# ---------------------------------------------------------------------------
# Tuple reductions
#
# Recursive tuple dispatch avoids hard-coding Float64 accumulator values.
# It also allows Julia to specialize the reduction on the concrete tuple
# of Hamiltonian terms stored by SpinModel.
# ---------------------------------------------------------------------------

@inline _sum_pair_strengths(
    ::AbstractGeometry,
    ::Tuple{},
    ::Integer,
    ::Integer,
) = 0

@inline function _sum_pair_strengths(
    geometry::AbstractGeometry,
    model_terms::Tuple,
    i::Integer,
    j::Integer,
)
    first_term = first(model_terms)
    remaining_terms = Base.tail(model_terms)

    return _pair_strength(
        geometry,
        first_term,
        i,
        j,
    ) + _sum_pair_strengths(
        geometry,
        remaining_terms,
        i,
        j,
    )
end


@inline _sum_pair_components(
    ::AbstractGeometry,
    ::Tuple{},
    ::Integer,
    ::Integer,
) = (0, 0, 0)

@inline function _sum_pair_components(
    geometry::AbstractGeometry,
    model_terms::Tuple,
    i::Integer,
    j::Integer,
)
    x1, y1, z1 = _pair_components(
        geometry,
        first(model_terms),
        i,
        j,
    )

    x2, y2, z2 = _sum_pair_components(
        geometry,
        Base.tail(model_terms),
        i,
        j,
    )

    return (
        x1 + x2,
        y1 + y2,
        z1 + z2,
    )
end


# ---------------------------------------------------------------------------
# Nearest-neighbor XXZ
# ---------------------------------------------------------------------------

@inline function _coupling_strength(
    geometry::Chain,
    profile::NearestNeighbor,
    i::Integer,
    j::Integer,
)
    T = typeof(profile.interaction.J)

    i == j && return zero(T)

    return distance(geometry, i, j) == 1 ? one(T) : zero(T)
end


@inline function _coupling_components(
    geometry::Chain,
    profile::NearestNeighbor{<:XXZ},
    i::Integer,
    j::Integer,
)
    interaction = profile.interaction
    T = typeof(interaction.J)

    if i == j
        z = zero(T)
        return (z, z, z)
    end

    w = _coupling_strength(
        geometry,
        profile,
        i,
        j,
    )

    Jxy = interaction.J * w
    Jz = interaction.J * interaction.Δ * w

    return (Jxy, Jxy, Jz)
end


# ---------------------------------------------------------------------------
# Power-law XXZ
# ---------------------------------------------------------------------------

@inline function _coupling_strength(
    geometry::Chain,
    profile::PowerLaw,
    i::Integer,
    j::Integer,
)
    T = promote_type(
        typeof(profile.interaction.J),
        typeof(profile.α),
    )

    i == j && return zero(T)

    r = distance(geometry, i, j)

    return one(T) / (T(r)^profile.α)
end


@inline function _coupling_components(
    geometry::Chain,
    profile::PowerLaw{<:XXZ},
    i::Integer,
    j::Integer,
)
    interaction = profile.interaction

    if i == j
        T = promote_type(
            typeof(interaction.J),
            typeof(profile.α),
        )
        z = zero(T)
        return (z, z, z)
    end

    w = _coupling_strength(
        geometry,
        profile,
        i,
        j,
    )

    Jxy = interaction.J * w
    Jz = interaction.J * interaction.Δ * w

    return (Jxy, Jxy, Jz)
end
