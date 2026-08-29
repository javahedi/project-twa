"""
    AbstractState

Abstract supertype for physical initial states.

State objects describe the quantum state itself, independently of the
approximation method used later for phase-space sampling. DTWA and CTWA
will therefore consume the same state description but generate different
sampled phase-space representations.
"""
abstract type AbstractState end


"""
    Polarized(axis; sign=1)

Uniform product state polarized along one Cartesian spin axis.

# Arguments

- `axis`: one of `:x`, `:y`, or `:z`.
- `sign`: either `+1` or `-1`.

# Examples

    Polarized(:z)
    Polarized(:x; sign=-1)

For spin-1/2 systems, this represents a product state in which every site
is aligned with the chosen Cartesian axis and sign.
"""
struct Polarized <: AbstractState
    axis::Symbol
    sign::Int

    function Polarized(
        axis::Symbol;
        sign::Integer=1,
    )
        _check_state_axis(axis)
        _check_state_sign(sign)

        new(axis, Int(sign))
    end
end


"""
    Up()

Convenience constructor for a product state polarized along `+z`.
"""
Up() = Polarized(:z; sign=1)


"""
    Down()

Convenience constructor for a product state polarized along `-z`.
"""
Down() = Polarized(:z; sign=-1)


"""
    DomainWall(; axis=:z)

Product state with opposite polarization on the two halves of the system.

The first half is polarized along the positive direction of `axis`; the
second half is polarized along the negative direction.

For an odd number of sites, the central site belongs to the first half.

# Example

    DomainWall()
    DomainWall(axis=:x)

The geometry size is intentionally not stored in the state object. The
site-dependent direction is resolved later using the model geometry.
"""
struct DomainWall <: AbstractState
    axis::Symbol

    function DomainWall(; axis::Symbol=:z)
        _check_state_axis(axis)
        new(axis)
    end
end


"""
    state_direction(state, site, nsites)

Return the Cartesian unit direction `(sx, sy, sz)` associated with a
physical product state at one site.

The values are integers `-1`, `0`, or `+1`. They describe polarization
direction only; DTWA/CTWA sampling conventions are applied later.
"""
@inline function state_direction(
    state::Polarized,
    site::Integer,
    nsites::Integer,
)
    _check_state_site(site, nsites)

    s = state.sign

    if state.axis === :x
        return (s, 0, 0)
    elseif state.axis === :y
        return (0, s, 0)
    else
        return (0, 0, s)
    end
end


@inline function state_direction(
    state::DomainWall,
    site::Integer,
    nsites::Integer,
)
    _check_state_site(site, nsites)

    midpoint = cld(nsites, 2)
    s = site <= midpoint ? 1 : -1

    if state.axis === :x
        return (s, 0, 0)
    elseif state.axis === :y
        return (0, s, 0)
    else
        return (0, 0, s)
    end
end


"""
    _check_state_axis(axis)

Internal validation helper for Cartesian state axes.
"""
@inline function _check_state_axis(axis::Symbol)
    axis in (:x, :y, :z) ||
        throw(ArgumentError(
            "axis must be one of :x, :y, or :z; got :$axis",
        ))

    return nothing
end


"""
    _check_state_sign(sign)

Internal validation helper for polarization sign.
"""
@inline function _check_state_sign(sign::Integer)
    sign in (-1, 1) ||
        throw(ArgumentError(
            "sign must be either -1 or +1; got $sign",
        ))

    return nothing
end


"""
    _check_state_site(site, nsites)

Internal validation helper for site-dependent state queries.
"""
@inline function _check_state_site(
    site::Integer,
    nsites::Integer,
)
    nsites > 0 ||
        throw(ArgumentError("nsites must be positive"))

    1 <= site <= nsites ||
        throw(BoundsError(1:nsites, site))

    return nothing
end
