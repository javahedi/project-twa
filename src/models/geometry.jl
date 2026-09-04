

"""
    AbstractGeometry

Abstract supertype for all lattice/geometry descriptions used by TWA.jl.

Concrete geometries should encode only spatial/topological information
(e.g. site count, boundary conditions, connectivity, and distances).
They should not contain Hamiltonian parameters such as coupling strengths.
"""
abstract type AbstractGeometry end


"""
    AbstractBoundaryCondition

Abstract supertype for boundary conditions associated with a geometry.
"""
abstract type AbstractBoundaryCondition end


"""
    OpenBoundary()

Open boundary conditions.

For a one-dimensional chain with `N` sites, nearest-neighbor bonds are

`(1,2), (2,3), ..., (N-1,N)`.
"""
struct OpenBoundary <: AbstractBoundaryCondition end


"""
    PeriodicBoundary()

Periodic boundary conditions.

For a one-dimensional chain with `N > 2` sites, the nearest-neighbor bonds
include the closing bond `(N,1)` in addition to the open-chain bonds.

Each undirected physical bond is represented exactly once. Therefore,
a two-site periodic chain contains only the bond `(1,2)`, not both
`(1,2)` and `(2,1)`.
"""
struct PeriodicBoundary <: AbstractBoundaryCondition end





"""
    nsites(geometry)

Return the number of physical sites in a geometry.

All concrete geometries used by `TWA.jl` should implement this interface.
"""
nsites(geometry::AbstractGeometry) =
    throw(MethodError(nsites, (geometry,)))







"""
    Chain(nsites)
    Chain(nsites, boundary)

One-dimensional chain geometry.

# Arguments
- `nsites`: Number of physical sites. Must be positive.
- `boundary`: Boundary condition. Defaults to [`OpenBoundary`](@ref).

# Examples
```julia
Chain(50)
Chain(50, PeriodicBoundary())
```

The boundary-condition type is part of the concrete `Chain` type, allowing
Julia to specialize geometry operations such as bond iteration without
runtime boolean checks.
"""
struct Chain{B<:AbstractBoundaryCondition} <: AbstractGeometry
    nsites::Int
    boundary::B

    function Chain(
        nsites::Integer,
        boundary::B = OpenBoundary(),
    ) where {B<:AbstractBoundaryCondition}

        nsites > 0 || throw(ArgumentError("nsites must be positive"))

        new{B}(Int(nsites), boundary)
    end
end





"""
    nsites(chain::Chain)

Return the number of sites in a one-dimensional chain.
"""
nsites(chain::Chain) = chain.nsites




"""
    SquareLattice(Lx, Ly)
    SquareLattice(Lx, Ly, boundary)

Two-dimensional square lattice with dimensions `Lx × Ly`.

Sites are numbered in row-major order,

```text
1        2        ...      Lx
Lx + 1   Lx + 2   ...      2Lx
...

(Ly - 1)Lx + 1  ...      LxLy
```

The total number of physical sites is Lx * Ly.
Arguments
Lx: Number of sites along the x direction. Must be positive.
Ly: Number of sites along the y direction. Must be positive.
boundary: Boundary condition. Defaults to OpenBoundary.


Examples
SquareLattice(4, 3)
SquareLattice(4, 3, PeriodicBoundary())


"""
struct SquareLattice{B<:AbstractBoundaryCondition} <: AbstractGeometry
    Lx::Int
    Ly::Int
    boundary::B
    function SquareLattice(
        Lx::Integer,
        Ly::Integer,
        boundary::B = OpenBoundary(),
    ) where {B<:AbstractBoundaryCondition}

        Lx > 0 || throw(ArgumentError("Lx must be positive"))
        Ly > 0 || throw(ArgumentError("Ly must be positive"))

        new{B}(Int(Lx), Int(Ly), boundary)
    end
end
"""
nsites(lattice::SquareLattice)
Return the total number of physical sites in a square lattice.
"""
nsites(lattice::SquareLattice) = lattice.Lx * lattice.Ly

"""
coordinates(lattice::SquareLattice, site)
Return the integer lattice coordinates (x, y) associated with site.
Sites are numbered in row-major order.
"""
function coordinates(
    lattice::SquareLattice,
    site::Integer,
    )
    _check_site(lattice, site)
    s = Int(site) - 1

    x = mod(s, lattice.Lx) + 1
    y = div(s, lattice.Lx) + 1

    return (x, y)
end



@inline function _check_site(
    lattice::SquareLattice,
    site::Integer,
)
    N = nsites(lattice)

    1 <= site <= N ||
        throw(BoundsError(1:N, site))

    return nothing
end



"""
    ChainBonds(chain)

Lightweight iterator over nearest-neighbor bonds of a [`Chain`](@ref).

This iterator avoids constructing and storing a `Vector` of bond tuples.
That becomes useful when topology is traversed repeatedly during model
construction, CTWA cluster mapping, and other performance-sensitive code.

Users will normally call [`bonds`](@ref) rather than construct
`ChainBonds` directly.
"""
struct ChainBonds{B<:AbstractBoundaryCondition}
    chain::Chain{B}
end


"""
    bonds(chain::Chain)

Return a lightweight iterator over the nearest-neighbor bonds of `chain`.

Each undirected physical bond appears exactly once.

# Examples
```julia
collect(bonds(Chain(4)))
# [(1, 2), (2, 3), (3, 4)]

collect(bonds(Chain(4, PeriodicBoundary())))
# [(1, 2), (2, 3), (3, 4), (4, 1)]
```

Package internals should generally iterate directly,

```julia
for (i, j) in bonds(chain)
    # work with bond (i, j)
end
```

rather than calling `collect`, so that no intermediate bond array is needed.
"""
bonds(chain::Chain) = ChainBonds(chain)


# Tell Julia that ChainBonds has a known length and tuple element type.
# These traits improve iterator behavior and make the intended API explicit.
Base.IteratorSize(::Type{<:ChainBonds}) = Base.HasLength()
Base.eltype(::Type{<:ChainBonds}) = Tuple{Int,Int}


"""
    length(bonds(chain))

Return the number of unique nearest-neighbor bonds in an open chain.
"""
Base.length(b::ChainBonds{OpenBoundary}) =
    max(b.chain.nsites - 1, 0)


"""
    length(bonds(chain))

Return the number of unique nearest-neighbor bonds in a periodic chain.

For `N = 1`, there are no bonds.
For `N = 2`, the single undirected bond `(1,2)` is counted once.
For `N > 2`, there are `N` bonds.
"""
Base.length(b::ChainBonds{PeriodicBoundary}) = begin
    n = b.chain.nsites

    if n <= 1
        0
    elseif n == 2
        1
    else
        n
    end
end


# Open-chain iteration:
# state `i` denotes the left site of the next bond `(i, i+1)`.
function Base.iterate(
    b::ChainBonds{OpenBoundary},
    i::Int = 1,
)
    n = b.chain.nsites

    i >= n && return nothing

    return ((i, i + 1), i + 1)
end


# Periodic-chain iteration:
# first emit the open-chain bonds, then emit the closing `(N,1)` bond
# only when N > 2. This prevents double-counting the N = 2 bond.
function Base.iterate(
    b::ChainBonds{PeriodicBoundary},
    i::Int = 1,
)
    n = b.chain.nsites

    n <= 1 && return nothing

    if i < n
        return ((i, i + 1), i + 1)
    elseif i == n && n > 2
        return ((n, 1), i + 1)
    else
        return nothing
    end
end




"""
    distance(chain::Chain, i::Integer, j::Integer)

Return the shortest lattice distance between sites `i` and `j`.

For an open chain, the distance is d(i,j) = |i-j|.
For a periodic chain, the shorter path around the ring is used,
    d(i,j) = min(|i-j|, N-|i-j|).
Both site indices must lie in 1:chain.nsites.
   
# 4
```julia
    distance(Chain(6), 1, 5)
    # 4

    distance(Chain(6, PeriodicBoundary()), 1, 5)
    # 2
```
"""
function distance(chain::Chain{OpenBoundary},
    i::Integer,
    j::Integer)

    _check_site(chain, i)
    _check_site(chain, j)
    return abs(i - j)
end

function distance(chain::Chain{PeriodicBoundary},
    i::Integer,
    j::Integer)

    _check_site(chain, i)
    _check_site(chain, j)
    d = abs(i - j)

    return min(d, chain.nsites - d)
end

"""
_check_site(chain, i)
Internal helper that verifies that i is a valid site index.
This function is intentionally kept separate from distance so that
future geometry operations can reuse the same bounds checking.
"""
@inline function _check_site(chain::Chain,i::Integer)
    1 <= i <= chain.nsites ||
    throw(BoundsError(1:chain.nsites, i))
    return nothing
end






struct SquareLatticeBonds{B<:AbstractBoundaryCondition}
    lattice::SquareLattice{B}
end

bonds(lattice::SquareLattice) = SquareLatticeBonds(lattice)

Base.IteratorSize(::Type{<:SquareLatticeBonds}) = Base.HasLength()
Base.eltype(::Type{<:SquareLatticeBonds}) = Tuple{Int,Int}


@inline function _site_index(
    lattice::SquareLattice,
    x::Integer,
    y::Integer,
)
    return (Int(y) - 1) * lattice.Lx + Int(x)
end


function bonds(lattice::SquareLattice{OpenBoundary})
    Lx = lattice.Lx
    Ly = lattice.Ly

    result = Tuple{Int,Int}[]

    sizehint!(
        result,
        max(Lx - 1, 0) * Ly +
        Lx * max(Ly - 1, 0),
    )

    for y in 1:Ly
        for x in 1:Lx
            i = _site_index(lattice, x, y)

            if x < Lx
                push!(
                    result,
                    (i, _site_index(lattice, x + 1, y)),
                )
            end

            if y < Ly
                push!(
                    result,
                    (i, _site_index(lattice, x, y + 1)),
                )
            end
        end
    end

    return result
end


function bonds(lattice::SquareLattice{PeriodicBoundary})
    Lx = lattice.Lx
    Ly = lattice.Ly

    result = Tuple{Int,Int}[]

    for y in 1:Ly
        for x in 1:Lx
            i = _site_index(lattice, x, y)

            # +x direction
            if Lx > 1
                xp = x == Lx ? 1 : x + 1
                j = _site_index(lattice, xp, y)

                if i < j || Lx > 2
                    push!(result, (i, j))
                end
            end

            # +y direction
            if Ly > 1
                yp = y == Ly ? 1 : y + 1
                j = _site_index(lattice, x, yp)

                if i < j || Ly > 2
                    push!(result, (i, j))
                end
            end
        end
    end

    return result
end



"""
    distance(lattice::SquareLattice, i::Integer, j::Integer)

Return the Euclidean lattice distance between sites `i` and `j`.

For periodic boundaries, the minimum-image displacement is used
independently along the x and y directions.
"""
function distance(
    lattice::SquareLattice{OpenBoundary},
    i::Integer,
    j::Integer,
)
    _check_site(lattice, i)
    _check_site(lattice, j)

    xi, yi = coordinates(lattice, i)
    xj, yj = coordinates(lattice, j)

    dx = abs(xi - xj)
    dy = abs(yi - yj)

    return hypot(dx, dy)
end


function distance(
    lattice::SquareLattice{PeriodicBoundary},
    i::Integer,
    j::Integer,
)
    _check_site(lattice, i)
    _check_site(lattice, j)

    xi, yi = coordinates(lattice, i)
    xj, yj = coordinates(lattice, j)

    dx = abs(xi - xj)
    dy = abs(yi - yj)

    dx = min(dx, lattice.Lx - dx)
    dy = min(dy, lattice.Ly - dy)

    return hypot(dx, dy)
end