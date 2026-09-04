

"""
    Clustering(nsites, cluster_size)
    Clustering(nsites, sites)

Equal-size partition of a spin system into CTWA clusters.

A clustering stores an ordered collection of physical sites for each cluster.
The physical site indices do not need to be consecutive.

`Clustering(nsites, cluster_size)` preserves the traditional contiguous
site-index partition.

More general cluster layouts can be specified by supplying a matrix whose
columns are clusters and whose rows are local positions within each cluster.

Every physical site in `1:nsites` must occur exactly once, and all clusters
must have equal size.
"""
struct Clustering
    nsites::Int
    cluster_size::Int
    nclusters::Int

    sites::Matrix{Int}
    site_to_cluster::Vector{Int}
    site_to_position::Vector{Int}

    function Clustering(
        nsites::Integer,
        sites::AbstractMatrix{<:Integer},
    )
        nsites > 0 ||
            throw(ArgumentError(
                "nsites must be positive; got $nsites",
            ))

        n = Int(nsites)
        k = size(sites, 1)
        nc = size(sites, 2)

        k > 0 ||
            throw(ArgumentError(
                "clusters must contain at least one site",
            ))

        nc > 0 ||
            throw(ArgumentError(
                "at least one cluster is required",
            ))

        length(sites) == n ||
            throw(ArgumentError(
                "clustering must contain exactly nsites=$n physical sites; " *
                "got $(length(sites)) entries",
            ))

        stored_sites = Matrix{Int}(sites)
        site_to_cluster = zeros(Int, n)
        site_to_position = zeros(Int, n)
        seen = falses(n)

        for cluster in 1:nc
            for position in 1:k
                site = stored_sites[position, cluster]

                1 <= site <= n ||
                    throw(ArgumentError(
                        "physical site $site lies outside 1:$n",
                    ))

                !seen[site] ||
                    throw(ArgumentError(
                        "physical site $site appears more than once",
                    ))

                seen[site] = true
                site_to_cluster[site] = cluster
                site_to_position[site] = position
            end
        end

        all(seen) ||
            throw(ArgumentError(
                "every physical site in 1:$n must appear exactly once",
            ))

        new(
            n,
            k,
            nc,
            stored_sites,
            site_to_cluster,
            site_to_position,
        )
    end
end


function Clustering(
    nsites::Integer,
    cluster_size::Integer,
)
    nsites > 0 ||
        throw(ArgumentError(
            "nsites must be positive; got $nsites",
        ))

    cluster_size > 0 ||
        throw(ArgumentError(
            "cluster_size must be positive; got $cluster_size",
        ))

    cluster_size <= nsites ||
        throw(ArgumentError(
            "cluster_size must not exceed nsites; got " *
            "cluster_size=$cluster_size, nsites=$nsites",
        ))

    nsites % cluster_size == 0 ||
        throw(ArgumentError(
            "nsites must be divisible by cluster_size for equal-size " *
            "clustering; got nsites=$nsites, cluster_size=$cluster_size",
        ))

    n = Int(nsites)
    k = Int(cluster_size)

    sites = reshape(
        collect(1:n),
        k,
        n ÷ k,
    )

    return Clustering(n, sites)
end


Clustering(
    geometry::AbstractGeometry,
    cluster_size::Integer,
) = Clustering(
    nsites(geometry),
    cluster_size,
)




"""
    Clustering(lattice::SquareLattice, block_shape)

Construct a spatially compact rectangular CTWA clustering of a square lattice.

`block_shape = (bx, by)` specifies the number of lattice sites along the
x and y directions within each cluster.

Both lattice dimensions must be divisible by the corresponding block
dimensions.

For example,

    Clustering(SquareLattice(4, 4), (2, 2))

produces the clusters

    (1, 2, 5, 6)
    (3, 4, 7, 8)
    (9, 10, 13, 14)
    (11, 12, 15, 16)

for the row-major site ordering of `SquareLattice`.
"""
function Clustering(
    lattice::SquareLattice,
    block_shape::Tuple{<:Integer,<:Integer},
)
    bx = Int(block_shape[1])
    by = Int(block_shape[2])

    bx > 0 ||
        throw(ArgumentError(
            "block size along x must be positive; got $bx",
        ))

    by > 0 ||
        throw(ArgumentError(
            "block size along y must be positive; got $by",
        ))

    lattice.Lx % bx == 0 ||
        throw(ArgumentError(
            "lattice size Lx=$(lattice.Lx) must be divisible by bx=$bx",
        ))

    lattice.Ly % by == 0 ||
        throw(ArgumentError(
            "lattice size Ly=$(lattice.Ly) must be divisible by by=$by",
        ))

    cluster_size = bx * by

    nclusters_x = lattice.Lx ÷ bx
    nclusters_y = lattice.Ly ÷ by
    nclusters = nclusters_x * nclusters_y

    sites =
        Matrix{Int}(
            undef,
            cluster_size,
            nclusters,
        )

    cluster = 0

    for block_y in 0:(nclusters_y - 1)
        for block_x in 0:(nclusters_x - 1)
            cluster += 1
            position = 0

            for local_y in 1:by
                for local_x in 1:bx
                    position += 1

                    x =
                        block_x * bx +
                        local_x

                    y =
                        block_y * by +
                        local_y

                    site =
                        (y - 1) * lattice.Lx +
                        x

                    sites[
                        position,
                        cluster,
                    ] = site
                end
            end
        end
    end

    return Clustering(
        nsites(lattice),
        sites,
    )
end


"""
    cluster_count(clustering)

Number of clusters.
"""
cluster_count(
    clustering::Clustering,
) = clustering.nclusters


Base.length(
    clustering::Clustering,
) = clustering.nclusters



function cluster_sites(
    clustering::Clustering,
    cluster::Integer,
)
    _check_cluster_index(
        clustering,
        cluster,
    )

    return @view clustering.sites[:, Int(cluster)]
end


function site_cluster(
    clustering::Clustering,
    site::Integer,
)
    _check_clustering_site(
        clustering,
        site,
    )

    return clustering.site_to_cluster[Int(site)]
end


function site_position(
    clustering::Clustering,
    site::Integer,
)
    _check_clustering_site(
        clustering,
        site,
    )

    return clustering.site_to_position[Int(site)]
end


function site_cluster_position(
    clustering::Clustering,
    site::Integer,
)
    _check_clustering_site(
        clustering,
        site,
    )

    s = Int(site)

    return (
        clustering.site_to_cluster[s],
        clustering.site_to_position[s],
    )
end


"""
    same_cluster(clustering, i, j)

Return `true` when physical sites `i` and `j` belong to the same cluster.
"""
function same_cluster(
    clustering::Clustering,
    i::Integer,
    j::Integer,
)
    _check_clustering_site(
        clustering,
        i,
    )

    _check_clustering_site(
        clustering,
        j,
    )

    return site_cluster(
        clustering,
        i,
    ) == site_cluster(
        clustering,
        j,
    )
end


"""
    local_pauli_digits(clustering, site, pauli_digit)

Return the cluster-local Pauli string corresponding to a one-site physical
Pauli operator.

For example, with two-site clusters,

    site 1, X -> (X, I) -> (1, 0)
    site 2, X -> (I, X) -> (0, 1)

and the pattern repeats independently in each cluster.

This is the bridge from physical-site operators to the canonical
`PauliStringBasis` ordering.
"""
function local_pauli_digits(
    clustering::Clustering,
    site::Integer,
    pauli_digit::Integer,
)
    _check_clustering_site(
        clustering,
        site,
    )

    _check_pauli_digit(
        pauli_digit,
    )

    pauli_digit != 0 ||
        throw(ArgumentError(
            "local physical operator must be X, Y, or Z; got identity",
        ))

    position = site_position(
        clustering,
        site,
    )

    n = clustering.cluster_size

    return ntuple(n) do local_site
        local_site == position ?
            Int(pauli_digit) :
            0
    end
end


"""
    local_pauli_index(basis, clustering, site, pauli_digit)

Map a one-site physical Pauli operator directly to its canonical cluster-basis
index.

`basis.cluster_size` must agree with `clustering.cluster_size`.
"""
function local_pauli_index(
    basis::PauliStringBasis,
    clustering::Clustering,
    site::Integer,
    pauli_digit::Integer,
)
    _check_basis_clustering_compatibility(
        basis,
        clustering,
    )

    digits = local_pauli_digits(
        clustering,
        site,
        pauli_digit,
    )

    return pauli_index(
        basis,
        digits,
    )
end


@inline function _check_clustering_site(
    clustering::Clustering,
    site::Integer,
)
    1 <= site <= clustering.nsites ||
        throw(BoundsError(
            1:clustering.nsites,
            site,
        ))

    return nothing
end


@inline function _check_cluster_index(
    clustering::Clustering,
    cluster::Integer,
)
    1 <= cluster <= clustering.nclusters ||
        throw(BoundsError(
            1:clustering.nclusters,
            cluster,
        ))

    return nothing
end


@inline function _check_basis_clustering_compatibility(
    basis::PauliStringBasis,
    clustering::Clustering,
)
    basis.cluster_size ==
        clustering.cluster_size ||
        throw(DimensionMismatch(
            "basis cluster size $(basis.cluster_size) does not match " *
            "clustering cluster size $(clustering.cluster_size)",
        ))

    return nothing
end
