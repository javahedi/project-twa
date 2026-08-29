"""
    Clustering(nsites, cluster_size)

Contiguous, equal-size partition of a one-dimensional spin system into CTWA
clusters.

The first implementation is intentionally strict:

- clusters are contiguous;
- every cluster has exactly `cluster_size` sites;
- `nsites` must be divisible by `cluster_size`.

For example,

    Clustering(6, 2)

defines

    C₁ = (1, 2)
    C₂ = (3, 4)
    C₃ = (5, 6).

This object is purely geometric/bookkeeping information. It does not contain a
Pauli basis or Hamiltonian coefficients.
"""
struct Clustering
    nsites::Int
    cluster_size::Int
    nclusters::Int

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

        new(
            n,
            k,
            n ÷ k,
        )
    end
end


"""
    Clustering(chain, cluster_size)

Construct a contiguous equal-size clustering from a `Chain`.
"""
Clustering(
    chain::Chain,
    cluster_size::Integer,
) = Clustering(
    chain.nsites,
    cluster_size,
)


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


"""
    cluster_sites(clustering, cluster)

Return the physical site range belonging to one cluster.

For example, with `Clustering(6, 2)`:

    cluster_sites(clustering, 1) == 1:2
    cluster_sites(clustering, 2) == 3:4
    cluster_sites(clustering, 3) == 5:6

The returned `UnitRange` is allocation-free and can be iterated directly.
"""
function cluster_sites(
    clustering::Clustering,
    cluster::Integer,
)
    _check_cluster_index(
        clustering,
        cluster,
    )

    first_site =
        (Int(cluster) - 1) *
        clustering.cluster_size +
        1

    last_site =
        first_site +
        clustering.cluster_size -
        1

    return first_site:last_site
end


"""
    site_cluster(clustering, site)

Return the cluster index containing physical `site`.
"""
function site_cluster(
    clustering::Clustering,
    site::Integer,
)
    _check_clustering_site(
        clustering,
        site,
    )

    return fld(
        Int(site) - 1,
        clustering.cluster_size,
    ) + 1
end


"""
    site_position(clustering, site)

Return the local position of physical `site` inside its cluster.

The position is in `1:cluster_size`.

For `Clustering(6, 2)`:

    site_position(clustering, 1) == 1
    site_position(clustering, 2) == 2
    site_position(clustering, 3) == 1
    site_position(clustering, 4) == 2
"""
function site_position(
    clustering::Clustering,
    site::Integer,
)
    _check_clustering_site(
        clustering,
        site,
    )

    return mod(
        Int(site) - 1,
        clustering.cluster_size,
    ) + 1
end


"""
    site_cluster_position(clustering, site)

Return `(cluster, position)` for a physical site.
"""
function site_cluster_position(
    clustering::Clustering,
    site::Integer,
)
    _check_clustering_site(
        clustering,
        site,
    )

    zero_based_site = Int(site) - 1

    cluster =
        fld(
            zero_based_site,
            clustering.cluster_size,
        ) + 1

    position =
        mod(
            zero_based_site,
            clustering.cluster_size,
        ) + 1

    return (
        cluster,
        position,
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
