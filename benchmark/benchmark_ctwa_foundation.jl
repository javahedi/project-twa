using TWA
using Printf

const NSITES = 12
const CLUSTER_SIZES = (2, 3, 4, 6)

function bytes_string(bytes::Integer)
    if bytes < 1024
        return "$(bytes) B"
    elseif bytes < 1024^2
        return @sprintf("%.2f KiB", bytes / 1024)
    elseif bytes < 1024^3
        return @sprintf("%.2f MiB", bytes / 1024^2)
    else
        return @sprintf("%.2f GiB", bytes / 1024^3)
    end
end

function benchmark_call(f; samples=10_000)
    # Warm-up to remove compilation from the measurement.
    f()

    elapsed = @elapsed begin
        for _ in 1:samples
            f()
        end
    end

    allocated = @allocated begin
        for _ in 1:samples
            f()
        end
    end

    return (
        time_ns = elapsed / samples * 1e9,
        bytes = allocated ÷ samples,
    )
end

function benchmark_basis_construction(k)
    # Construction is setup work, so fewer repeats are enough.
    return benchmark_call(
        () -> PauliStringBasis(k);
        samples=10_000,
    )
end

function benchmark_clustering_construction(k)
    return benchmark_call(
        () -> Clustering(NSITES, k);
        samples=10_000,
    )
end

function benchmark_local_mapping(k)
    basis = PauliStringBasis(k)
    clustering = Clustering(NSITES, k)

    site = min(k, NSITES)

    return benchmark_call(
        () -> local_pauli_index(
            basis,
            clustering,
            site,
            1,
        );
        samples=100_000,
    )
end

function benchmark_pauli_product(k)
    basis = PauliStringBasis(k)

    # Pick two deterministic, nontrivial basis elements spread through
    # the basis. The exact pair is not physically important here; we are
    # measuring the current arithmetic/decoding implementation.
    left = max(1, basis_size(basis) ÷ 3)
    right = max(1, 2 * basis_size(basis) ÷ 3)

    return benchmark_call(
        () -> pauli_product(
            basis,
            left,
            right,
        );
        samples=100_000,
    )
end

function benchmark_pauli_commutator(k)
    basis = PauliStringBasis(k)

    left = max(1, basis_size(basis) ÷ 3)
    right = max(1, 2 * basis_size(basis) ÷ 3)

    return benchmark_call(
        () -> pauli_commutator(
            basis,
            left,
            right,
        );
        samples=100_000,
    )
end

function algebra_pair_count(k)
    d = 4^k - 1
    return d^2
end

println()
println("CTWA foundation benchmark")
println("=========================")
println("physical sites: ", NSITES)
println()

header = @sprintf(
    "%7s %10s %10s %16s",
    "cluster",
    "clusters",
    "basis dim",
    "basis pairs",
)
println(header)
println("-"^length(header))

for k in CLUSTER_SIZES
    nclusters = NSITES ÷ k
    dim = 4^k - 1
    pairs = algebra_pair_count(k)

    @printf(
        "%7d %10d %10d %16d\n",
        k,
        nclusters,
        dim,
        pairs,
    )
end

println()
println("Per-call timings after warm-up")
println("==============================")

for k in CLUSTER_SIZES
    println()
    println("cluster size = $k")
    println("----------------")

    basis = benchmark_basis_construction(k)
    clustering = benchmark_clustering_construction(k)
    mapping = benchmark_local_mapping(k)
    product = benchmark_pauli_product(k)
    commutator = benchmark_pauli_commutator(k)

    @printf(
        "PauliStringBasis construction : %10.2f ns   %s/call\n",
        basis.time_ns,
        bytes_string(basis.bytes),
    )

    @printf(
        "Clustering construction       : %10.2f ns   %s/call\n",
        clustering.time_ns,
        bytes_string(clustering.bytes),
    )

    @printf(
        "local_pauli_index             : %10.2f ns   %s/call\n",
        mapping.time_ns,
        bytes_string(mapping.bytes),
    )

    @printf(
        "pauli_product                 : %10.2f ns   %s/call\n",
        product.time_ns,
        bytes_string(product.bytes),
    )

    @printf(
        "pauli_commutator              : %10.2f ns   %s/call\n",
        commutator.time_ns,
        bytes_string(commutator.bytes),
    )
end

println()
println("Important:")
println("  * These are microbenchmarks of the foundation layer only.")
println("  * Full algebra materialization would scale as (4^k - 1)^2 pairs.")
println("  * cluster_size = 6 already means 4095 generators and 16,769,025 pairs.")
println("  * We should benchmark Hamiltonian compilation and RHS separately later.")


#julia --project=. benchmark/benchmark_ctwa_foundation.jl