using TWA
using Printf
using Random

const NSITES = 12
const CLUSTER_SIZES = (1, 2, 3, 4, 6)

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

function benchmark_call(f; samples)
    # Warm-up first so compilation is excluded as much as possible.
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

function compile_measurement(model, clustering)
    # Warm-up.
    compile_ctwa_hamiltonian(
        model,
        clustering,
    )

    elapsed = @elapsed begin
        compiled =
            compile_ctwa_hamiltonian(
                model,
                clustering,
            )
    end

    allocated = @allocated begin
        compiled_alloc =
            compile_ctwa_hamiltonian(
                model,
                clustering,
            )
    end

    return (
        time_ns = elapsed * 1e9,
        bytes = allocated,
    )
end

function rhs_measurement(compiled, rng)
    dimension =
        basis_size(
            compiled.basis,
        )

    nclusters =
        cluster_count(
            compiled.clustering,
        )

    x = randn(
        rng,
        Float64,
        dimension,
        nclusters,
    )

    dx = similar(x)

    f = () -> ctwa_rhs!(
        dx,
        x,
        compiled,
        0.0,
    )

    # Keep the number of repetitions moderate for large clusters.
    samples =
        dimension <= 15  ? 20_000 :
        dimension <= 63  ? 5_000  :
        dimension <= 255 ? 1_000  :
                           50

    measurement =
        benchmark_call(
            f;
            samples=samples,
        )

    return (
        time_ns = measurement.time_ns,
        bytes = measurement.bytes,
        samples = samples,
    )
end

model = SpinModel(
    Chain(NSITES),
    XXZ(
        J=1.0,
        Δ=0.5,
    ) +
    Field(
        :z,
        0.3,
    ),
)

rng = Xoshiro(12345)

println()
println("CTWA RHS benchmark")
println("==================")
println("physical sites: ", NSITES)
println("Hamiltonian: nearest-neighbor XXZ(J=1, Δ=0.5) + z field 0.3")
println()

header = @sprintf(
    "%7s %9s %10s %11s %11s",
    "cluster",
    "clusters",
    "basis dim",
    "local terms",
    "inter terms",
)

println(header)
println("-"^length(header))

compiled_models = Dict{Int,Any}()

for k in CLUSTER_SIZES
    clustering =
        Clustering(
            NSITES,
            k,
        )

    compiled =
        compile_ctwa_hamiltonian(
            model,
            clustering,
        )

    compiled_models[k] =
        compiled

    @printf(
        "%7d %9d %10d %11d %11d\n",
        k,
        cluster_count(clustering),
        basis_size(compiled.basis),
        length(local_terms(compiled)),
        length(intercluster_terms(compiled)),
    )
end

println()
println("Compilation cost")
println("================")

for k in CLUSTER_SIZES
    clustering =
        Clustering(
            NSITES,
            k,
        )

    result =
        compile_measurement(
            model,
            clustering,
        )

    @printf(
        "cluster %d : %12.2f μs   %s\n",
        k,
        result.time_ns / 1e3,
        bytes_string(result.bytes),
    )
end

println()
println("Single RHS evaluation")
println("=====================")

for k in CLUSTER_SIZES
    compiled =
        compiled_models[k]

    result =
        rhs_measurement(
            compiled,
            rng,
        )

    @printf(
        "cluster %d : %12.2f μs   %10s/call   samples=%d\n",
        k,
        result.time_ns / 1e3,
        bytes_string(result.bytes),
        result.samples,
    )
end

println()
println("Approximate state size per trajectory")
println("=====================================")

for k in CLUSTER_SIZES
    compiled =
        compiled_models[k]

    dimension =
        basis_size(
            compiled.basis,
        )

    nclusters =
        cluster_count(
            compiled.clustering,
        )

    values =
        dimension *
        nclusters

    bytes =
        values *
        sizeof(Float64)

    @printf(
        "cluster %d : %8d variables   %s\n",
        k,
        values,
        bytes_string(bytes),
    )
end

println()
println("Notes:")
println("  * RHS timing measures the current sparse, matrix-free CTWA dynamics.")
println("  * No dense structure-constant tensor is constructed.")
println("  * Compilation allocations are setup-time only.")
println("  * RHS allocations should ideally be zero.")
println("  * The important scaling variable is basis dimension 4^k - 1.")
