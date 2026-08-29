using TWA
using Printf
using Random

const NSITES = 16
const CLUSTER_SIZES = (1, 2,  4, 8)

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
println("CTWA cached RHS benchmark")
println("=========================")
println("physical sites: ", NSITES)
println()

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

    algebra =
        compile_ctwa_algebra(
            compiled,
        )

    dimension =
        basis_size(
            compiled.basis,
        )

    nclusters =
        cluster_count(
            clustering,
        )

    x = randn(
        rng,
        Float64,
        dimension,
        nclusters,
    )

    dx_reference = similar(x)
    dx_cached = similar(x)

    samples =
        dimension <= 15  ? 20_000 :
        dimension <= 63  ? 5_000  :
        dimension <= 255 ? 1_000  :
                           100

    reference =
        benchmark_call(
            () -> ctwa_rhs!(
                dx_reference,
                x,
                compiled,
                0.0,
            );
            samples=samples,
        )

    cached =
        benchmark_call(
            () -> ctwa_rhs!(
                dx_cached,
                x,
                compiled,
                algebra,
                0.0,
            );
            samples=samples,
        )

    used_actions =
        count(
            action -> action !== nothing,
            algebra.actions,
        )

    cached_entries =
        sum(
            action === nothing ?
            0 :
            length(action.entries)
            for action in algebra.actions
        )

    cache_bytes =
        Base.summarysize(
            algebra,
        )

    speedup =
        reference.time_ns /
        cached.time_ns

    println("cluster size = $k")
    @printf(
        "  basis dimension       : %d\n",
        dimension,
    )
    @printf(
        "  cached generators     : %d\n",
        used_actions,
    )
    @printf(
        "  cached action entries : %d\n",
        cached_entries,
    )
    @printf(
        "  cache size            : %s\n",
        bytes_string(cache_bytes),
    )
    @printf(
        "  reference RHS         : %10.2f μs   %s/call\n",
        reference.time_ns / 1e3,
        bytes_string(reference.bytes),
    )
    @printf(
        "  cached RHS            : %10.2f μs   %s/call\n",
        cached.time_ns / 1e3,
        bytes_string(cached.bytes),
    )
    @printf(
        "  speedup               : %10.2fx\n",
        speedup,
    )
    println()
end
