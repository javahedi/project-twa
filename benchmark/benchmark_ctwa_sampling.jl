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

println()
println("CTWA sampling benchmark")
println("=======================")

for k in CLUSTER_SIZES
    clustering =
        Clustering(
            NSITES,
            k,
        )

    plan =
        compile_ctwa_sampling(
            Up(),
            clustering,
        )

    rng =
        Xoshiro(1234)

    x =
        similar(
            plan.means,
        )

    # Warmup.
    sample_ctwa!(
        x,
        plan;
        rng=rng,
    )

    samples =
        k <= 2 ? 10_000 :
        k <= 4 ? 2_000 :
                 500

    elapsed =
        @elapsed begin
            for _ in 1:samples
                sample_ctwa!(
                    x,
                    plan;
                    rng=rng,
                )
            end
        end

    allocated =
        @allocated begin
            for _ in 1:samples
                sample_ctwa!(
                    x,
                    plan;
                    rng=rng,
                )
            end
        end

    dimension =
        basis_size(
            plan.basis,
        )

    nclusters =
        cluster_count(
            clustering,
        )

    nmasks =
        2^k - 1

    @printf(
        "cluster %d : dim=%5d clusters=%2d masks/cluster=%3d  %9.2f μs  %s/call\n",
        k,
        dimension,
        nclusters,
        nmasks,
        elapsed / samples * 1e6,
        bytes_string(allocated ÷ samples),
    )
end
