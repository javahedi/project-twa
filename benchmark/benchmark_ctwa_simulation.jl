using TWA
using Random
using Printf

println()
println("CTWA simulation benchmark")
println("=========================")

const NSITES = 12
const CLUSTER_SIZES = (2, 3, 4, 6)

model =
    SpinModel(
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

for k in CLUSTER_SIZES
    method =
        CTWA(
            cluster_size=k,
            trajectories=10,
        )

    # Warm up compiler + solver path.
    simulate(
        model,
        Up(),
        CTWA(
            cluster_size=k,
            trajectories=1,
        );
        tspan=(0.0, 0.05),
        saveat=[
            0.0,
            0.05,
        ],
        rng=Xoshiro(1),
    )

    elapsed =
        @elapsed result =
            simulate(
                model,
                Up(),
                method;
                tspan=(0.0, 0.5),
                saveat=range(
                    0.0,
                    0.5;
                    length=11,
                ),
                rng=Xoshiro(1234),
            )

    nvariables =
        basis_size(
            result.basis,
        ) *
        cluster_count(
            result.clustering,
        )

    @printf(
        "cluster %d : variables/trajectory=%5d  %9.3f ms total  %8.3f ms/trajectory\n",
        k,
        nvariables,
        elapsed * 1e3,
        elapsed / method.trajectories * 1e3,
    )
end
