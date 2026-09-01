using TWA
using CairoMakie
using Random

# Long-range transverse Ising benchmark:
#
# H = sum_{i<j} J / r_ij^α σᵢˣ σⱼˣ
#
# Initial state: |↑z ↑z ... ↑z>
#
# Compare:
#   Exact
#   traditional TWA = CTWA(k=1, GaussianSampling())
#   DTWA
#   dcTWA with k=2
#   dcTWA with k=4

const L = 100
const J = 1.0
const ALPHAS = (0.0, 3.0)

const TMAX = 4.0
const DT = 0.02

# Use a finer time grid only for plotting the analytical exact solution.
# This does NOT change the simulation save times or simulation cost.
const EXACT_DT = 0.02

const NTRAJ = 1000

const OUTPUT_TIMES = collect(0.0:DT:TMAX)
const EXACT_TIMES  = collect(0.0:EXACT_DT:TMAX)

# Note the trailing comma: (2,) is a one-element tuple.
const DC_CLUSTER_SIZES = (2,) # use (2, 4) when both are wanted

function exact_magnetization_z(model::SpinModel, ts)
    nsites = geometry(model).nsites
    values = zeros(Float64, length(ts))

    for (time_index, t) in pairs(ts)
        total = 0.0

        for i in 1:nsites
            local_magnetization = 1.0

            for j in 1:nsites
                i == j && continue

                Jx, _, _ = coupling_components(model, i, j)
                local_magnetization *= cos(2 * Jx * t)
            end

            total += local_magnetization
        end

        values[time_index] = total / nsites
    end

    return values
end

function run_benchmark(α::Real)
    println()
    println("Running α = $α")

    model = SpinModel(
        Chain(L),
        PowerLaw(
            Ising(:x; J=J);
            α=α,
        ),
    )

    state = Up()

    # Exact solution is evaluated on its own finer plotting grid.
    exact = exact_magnetization_z(model, EXACT_TIMES)

    println("  traditional TWA (Gaussian, k=1) ...")
    twa_result = simulate(
        model,
        state,
        CTWA(
            cluster_size=1,
            trajectories=NTRAJ,
            sampling=GaussianSampling(),
        );
        tspan=(0.0, TMAX),
        saveat=OUTPUT_TIMES,
        rng=Xoshiro(100 + round(Int, 10α)),
    )

    twa = collect(
        expectation(
            twa_result,
            Magnetization(:z),
        ),
    )

    println("  DTWA ...")
    dtwa_result = simulate(
        model,
        state,
        DTWA(
            trajectories=NTRAJ,
        );
        tspan=(0.0, TMAX),
        saveat=OUTPUT_TIMES,
        rng=Xoshiro(200 + round(Int, 10α)),
    )

    dtwa = collect(
        expectation(
            dtwa_result,
            Magnetization(:z),
        ),
    )

    dctwa = Dict{Int,Vector{Float64}}()
    gctwa = Dict{Int,Vector{Float64}}()

    

    for cluster_size in DC_CLUSTER_SIZES
        println("  dcTWA k=$cluster_size ...")

        result = simulate(
            model,
            state,
            CTWA(
                cluster_size=cluster_size,
                trajectories=NTRAJ,
                sampling=DiscreteSampling(),
            );
            tspan=(0.0, TMAX),
            saveat=OUTPUT_TIMES,
            rng=Xoshiro(
                1000 +
                100 * cluster_size +
                round(Int, 10α),
            ),
        )

        dctwa[cluster_size] = collect(
            expectation(
                result,
                Magnetization(:z),
            ),
        )



        println("  gcTWA k=2 ...")

        gctwa_result = simulate(
            model,
            state,
            CTWA(
                cluster_size=2,
                trajectories=NTRAJ,
                sampling=GaussianSampling(),
            );
            tspan=(0.0, TMAX),
            saveat=OUTPUT_TIMES,
            rng=Xoshiro(
                500 +
                round(Int, 10α),
            ),
        )

        gctwa[cluster_size]  = collect(
            expectation(
                gctwa_result,
                Magnetization(:z),
            ),
        )
    end

    return (
        exact=exact,
        twa=twa,
        dtwa=dtwa,
        gctwa=gctwa,
        dctwa=dctwa,
    )
end

results = Dict(
    α => run_benchmark(α)
    for α in ALPHAS
)


fig = Figure(size=(700, 300))

axes = Axis[]

for (column, α) in enumerate(ALPHAS)
    ax = Axis(
        fig[1, column],
        xlabel="tJ",
        ylabel=column == 1 ? "⟨Sᶻ⟩ / L" : "",
        title="α = $α",
        limits=(
            0.0,
            4.0,
            column == 1 ? -1.0 : 0.0,
            1.0,
        ),
    )

    push!(axes, ax)

    result = results[α]

    # Exact result: finer plotting grid than the numerical simulations.
    scatter!(
        ax,
        EXACT_TIMES,
        result.exact;
        label="Exact",
        markersize=6,
        color="black",
    )

    lines!(
        ax,
        OUTPUT_TIMES,
        result.twa;
        label="TWA (Gaussian)",
        linestyle=:dash,
        linewidth=1.5,
        color="red",
    )

    lines!(
        ax,
        OUTPUT_TIMES,
        result.dtwa;
        label="DTWA",
        linewidth=2.5,
        color="blue",
    )

    for cluster_size in DC_CLUSTER_SIZES
        lines!(
            ax,
            OUTPUT_TIMES,
            result.dctwa[cluster_size];
            label="dcTWA k=$cluster_size",
            linewidth=2,
            linestyle=:dot,
            color="green",
        )

         lines!(
            ax,
            OUTPUT_TIMES,
            result.gctwa[cluster_size];
            label="gcTWA k=$cluster_size",
            linewidth=2,
            color="yellow",
        )
    end
end

# Same x range is already imposed explicitly on both axes.
# Do not link y-axes because the requested y-ranges are different.

axislegend(
    axes[1];
    position=:rb,
    labelsize=8,
    patchsize=(20, 10),
    rowgap=2,
)

Label(
    fig[0, :],
    "Long-range Ising dynamics — L = $L",
    fontsize=12,
)

output = joinpath(
    @__DIR__,
    "..",
    "docs",
    "src",
    "assets",
    "ising_twa_dtwa_dctwa_comparison.png",
)

mkpath(dirname(output))
save(output, fig)

println()
println("Saved figure to:")
println(output)

fig
