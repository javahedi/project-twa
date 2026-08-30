using TWA
using CairoMakie
using Random

# ---------------------------------------------------------------------------
# Long-range transverse Ising benchmark
#
# H = sum_{i<j} J / r_ij^α σᵢˣ σⱼˣ
#
# Initial state:
#     |↑z ↑z ... ↑z>
#
# For this commuting Ising Hamiltonian,
#
# <σᵢᶻ(t)> = ∏_{j ≠ i} cos(2 J_ij t).
#
# We compare that exact result against DTWA and CTWA with several
# cluster sizes.
# ---------------------------------------------------------------------------

const L = 24
const J = 1.0
const α = 3.0

const TMAX = 4.0
const DT = 0.05

# Keep this deliberately modest while developing the docs figure.
# Increase later for the final committed asset.
const NTRAJ = 32

chain = Chain(L)

model = SpinModel(
    chain,
    PowerLaw(
        Ising(:x; J=J);
        α=α,
    ),
)

state = Up()

output_times = collect(0.0:DT:TMAX)

# ---------------------------------------------------------------------------
# Exact solution
# ---------------------------------------------------------------------------

function exact_magnetization_z(
    model::SpinModel,
    ts,
)
    nsites = geometry(model).nsites

    values = zeros(
        Float64,
        length(ts),
    )

    for (time_index, t) in pairs(ts)
        total = 0.0

        for i in 1:nsites
            local_magnetization = 1.0

            for j in 1:nsites
                i == j && continue

                Jx, _, _ =
                    coupling_components(
                        model,
                        i,
                        j,
                    )

                local_magnetization *=
                    cos(
                        2 *
                        Jx *
                        t,
                    )
            end

            total += local_magnetization
        end

        values[time_index] =
            total / nsites
    end

    return values
end

exact =
    exact_magnetization_z(
        model,
        output_times,
    )

# ---------------------------------------------------------------------------
# DTWA
# ---------------------------------------------------------------------------

println("Running DTWA ...")

dtwa_result = simulate(
    model,
    state,
    DTWA(
        trajectories=NTRAJ,
    );
    tspan=(0.0, TMAX),
    saveat=output_times,
    rng=Xoshiro(100),
)

dtwa =
    expectation(
        dtwa_result,
        Magnetization(:z),
    )

# ---------------------------------------------------------------------------
# CTWA
# ---------------------------------------------------------------------------

cluster_sizes =
    (2, 3, 4, 6)

ctwa =
    Dict{Int,Vector{Float64}}()

for cluster_size in cluster_sizes
    println(
        "Running CTWA cluster size $cluster_size ...",
    )

    result = simulate(
        model,
        state,
        CTWA(
            cluster_size=cluster_size,
            trajectories=NTRAJ,
        );
        tspan=(0.0, TMAX),
        saveat=output_times,
        rng=Xoshiro(
            1000 +
            cluster_size,
        ),
    )

    ctwa[cluster_size] =
        collect(
            expectation(
                result,
                Magnetization(:z),
            ),
        )
end

# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

fig = Figure(
    size=(820, 700),
)

ax1 = Axis(
    fig[1, 1],
    xlabel="tJ",
    ylabel="⟨Sᶻ⟩ / L",
    title="Long-range Ising dynamics — L = $L, α = $α",
)

lines!(
    ax1,
    output_times,
    exact;
    label="Exact",
    linewidth=3,
)

lines!(
    ax1,
    output_times,
    dtwa;
    label="DTWA",
)

for cluster_size in cluster_sizes
    lines!(
        ax1,
        output_times,
        ctwa[cluster_size];
        label="CTWA k=$cluster_size",
    )
end

axislegend(
    ax1;
    position=:rt,
)

ax2 = Axis(
    fig[2, 1],
    xlabel="tJ",
    ylabel="|mᶻ - mᶻexact|",
    yscale=log10,
)

lines!(
    ax2,
    output_times,
    abs.(dtwa .- exact) .+ eps();
    label="DTWA",
)

for cluster_size in cluster_sizes
    lines!(
        ax2,
        output_times,
        abs.(
            ctwa[cluster_size] .-
            exact
        ) .+
        eps();
        label="CTWA k=$cluster_size",
    )
end

linkxaxes!(
    ax1,
    ax2,
)

output =
    joinpath(
        @__DIR__,
        "..",
        "docs",
        "src",
        "assets",
        "long_range_ising_comparison.png",
    )

mkpath(
    dirname(output),
)

save(
    output,
    fig,
)

println("Saved figure to:")
println(output)

fig
