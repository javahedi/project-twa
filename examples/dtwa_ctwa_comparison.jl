using TWA
using CairoMakie
using Random

Random.seed!(1234)

model = SpinModel(
    Chain(8),
    XXZ(J=1.0, Δ=0.5),
)

state = DomainWall()

dtwa = simulate(
    model,
    state,
    DTWA(trajectories=500);
    tspan=(0.0, 6.0),
    saveat=0.1,
)

ctwa = simulate(
    model,
    state,
    CTWA(
        cluster_size=2,
        trajectories=500,
    );
    tspan=(0.0, 6.0),
    saveat=0.1,
)

mz_dtwa = expectation(
    dtwa,
    Magnetization(:z),
)

mz_ctwa = expectation(
    ctwa,
    Magnetization(:z),
)

fig = Figure(size=(760, 480))

ax = Axis(
    fig[1, 1],
    xlabel="Time",
    ylabel="⟨Sᶻ⟩ / N",
    title="DTWA and CTWA dynamics",
)

lines!(
    ax,
    times(dtwa),
    mz_dtwa,
    label="DTWA",
)

lines!(
    ax,
    times(ctwa),
    mz_ctwa,
    label="CTWA, cluster size 2",
)

axislegend(ax)

save(
    joinpath(
        @__DIR__,
        "..",
        "docs",
        "src",
        "assets",
        "dtwa_ctwa_magnetization.png",
    ),
    fig,
)

fig