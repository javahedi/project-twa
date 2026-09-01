# Quick Start

This page introduces the basic `TWA.jl` workflow.

A simulation is built from three main ingredients:

1. a physical model,
2. an initial state,
3. an approximation method.

The resulting trajectories can then be analyzed through observables.

## 1. Define a model

A spin model combines a geometry with a Hamiltonian.

For example, a nearest-neighbor XXZ chain can be written as

```julia
using TWA

model = SpinModel(
    Chain(50),
    XXZ(
        J=1.0,
        Δ=0.5,
    ),
)
```

A long-range Ising model can instead be constructed as

```julia
model = SpinModel(
    Chain(100),
    PowerLaw(
        Ising(:x; J=1.0);
        α=3.0,
    ),
)
```

The model definition is independent of the approximation method used
later.

## 2. Choose an initial state

Initial product states are represented separately from the Hamiltonian.

For example,

```julia
state = Up()
```

creates a fully polarized state along the positive ``z`` direction.

Other available product-state descriptions include

```julia
Down()
DomainWall()
```

and more general polarized states.

## 3. Choose an approximation

### DTWA

The discrete truncated-Wigner approximation is selected with

```julia
method = DTWA(
    trajectories=1000,
)
```

The `trajectories` parameter controls the number of Monte Carlo
trajectories used for phase-space averaging.

### CTWA

Cluster TWA groups physical spins into clusters and evolves the
corresponding cluster phase-space variables.

For example,

```julia
method = CTWA(
    cluster_size=2,
    trajectories=1000,
    sampling=DiscreteSampling(),
)
```

The sampling prescription can be changed independently:

```julia
method = CTWA(
    cluster_size=2,
    trajectories=1000,
    sampling=GaussianSampling(),
)
```

Thus, for the same cluster size, one can compare Gaussian and discrete
initial sampling.

## 4. Run the simulation

Once the model, state, and approximation are defined, call `simulate`:

```julia
result = simulate(
    model,
    state,
    method;
    tspan=(0.0, 4.0),
    saveat=0.05,
)
```

The same simulation interface is used for both DTWA and CTWA.

For example, changing

```julia
method = DTWA(
    trajectories=1000,
)
```

to

```julia
method = CTWA(
    cluster_size=2,
    trajectories=1000,
    sampling=DiscreteSampling(),
)
```

does not require any change to the model or the initial state.

## 5. Measure observables

Observables are evaluated from the simulation result.

For example, the magnetization density along the ``z`` direction is

```julia
mz = expectation(
    result,
    Magnetization(:z),
)
```

The corresponding simulation times are available through

```julia
ts = times(result)
```

A simple plot can then be produced with any Julia plotting package.

For example, with CairoMakie:

```julia
using CairoMakie

fig = Figure()

ax = Axis(
    fig[1, 1],
    xlabel="t",
    ylabel="⟨Sᶻ⟩ / N",
)

lines!(
    ax,
    ts,
    mz,
)

fig
```

## Comparing approximation methods

One of the main goals of the package is to make approximation
comparisons simple.

Consider the same model and initial state:

```julia
model = SpinModel(
    Chain(100),
    PowerLaw(
        Ising(:x; J=1.0);
        α=3.0,
    ),
)

state = Up()
```

Single-spin Gaussian TWA can be represented as CTWA with cluster size one:

```julia
twa = CTWA(
    cluster_size=1,
    trajectories=1000,
    sampling=GaussianSampling(),
)
```

Discrete TWA is

```julia
dtwa = DTWA(
    trajectories=1000,
)
```

Gaussian cluster TWA with two-spin clusters is

```julia
gctwa = CTWA(
    cluster_size=2,
    trajectories=1000,
    sampling=GaussianSampling(),
)
```

and discrete cluster TWA is

```julia
dctwa = CTWA(
    cluster_size=2,
    trajectories=1000,
    sampling=DiscreteSampling(),
)
```

Conceptually, these methods can be organized as

| Sampling | Single-spin | Cluster |
|---|---|---|
| Gaussian | TWA | gcTWA |
| Discrete | DTWA | dcTWA |

This separation between cluster size and sampling strategy is central to
the `TWA.jl` API.

## Next steps

The following manual pages provide more detail:

- [Models](manual/models.md)
- [Initial States](manual/states.md)
- [DTWA](manual/dtwa.md)
- [CTWA](manual/ctwa.md)
- [Observables](manual/observables.md)

The [Long-range Ising benchmark](examples/long_range_ising.md) shows how
the different approximation schemes can be compared against an exact
solution.