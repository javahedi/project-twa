# TWA.jl

`TWA.jl` is a Julia package for semiclassical simulation of interacting
quantum spin systems using truncated-Wigner methods.

The package provides a unified interface for discrete TWA (DTWA) and
cluster TWA (CTWA), while keeping the physical model, initial state,
approximation method, and observables separate.

## At a glance

```julia
using TWA

model = SpinModel(
    Chain(50),
    XXZ(J=1.0, Δ=0.5),
)

state = DomainWall()

method = DTWA(
    trajectories=1000,
)

result = simulate(
    model,
    state,
    method;
    tspan=(0.0, 12.0),
)

mz = expectation(
    result,
    Magnetization(:z),
)
```

The same model and initial state can be simulated with another
approximation without changing the physical problem:

```julia
result = simulate(
    model,
    state,
    CTWA(
        cluster_size=2,
        trajectories=1000,
        sampling=DiscreteSampling(),
    );
    tspan=(0.0, 12.0),
)
```

## Package design

A simulation in `TWA.jl` is assembled from independent components:

```text
Model × State × Approximation
            │
            ▼
         simulate
            │
            ▼
          Result
            │
            ▼
       Observables
```

This structure makes it straightforward to compare approximation schemes
while keeping the Hamiltonian and initial state fixed.

## Methods

`TWA.jl` currently supports:

- **DTWA** — discrete phase-space sampling of individual spins.
- **CTWA** — cluster truncated-Wigner dynamics with configurable cluster
  size.
- **Gaussian CTWA (gcTWA)** — CTWA with Gaussian initial sampling.
- **Discrete CTWA (dcTWA)** — CTWA with discrete cluster sampling.

For CTWA, the cluster approximation and the initial sampling prescription
are independent choices:

```julia
gaussian = CTWA(
    cluster_size=2,
    trajectories=1000,
    sampling=GaussianSampling(),
)

discrete = CTWA(
    cluster_size=2,
    trajectories=1000,
    sampling=DiscreteSampling(),
)
```

This makes it possible to study separately the effect of clustering and
the effect of phase-space sampling.

## Start here

If you are new to the package, continue with the
[Quick Start](quickstart.md).

For more detail, see the manual pages on models, initial states, DTWA,
CTWA, and observables.