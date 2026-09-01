# TWA.jl

**TWA.jl** is a Julia package for semiclassical phase-space simulations of
interacting quantum spin systems using truncated-Wigner methods.

The package provides a unified interface for

- discrete truncated-Wigner approximation (**DTWA**),
- cluster truncated-Wigner approximation (**CTWA**),
- Gaussian cluster sampling (**gcTWA**), and
- discrete cluster sampling (**dcTWA**).

Models, initial states, approximation methods, and observables are kept
independent, making it possible to change the approximation without
redefining the physical problem.

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

## Quick start

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
    saveat=0.05,
)

mz = expectation(
    result,
    Magnetization(:z),
)
```

The same model and initial state can be simulated with CTWA simply by
changing the approximation method:

```julia
method = CTWA(
    cluster_size=2,
    trajectories=1000,
)

result = simulate(
    model,
    state,
    method;
    tspan=(0.0, 12.0),
    saveat=0.05,
)
```

## CTWA sampling

Cluster size and initial phase-space sampling are independent choices.

### Gaussian cluster sampling

```julia
method = CTWA(
    cluster_size=2,
    trajectories=1000,
    sampling=GaussianSampling(),
)
```

We refer to this as **gcTWA**.

### Discrete cluster sampling

```julia
method = CTWA(
    cluster_size=2,
    trajectories=1000,
    sampling=DiscreteSampling(),
)
```

We refer to this as **dcTWA**.

`DiscreteSampling()` is the default CTWA sampling strategy.

The resulting methods can be viewed as

| | Single spin | Cluster |
|---|---|---|
| Gaussian sampling | TWA | gcTWA |
| Discrete sampling | DTWA | dcTWA |

## Models

TWA.jl separates geometry from Hamiltonian terms.

For example, a nearest-neighbor XXZ chain is

```julia
model = SpinModel(
    Chain(100),
    XXZ(J=1.0, Δ=0.5),
)
```

while a long-range Ising model can be constructed as

```julia
model = SpinModel(
    Chain(100),
    PowerLaw(
        Ising(:x; J=1.0);
        α=3.0,
    ),
)
```

Open and periodic chains are supported through

```julia
Chain(100, OpenBoundary())
Chain(100, PeriodicBoundary())
```

The Hamiltonian convention uses Pauli matrices directly:

$$
H =\sum_{i<j}
\left( J^x_{ij}\sigma_i^x\sigma_j^x +
J^y_{ij}\sigma_i^y\sigma_j^y +
J^z_{ij}\sigma_i^z\sigma_j^z
\right)
+ \sum_i \mathbf h_i\cdot\boldsymbol{\sigma}_i.
$$

In particular, TWA.jl uses `σ`, not the convention
`S = σ/2`.

## Initial states

Available product-state descriptions include

```julia
Up()
Down()
DomainWall()

Polarized(:x)
Polarized(:y; sign=-1)
```

Initial states describe the physical quantum state independently of the
phase-space approximation used to simulate it.

## Observables

Physical observables are evaluated through the common `expectation`
interface:

```julia
mz = expectation(
    result,
    Magnetization(:z),
)
```

Available observable functionality includes

- local magnetization,
- collective magnetization,
- two-point correlations,
- connected correlations,
- collective second moments,
- collective variances,
- distance-resolved correlation profiles, and
- static structure factors.

The saved simulation times are obtained with

```julia
t = times(result)
```

The observable API works with both DTWA and CTWA results, so physical
measurements do not depend on the internal phase-space representation.

## Installation

TWA.jl is currently under development.

Install directly from the repository with Julia's package manager:

```julia
using Pkg

Pkg.add(
    url="https://github.com/javahedi/project-twa"
)
```

For development,

```julia
using Pkg

Pkg.develop(
    url="https://github.com/javahedi/project-twa"
)
```

## Documentation

The documentation contains

- a Quick Start,
- model and initial-state guides,
- DTWA and CTWA introductions,
- observable definitions,
- a long-range Ising example,
- an API reference, and
- methodological references and citation information.

To build the documentation locally:

```bash
julia --project=docs docs/make.jl
```

Then open

```text
docs/build/index.html
```

in a browser.

## References

### Discrete truncated-Wigner approximation

J. Schachenmayer, A. Pikovski, and A. M. Rey,  
*Many-Body Quantum Spin Dynamics with Monte Carlo Trajectories on a
Discrete Phase Space*,  
**Physical Review X 5**, 011022 (2015).  
https://doi.org/10.1103/PhysRevX.5.011022

### Cluster truncated-Wigner approximation

J. Wurtz, A. Polkovnikov, and D. Sels,  
*Cluster truncated Wigner approximation in strongly interacting systems*,  
**Annals of Physics 395**, 341–365 (2018).  
https://doi.org/10.1016/j.aop.2018.06.001

### Discrete cluster sampling

A. Braemer, J. Vahedi, and M. Gärttner,  
*Cluster truncated Wigner approximation for bond-disordered Heisenberg
spin models*,  
**Physical Review B 110**, 054204 (2024).  
https://doi.org/10.1103/PhysRevB.110.054204

## Citation

If you use TWA.jl in scientific work, please cite the methodological
references relevant to the approximation used in your calculation.

A permanent software citation for TWA.jl will be added when an archival
release is available.

## License

See [`LICENSE`](LICENSE) for licensing information.