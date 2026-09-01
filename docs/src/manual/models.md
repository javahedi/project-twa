# Models

A physical spin model in `TWA.jl` is represented by a `SpinModel`.

The model describes the Hamiltonian and the geometry of the spin system.
It is independent of the initial state and of the approximation method
used to simulate the dynamics.

A typical model is constructed as

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

The same `model` can subsequently be simulated using DTWA, Gaussian
CTWA, or discrete CTWA.

## Hamiltonian convention

`TWA.jl` uses Pauli matrices as the fundamental spin variables.

For pair interactions and local fields, the Hamiltonian convention is

```math
H =
\sum_{i<j}
\left(
J_{ij}^{x}\sigma_i^x\sigma_j^x
+
J_{ij}^{y}\sigma_i^y\sigma_j^y
+
J_{ij}^{z}\sigma_i^z\sigma_j^z
\right)
+
\sum_i
\left(
h_i^x\sigma_i^x
+
h_i^y\sigma_i^y
+
h_i^z\sigma_i^z
\right).
```

Here ``\sigma_i^\alpha`` denotes a Pauli matrix acting on site ``i``.

!!! note "Pauli-matrix convention"
    The package uses ``\sigma^\alpha`` rather than
    ``S^\alpha=\sigma^\alpha/2`` as the microscopic spin variables.
    This convention determines the numerical values of coupling constants,
    fields, equations of motion, and observables.

## Geometry

The geometry specifies the physical sites and their spatial relationships.

### One-dimensional chains

A chain of `L` spins is constructed with

```julia
geometry = Chain(L)
```

For example,

```julia
geometry = Chain(100)
```

creates a one-dimensional chain containing 100 spins.

The geometry is passed as the first argument of `SpinModel`:

```julia
model = SpinModel(
    Chain(100),
    Ising(:x; J=1.0),
)
```

The geometry is responsible for defining site relationships used by
spatially dependent interactions.

## XXZ interactions

The XXZ interaction can be constructed with

```julia
XXZ(
    J=1.0,
    Δ=0.5,
)
```

It corresponds to the spin interaction

```math
H_{\mathrm{XXZ}}
=
J
\sum_{i<j}
\left(
\sigma_i^x\sigma_j^x
+
\sigma_i^y\sigma_j^y
+
\Delta\,\sigma_i^z\sigma_j^z
\right),
```

with the set of interacting pairs determined by the spatial interaction
used in the model.

For example,

```julia
model = SpinModel(
    Chain(50),
    XXZ(
        J=1.0,
        Δ=0.5,
    ),
)
```

constructs the corresponding short-range chain model.

The isotropic point is obtained with

```julia
XXZ(
    J=1.0,
    Δ=1.0,
)
```

while other values of `Δ` introduce exchange anisotropy.

## Ising interactions

An Ising interaction is specified by its spin axis and coupling strength.

For example,

```julia
Ising(
    :x;
    J=1.0,
)
```

describes an interaction proportional to

```math
\sigma_i^x\sigma_j^x.
```

The interaction axis can be selected using `:x`, `:y`, or `:z`:

```julia
Ising(:x; J=1.0)
Ising(:y; J=1.0)
Ising(:z; J=1.0)
```

For an ``x``-Ising model,

```math
H_{\mathrm{Ising}}
=
\sum_{i<j}
J_{ij}
\sigma_i^x\sigma_j^x.
```

## Power-law interactions

Spatially long-range interactions can be introduced with `PowerLaw`.

For example,

```julia
model = SpinModel(
    Chain(100),
    PowerLaw(
        Ising(:x; J=1.0);
        α=3.0,
    ),
)
```

produces couplings that decay with distance as

```math
J_{ij}
=
\frac{J}{r_{ij}^{\alpha}},
```

where ``r_{ij}`` is the distance between sites ``i`` and ``j``.

The exponent `α` controls the interaction range.

For example,

```julia
α = 0.0
```

gives distance-independent all-to-all couplings, while

```julia
α = 3.0
```

gives a dipolar power-law decay in one dimension.

The same spatial wrapper can be applied to other supported pair
interactions. For example,

```julia
model = SpinModel(
    Chain(50),
    PowerLaw(
        XXZ(
            J=1.0,
            Δ=0.5,
        );
        α=3.0,
    ),
)
```

defines a long-range XXZ chain.

## Local fields

Local single-spin terms are represented with `Field`.

They contribute terms of the form

```math
H_{\mathrm{field}}
=
\sum_i
\mathbf{h}_i\cdot\boldsymbol{\sigma}_i.
```

For example, a field along the ``z`` direction can be specified with

```julia
Field(:z, 0.5)
```

and similarly for the other Cartesian directions.

Fields and pair interactions can be combined to construct more general
Hamiltonians.

## Building a spin model

The basic pattern is

```julia
model = SpinModel(
    geometry,
    hamiltonian,
)
```

For example, a long-range transverse Ising interaction is

```julia
model = SpinModel(
    Chain(100),
    PowerLaw(
        Ising(:x; J=1.0);
        α=3.0,
    ),
)
```

Once constructed, the model does not depend on the approximation used
to evolve it.

For example, the same model can be used with DTWA,

```julia
result = simulate(
    model,
    Up(),
    DTWA(
        trajectories=1000,
    );
    tspan=(0.0, 4.0),
)
```

or with CTWA,

```julia
result = simulate(
    model,
    Up(),
    CTWA(
        cluster_size=2,
        trajectories=1000,
        sampling=DiscreteSampling(),
    );
    tspan=(0.0, 4.0),
)
```

This separation between the physical model and the approximation method
is one of the central design principles of `TWA.jl`.

## Model, state, and method are independent

It is useful to think of a simulation as

```text
                 SpinModel
              ┌─────────────┐
              │  Geometry   │
              │      +      │
              │ Hamiltonian │
              └─────────────┘
                     │
                     │
       State ────────┼──────── Approximation
                     │
                     ▼
                  simulate
                     │
                     ▼
                   Result
```

The model answers:

> **What physical system is being simulated?**

The state answers:

> **How is the system initialized?**

The approximation answers:

> **How is the quantum dynamics approximated?**

Keeping these concepts separate makes it possible to reuse a physical
model across different initial states and approximation methods.

## Next steps

Continue with [Initial States](states.md) to define the starting quantum
state, or see [DTWA](dtwa.md) and [CTWA](ctwa.md) for the available
approximation methods.