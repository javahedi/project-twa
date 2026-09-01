# Initial States

Initial states in `TWA.jl` describe the quantum state of the spin system
at ``t=0``.

The initial state is specified independently of both the Hamiltonian and
the approximation method:

```julia
using TWA

model = SpinModel(
    Chain(50),
    XXZ(
        J=1.0,
        Δ=0.5,
    ),
)

state = Up()

method = DTWA(
    trajectories=1000,
)
```

The same `state` can therefore be used with DTWA, Gaussian CTWA, or
discrete CTWA.

## Fully polarized states

### Up

```julia
state = Up()
```

represents the product state in which every spin is polarized along the
positive ``z`` direction,

```math
|\psi_0\rangle
=
|\uparrow_z\rangle^{\otimes N}.
```

For every site,

```math
\langle \sigma^z_i\rangle = 1,
\qquad
\langle \sigma^x_i\rangle
=
\langle \sigma^y_i\rangle
=
0.
```

For example:

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

### Down

```julia
state = Down()
```

represents the product state

```math
|\psi_0\rangle
=
|\downarrow_z\rangle^{\otimes N},
```

with

```math
\langle \sigma^z_i\rangle = -1.
```

It can be used in exactly the same way:

```julia
result = simulate(
    model,
    Down(),
    DTWA(
        trajectories=1000,
    );
    tspan=(0.0, 4.0),
)
```

## Polarized states

Uniform product states polarized along an arbitrary Cartesian axis can
be constructed with

```julia
Polarized(axis; sign=1)
```

where `axis` is one of

```julia
:x
:y
:z
```

and `sign` is either `+1` or `-1`.

For example,

```julia
Polarized(:x)
```

represents polarization along ``+x``, while

```julia
Polarized(:y; sign=-1)
```

represents polarization along ``-y``.

The convenience constructors

```julia
Up()
Down()
```

are equivalent to

```julia
Polarized(:z; sign=1)
Polarized(:z; sign=-1)
```

respectively.

## Domain-wall state

A spatially inhomogeneous product state can be constructed with

```julia
state = DomainWall()
```

The domain-wall state divides the chain into oppositely polarized
regions. It is useful for studying transport, relaxation, and the
spreading of an initially inhomogeneous spin profile.

Schematically,

```text
↑ ↑ ↑ ↑ ↑ │ ↓ ↓ ↓ ↓ ↓
```

For example:

```julia
model = SpinModel(
    Chain(50),
    XXZ(
        J=1.0,
        Δ=0.5,
    ),
)

result = simulate(
    model,
    DomainWall(),
    DTWA(
        trajectories=1000,
    );
    tspan=(0.0, 12.0),
)
```

The same state can be simulated with CTWA:

```julia
result = simulate(
    model,
    DomainWall(),
    CTWA(
        cluster_size=2,
        trajectories=1000,
        sampling=DiscreteSampling(),
    );
    tspan=(0.0, 12.0),
)
```

## States and phase-space sampling

A state specifies the **physical quantum initial condition**.

It does not specify how that state is represented in phase space.

For example,

```julia
state = Up()
```

can be combined with

```julia
DTWA(
    trajectories=1000,
)
```

or with Gaussian cluster sampling,

```julia
CTWA(
    cluster_size=2,
    trajectories=1000,
    sampling=GaussianSampling(),
)
```

or with discrete cluster sampling,

```julia
CTWA(
    cluster_size=2,
    trajectories=1000,
    sampling=DiscreteSampling(),
)
```

In all three cases the physical initial state is the same. What changes
is its phase-space representation and the approximation used for the
subsequent dynamics.

This distinction is important:

```text
Physical initial state
         │
         ▼
       State
         │
         ├──────── Gaussian sampling
         │
         └──────── Discrete sampling
                     │
                     ▼
             Initial trajectories
```

The sampling procedure is handled by the approximation method rather
than by the state itself.

## Example: polarized Ising dynamics

The long-range Ising benchmark uses

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

The initial spins point along ``z``, while the interaction acts along
``x``.

The same physical problem can then be studied using different
phase-space approximations:

```julia
twa = CTWA(
    cluster_size=1,
    trajectories=1000,
    sampling=GaussianSampling(),
)

dtwa = DTWA(
    trajectories=1000,
)

gctwa = CTWA(
    cluster_size=2,
    trajectories=1000,
    sampling=GaussianSampling(),
)

dctwa = CTWA(
    cluster_size=2,
    trajectories=1000,
    sampling=DiscreteSampling(),
)
```

This is the comparison developed in the
[Long-range Ising benchmark](../examples/long_range_ising.md).

## Next steps

Continue with [DTWA](dtwa.md) to learn how discrete single-spin
phase-space sampling works, or [CTWA](ctwa.md) for cluster
phase-space dynamics and Gaussian versus discrete cluster sampling.