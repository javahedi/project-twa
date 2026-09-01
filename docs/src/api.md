# API Reference

This page documents the main public interfaces of `TWA.jl`.

For a guided introduction, start with the [Quick Start](quickstart.md).
For the physical and numerical meaning of the individual components,
see the corresponding pages in the Manual.

## Models

### Geometry

```@docs
Chain
OpenBoundary
PeriodicBoundary
bonds
distance
```

### Interactions

```@docs
XXZ
Ising
NearestNeighbor
PowerLaw
```

A bare pair interaction passed to `SpinModel` is interpreted as a
nearest-neighbor interaction. For example,

```julia
SpinModel(
    Chain(50),
    XXZ(J=1.0, Δ=0.5),
)
```

is equivalent to using an explicit nearest-neighbor coupling profile.

Use `PowerLaw` when the spatial dependence should instead be long
ranged.

### Fields and Hamiltonians

```@docs
Field
field_components
Hamiltonian
```

### Spin models

```@docs
SpinModel
geometry
hamiltonian
terms
coupling_strength
coupling_components
```

See [Models](manual/models.md) for examples of constructing
Hamiltonians and geometries.

## Initial States

```@docs
AbstractState
Polarized
Up
Down
DomainWall
state_direction
```

Uniform polarization along an arbitrary Cartesian direction can be
specified with

```julia
Polarized(:x)
Polarized(:y; sign=-1)
Polarized(:z; sign=1)
```

The convenience constructors

```julia
Up()
Down()
```

correspond to polarization along `+z` and `-z`, respectively.

See [Initial States](manual/states.md).

## Approximation Methods

```@docs
AbstractApproximation
DTWA
CTWA
```

### CTWA sampling

CTWA separates the cluster approximation from the representation used
to sample the initial quantum state.

```@docs
AbstractCTWASampling
GaussianSampling
DiscreteSampling
```

For example,

```julia
CTWA(
    cluster_size=2,
    trajectories=1000,
    sampling=GaussianSampling(),
)
```

selects Gaussian cluster sampling, while

```julia
CTWA(
    cluster_size=2,
    trajectories=1000,
    sampling=DiscreteSampling(),
)
```

selects discrete cluster sampling.

See [DTWA](manual/dtwa.md) and [CTWA](manual/ctwa.md) for the
methodological background.

## Simulation

The main simulation interface is

```@docs
simulate
```

Both DTWA and CTWA use the same high-level pattern:

```julia
result = simulate(
    model,
    state,
    method;
    tspan=(0.0, 10.0),
    saveat=0.05,
)
```

### Results

```@docs
AbstractSimulationResult
DTWAResult
CTWAResult
times
```

DTWA additionally provides direct trajectory inspection through

```@docs
ntrajectories
trajectory
```

!!! note
    The internal storage layouts of `DTWAResult` and `CTWAResult` are
    different. Physical observables should normally be accessed through
    `expectation` rather than by depending directly on the trajectory
    array representation.

## Observables

All physical observables derive from

```@docs
AbstractObservable
```

and are evaluated using

```@docs
expectation
```

### Magnetization

```@docs
LocalMagnetization
Magnetization
```

For example,

```julia
mz = expectation(
    result,
    Magnetization(:z),
)
```

returns the magnetization density

```math
\frac{\langle S^z\rangle}{N}.
```

### Two-point correlations

```@docs
Correlation
ConnectedCorrelation
```

For example,

```julia
Correlation(
    :z, 1,
    :z, 2,
)
```

represents

```math
\langle
\sigma_1^z
\sigma_2^z
\rangle.
```

### Collective fluctuations

```@docs
SpinSecondMoment
SpinVariance
```

These use the normalization

```math
\frac{\langle(S^\alpha)^2\rangle}{N^2}
```

and

```math
\frac{
\langle(S^\alpha)^2\rangle
-
\langle S^\alpha\rangle^2
}{N^2},
```

respectively.

### Correlation profiles

```@docs
CorrelationProfile
correlation_profile
```

### Structure factors

```@docs
StructureFactor
structure_factor
```

See [Observables](manual/observables.md) for operator conventions,
normalizations, and examples.

# Advanced API

The interfaces below expose lower-level components of the DTWA and CTWA
implementations.

They are useful for developing new sampling schemes, inspecting the
cluster representation, implementing specialized calculations, or
working directly with the compiled dynamics.

Most users do not need these interfaces.

## Pauli-string basis

CTWA represents an `n`-spin cluster using the `4^n - 1` nontrivial
Pauli strings.

The basis representation can be inspected with

```julia
basis = PauliStringBasis(2)
```

The associated exported utilities are

```text
PauliStringBasis
basis_size

pauli_code
pauli_index
pauli_digits

pauli_symbol
pauli_symbols
```

## Pauli algebra

The sparse Pauli-string algebra is exposed through

```text
pauli_product_digit
pauli_product
pauli_commutator
pauli_anticommutes
```

These functions underlie the CTWA equations of motion.

## Clustering

The cluster decomposition is represented by `Clustering`.

Related inspection functions include

```text
cluster_count
cluster_sites

site_cluster
site_position
site_cluster_position

same_cluster

local_pauli_digits
local_pauli_index
```

## Compilation

The package exposes compiled representations for both DTWA and CTWA.

Relevant advanced types and functions include

```text
CompiledDTWA

CompiledCTWA
CompiledCTWAHamiltonian
CompiledCTWAAlgebra

compile
compile_ctwa_hamiltonian
compile_ctwa_algebra
```

These interfaces are primarily intended for advanced use and package
development.

## Sampling

Low-level sampling interfaces include

```text
sample_dtwa
sample_initial

compile_ctwa_sampling
sample_ctwa
sample_ctwa!
```

Normal simulations perform sampling automatically through `simulate`.

## Low-level dynamics

The raw equations of motion are available through

```text
dtwa_rhs!
ctwa_rhs!
```

These functions are intended primarily for advanced use and internal
development.