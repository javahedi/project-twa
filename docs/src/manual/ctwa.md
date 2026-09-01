# Cluster Truncated-Wigner Approximation

The cluster truncated-Wigner approximation (CTWA) extends the
truncated-Wigner approach by grouping physical spins into clusters.

Instead of representing every spin only by three classical spin
components, CTWA represents a cluster of ``n`` spin-``1/2`` degrees of
freedom using the complete set of nontrivial Pauli strings acting within
that cluster.

This allows quantum correlations **inside each cluster** to be included
directly in the cluster phase-space description, while interactions
between different clusters are treated semiclassically.

In `TWA.jl`, a CTWA method is constructed with

```julia
method = CTWA(
    cluster_size=2,
    trajectories=1000,
    sampling=DiscreteSampling(),
)
```

and uses the same `simulate` interface as DTWA:

```julia
result = simulate(
    model,
    state,
    method;
    tspan=(0.0, 10.0),
    saveat=0.05,
)
```

## Clustering the spin system

Consider a chain of ``N`` physical spins.

For cluster size

```julia
cluster_size = 2
```

the chain is partitioned schematically as

```text
spin:      1   2 | 3   4 | 5   6 | 7   8
cluster:     1   |   2   |   3   |   4
```

For cluster size four,

```text
spin:      1   2   3   4 | 5   6   7   8
cluster:          1       |         2
```

The current contiguous clustering scheme requires the number of physical
sites to be divisible by the requested cluster size.

For example, a 100-site chain can be simulated with

```julia
CTWA(cluster_size=2, trajectories=1000)
```

or

```julia
CTWA(cluster_size=4, trajectories=1000)
```

because both cluster sizes divide 100.

## Cluster phase space

A cluster containing ``n`` spin-``1/2`` particles has Hilbert-space
dimension

```math
D = 2^n.
```

A complete operator basis consists of the identity together with all
nontrivial Pauli strings.

The number of nontrivial cluster generators is therefore

```math
D^2-1
=
4^n-1.
```

For example:

| Cluster size ``n`` | Hilbert dimension ``2^n`` | Phase-space coordinates ``4^n-1`` |
|---:|---:|---:|
| 1 | 2 | 3 |
| 2 | 4 | 15 |
| 3 | 8 | 63 |
| 4 | 16 | 255 |
| 5 | 32 | 1023 |
| 6 | 64 | 4095 |

This exponential growth is the central computational tradeoff of CTWA.

Larger clusters retain a richer description of intra-cluster quantum
correlations, but the number of dynamical variables per cluster grows
as ``4^n-1``.

## Pauli-string basis

For a two-spin cluster, the nontrivial operator basis contains

```text
IX  IY  IZ
XI  XX  XY  XZ
YI  YX  YY  YZ
ZI  ZX  ZY  ZZ
```

where, for example,

```math
XY =
\sigma_1^x\sigma_2^y.
```

Each nontrivial Pauli string is associated with one classical CTWA
phase-space coordinate.

For an ``n``-spin cluster, the cluster variables can therefore be
written as

```math
x^\mu,
\qquad
\mu=1,\ldots,4^n-1.
```

The single-spin case ``n=1`` contains only

```math
X,\quad Y,\quad Z,
```

and therefore reduces structurally to the familiar three-component
classical-spin phase space.

## CTWA dynamics

Let ``X_\mu`` denote the cluster Pauli-string generators.

Their commutators define the cluster algebra,

```math
[X_\mu,X_\nu]
=
2i
\sum_\rho
f_{\mu\nu\rho}X_\rho,
```

where ``f_{\mu\nu\rho}`` are the structure constants.

The corresponding classical CTWA equations take the form

```math
\dot{x}_c^\mu
=
2
\sum_{\nu,\rho}
f_{\mu\nu\rho}
\frac{\partial H_W}{\partial x_c^\nu}
x_c^\rho ,
```

where ``c`` labels the cluster and ``H_W`` is the phase-space
Hamiltonian.

Interactions entirely inside a cluster contribute to the cluster
Hamiltonian directly, whereas interactions connecting different
clusters produce couplings between their phase-space variables.

The construction reduces to the ordinary single-spin classical
equations when

```julia
cluster_size=1
```

so the single-spin and cluster approximations share the same underlying
framework.

## Sampling is independent of clustering

An important design choice in `TWA.jl` is that

> **cluster size and initial phase-space sampling are independent
> approximation choices.**

The cluster dynamics are selected by `CTWA`, while the initial
phase-space representation is selected separately through the
`sampling` argument.

Two sampling strategies are currently available:

```julia
GaussianSampling()
```

and

```julia
DiscreteSampling()
```

Thus,

```julia
CTWA(
    cluster_size=2,
    trajectories=1000,
    sampling=GaussianSampling(),
)
```

and

```julia
CTWA(
    cluster_size=2,
    trajectories=1000,
    sampling=DiscreteSampling(),
)
```

have the **same cluster size and cluster equations of motion** but
different representations of the initial quantum state.

This makes direct comparisons between Gaussian and discrete cluster
sampling possible.

## Gaussian CTWA

Gaussian CTWA, abbreviated here as **gcTWA**, approximates the initial
cluster phase-space distribution using Gaussian random variables chosen
to reproduce the required initial moments.

It is selected with

```julia
method = CTWA(
    cluster_size=2,
    trajectories=1000,
    sampling=GaussianSampling(),
)
```

For example:

```julia
result = simulate(
    model,
    state,
    method;
    tspan=(0.0, 4.0),
    saveat=0.05,
)
```

The Gaussian construction provides a continuous approximation to the
initial cluster Wigner distribution.

## Discrete CTWA

Discrete CTWA, abbreviated here as **dcTWA**, uses discrete
phase-space samples to initialize the cluster coordinates.

It is selected with

```julia
method = CTWA(
    cluster_size=2,
    trajectories=1000,
    sampling=DiscreteSampling(),
)
```

The key point is that the ``4^n-1`` cluster coordinates are **not
sampled independently**.

Instead, the physical spins inside the cluster are first sampled from
their discrete single-spin phase spaces.

For an ``n``-spin product state, let

```math
s_{j,\alpha}
```

denote the sampled component ``\alpha`` of physical spin ``j``.

A cluster Pauli string

```math
P_\mu
=
\sigma_1^{a_1}
\sigma_2^{a_2}
\cdots
\sigma_n^{a_n}
```

is then assigned the phase-space value

```math
x^\mu(0)
=
\prod_{j:a_j\neq I}
s_{j,a_j}.
```

For example, in a two-spin cluster,

```math
x^{XI} = s_{1,x},
\qquad
x^{IY} = s_{2,y},
```

and consequently

```math
x^{XY}
=
s_{1,x}s_{2,y}
=
x^{XI}x^{IY}.
```

Similarly,

```math
x^{ZZ}
=
s_{1,z}s_{2,z}.
```

The multi-spin cluster coordinates are therefore constructed
consistently from the **same sampled physical-spin configuration**.

This preserves the product structure of the discrete phase-space
representation instead of treating every cluster generator as an
independent random variable.

## Example: an initially polarized cluster

For

```julia
state = Up()
```

each physical spin is polarized along ``+z``.

In the discrete single-spin representation,

```math
s_z=1,
\qquad
s_x=\pm1,
\qquad
s_y=\pm1.
```

For a two-spin cluster this immediately implies

```math
x^{ZI}=1,
\qquad
x^{IZ}=1,
\qquad
x^{ZZ}=1,
```

while coordinates containing transverse components fluctuate between
trajectories.

For example,

```math
x^{XX}
=
x^{XI}x^{IX}.
```

This illustrates why dcTWA samples physical spins first and constructs
the complete cluster phase-space point afterward.

## TWA, DTWA, gcTWA, and dcTWA

The relationship between the methods can be summarized as

| | Single-spin ``n=1`` | Cluster ``n>1`` |
|---|---|---|
| **Gaussian sampling** | TWA | gcTWA |
| **Discrete sampling** | DTWA | dcTWA |

In `TWA.jl`, these cases can be written as follows.

Traditional single-spin Gaussian TWA:

```julia
twa = CTWA(
    cluster_size=1,
    trajectories=1000,
    sampling=GaussianSampling(),
)
```

Discrete TWA:

```julia
dtwa = DTWA(
    trajectories=1000,
)
```

Gaussian CTWA:

```julia
gctwa = CTWA(
    cluster_size=2,
    trajectories=1000,
    sampling=GaussianSampling(),
)
```

Discrete CTWA:

```julia
dctwa = CTWA(
    cluster_size=2,
    trajectories=1000,
    sampling=DiscreteSampling(),
)
```

This organization is intentional: Gaussian and discrete sampling are
properties of the initial phase-space representation rather than
separate dynamical approximations.

## Choosing a cluster size

Increasing the cluster size increases the amount of quantum information
represented explicitly inside each cluster.

However, the phase-space dimension grows exponentially:

```math
4^n-1.
```

Consequently, cluster size should be treated as a convergence and
accuracy parameter subject to computational constraints.

A practical study might compare

```julia
CTWA(
    cluster_size=1,
    trajectories=1000,
    sampling=DiscreteSampling(),
)

CTWA(
    cluster_size=2,
    trajectories=1000,
    sampling=DiscreteSampling(),
)

CTWA(
    cluster_size=4,
    trajectories=1000,
    sampling=DiscreteSampling(),
)
```

when the system size is compatible with all three clusterings.

!!! note
    Increasing the cluster size does not by itself guarantee a specific
    error reduction for every observable, model, and time scale.
    Results should be checked for convergence in both cluster size and
    trajectory number whenever possible.

## Monte Carlo convergence

Like DTWA, CTWA obtains observables from an ensemble of trajectories.

The parameter

```julia
trajectories=1000
```

controls the number of sampled initial conditions.

Convergence with trajectory number should be checked separately from
convergence with cluster size.

These are two different numerical questions:

```text
trajectories
    │
    └── controls Monte Carlo sampling noise

cluster_size
    │
    └── controls the cluster approximation
```

For gcTWA and dcTWA, the sampling strategy introduces a third
independent choice:

```text
CTWA calculation
     │
     ├── cluster size
     ├── sampling strategy
     └── number of trajectories
```

## Example

Consider the long-range Ising model

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

A two-spin discrete CTWA calculation is

```julia
result = simulate(
    model,
    state,
    CTWA(
        cluster_size=2,
        trajectories=1000,
        sampling=DiscreteSampling(),
    );
    tspan=(0.0, 4.0),
    saveat=0.05,
)
```

The corresponding Gaussian calculation requires only changing the
sampling strategy:

```julia
result = simulate(
    model,
    state,
    CTWA(
        cluster_size=2,
        trajectories=1000,
        sampling=GaussianSampling(),
    );
    tspan=(0.0, 4.0),
    saveat=0.05,
)
```

The model, state, cluster size, trajectory number, and simulation
interval remain unchanged.

## Reference

For cluster TWA and the discrete cluster sampling formulation used here,
see:


A. Braemer, J. Vahedi, and M. Gärttner,  
*Cluster truncated Wigner approximation for bond-disordered Heisenberg
spin models*,  
**Physical Review B 110**, 054204 (2024).  
DOI: `10.1103/PhysRevB.110.054204`

See [References and Citation](../references.md) for the complete
citation and BibTeX entry.

## Next steps

See [Observables](observables.md) for extracting physical quantities
from DTWA and CTWA results.

The [Long-range Ising benchmark](../examples/long_range_ising.md)
provides a direct comparison of the exact solution, TWA, DTWA, gcTWA,
and dcTWA.