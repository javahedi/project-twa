# Observables

Simulation results in `TWA.jl` contain an ensemble of phase-space
trajectories. Physical observables are obtained by averaging the
corresponding quantities over this ensemble.

The main user interface is

```julia
expectation(
    result,
    observable,
)
```

The same observable interface is used for DTWA and CTWA results.

## Conventions

`TWA.jl` uses Pauli matrices as the microscopic spin variables,

```math
\sigma_i^x,\qquad
\sigma_i^y,\qquad
\sigma_i^z.
```

The collective Pauli operator is defined as

```math
S^\alpha
=
\sum_{i=1}^{N}
\sigma_i^\alpha.
```

Consequently, the magnetization density is

```math
m_\alpha
=
\frac{\langle S^\alpha\rangle}{N}
=
\frac{1}{N}
\sum_i
\langle\sigma_i^\alpha\rangle.
```

For a fully polarized `Up()` state,

```math
m_z(0)=1.
```

!!! warning "Pauli versus spin-1/2 convention"
    Some literature defines physical spin operators as

    ```math
    \hat{s}_i^\alpha=\frac{1}{2}\sigma_i^\alpha.
    ```

    In that convention the fully polarized magnetization per spin is
    ``1/2`` rather than `1`.

    Keep this factor of two in mind when comparing `TWA.jl` results with
    formulas or figures using spin-``1/2`` operators.

## Simulation times

The time points stored in a simulation result are obtained with

```julia
ts = times(result)
```

For example,

```julia
result = simulate(
    model,
    state,
    DTWA(trajectories=1000);
    tspan=(0.0, 4.0),
    saveat=0.05,
)

ts = times(result)
```

The returned observable arrays use these saved time points.

## Magnetization

The magnetization density along a Cartesian direction is represented by
`Magnetization`.

For example,

```julia
mz = expectation(
    result,
    Magnetization(:z),
)
```

returns

```math
m_z(t)
=
\frac{1}{N}
\sum_i
\langle\sigma_i^z(t)\rangle.
```

Similarly,

```julia
mx = expectation(
    result,
    Magnetization(:x),
)

my = expectation(
    result,
    Magnetization(:y),
)
```

give the ``x`` and ``y`` magnetization densities.

A typical plotting workflow is

```julia
using CairoMakie

ts = times(result)

mz = expectation(
    result,
    Magnetization(:z),
)

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

## Local observables

Global magnetization can hide spatial structure.

Local spin expectation values instead resolve individual physical sites,

```math
\langle\sigma_i^\alpha(t)\rangle.
```

This is particularly useful for inhomogeneous initial states such as

```julia
DomainWall()
```

where transport and relaxation can be studied through the spatial
magnetization profile.

The observable layer operates on **physical spin sites**, even for CTWA
simulations. Users therefore do not need to convert cluster generators
back into physical spin variables manually.

## Two-point correlations

Two-point spin correlations have the form

```math
C_{ij}^{\alpha\beta}(t)
=
\left\langle
\sigma_i^\alpha(t)
\sigma_j^\beta(t)
\right\rangle.
```

They provide information that is not contained in single-spin
magnetization alone.

For equal components,

```math
C_{ij}^{\alpha\alpha}
=
\langle
\sigma_i^\alpha
\sigma_j^\alpha
\rangle.
```

For different components,

```math
C_{ij}^{\alpha\beta}
=
\langle
\sigma_i^\alpha
\sigma_j^\beta
\rangle.
```

The same physical definition is used whether sites ``i`` and ``j`` lie
inside the same CTWA cluster or in different clusters.

This is important because the clustering is an approximation detail,
not part of the physical observable definition.

## Connected correlations

The connected two-point correlation subtracts the product of the
one-point expectation values,

```math
C_{ij,\mathrm{conn}}^{\alpha\beta}
=
\langle
\sigma_i^\alpha\sigma_j^\beta
\rangle
-
\langle\sigma_i^\alpha\rangle
\langle\sigma_j^\beta\rangle.
```

Connected correlations isolate correlations beyond the product of the
local mean values.

For an initially uncorrelated product state, connected correlations
typically vanish initially between different sites and can develop as
the system evolves.

## Same-site Pauli products

Pauli matrices satisfy

```math
(\sigma_i^\alpha)^2 = I.
```

Therefore,

```math
\left\langle
(\sigma_i^\alpha)^2
\right\rangle
=
1
```

exactly.

`TWA.jl` respects this operator identity when evaluating same-site,
same-axis second moments.

This distinction matters because simply squaring a classical
phase-space coordinate is not, in general, the correct quantum
estimator for a same-site Pauli product.

For different Cartesian components on the same site,

```math
\sigma_i^\alpha\sigma_i^\beta
=
i\epsilon_{\alpha\beta\gamma}
\sigma_i^\gamma,
\qquad
\alpha\neq\beta,
```

so operator ordering becomes relevant.

The observable API therefore does not silently interpret an ordered
same-site, different-axis product as an ordinary product of commuting
classical numbers.

## Collective observables

Define the collective Pauli operators

```math
S^\alpha
=
\sum_i\sigma_i^\alpha.
```

Their first moment is

```math
\langle S^\alpha\rangle
=
\sum_i
\langle\sigma_i^\alpha\rangle.
```

The normalized magnetization used by `Magnetization` is

```math
\frac{\langle S^\alpha\rangle}{N}.
```

Collective second moments contain both local and nonlocal contributions,

```math
\left\langle
(S^\alpha)^2
\right\rangle
=
\sum_{i,j}
\left\langle
\sigma_i^\alpha\sigma_j^\alpha
\right\rangle.
```

Using ``(\sigma_i^\alpha)^2=I``,

```math
\left\langle
(S^\alpha)^2
\right\rangle
=
N
+
\sum_{i\neq j}
\left\langle
\sigma_i^\alpha\sigma_j^\alpha
\right\rangle.
```

Thus collective fluctuations depend directly on two-point correlations.

## Collective variance

The variance of a collective component is

```math
\operatorname{Var}(S^\alpha)
=
\left\langle
(S^\alpha)^2
\right\rangle
-
\left\langle
S^\alpha
\right\rangle^2.
```

Depending on the physical comparison, one may normalize this quantity by
``N`` or ``N^2``.

For example,

```math
\frac{\operatorname{Var}(S^\alpha)}{N^2}
```

measures fluctuations relative to the square of the total system size.

When comparing against a paper or exact calculation, always check both
the operator convention and the normalization.

## Mixed collective correlations

Mixed collective moments have the form

```math
\langle
S^\alpha S^\beta
\rangle
=
\sum_{i,j}
\langle
\sigma_i^\alpha\sigma_j^\beta
\rangle.
```

For ``\alpha\neq\beta``, these quantities can be complex because the
corresponding quantum operators need not commute.

A frequently useful real quantity is the symmetrized correlation,

```math
\frac{1}{2}
\left\langle
S^\alpha S^\beta
+
S^\beta S^\alpha
\right\rangle,
```

which is equivalent to

```math
\operatorname{Re}
\langle
S^\alpha S^\beta
\rangle
```

for Hermitian ``S^\alpha`` and ``S^\beta``.

Such quantities appear naturally in collective-spin fluctuation and
correlation benchmarks.

## Correlation profiles

For translational or approximately translational systems, it is often
useful to organize two-point correlations by physical separation rather
than by individual site pairs.

A correlation profile groups

```math
C_{ij}^{\alpha\beta}
```

according to the distance

```math
r = r_{ij}.
```

This provides a spatial picture of how correlations spread through the
system.

The separation is determined from the physical geometry and therefore
remains meaningful independently of the CTWA clustering.

## Structure factors

Spatial correlations can also be studied in momentum space through
structure factors.

A generic spin structure factor has the form

```math
S^{\alpha\beta}(q)
=
\frac{1}{N}
\sum_{i,j}
e^{iq(r_i-r_j)}
\left\langle
\sigma_i^\alpha
\sigma_j^\beta
\right\rangle.
```

Connected structure factors can similarly be constructed from connected
correlations.

Structure factors are useful for identifying characteristic spatial
wavevectors and ordering patterns that may not be apparent from global
magnetization alone.

## Physical sites versus CTWA clusters

One of the goals of the observable API is to hide the internal cluster
representation.

Suppose a CTWA calculation uses

```julia
CTWA(
    cluster_size=4,
    trajectories=1000,
    sampling=DiscreteSampling(),
)
```

Internally, each cluster evolves ``4^4-1=255`` phase-space coordinates.

Nevertheless, physical observables are still expressed in terms of

```math
\sigma_i^\alpha
```

on the original physical sites.

Conceptually,

```text
CTWA cluster coordinates
          │
          ▼
 physical Pauli operators
          │
          ▼
   correlations / moments
          │
          ▼
    user observables
```

The cluster basis is therefore an implementation of the approximation,
not a replacement for the physical observable language.

## Statistical convergence

DTWA and CTWA observables are Monte Carlo estimates.

For an observable ``O``,

```math
\langle O(t)\rangle
\approx
\frac{1}{M}
\sum_{m=1}^{M}
O_m(t),
```

where ``M`` is the number of trajectories.

Increasing

```julia
trajectories=M
```

reduces statistical sampling noise.

Observable convergence should therefore be checked with respect to
trajectory number, especially for correlation functions and higher
moments, which can require more samples than simple magnetization
observables.

For CTWA, trajectory convergence and cluster-size convergence are
separate questions.

## Example workflow

A typical analysis has the form

```julia
using TWA

model = SpinModel(
    Chain(100),
    PowerLaw(
        Ising(:x; J=1.0);
        α=3.0,
    ),
)

result = simulate(
    model,
    Up(),
    DTWA(
        trajectories=1000,
    );
    tspan=(0.0, 4.0),
    saveat=0.05,
)

ts = times(result)

mz = expectation(
    result,
    Magnetization(:z),
)
```

The same observable call can be used after changing the approximation:

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
    saveat=0.05,
)

mz = expectation(
    result,
    Magnetization(:z),
)
```

This makes comparisons between DTWA and CTWA possible without changing
the physical observable definition.

## Next steps

See the [Long-range Ising benchmark](../examples/long_range_ising.md)
for a comparison of exact dynamics, TWA, DTWA, gcTWA, and dcTWA.

For constructors and the complete public observable interface, see the
[API Reference](../api.md).