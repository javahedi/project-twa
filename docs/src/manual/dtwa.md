# Discrete Truncated-Wigner Approximation

The discrete truncated-Wigner approximation (DTWA) is a semiclassical
phase-space method for simulating the dynamics of interacting quantum
spin systems.

Instead of evolving the exponentially large many-body wavefunction,
DTWA samples an ensemble of classical spin configurations from a
discrete phase-space representation of the initial quantum state.
Each configuration is evolved independently using classical equations
of motion, and observables are obtained by averaging over the resulting
trajectories.

In `TWA.jl`, DTWA is selected with

```julia
method = DTWA(
    trajectories=1000,
)
```

and used through the common simulation interface:

```julia
result = simulate(
    model,
    state,
    method;
    tspan=(0.0, 10.0),
    saveat=0.05,
)
```

## Basic idea

Consider a system of ``N`` spin-``1/2`` degrees of freedom.

`TWA.jl` uses Pauli variables

```math
\mathbf{s}_i =
\left(
s_i^x,
s_i^y,
s_i^z
\right)
```

for every physical site ``i``.

The DTWA calculation consists schematically of

```text
Quantum initial state
        │
        ▼
Discrete phase-space sampling
        │
        ├── trajectory 1
        ├── trajectory 2
        ├── trajectory 3
        │       ⋮
        └── trajectory M
                │
                ▼
       Classical evolution
                │
                ▼
       Ensemble averaging
                │
                ▼
          Observables
```

The number of sampled trajectories is controlled by

```julia
DTWA(
    trajectories=M,
)
```

Increasing `M` reduces Monte Carlo sampling noise but increases the
computational cost.

## Discrete initial sampling

The defining feature of DTWA is the use of a **discrete spin phase
space** for the initial conditions.

For example, consider a spin polarized along ``+z``,

```math
|\uparrow_z\rangle .
```

Its exact one-spin expectation values are

```math
\langle \sigma^x\rangle = 0,
\qquad
\langle \sigma^y\rangle = 0,
\qquad
\langle \sigma^z\rangle = 1.
```

In the discrete representation used by DTWA, an individual trajectory
can instead start with

```math
s^z = 1,
\qquad
s^x = \pm 1,
\qquad
s^y = \pm 1,
```

where the transverse components are sampled from the appropriate
discrete distribution.

Their ensemble averages reproduce the initial spin expectation values:

```math
\overline{s^x}=0,
\qquad
\overline{s^y}=0,
\qquad
\overline{s^z}=1.
```

The fluctuations of the individual trajectories encode information
about the initial quantum fluctuations.

The user normally does not perform this sampling explicitly. Given

```julia
state = Up()
```

and

```julia
method = DTWA(
    trajectories=1000,
)
```

`simulate` constructs the initial discrete trajectories automatically.

## Equations of motion

For a Hamiltonian written in the Pauli convention

```math
H =
\sum_{i<j}
\sum_{\alpha=x,y,z}
J_{ij}^{\alpha}
\sigma_i^\alpha\sigma_j^\alpha
+
\sum_i
\sum_{\alpha=x,y,z}
h_i^\alpha\sigma_i^\alpha,
```

the classical spin variables evolve according to precession in an
effective field.

In the convention used by `TWA.jl`,

```math
\dot{\mathbf{s}}_i
=
2\,
\mathbf{h}_{\mathrm{eff},i}
\times
\mathbf{s}_i ,
```

with

```math
h_{\mathrm{eff},i}^{\alpha}
=
h_i^\alpha
+
\sum_{j\neq i}
J_{ij}^{\alpha}s_j^\alpha .
```

!!! note "Spin convention"
    `TWA.jl` uses Pauli matrices ``\sigma^\alpha`` as the microscopic
    variables. This is why the equations contain the factor of `2`
    shown above.

Each sampled initial condition is evolved independently using these
equations.

## Example

Consider a long-range Ising chain,

```julia
using TWA

model = SpinModel(
    Chain(100),
    PowerLaw(
        Ising(:x; J=1.0);
        α=3.0,
    ),
)

state = Up()
```

A DTWA calculation can be performed with

```julia
result = simulate(
    model,
    state,
    DTWA(
        trajectories=1000,
    );
    tspan=(0.0, 4.0),
    saveat=0.05,
)
```

The longitudinal magnetization is then obtained with

```julia
mz = expectation(
    result,
    Magnetization(:z),
)
```

and the corresponding times with

```julia
ts = times(result)
```

## Monte Carlo convergence

DTWA observables are estimated from a finite ensemble of trajectories.

Consequently, numerical results contain statistical sampling noise.

A calculation should therefore be checked for convergence with respect
to the number of trajectories.

For example,

```julia
DTWA(trajectories=100)
DTWA(trajectories=1000)
DTWA(trajectories=10000)
```

represent increasingly large Monte Carlo ensembles.

The appropriate number of trajectories depends on the observable,
system size, model, and desired statistical accuracy.

## DTWA and TWA

The distinction between traditional TWA and DTWA is particularly clear
for spin-``1/2`` systems.

Both approaches evolve semiclassical spin variables. The important
difference is the representation of the initial quantum fluctuations.

In the terminology used throughout these documentation examples:

| Method | Cluster size | Initial sampling |
|---|---:|---|
| TWA | 1 | Gaussian |
| DTWA | 1 | Discrete |

Within the unified CTWA machinery, the Gaussian single-spin case can
therefore be represented as

```julia
CTWA(
    cluster_size=1,
    trajectories=1000,
    sampling=GaussianSampling(),
)
```

whereas discrete single-spin TWA is exposed directly through

```julia
DTWA(
    trajectories=1000,
)
```

This distinction is useful when comparing the effect of Gaussian and
discrete phase-space sampling while keeping the underlying classical
single-spin dynamics fixed.

## Reference

The DTWA approach for interacting spin systems is based on the discrete
phase-space formulation developed in:

J. Schachenmayer, A. Pikovski, and A. M. Rey,  
*Many-Body Quantum Spin Dynamics with Monte Carlo Trajectories on a
Discrete Phase Space*,  
**Physical Review X 5**, 011022 (2015).  
DOI: `10.1103/PhysRevX.5.011022`

See [References and Citation](../references.md) for the full citation
and BibTeX entry.

## Next steps

See [CTWA](ctwa.md) for the cluster generalization and the distinction
between Gaussian and discrete cluster sampling.

The [Long-range Ising benchmark](../examples/long_range_ising.md)
compares TWA, DTWA, Gaussian CTWA, and discrete CTWA against an exact
solution.