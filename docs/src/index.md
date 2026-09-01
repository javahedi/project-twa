# TWA.jl

`TWA.jl` provides a unified interface for truncated-Wigner approximations
for interacting spin systems.

The package separates the physical model, initial state, approximation
method, and simulation:

```julia
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
    tspan=(0, 12),
)
```

The same model and initial state can be evolved using different
truncated-Wigner approximations without changing the model definition.

## DTWA and CTWA

`TWA.jl` currently supports both the discrete truncated-Wigner
approximation (DTWA) and the cluster truncated-Wigner approximation
(CTWA).

For CTWA, the initial phase-space sampling can be selected independently
of the cluster dynamics:

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

This separation makes it possible to compare Gaussian and discrete
phase-space sampling at the same cluster size.

## Long-range Ising benchmark

As a benchmark, consider a one-dimensional spin chain with long-range
Ising interactions,

```math
H =
\sum_{i<j}
\frac{J}{r_{ij}^{\alpha}}
\sigma_i^x \sigma_j^x ,
```

starting from the fully polarized state

```math
|\psi_0\rangle =
|\uparrow_z \uparrow_z \cdots \uparrow_z\rangle .
```

For this commuting Ising Hamiltonian, the longitudinal magnetization can
be calculated exactly. For an individual spin,

```math
\langle \sigma_i^z(t) \rangle
=
\prod_{j\neq i}
\cos\left(2J_{ij}t\right),
```

and therefore

```math
m_z(t)
=
\frac{1}{L}
\sum_i
\langle \sigma_i^z(t) \rangle .
```

The figure below compares this exact solution with four
phase-space approximations for a chain of `L = 100` spins:

- **TWA** — single-spin Gaussian sampling (`cluster_size=1`).
- **DTWA** — single-spin discrete phase-space sampling.
- **gcTWA** — cluster TWA with Gaussian sampling and `cluster_size=2`.
- **dcTWA** — cluster TWA with discrete sampling and `cluster_size=2`.

The two panels show all-to-all interactions (`α = 0`) and dipolar
power-law interactions (`α = 3`).

![Comparison of exact Ising dynamics, TWA, DTWA, gcTWA, and dcTWA](assets/ising_twa_dtwa_dctwa_comparison.png)

This comparison illustrates two independent approximation choices in
`TWA.jl`: the **phase-space sampling prescription** and the
**cluster size**.

At cluster size one, Gaussian and discrete sampling correspond to the
traditional TWA and DTWA descriptions, respectively. Increasing the
cluster size incorporates intra-cluster quantum correlations directly
into the cluster phase space. CTWA can then be combined with either
Gaussian (`gcTWA`) or discrete (`dcTWA`) initial sampling.

For example:

```julia
model = SpinModel(
    Chain(100),
    PowerLaw(
        Ising(:x; J=1.0);
        α=3.0,
    ),
)

state = Up()

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

The approximation method can therefore be changed without modifying
either the Hamiltonian or the initial-state specification.