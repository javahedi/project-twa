# References and Citation

This page collects the main methodological references underlying the
approximations implemented in `TWA.jl`.

## Discrete truncated-Wigner approximation

The discrete truncated-Wigner approximation for interacting spin
systems is described in:

J. Schachenmayer, A. Pikovski, and A. M. Rey,  
*Many-Body Quantum Spin Dynamics with Monte Carlo Trajectories on a
Discrete Phase Space*,  
**Physical Review X 5**, 011022 (2015).  
DOI: [`10.1103/PhysRevX.5.011022`](https://doi.org/10.1103/PhysRevX.5.011022)

This work introduced a discrete phase-space Monte Carlo formulation for
the semiclassical dynamics of interacting quantum spin systems and is a
primary reference for DTWA.

### BibTeX

```bibtex
@article{PhysRevX.5.011022,
  title = {Many-Body Quantum Spin Dynamics with Monte Carlo Trajectories on a Discrete Phase Space},
  author = {Schachenmayer, J. and Pikovski, A. and Rey, A. M.},
  journal = {Phys. Rev. X},
  volume = {5},
  issue = {1},
  pages = {011022},
  numpages = {10},
  year = {2015},
  month = {Feb},
  publisher = {American Physical Society},
  doi = {10.1103/PhysRevX.5.011022},
  url = {https://link.aps.org/doi/10.1103/PhysRevX.5.011022}
}
```

## Cluster truncated-Wigner approximation

The cluster truncated-Wigner approximation was introduced by:

J. Wurtz, A. Polkovnikov, and D. Sels,  
*Cluster truncated Wigner approximation in strongly interacting systems*,  
**Annals of Physics 395**, 341–365 (2018).  
DOI: [`10.1016/j.aop.2018.06.001`](https://doi.org/10.1016/j.aop.2018.06.001)

This work introduced the cluster truncated-Wigner approximation as a
phase-space approach for treating quantum dynamics in strongly
interacting systems. The central idea is to group degrees of freedom
into clusters, treat the dynamics within each cluster using an enlarged
phase space, and approximate the coupling between clusters
semiclassically.

### BibTeX

```bibtex
@article{WURTZ2018341,
  title = {Cluster truncated Wigner approximation in strongly interacting systems},
  journal = {Annals of Physics},
  volume = {395},
  pages = {341--365},
  year = {2018},
  issn = {0003-4916},
  doi = {10.1016/j.aop.2018.06.001},
  url = {https://www.sciencedirect.com/science/article/pii/S0003491618301647},
  author = {Jonathan Wurtz and Anatoli Polkovnikov and Dries Sels},
  keywords = {Semiclassical dynamics, Phase space, Wigner functions, Strongly interacting quantum systems, Quantum dynamics}
}
```

## Discrete cluster sampling

For the development and application of discrete sampling within the
cluster truncated-Wigner framework, see:

A. Braemer, J. Vahedi, and M. Gärttner,  
*Cluster truncated Wigner approximation for bond-disordered Heisenberg
spin models*,  
**Physical Review B 110**, 054204 (2024).  
DOI: [`10.1103/PhysRevB.110.054204`](https://doi.org/10.1103/PhysRevB.110.054204)

This work applies CTWA to interacting disordered spin systems and
develops a discrete sampling scheme for the initial cluster phase space
as an alternative to Gaussian sampling. In the terminology used
throughout the `TWA.jl` documentation, these two sampling choices are
referred to as gcTWA and dcTWA.

### BibTeX

```bibtex
@article{PhysRevB.110.054204,
  title = {Cluster truncated Wigner approximation for bond-disordered Heisenberg spin models},
  author = {Braemer, Adrian and Vahedi, Javad and G\"arttner, Martin},
  journal = {Phys. Rev. B},
  volume = {110},
  issue = {5},
  pages = {054204},
  numpages = {14},
  year = {2024},
  month = {Aug},
  publisher = {American Physical Society},
  doi = {10.1103/PhysRevB.110.054204},
  url = {https://link.aps.org/doi/10.1103/PhysRevB.110.054204}
}
```

## Discrete cluster sampling

For the development and application of discrete sampling within the
cluster truncated-Wigner framework, see:

A. Braemer, J. Vahedi, and M. Gärttner,  
*Cluster truncated Wigner approximation for bond-disordered Heisenberg
spin models*,  
**Physical Review B 110**, 054204 (2024).  
DOI: [`10.1103/PhysRevB.110.054204`](https://doi.org/10.1103/PhysRevB.110.054204)

This work applies CTWA to interacting disordered spin systems and
develops a discrete sampling scheme for the initial cluster phase space
as an alternative to Gaussian sampling. In the terminology used
throughout the `TWA.jl` documentation, these two sampling choices are
referred to as gcTWA and dcTWA.

### BibTeX

```bibtex
@article{PhysRevB.110.054204,
  title = {Cluster truncated Wigner approximation for bond-disordered Heisenberg spin models},
  author = {Braemer, Adrian and Vahedi, Javad and G\"arttner, Martin},
  journal = {Phys. Rev. B},
  volume = {110},
  issue = {5},
  pages = {054204},
  numpages = {14},
  year = {2024},
  month = {Aug},
  publisher = {American Physical Society},
  doi = {10.1103/PhysRevB.110.054204},
  url = {https://link.aps.org/doi/10.1103/PhysRevB.110.054204}
}
```

## Which reference should I cite?

If you use **DTWA**, please cite the discrete phase-space work by
Schachenmayer, Pikovski, and Rey (2015).

If you use **CTWA**, please cite the original CTWA work by Wurtz,
Polkovnikov, and Sels (2018).

If you use **discrete cluster sampling (dcTWA)**, please additionally
cite Braemer, Vahedi, and Gärttner (2024), which develops the discrete
sampling formulation used here.

When appropriate, please also cite `TWA.jl` itself once a permanent
software citation is available.


## Method overview

The relationship between the terminology used in the documentation and
the corresponding references is summarized below.

| Method | Phase space | Sampling | Primary methodological reference |
|---|---|---|---|
| TWA | single spin | Gaussian | traditional truncated-Wigner framework |
| DTWA | single spin | discrete | Schachenmayer, Pikovski & Rey (2015) |
| gcTWA | cluster | Gaussian | Wurtz, Polkovnikov & Sels (2018) |
| dcTWA | cluster | discrete | Wurtz, Polkovnikov & Sels (2018); Braemer, Vahedi & Gärttner (2024) |


In the `TWA.jl` API, gcTWA and dcTWA are not separate approximation
types. They are selected by choosing the CTWA sampling strategy:

```julia
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

This reflects the distinction between the cluster dynamical
approximation and the representation used to sample the initial quantum
state.

## Further reading

See the following documentation pages for the implementation-facing
description of these methods:

- [Discrete Truncated-Wigner Approximation](manual/dtwa.md)
- [Cluster Truncated-Wigner Approximation](manual/ctwa.md)
- [Long-range Ising Benchmark](examples/long_range_ising.md)



