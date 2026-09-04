"""
    CompiledCTWA

Compiled closed-system CTWA problem data.

This object owns the setup-time structures needed repeatedly during trajectory
sampling and evolution:

* sparse cluster Hamiltonian;
* sparse cached commutator algebra;
* matrix-free Gaussian initial-state sampling plan.

It deliberately does not own a solver, RNG, trajectory result array, or open
system environment.
"""
struct CompiledCTWA{
    H<:CompiledCTWAHamiltonian,
    A<:CompiledCTWAAlgebra,
    S,
}
    hamiltonian::H
    algebra::A
    sampling::S
end


"""
    compile(model, state, method::CTWA; T=Float64)

CTWA currently requires equal-size clusters. When no explicit clustering is
provided, clusters are formed from the canonical site ordering of the model
geometry. An explicit `Clustering` may be supplied through `CTWA` for
geometry-aware spatial cluster layouts.

"""
function compile(
    model::SpinModel,
    state::AbstractState,
    method::CTWA;
    T::Type{<:AbstractFloat}=Float64,
)
    clustering =
        if method.clustering === nothing
            Clustering(
                geometry(model),
                method.cluster_size,
            )
        else
            method.clustering
        end

    clustering.nsites == nsites(geometry(model)) ||
        throw(ArgumentError(
            "explicit clustering contains $(clustering.nsites) sites, " *
            "but model geometry contains $(nsites(geometry(model))) sites",
        ))

    hamiltonian =
        compile_ctwa_hamiltonian(
            model,
            clustering;
            T=T,
        )

    algebra =
        compile_ctwa_algebra(
            hamiltonian,
        )

    sampling =
        compile_ctwa_sampling(
            state,
            clustering,
            method.sampling;
            T=T,
        )

    return CompiledCTWA(
        hamiltonian,
        algebra,
        sampling,
    )
end


"""
    clustering(compiled::CompiledCTWA)

Return the CTWA clustering used by the compiled problem.
"""
clustering(
    compiled::CompiledCTWA,
) =
    compiled.hamiltonian.clustering


"""
    basis(compiled::CompiledCTWA)

Return the canonical cluster Pauli-string basis.
"""
basis(
    compiled::CompiledCTWA,
) =
    compiled.hamiltonian.basis


"""
    sample_ctwa(compiled::CompiledCTWA; rng=nothing)

Draw one initial CTWA phase-space state from the compiled sampling plan.
"""
sample_ctwa(
    compiled::CompiledCTWA;
    rng=nothing,
) =
    sample_ctwa(
        compiled.sampling;
        rng=rng,
    )


"""
    sample_ctwa!(x, compiled::CompiledCTWA; rng=nothing)

Fill a preallocated CTWA phase-space state from the compiled sampling plan.
"""
sample_ctwa!(
    x::AbstractMatrix,
    compiled::CompiledCTWA;
    rng=nothing,
) =
    sample_ctwa!(
        x,
        compiled.sampling;
        rng=rng,
    )


"""
    ctwa_rhs!(dx, x, compiled::CompiledCTWA, t)

Production CTWA deterministic RHS using the cached sparse commutator algebra.
"""
ctwa_rhs!(
    dx::AbstractMatrix,
    x::AbstractMatrix,
    compiled::CompiledCTWA,
    t,
) =
    ctwa_rhs!(
        dx,
        x,
        compiled.hamiltonian,
        compiled.algebra,
        t,
    )
