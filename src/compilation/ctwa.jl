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

Compile a model, product initial state, and CTWA approximation into the sparse
representation used by CTWA trajectories.

The current clustering implementation requires equal-size contiguous clusters,
so the number of physical sites must be divisible by `method.cluster_size`.
"""
function compile(
    model::SpinModel,
    state::AbstractState,
    method::CTWA;
    T::Type{<:AbstractFloat}=Float64,
)
    clustering =
        Clustering(
            geometry(model),
            method.cluster_size,
        )

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
