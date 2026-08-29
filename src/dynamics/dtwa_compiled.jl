"""
    CompiledDTWA

Precomputed closed-system data used by the optimized DTWA dynamics kernel.

The compiled representation stores:

- the total uniform Cartesian field;
- only nonzero pair couplings;
- the number of lattice sites.

Each pair is stored as a compact immutable tuple

    (i, j, Jx, Jy, Jz)

so the RHS can iterate directly over interacting pairs without repeatedly
querying the high-level Hamiltonian representation.
"""
struct CompiledDTWA{T<:Real,P<:Tuple}
    nsites::Int
    field::NTuple{3,T}
    pairs::P
end


"""
    compile(model, method::DTWA)

Compile a `SpinModel` into data specialized for closed-system DTWA
dynamics.

This step interprets the model once and extracts the nonzero couplings
needed by the time-evolution kernel.

The `DTWA` method is accepted so the compilation API follows the package
architecture

    compile(model, method)

and can later be extended independently for CTWA and open-system
environments.
"""
function compile(
    model::SpinModel,
    ::DTWA,
)
    nsites = model.geometry.nsites

    field = _promoted_model_field_components(model)
    pairs = _compile_dtwa_pairs(model)

    T = _compiled_dtwa_numeric_type(
        field,
        pairs,
    )

    typed_field = (
        T(field[1]),
        T(field[2]),
        T(field[3]),
    )

    typed_pairs = map(pairs) do pair
        (
            pair[1],
            pair[2],
            T(pair[3]),
            T(pair[4]),
            T(pair[5]),
        )
    end

    return CompiledDTWA(
        nsites,
        typed_field,
        typed_pairs,
    )
end


"""
    dtwa_rhs!(du, u, compiled::CompiledDTWA, t)

Evaluate the closed-system DTWA equations using precompiled model data.

The layout and Pauli-variable convention are identical to the reference
model-based `dtwa_rhs!` method.
"""
function dtwa_rhs!(
    du::AbstractMatrix,
    u::AbstractMatrix,
    compiled::CompiledDTWA,
    t,
)
    nsites = compiled.nsites

    _check_dtwa_dynamics_layout(
        du,
        u,
        nsites,
    )

    fill!(du, zero(eltype(du)))

    hx, hy, hz = compiled.field

    @inbounds for i in 1:nsites
        sx = u[i, 1]
        sy = u[i, 2]
        sz = u[i, 3]

        du[i, 1] += 2 * (hy * sz - hz * sy)
        du[i, 2] += 2 * (hz * sx - hx * sz)
        du[i, 3] += 2 * (hx * sy - hy * sx)
    end

    @inbounds for pair in compiled.pairs
        i, j, Jx, Jy, Jz = pair

        six = u[i, 1]
        siy = u[i, 2]
        siz = u[i, 3]

        sjx = u[j, 1]
        sjy = u[j, 2]
        sjz = u[j, 3]

        hix = Jx * sjx
        hiy = Jy * sjy
        hiz = Jz * sjz

        du[i, 1] += 2 * (hiy * siz - hiz * siy)
        du[i, 2] += 2 * (hiz * six - hix * siz)
        du[i, 3] += 2 * (hix * siy - hiy * six)

        hjx = Jx * six
        hjy = Jy * siy
        hjz = Jz * siz

        du[j, 1] += 2 * (hjy * sjz - hjz * sjy)
        du[j, 2] += 2 * (hjz * sjx - hjx * sjz)
        du[j, 3] += 2 * (hjx * sjy - hjy * sjx)
    end

    return nothing
end


"""
    _compile_dtwa_pairs(model)

Return an immutable tuple containing only nonzero pair couplings.
"""
function _compile_dtwa_pairs(
    model::SpinModel,
)
    nsites = model.geometry.nsites

    pairs = Tuple[]

    for i in 1:(nsites - 1)
        for j in (i + 1):nsites
            Jx, Jy, Jz = coupling_components(
                model,
                i,
                j,
            )

            iszero(Jx) && iszero(Jy) && iszero(Jz) &&
                continue

            push!(
                pairs,
                (i, j, Jx, Jy, Jz),
            )
        end
    end

    return Tuple(pairs)
end


"""
    _promoted_model_field_components(model)

Return the total field with all components promoted to one common type.
"""
@inline function _promoted_model_field_components(
    model::SpinModel,
)
    hx, hy, hz = _model_field_components(model)
    hx, hy, hz = promote(hx, hy, hz)

    return (hx, hy, hz)
end


"""
    _compiled_dtwa_numeric_type(field, pairs)

Determine a common numeric type for all compiled field and pair data.
"""
function _compiled_dtwa_numeric_type(
    field::Tuple,
    pairs::Tuple,
)
    T = promote_type(
        typeof(field[1]),
        typeof(field[2]),
        typeof(field[3]),
    )

    for pair in pairs
        T = promote_type(
            T,
            typeof(pair[3]),
            typeof(pair[4]),
            typeof(pair[5]),
        )
    end

    return T
end
