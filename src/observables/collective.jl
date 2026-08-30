"""
    SpinSecondMoment(axis)

Normalized collective-spin second moment

    ⟨(Sᵅ)²⟩ / N²,

with

    Sᵅ = Σᵢ σᵢᵅ.

For DTWA results this is evaluated using the operator identity

    (σᵢᵅ)² = I

for on-site terms and ensemble-averaged phase-space products for distinct
sites:

    ⟨(Sᵅ)²⟩
    = N + 2 Σᵢ<ⱼ ⟨σᵢᵅ σⱼᵅ⟩.

The explicit on-site identity is important: squaring the classical collective
coordinate trajectory by trajectory would not, in general, be the correct
Weyl-symbol prescription.
"""
struct SpinSecondMoment <: AbstractObservable
    axis::Symbol

    function SpinSecondMoment(axis::Symbol)
        _check_observable_axis(axis)
        new(axis)
    end
end


"""
    SpinVariance(axis)

Normalized collective-spin variance

    (ΔSᵅ)² / N²
    = [⟨(Sᵅ)²⟩ - ⟨Sᵅ⟩²] / N².

Because `Magnetization(axis)` already returns `⟨Sᵅ⟩ / N`, this observable is
computed as

    expectation(result, SpinSecondMoment(axis))
    - expectation(result, Magnetization(axis)).^2.
"""
struct SpinVariance <: AbstractObservable
    axis::Symbol

    function SpinVariance(axis::Symbol)
        _check_observable_axis(axis)
        new(axis)
    end
end


function expectation(
    result::DTWAResult,
    observable::SpinSecondMoment,
)
    axis = _axis_index(observable.axis)

    nsites = size(result.trajectories, 1)
    ntimes = length(result.t)
    ntraj = ntrajectories(result)

    Tout = promote_type(
        eltype(result.trajectories),
        Float64,
    )

    values = Vector{Tout}(undef, ntimes)

    normalization = nsites * nsites

    @inbounds for time_index in 1:ntimes
        # On-site contribution:
        #
        #     Σᵢ ⟨(σᵢᵅ)²⟩ = N.
        #
        total = Tout(nsites)

        # Distinct-site contribution:
        #
        #     2 Σᵢ<ⱼ ⟨σᵢᵅ σⱼᵅ⟩.
        #
        for i in 1:(nsites - 1)
            for j in (i + 1):nsites
                pair_total = zero(Tout)

                for trajectory_index in 1:ntraj
                    pair_total += (
                        result.trajectories[
                            i,
                            axis,
                            time_index,
                            trajectory_index,
                        ] *
                        result.trajectories[
                            j,
                            axis,
                            time_index,
                            trajectory_index,
                        ]
                    )
                end

                total += 2 * pair_total / ntraj
            end
        end

        values[time_index] =
            total / normalization
    end

    return values
end


function expectation(
    result::DTWAResult,
    observable::SpinVariance,
)
    second_moment = expectation(
        result,
        SpinSecondMoment(observable.axis),
    )

    magnetization = expectation(
        result,
        Magnetization(observable.axis),
    )

    return second_moment .- magnetization .* magnetization
end
"""
Collective observables for `CTWAResult`.

These methods deliberately build on the already validated primitive
`expectation(result, Correlation(...))` and
`expectation(result, Magnetization(...))` paths.

With the package's Pauli convention

    S_a = sum_i sigma_i^a,

the normalized second moment is

    <S_a^2>/N^2
      = [N + 2 sum_{i<j} <sigma_i^a sigma_j^a>] / N^2,

because (sigma_i^a)^2 = I exactly.

The normalized variance is

    (Delta S_a)^2 / N^2
      = <S_a^2>/N^2 - (<S_a>/N)^2.

Using physical correlations here is important for CTWA because
same-cluster pairs must be represented by their joint Pauli-string
coordinate, not by multiplying two one-site coordinates.
"""


"""
    expectation(result::CTWAResult, observable::SpinSecondMoment)

Return the normalized collective second moment

    <S_a^2> / N^2

as a function of time.
"""
function expectation(
    result::CTWAResult,
    observable::SpinSecondMoment,
)
    axis =
        observable.axis

    nsites =
        result.clustering.nsites

    ntimes =
        length(
            result.t,
        )

    T =
        eltype(
            result.trajectories,
        )

    values =
        fill(
            T(nsites),
            ntimes,
        )

    for site_i in 1:(nsites - 1)
        for site_j in (site_i + 1):nsites
            values .+=
                T(2) .*
                expectation(
                    result,
                    Correlation(
                        axis,
                        site_i,
                        axis,
                        site_j,
                    ),
                )
        end
    end

    values ./=
        T(nsites)^2

    return values
end


"""
    expectation(result::CTWAResult, observable::SpinVariance)

Return the normalized collective variance

    (Delta S_a)^2 / N^2

as a function of time.
"""
function expectation(
    result::CTWAResult,
    observable::SpinVariance,
)
    axis =
        observable.axis

    second_moment =
        expectation(
            result,
            SpinSecondMoment(
                axis,
            ),
        )

    mean =
        expectation(
            result,
            Magnetization(
                axis,
            ),
        )

    return second_moment .-
           mean .^ 2
end
