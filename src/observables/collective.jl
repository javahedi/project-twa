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
