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
