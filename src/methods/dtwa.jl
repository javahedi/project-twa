"""
    DTWA(; trajectories=1000)

Discrete truncated-Wigner approximation.

# Keyword arguments

- `trajectories`: number of stochastic initial phase-space trajectories
  used for ensemble averaging. Must be positive.

# Example

    method = DTWA(trajectories=1000)

The method object stores approximation-specific configuration only.
Random-number generators, floating-point precision, solver settings, and
open-system environments are configured elsewhere.
"""
struct DTWA <: AbstractApproximation
    trajectories::Int

    function DTWA(; trajectories::Integer=1000)
        trajectories > 0 ||
            throw(ArgumentError(
                "trajectories must be positive; got $trajectories",
            ))

        new(Int(trajectories))
    end
end
