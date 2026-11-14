# Keep the ODE solver options on the coherent path. The dissipative path
# integrates the Hamiltonian and local Itô terms with a fixed maximum step.
function _validate_open_options(rates, solver, kwargs)
    iszero(rates) && return nothing
    isnothing(solver) ||
        throw(ArgumentError("local dissipation uses the built-in integrator; omit solver"))
    isempty(kwargs) ||
        throw(ArgumentError("unsupported open-system solver keywords: $(keys(kwargs))"))
    return nothing
end

function _open_trajectory!(output, u, rhs!, compiled, rates, t0, times, dt, rng)
    T = eltype(u)
    step = T(dt)
    isfinite(step) && step > 0 ||
        throw(ArgumentError("dt must be finite and positive"))

    next = similar(u)
    stage = similar(u)
    derivative = similar(u)
    t = t0
    for (index, target) in enumerate(times)
        while t < target
            new_time = min(t + step, target)
            new_time > t || throw(ArgumentError("dt is too small to advance time"))
            h = new_time - t
            _coherent_step!(next, stage, derivative, u, rhs!, compiled, t, h)
            _dissipative_step!(next, u, compiled, rates, h, rng)
            u, next = next, u
            t = new_time
        end
        all(isfinite, u) ||
            throw(ErrorException("open-system trajectory became nonfinite at t=$target; reduce dt"))
        output[:, :, index] .= u
    end
    return nothing
end

# Integrate the Hamiltonian part with RK4. Euler updates of the precession
# term can produce runaway spin lengths on long open-DTWA runs.
function _coherent_step!(next, stage, derivative, u, rhs!, compiled, t, h)
    half_step = h / 2
    rhs!(derivative, u, compiled, t)
    @. next = u + (h / 6) * derivative
    @. stage = u + half_step * derivative

    rhs!(derivative, stage, compiled, t + half_step)
    @. next += (h / 3) * derivative
    @. stage = u + half_step * derivative

    rhs!(derivative, stage, compiled, t + half_step)
    @. next += (h / 3) * derivative
    @. stage = u + h * derivative

    rhs!(derivative, stage, compiled, t + h)
    @. next += (h / 6) * derivative
    return nothing
end
