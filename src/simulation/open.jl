# Keep the ODE solver options on the coherent path. The dissipative path
# directly evaluates the legacy Itô Euler–Maruyama update.
function _validate_open_options(rates, solver, kwargs)
    iszero(rates) && return nothing
    isnothing(solver) ||
        throw(ArgumentError("local dissipation uses Euler–Maruyama; omit solver"))
    isempty(kwargs) ||
        throw(ArgumentError("unsupported Euler–Maruyama keywords: $(keys(kwargs))"))
    return nothing
end

function _open_trajectory!(output, u, rhs!, compiled, rates, t0, times, dt, rng)
    T = eltype(u)
    step = T(dt)
    isfinite(step) && step > 0 ||
        throw(ArgumentError("dt must be finite and positive"))

    next = similar(u)
    t = t0
    for (index, target) in enumerate(times)
        while t < target
            new_time = min(t + step, target)
            new_time > t || throw(ArgumentError("dt is too small to advance time"))
            h = new_time - t
            rhs!(next, u, compiled, t)
            @. next = u + h * next
            _dissipative_step!(next, u, compiled, rates, h, rng)
            u, next = next, u
            t = new_time
        end
        output[:, :, index] .= u
    end
    return nothing
end
