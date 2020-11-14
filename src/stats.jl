# Statistical utilities.

# Ambegaokar, V., & Troyer, M. (2010). Estimating errors reliably in Monte
# Carlo simulations of the Ehrenfest model. American Journal of Physics, 78(2),
# 150-157.

"""
    stderr_bin(dt::Float64, xs::AbstractVector)

Compute the standard error and autocorrelation time of `xs` by binning.

Stops when there are too few bins left. Returns the largest error attained and
the corresponding autocorrelation time in units of `dt`.
"""
function stderr_bin(dt::Float64, xs::AbstractVector)
    stderr_full = std(xs) / sqrt(length(xs))
    stderr_max = stderr_full

    while length(xs) >= 60
        xs = [0.5 * (xs[i] + xs[i+1]) for i in (1+length(xs)%2):2:length(xs)]
        stderr_max = max(stderr_max, std(xs) / sqrt(length(xs)))
    end

    autocorrelation_steps = 0.5 * ((stderr_max / stderr_full)^2 - 1)

    stderr_max, autocorrelation_steps * dt
end

"""
    aggregate(dt::Float64, xs::AbstractVector)

Perform a statistical analysis of `xs` for time step `dt`.
"""
function aggregate(dt::Float64, xs::AbstractVector)
    mean(xs), stderr_bin(dt, xs)...
end
