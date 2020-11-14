# Weighted histogram analysis method.

# Kumar, S., Rosenberg, J. M., Bouzida, D., Swendsen, R. H., & Kollman, P. A.
# (1992). The weighted histogram analysis method for free‐energy calculations
# on biomolecules. I. The method. Journal of Computational Chemistry, 13(8),
# 1011-1021.

function unbias_histograms(T::Float64, bias_force_constant::Float64,
                           bin_size::Float64,
                           trajectories::Dict{Float64,Vector{Float64}},
                           num_iter::Int; aux_output::IO=stderr)
    # Most extreme positions of restrained coordinate.
    min_xi = Inf
    max_xi = -Inf

    for trajectory in values(trajectories)
        min_xi = min(min_xi, minimum(trajectory))
        max_xi = max(max_xi, maximum(trajectory))
    end

    # Round to one decimal place (0.1 nm).
    min_xi = round(min_xi, RoundDown; digits=1)
    max_xi = round(max_xi, RoundUp; digits=1)

    bins = min_xi:bin_size:max_xi
    num_edges = length(bins)+1
    edges = range(min_xi - bin_size / 2, max_xi + bin_size / 2;
                  length=num_edges)

    total = zeros(Float64, size(bins))
    windows = Float64[]
    sizes = Float64[]

    for (window, trajectory) in sort(collect(trajectories))
        _, autocorrelation_steps = stderr_bin(1.0, trajectory)
        inefficiency = 1 + 2 * autocorrelation_steps
        hist = fit(Histogram, trajectory, edges).weights / inefficiency
        total .+= hist
        push!(windows, window)
        push!(sizes, sum(hist))
    end

    offsets = bins .- windows'
    exppot = exp.(-beta(T) .* (0.5 .* bias_force_constant .* offsets.^2))

    if size(exppot) != (length(total), length(sizes))
        error("Mismatched dimensions.")
    end

    # Apply WHAM equations iteratively.
    prob = nothing
    weights = ones(Float64, size(sizes))

    meter = Progress(num_iter, output=aux_output)

    for _ in ProgressWrapper(1:num_iter, meter)
        factors = dropdims(sum(sizes' .* exppot ./ weights'; dims=2); dims=2)
        prob = total ./ factors
        # Normalize.
        prob ./= sum(prob)
        weights = dropdims(sum(prob .* exppot; dims=1); dims=1)
    end

    pmf = -log.(prob) ./ beta(T)

    bins, pmf
end
