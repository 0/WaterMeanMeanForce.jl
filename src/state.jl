# Atomic state.

abstract type AbstractState end


"""
    State{P,D,A,N}

State of `N` water molecules with `A` atoms each in `D` dimensions copied over
`P` beads, in Cartesian coordinates.
"""
struct State{P,D,A,N} <: AbstractState
    "Absolute spatial positions (nm)."
    qs::Array{Float64,4}
    "Conjugate momenta (g nm / mol ps)."
    ps::Array{Float64,4}
end

function State(P::Int, D::Int, A::Int, N::Int; R::Maybe{Float64}=nothing,
               qs_init::Maybe{Vector{Float64}}=nothing)
    N >= max(CN1, CN2) || error("Too few molecules.")
    A == length(MS) || error("Water: H-O-H.")
    D >= 2 || error("At least 2 dimensions.")
    P >= 2 || error("At least 2 beads.")

    if isnothing(R) == isnothing(qs_init)
        error("Expected either R or qs_init.")
    end

    qs = zeros(Float64, P, D, A, N)
    ps = zeros(Float64, P, D, A, N)

    if !isnothing(qs_init)
        qs .= reshape(qs_init, P, D, A, N)
    else
        for n in 1:N
            # Evenly spaced in one direction.
            offset = (n-1) * R # nm

            for j in 1:P
                qs[j, 1, AO, n] += offset
                qs[j, 1, AH1, n] += offset
                qs[j, 1, AH2, n] += offset

                # Roughly watery shape.
                qs[j, 2, AH1, n] += 0.1
                qs[j, 1, AH2, n] += 0.1 * sin(112.0 * deg_to_rad)
                qs[j, 2, AH2, n] += 0.1 * cos(112.0 * deg_to_rad)
            end
        end
    end

    State{P,D,A,N}(qs, ps)
end


"""
    NormalModeState{P,D,A,N}

State of `N` water molecules with `A` atoms each in `D` dimensions copied over
`P` beads, in path normal mode coordinates.
"""
struct NormalModeState{P,D,A,N} <: AbstractState
    "Normal mode positions (nm)."
    Qs::Array{Float64,4}
    "Normal mode momenta (g nm / mol ps)."
    Ps::Array{Float64,4}
end

function NormalModeState(s::State{P,D,A,N}) where {P,D,A,N}
    Qs = zeros(Float64, P, D, A, N)
    Ps = zeros(Float64, P, D, A, N)

    NormalModeState{P,D,A,N}(Qs, Ps)
end
