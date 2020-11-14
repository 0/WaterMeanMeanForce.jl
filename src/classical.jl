# Utilities for classical mechanics.

"""
    kinetic_energy(ms::AbstractVector{Float64}, P::Int, D::Int, A::Int, N::Int,
                   ps::Array{Float64,4})

Compute the classical kinetic energy for `ps` (of shape `P`, `D`, `A`, `N`)
given masses `ms` (of length `A`).
"""
function kinetic_energy(ms::AbstractVector{Float64}, P::Int, D::Int, A::Int,
                        N::Int, ps::Array{Float64,4})
    result = 0.0 # kJ / mol

    for n in 1:N
        for a in 1:A
            for d in 1:D
                for j in 1:P
                    result += ps[j, d, a, n]^2 / ms[a]
                end
            end
        end
    end

    0.5 * result
end

"""
    kinetic_energy(ms::AbstractVector{Float64}, s::State)

Compute the classical kinetic energy for `s` given masses `ms`.
"""
function kinetic_energy(ms::AbstractVector{Float64},
                        s::State{P,D,A,N}) where {P,D,A,N}
    kinetic_energy(ms, P, D, A, N, s.ps) # kJ / mol
end

"""
    kinetic_energy(ms::AbstractVector{Float64}, nms::NormalModeState)

Compute the classical kinetic energy for `nms` given masses `ms`.
"""
function kinetic_energy(ms::AbstractVector{Float64},
                        nms::NormalModeState{P,D,A,N}) where {P,D,A,N}
    kinetic_energy(ms, P, D, A, N, nms.Ps) # kJ / mol
end

"""
    temperature(ms::AbstractVector{Float64}, N_dof::Int, as::AbstractState)

Compute the instantaneous classical temperature for `as` with `N_dof` degrees
of freedom, given masses `ms`.
"""
function temperature(ms::AbstractVector{Float64}, N_dof::Int,
                     as::AbstractState)
    2.0 * kinetic_energy(ms, as) / (kB * N_dof) # K
end

"""
    V_springs(beta::Float64, s::State)

Compute the path integral spring potential energy for `s` at `beta`.
"""
function V_springs(beta::Float64, s::State{P,D,A,N}) where {P,D,A,N}
    result = 0.0 # kJ / mol

    for n in 1:N
        for a in 1:A
            k = MS[a] * P / (hbar^2 * beta^2) # kJ / nm^2 mol

            for d in 1:D
                for j in 1:P
                    j_p = bead_next(P, j)

                    result += k * (s.qs[j, d, a, n] - s.qs[j_p, d, a, n])^2
                end
            end
        end
    end

    0.5 * result
end
