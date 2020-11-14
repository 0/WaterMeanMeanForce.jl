# Harmonic "water" model.

"""
    HarmonicWater

Purely harmonic "water" model.
"""
struct HarmonicWater <: WaterModel
    "Potential terms (returning kJ / mol)."
    potential::TermContainer
    "Force terms (returning kJ / nm mol)."
    force::TermContainer
end

function HarmonicWater(N::Int)
    build_model(HarmonicWater,
                harmonic_potential_intra, harmonic_potential_inter,
                harmonic_force_intra, harmonic_force_inter,
                N)
end


const harmonic_k_intra = 1.0e5 # kJ / nm^2 mol
const harmonic_k_inter = 1.0e2 # kJ / nm^2 mol

function harmonic_potential_intra(n::Int)
    function term_f(j::Int, s::State{P,D,A,N}) where {P,D,A,N}
        result = 0.0 # kJ / mol

        for a1 in 1:A
            for a2 in (a1+1):A
                R = norm(s.qs[j, :, a1, n] .- s.qs[j, :, a2, n]) # nm
                result += harmonic_k_intra * R^2
            end
        end

        0.5 * result
    end
end

function harmonic_potential_inter(n1::Int, n2::Int)
    function term_f(j::Int, s::State{P,D,A,N}) where {P,D,A,N}
        result = 0.0 # kJ / mol

        for a1 in 1:A
            for a2 in 1:A
                R = norm(s.qs[j, :, a1, n1] .- s.qs[j, :, a2, n2]) # nm
                result += harmonic_k_inter * R^2
            end
        end

        0.5 * result
    end
end

function harmonic_force_intra(n::Int)
    function term_f!(result::ForceContainer, j::Int,
                     s::State{P,D,A,N}) where {P,D,A,N}
        Rvec = Array{Float64}(undef, 3) # nm
        F = Array{Float64}(undef, 3) # kJ / nm mol

        for a1 in 1:A
            for a2 in (a1+1):A
                @. Rvec = s.qs[j, :, a1, n] - s.qs[j, :, a2, n]
                @. F = -harmonic_k_intra * Rvec
                @. result.f[j, :, a1, n] += F
                @. result.f[j, :, a2, n] -= F
            end
        end

        nothing
    end
end

function harmonic_force_inter(n1, n2)
    function term_f!(result::ForceContainer, j::Int,
                     s::State{P,D,A,N}) where {P,D,A,N}
        Rvec = Array{Float64}(undef, 3) # nm
        F = Array{Float64}(undef, 3) # kJ / nm mol

        for a1 in 1:A
            for a2 in 1:A
                @. Rvec = s.qs[j, :, a1, n1] - s.qs[j, :, a2, n2]
                @. F = -harmonic_k_inter * Rvec
                @. result.f[j, :, a1, n1] += F
                @. result.f[j, :, a2, n2] -= F
            end
        end

        nothing
    end
end


registered_models["harmonic"] = HarmonicWater


function harmonic_exact_slope(N::Int, T::Float64)
    # Force constant matrix (kJ / nm^2 mol).
    K = zeros(Float64, 3N, 3N)

    idx1 = 0

    for n1 in 1:N
        for a1 in 1:3
            idx1 += 1
            idx2 = 0

            for n2 in 1:N
                for a2 in 1:3
                    idx2 += 1

                    idx1 == idx2 && continue

                    if n1 == n2
                        K[idx2, idx1] += harmonic_k_intra
                    else
                        K[idx2, idx1] += harmonic_k_inter
                    end
                end
            end
        end
    end

    A = diagm(0 => sum(K; dims=2)[:, 1]) - K # kJ / nm^2 mol
    sqrtM = diagm(0 => repeat(sqrt.(MS), N)) # g^1/2 / mol^1/2
    sqrtMinv = diagm(0 => repeat(1 ./ sqrt.(MS), N)) # mol^1/2 / g^1/2

    eigvals, eigvecs = eigen(sqrtMinv * A * sqrtMinv) # 1 / ps^2, 1
    A_norm = norm(sqrtMinv * A * sqrtMinv) # 1 / ps^2
    # Indices of zero eigenvalues.
    zidxs = Set{Int}()

    for (i, eigval) in enumerate(eigvals)
        if abs(eigval) / A_norm < 1e-10
            push!(zidxs, i)
        end
    end

    # Rank of A.
    R = 3N - length(zidxs)

    omega = zeros(Float64, R, R) # mol / g nm^2
    omega_inv = zeros(Float64, R, R) # g nm^2 / mol
    idx = 0

    for i in 1:R
        idx += 1

        while idx in zidxs
            idx += 1
        end

        freq = sqrt(eigvals[idx]) # 1 / ps
        omega[i, i] += 2 * freq * tanh(0.5 * beta(T) * hbar * freq) / hbar
        omega_inv[i, i] += 1 / omega[i, i]
    end

    # g^1/2 / mol^1/2
    G = (eigvecs' * sqrtM)[setdiff(begin:end, zidxs), :]
    C_r = zeros(Float64, R) # g^1/2 / mol^1/2
    C_Q = zeros(Float64, R, 3N-1) # g^1/2 / mol^1/2

    for i in 1:3
        @. C_r += 0.5 * G[:, i]
        @. C_r -= 0.5 * G[:, i+3]
        @. C_Q[:, 1] += G[:, i]
        @. C_Q[:, 1] += G[:, i+3]
    end

    # Total mass of molecule (g / mol).
    M_mol = sum(MS)

    @. C_Q[:, 2] += (MS[2]/(MS[1] + MS[2]) * G[:, 1]
                     - MS[1]/(MS[1] + MS[2]) * G[:, 2])
    @. C_Q[:, 3] += (MS[3]/M_mol * G[:, 1]
                     + MS[3]/M_mol * G[:, 2]
                     - (MS[1] + MS[2])/M_mol * G[:, 3])
    @. C_Q[:, 4] += (MS[2]/(MS[1] + MS[2]) * G[:, 4]
                     - MS[1]/(MS[1] + MS[2]) * G[:, 5])
    @. C_Q[:, 5] += (MS[3]/M_mol * G[:, 4]
                     + MS[3]/M_mol * G[:, 5]
                     - (MS[1] + MS[2])/M_mol * G[:, 6])

    for i in 1:(3N-6)
        @. C_Q[:, i+5] += G[:, i+6]
    end

    v = omega * C_r # mol^1/2 / g^1/2 nm^2
    # 1 / nm^2
    slope = v' * (omega_inv - C_Q * pinv(C_Q' * omega * C_Q) * C_Q') * v

    beta(T), slope
end
