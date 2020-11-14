# q-SPC/Fw water model.

# Paesani, F., Zhang, W., Case, D. A., Cheatham III, T. E., & Voth, G. A.
# (2006). An accurate and simple quantum model for liquid water. The Journal of
# Chemical Physics, 125(18), 184507.

"""
    QSPCFW

Implementation of the q-SPC/Fw flexible water model.
"""
struct QSPCFW <: WaterModel
    "Potential terms (returning kJ / mol)."
    potential::TermContainer
    "Force terms (returning kJ / mol nm)."
    force::TermContainer
end

function QSPCFW(N::Int)
    build_model(QSPCFW,
                qspcfw_potential_intra, qspcfw_potential_inter,
                qspcfw_force_intra, qspcfw_force_inter,
                N)
end


const qspcfw_k_b = 1059.162 * kcal_to_kJ / A_to_nm^2 # kJ / mol nm^2
const qspcfw_r_OH_eq = 1.0 * A_to_nm # nm
const qspcfw_k_a = 75.90 * kcal_to_kJ # kJ / mol
const qspcfw_theta_HOH_eq = 112.0 * deg_to_rad

const qspcfw_epsilon_OO = 0.1554252 * kcal_to_kJ # kJ / mol
const qspcfw_sigma_OO = 3.165492 * A_to_nm # nm

const qspcfw_epsilon_0 = (8.8541878128e-12 / e_to_C^2 / N_A / to_kilo
                          / to_nano) # e^2 mol / kJ nm
const qspcfw_q_O = -0.84 # e
const qspcfw_q_H = 0.42 # e
const qspcfw_qs = Vector{Float64}(undef, length(MS))
qspcfw_qs[AO] = qspcfw_q_O
qspcfw_qs[AH1] = qspcfw_q_H
qspcfw_qs[AH2] = qspcfw_q_H

function qspcfw_potential_intra(n::Int)
    function term_f(j::Int, s::State)
        rvec_OH_1 = s.qs[j, :, AO, n] .- s.qs[j, :, AH1, n] # nm
        rvec_OH_2 = s.qs[j, :, AO, n] .- s.qs[j, :, AH2, n] # nm
        r_OH_1 = norm(rvec_OH_1) # nm
        r_OH_2 = norm(rvec_OH_2) # nm
        theta_HOH = acos_(dot(rvec_OH_1 ./ r_OH_1, rvec_OH_2 ./ r_OH_2))

        result = 0.0 # kJ / mol
        result += 0.5 * qspcfw_k_b * (r_OH_1 - qspcfw_r_OH_eq)^2
        result += 0.5 * qspcfw_k_b * (r_OH_2 - qspcfw_r_OH_eq)^2
        result += 0.5 * qspcfw_k_a * (theta_HOH - qspcfw_theta_HOH_eq)^2

        result
    end
end

function qspcfw_elec_pot(n1::Int, a1::Int, n2::Int, a2::Int, j::Int, s::State)
    R = norm(s.qs[j, :, a1, n1] .- s.qs[j, :, a2, n2]) # nm
    qspcfw_qs[a1] * qspcfw_qs[a2] / (4pi * qspcfw_epsilon_0 * R) # kJ / mol
end

function qspcfw_potential_inter(n1::Int, n2::Int)
    function term_f(j::Int, s::State)
        result = 0.0 # kJ / mol

        # O-O (LJ)
        R = norm(s.qs[j, :, AO, n1] .- s.qs[j, :, AO, n2]) # nm
        diff = (qspcfw_sigma_OO / R)^12 - (qspcfw_sigma_OO / R)^6
        result += 4 * qspcfw_epsilon_OO * diff

        # O-O (electrostatic)
        result += qspcfw_elec_pot(n1, AO, n2, AO, j, s)

        # O-H (electrostatic)
        result += qspcfw_elec_pot(n1, AO, n2, AH1, j, s)
        result += qspcfw_elec_pot(n1, AO, n2, AH2, j, s)
        result += qspcfw_elec_pot(n2, AO, n1, AH1, j, s)
        result += qspcfw_elec_pot(n2, AO, n1, AH2, j, s)

        # H-H (electrostatic)
        result += qspcfw_elec_pot(n1, AH1, n2, AH1, j, s)
        result += qspcfw_elec_pot(n1, AH1, n2, AH2, j, s)
        result += qspcfw_elec_pot(n1, AH2, n2, AH1, j, s)
        result += qspcfw_elec_pot(n1, AH2, n2, AH2, j, s)

        result
    end
end

function qspcfw_force_intra(n::Int)
    function term_f!(result::ForceContainer, j::Int, s::State)
        rvec_OH_1 = s.qs[j, :, AO, n] .- s.qs[j, :, AH1, n] # nm
        rvec_OH_2 = s.qs[j, :, AO, n] .- s.qs[j, :, AH2, n] # nm
        r_OH_1 = norm(rvec_OH_1) # nm
        r_OH_2 = norm(rvec_OH_2) # nm
        theta_HOH = acos_(dot(rvec_OH_1 ./ r_OH_1, rvec_OH_2 ./ r_OH_2))

        # kJ / nm mol
        F_OH_1 = -qspcfw_k_b * (1.0 - qspcfw_r_OH_eq / r_OH_1) .* rvec_OH_1
        result.f[j, :, AO, n] .+= F_OH_1
        result.f[j, :, AH1, n] .-= F_OH_1

        # kJ / nm mol
        F_OH_2 = -qspcfw_k_b * (1.0 - qspcfw_r_OH_eq / r_OH_2) .* rvec_OH_2
        result.f[j, :, AO, n] .+= F_OH_2
        result.f[j, :, AH2, n] .-= F_OH_2

        x = (dot(rvec_OH_1, rvec_OH_1) * dot(rvec_OH_2, rvec_OH_2)
             - dot(rvec_OH_1, rvec_OH_2)^2) # nm^4
        # kJ / nm^2 mol
        pre_F_HOH = qspcfw_k_a * (theta_HOH - qspcfw_theta_HOH_eq) / sqrt(x)
        x_HOH_1 = (dot(rvec_OH_1, rvec_OH_2) / dot(rvec_OH_1, rvec_OH_1)
                       .* rvec_OH_1
                   .- rvec_OH_2) # nm
        x_HOH_2 = (dot(rvec_OH_1, rvec_OH_2) / dot(rvec_OH_2, rvec_OH_2)
                       .* rvec_OH_2
                   .- rvec_OH_1) # nm
        @. result.f[j, :, AO, n] -= pre_F_HOH * (x_HOH_1 + x_HOH_2)
        @. result.f[j, :, AH1, n] += pre_F_HOH * x_HOH_1
        @. result.f[j, :, AH2, n] += pre_F_HOH * x_HOH_2

        nothing
    end
end

function qspcfw_elec_force!(result::ForceContainer, n1::Int, a1::Int,
                            n2::Int, a2::Int, j::Int, s::State)
    Rvec = s.qs[j, :, a1, n1] .- s.qs[j, :, a2, n2] # nm
    R = norm(Rvec) # nm
    # kJ / nm mol
    F = qspcfw_qs[a1] * qspcfw_qs[a2] / (4pi * qspcfw_epsilon_0 * R^3) .* Rvec
    result.f[j, :, a1, n1] .+= F
    result.f[j, :, a2, n2] .-= F

    nothing
end

function qspcfw_force_inter(n1, n2)
    function term_f!(result::ForceContainer, j::Int, s::State)
        # O-O (LJ)
        Rvec = s.qs[j, :, AO, n1] .- s.qs[j, :, AO, n2] # nm
        R = norm(Rvec) # nm
        # 1 / nm^2
        diff = 12 * qspcfw_sigma_OO^12/R^14 - 6 * qspcfw_sigma_OO^6/R^8
        F = 4 * qspcfw_epsilon_OO * diff .* Rvec # kJ / nm mol
        result.f[j, :, AO, n1] .+= F
        result.f[j, :, AO, n2] .-= F

        # O-O (electrostatic)
        qspcfw_elec_force!(result, n1, AO, n2, AO, j, s)

        # O-H (electrostatic)
        qspcfw_elec_force!(result, n1, AO, n2, AH1, j, s)
        qspcfw_elec_force!(result, n1, AO, n2, AH2, j, s)
        qspcfw_elec_force!(result, n2, AO, n1, AH1, j, s)
        qspcfw_elec_force!(result, n2, AO, n1, AH2, j, s)

        # H-H (electrostatic)
        qspcfw_elec_force!(result, n1, AH1, n2, AH1, j, s)
        qspcfw_elec_force!(result, n1, AH1, n2, AH2, j, s)
        qspcfw_elec_force!(result, n1, AH2, n2, AH1, j, s)
        qspcfw_elec_force!(result, n1, AH2, n2, AH2, j, s)

        nothing
    end
end


registered_models["q-SPC/Fw"] = QSPCFW
