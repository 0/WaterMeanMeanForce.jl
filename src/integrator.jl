# Path integral molecular dynamics integration.

function make_Cmat(P::Int)
    Cmat = Array{Float64}(undef, P, P)
    # This works for both even and odd P.
    M = P / 2

    for k_idx in 1:P
        k = k_idx - 1

        for j in 1:P
            if k == 0
                Cmat[j, k_idx] = 1.0
            elseif 1 <= k < M
                Cmat[j, k_idx] = sqrt(2.0) * cos(2pi * j * k / P)
            elseif k == M
                Cmat[j, k_idx] = (-1.0)^j
            elseif M < k <= P-1
                Cmat[j, k_idx] = sqrt(2.0) * sin(2pi * j * k / P)
            else
                error()
            end
        end
    end

    Cmat ./ sqrt(P)
end

MS_fict(P::Int) = MS ./ P # g / mol


struct NormalModeTransformer
    "Transformation matrix."
    Cmat::Matrix{Float64}
end

function to_normal_modes!(nms::NormalModeState{P,D,A,N}, s::State{P,D,A,N},
                          nmt::NormalModeTransformer) where {P,D,A,N}
    for n in 1:N
        for a in 1:A
            for d in 1:D
                nms.Qs[:, d, a, n] .= (s.qs[:, d, a, n]' * nmt.Cmat)'
            end
        end
    end

    for n in 1:N
        for a in 1:A
            for d in 1:D
                nms.Ps[:, d, a, n] .= (s.ps[:, d, a, n]' * nmt.Cmat)'
            end
        end
    end

    nothing
end

function from_normal_modes!(s::State{P,D,A,N}, nms::NormalModeState{P,D,A,N},
                            nmt::NormalModeTransformer) where {P,D,A,N}
    for n in 1:N
        for a in 1:A
            for d in 1:D
                s.qs[:, d, a, n] .= nmt.Cmat * nms.Qs[:, d, a, n]
            end
        end
    end

    for n in 1:N
        for a in 1:A
            for d in 1:D
                s.ps[:, d, a, n] .= nmt.Cmat * nms.Ps[:, d, a, n]
            end
        end
    end

    nothing
end

function in_normal_modes(f, nms::NormalModeState{P,D,A,N}, s::State{P,D,A,N},
                         nmt::NormalModeTransformer) where {P,D,A,N}
    to_normal_modes!(nms, s, nmt)
    f()
    from_normal_modes!(s, nms, nmt)
end


struct OscillatorPropagator
    "Propagation coefficients."
    coefs::Vector{Vector{Float64}}
end

function OscillatorPropagator(dt::Float64, omegas::AbstractVector{Float64})
    coef1 = dt .* sinc.(dt .* omegas ./ pi) # ps
    coef2 = cos.(dt .* omegas)
    coef3 = cos.(dt .* omegas)
    coef4 = -omegas .* sin.(dt .* omegas) # 1 / ps

    coefs = [coef1, coef2, coef3, coef4]

    OscillatorPropagator(coefs)
end

function propagate_oscillators!(nms::NormalModeState{P,D,A,N},
                                oscprop::OscillatorPropagator) where {P,D,A,N}
    for n in 1:N
        for a in 1:A
            m = MS_fict(P)[a] # g / mol

            for d in 1:D
                Qs_p = (oscprop.coefs[1] .* nms.Ps[:, d, a, n] ./ m
                        .+ oscprop.coefs[2] .* nms.Qs[:, d, a, n])
                Ps_p = (oscprop.coefs[3] .* nms.Ps[:, d, a, n]
                        .+ oscprop.coefs[4] .* nms.Qs[:, d, a, n] .* m)

                nms.Qs[:, d, a, n] .= Qs_p
                nms.Ps[:, d, a, n] .= Ps_p
            end
        end
    end

    nothing
end


struct Constrainer
    "Constraint distance (nm)."
    R::Float64
    "Propagation coefficients rotated to normal mode basis."
    coefs_rot::Vector{Matrix{Float64}}
end

function Constrainer(R::Float64, oscprop::OscillatorPropagator,
                     Cmat::Matrix{Float64})
    coefs_rot = [Cmat * diagm(0 => coef) * Cmat' for coef in oscprop.coefs]

    Constrainer(R, coefs_rot)
end

struct ConstraintQuantities
    "Center of mass displacement (nm)."
    disp::Vector{Float64}
    "Center of mass distance (nm)."
    dist::Float64
    "Gradient of center of mass distance."
    grad::Array{Float64,3}
end

function ConstraintQuantities(s::State)
    disp = com_disp(CN1, CN2, CJ, s)
    dist = norm(disp)
    grad = com_r_grad(CN1, CN2, CJ, s)

    ConstraintQuantities(disp, dist, grad)
end

function apply_constraint!(s::State{P,D,A,N}, cstr::Constrainer,
                           cq::ConstraintQuantities) where {P,D,A,N}
    # Reduced mass between the two molecules (g / mol).
    M_red = 0.5 * sum(MS_fict(P))

    ddisp = com_disp(CN1, CN2, CJ, s) - cq.disp # nm
    x = dot(cq.disp ./ cq.dist, ddisp) # nm
    discriminant = cstr.R^2 - (norm(ddisp)^2 - x^2) # nm^2

    discriminant >= 0.0 || error("Cannot enforce constraint.")

    lambda = -(cq.dist + x - sqrt(discriminant)) * M_red # g nm / mol

    for n in 1:N
        for a in 1:A
            for d in 1:D
                # g nm / ps mol
                kp = cq.grad[d, a, n] * lambda / cstr.coefs_rot[1][1, 1]
                kq = kp / MS_fict(P)[a] # nm / ps
                s.qs[:, d, a, n] .+= kq .* cstr.coefs_rot[1][:, 1]
                s.ps[:, d, a, n] .+= kp .* cstr.coefs_rot[3][:, 1]
            end
        end
    end

    nothing
end

function apply_velocity_constraint!(s::State{P,D,A,N}) where {P,D,A,N}
    # Reduced mass between the two molecules (g / mol).
    M_red = 0.5 * sum(MS_fict(P))

    grad = com_r_grad(CN1, CN2, CJ, s)

    lambda = 0.0 # nm / ps

    for n in 1:N
        for d in 1:D
            lambda -= dot(grad[d, :, n], s.ps[CJ, d, :, n] ./ MS_fict(P))
        end
    end

    s.ps[CJ, :, :, :] .+= M_red .* lambda .* grad

    nothing
end


struct ForceApplicator
    "Water model, containing force implementations."
    model::WaterModel
    "Additional force terms."
    extra_forces::Dict{String,Vector{AdHocForce}}
    "Time step (ps)."
    dt::Float64
    "Result of the previous force evaluation."
    current_force::ForceContainer
end

function eval_pot(forceapp::ForceApplicator, s::State{P,D,A,N}) where {P,D,A,N}
    result = Dict{String,Float64}()

    result["model"] = eval_pot(forceapp.model, s) / P

    for (label, terms) in forceapp.extra_forces
        haskey(result, label) && error("Invalid label: '$(label)'.")

        result[label] = sum(eval_pot(term, s) for term in terms)
    end

    result
end

function apply_force!(s::State{P,D,A,N}, forceapp::ForceApplicator;
                      fresh_force::Bool=true) where {P,D,A,N}
    if fresh_force || !forceapp.current_force.initialized
        fill!(forceapp.current_force.f, 0.0)
        forceapp.current_force.initialized = true

        eval_force!(forceapp.current_force, forceapp.model, s)

        forceapp.current_force.f ./= P

        for terms in values(forceapp.extra_forces)
            for term in terms
                eval_force!(forceapp.current_force, term, s)
            end
        end
    end

    s.ps .+= forceapp.dt .* forceapp.current_force.f

    nothing
end


struct Thermostat
    "Momentum dissipation coefficient (P)."
    coef1::Vector{Float64}
    "Momentum fluctuation coefficient (A, P; g nm / ps mol)."
    coef2::Matrix{Float64}
end

function Thermostat(beta::Float64, gamma0::Float64, dt::Float64,
                    omegas::AbstractVector{Float64})
    P = length(omegas)

    frictions = 2omegas # 1 / ps
    frictions[1] = gamma0

    coef1 = exp.(-dt .* frictions)
    coef2 = sqrt.((1.0 .- coef1.^2)' .* MS_fict(P) ./ beta)

    Thermostat(coef1, coef2)
end

function apply_thermostat!(nms::NormalModeState{P,D,A,N},
                           thermostat::Thermostat) where {P,D,A,N}
    for n in 1:N
        for a in 1:A
            for d in 1:D
                nms.Ps[:, d, a, n] .*= thermostat.coef1
                nms.Ps[:, d, a, n] .+= thermostat.coef2[a, :] .* randn(P)
            end
        end
    end

    nothing
end


abstract type Integrator end

const registered_integrators = Dict{String,Type{<:Integrator}}()

list_integrators() = sort(collect(keys(registered_integrators)))

get_integrator(name::String) = registered_integrators[name]

include("integrator_baoab.jl")
include("integrator_cbaoab.jl")
include("integrator_cobabo.jl")
include("integrator_obabo.jl")
