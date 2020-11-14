# BAOAB integrator.

# Liu, J., Li, D., & Liu, X. (2016). A simple and accurate algorithm for path
# integral molecular dynamics with the Langevin thermostat. The Journal of
# Chemical Physics, 145(2), 024103.

struct BAOAB{P,D,A,N} <: Integrator
    "Reciprocal temperature (mol / kJ)."
    beta::Float64
    nmt::NormalModeTransformer
    nms::NormalModeState
    oscprop::OscillatorPropagator
    forceapp::ForceApplicator
    thermostat::Maybe{Thermostat}
    "Instantaneous temperature for latest step (K)."
    temperature::Ref{Float64}
end

function BAOAB(
        model::WaterModel, dt::Float64, beta::Float64, gamma0::Float64,
        fc::ForceContainer, s::State{P,D,A,N};
        extra_forces::Dict{String,<:Vector{<:AdHocForce}}=Dict{String,Vector{AdHocForce}}(),
        use_thermostat::Bool=true) where {P,D,A,N}
    Cmat = make_Cmat(P)
    omegas = [2P * sin(pi * k / P) / (hbar * beta) for k in 0:(P-1)] # 1 / ps

    nmt = NormalModeTransformer(Cmat)
    nms = NormalModeState(s)
    # A
    oscprop = OscillatorPropagator(0.5 * dt, omegas)
    # B
    forceapp = ForceApplicator(model, extra_forces, 0.5 * dt, fc)
    # O
    if use_thermostat
        thermostat = Thermostat(beta, gamma0, dt, omegas)
    else
        thermostat = nothing
    end

    BAOAB{P,D,A,N}(beta, nmt, nms, oscprop, forceapp, thermostat, NaN)
end

is_constrained(::Type{BAOAB}) = false

function N_dof(::BAOAB{P,D,A,N}) where {P,D,A,N}
    N * A * D * P
end

thermostat_enabled(int::BAOAB) = !isnothing(int.thermostat)

function step!(s::State{P,D,A,N}, int::BAOAB{P,D,A,N}) where {P,D,A,N}
    apply_force!(s, int.forceapp; fresh_force=false)

    in_normal_modes(int.nms, s, int.nmt) do
        propagate_oscillators!(int.nms, int.oscprop)

        if thermostat_enabled(int)
            apply_thermostat!(int.nms, int.thermostat)
        end
        int.temperature[] = temperature(MS_fict(P), N_dof(int), int.nms)

        propagate_oscillators!(int.nms, int.oscprop)
    end

    apply_force!(s, int.forceapp)

    nothing
end

function energy(int::BAOAB{P,D,A,N}, s::State{P,D,A,N}) where {P,D,A,N}
    result = Dict{String,Float64}()

    result["kinetic"] = kinetic_energy(MS_fict(P), s)
    result["spring"] = V_springs(int.beta, s)

    mergewith!(result, eval_pot(int.forceapp, s)) do
        error("Unexpected energy contribution.")
    end

    result
end

function get_force(int::BAOAB)
    int.forceapp.current_force.initialized || error("Force not initialized.")

    int.forceapp.current_force.f
end

get_temperature(int::BAOAB) = int.temperature[]


registered_integrators["BAOAB"] = BAOAB
