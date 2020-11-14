# Constrained BAOAB integrator.

struct ConstrainedBAOAB{P,D,A,N} <: Integrator
    "Reciprocal temperature (mol / kJ)."
    beta::Float64
    nmt::NormalModeTransformer
    nms::Maybe{NormalModeState}
    oscprop::OscillatorPropagator
    cstr::Constrainer
    forceapp::ForceApplicator
    thermostat::Maybe{Thermostat}
    "Instantaneous temperature for latest step (K)."
    temperature::Ref{Float64}
end

function ConstrainedBAOAB(
        model::WaterModel, R::Float64, dt::Float64, beta::Float64,
        gamma0::Float64, fc::ForceContainer, s::State{P,D,A,N};
        extra_forces::Dict{String,<:Vector{<:AdHocForce}}=Dict{String,Vector{AdHocForce}}(),
        use_thermostat::Bool=true) where {P,D,A,N}
    Cmat = make_Cmat(P)
    omegas = [2P * sin(pi * k / P) / (hbar * beta) for k in 0:(P-1)] # 1 / ps

    nmt = NormalModeTransformer(Cmat)
    nms = NormalModeState(s)
    # A
    oscprop = OscillatorPropagator(0.5 * dt, omegas)
    cstr = Constrainer(R, oscprop, Cmat)
    # B
    forceapp = ForceApplicator(model, extra_forces, 0.5 * dt, fc)
    # O
    if use_thermostat
        thermostat = Thermostat(beta, gamma0, dt, omegas)
    else
        thermostat = nothing
    end

    ConstrainedBAOAB{P,D,A,N}(beta, nmt, nms, oscprop, cstr, forceapp,
                              thermostat, NaN)
end

is_constrained(::Type{ConstrainedBAOAB}) = true

function N_dof(::ConstrainedBAOAB{P,D,A,N}) where {P,D,A,N}
    N * A * D * P - 1
end

thermostat_enabled(int::ConstrainedBAOAB) = !isnothing(int.thermostat)

function step!(s::State{P,D,A,N},
               int::ConstrainedBAOAB{P,D,A,N}) where {P,D,A,N}
    apply_force!(s, int.forceapp; fresh_force=false)
    apply_velocity_constraint!(s)

    cq = ConstraintQuantities(s)
    in_normal_modes(int.nms, s, int.nmt) do
        propagate_oscillators!(int.nms, int.oscprop)
    end
    apply_constraint!(s, int.cstr, cq)
    apply_velocity_constraint!(s)

    if thermostat_enabled(int)
        in_normal_modes(int.nms, s, int.nmt) do
            apply_thermostat!(int.nms, int.thermostat)
        end
        apply_velocity_constraint!(s)
    end
    int.temperature[] = temperature(MS_fict(P), N_dof(int), s)

    cq = ConstraintQuantities(s)
    in_normal_modes(int.nms, s, int.nmt) do
        propagate_oscillators!(int.nms, int.oscprop)
    end
    apply_constraint!(s, int.cstr, cq)

    apply_force!(s, int.forceapp)
    apply_velocity_constraint!(s)

    nothing
end

function energy(int::ConstrainedBAOAB{P,D,A,N},
                s::State{P,D,A,N}) where {P,D,A,N}
    result = Dict{String,Float64}()

    result["kinetic"] = kinetic_energy(MS_fict(P), s)
    result["spring"] = V_springs(int.beta, s)

    mergewith!(result, eval_pot(int.forceapp, s)) do
        error("Unexpected energy contribution.")
    end

    result
end

function get_force(int::ConstrainedBAOAB)
    int.forceapp.current_force.initialized || error("Force not initialized.")

    int.forceapp.current_force.f
end

get_temperature(int::ConstrainedBAOAB) = int.temperature[]


registered_integrators["c-BAOAB"] = ConstrainedBAOAB
