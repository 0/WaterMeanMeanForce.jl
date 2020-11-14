# Driver functions for the integrators.

struct MolecularDynamicsParameters
    "Centroid friction (1 / ps)."
    gamma0::Float64
    "Time step (ps)."
    dt::Float64
    "Total duration (ps)."
    duration::Float64
end

"""
    num_steps(mdp::MolecularDynamicsParameters)

Number of molecular dynamics steps for `mdp`.
"""
num_steps(mdp::MolecularDynamicsParameters) = ceil(Int, mdp.duration / mdp.dt)


function run_constrained_equilibration(
        model_T::Type{<:WaterModel}, int_T::Type{<:Integrator}, N::Int,
        T::Float64, R::Float64, P::Int, mdp_equil::MolecularDynamicsParameters;
        qs_init::Maybe{Vector{Float64}}=nothing, output_skip::Int=100)
    is_constrained(int_T) || error("Integrator is not constrained.")

    # Number of atoms in each water molecule.
    A = 3
    # Number of spatial dimensions.
    D = 3

    if !isnothing(qs_init)
        s = State(P, D, A, N; qs_init)
    else
        s = State(P, D, A, N; R)
    end

    model = model_T(N)
    fc = ForceContainer(P, D, A, N)
    int_equil = int_T(model, R, mdp_equil.dt, beta(T), mdp_equil.gamma0, fc, s)

    for i in 1:num_steps(mdp_equil)
        step!(s, int_equil)

        i % output_skip == 0 || continue

        time = i * mdp_equil.dt
        E = energy(int_equil, s)
        E_K = pop!(E, "kinetic")
        E_Vspring = pop!(E, "spring")
        E_Vmodel = pop!(E, "model")
        isempty(E) || error("Unexpected energy contribution.")
        E_total = E_K + E_Vspring + E_Vmodel
        T = get_temperature(int_equil)

        println("$(time) $(E_K) $(E_Vspring) $(E_Vmodel) $(E_total) $(T)")
    end

    s
end

function run_constrained_production(
        model_T::Type{<:WaterModel}, int_T::Type{<:Integrator}, N::Int,
        T::Float64, R::Float64, P::Int, mdp_equil::MolecularDynamicsParameters,
        mdp_prod::MolecularDynamicsParameters;
        qs_init::Maybe{Vector{Float64}}=nothing, aux_output::IO=stderr)
    is_constrained(int_T) || error("Integrator is not constrained.")

    # Number of atoms in each water molecule.
    A = 3
    # Number of spatial dimensions.
    D = 3

    # Number of molecular dynamics steps for equilibration.
    N_equil = num_steps(mdp_equil)
    # Number of molecular dynamics steps for production.
    N_prod = num_steps(mdp_prod)

    if !isnothing(qs_init)
        s = State(P, D, A, N; qs_init)
    else
        s = State(P, D, A, N; R)
    end

    model = model_T(N)
    fc = ForceContainer(P, D, A, N)
    int_equil = int_T(model, R, mdp_equil.dt, beta(T), mdp_equil.gamma0, fc, s)
    int_prod = int_T(model, R, mdp_prod.dt, beta(T), mdp_prod.gamma0, fc, s)

    meter_equil = Progress(N_equil, output=aux_output)

    for _ in ProgressWrapper(1:N_equil, meter_equil)
        step!(s, int_equil)
    end

    constraints = Array{Float64}(undef, N_prod)
    E_Ks = Array{Float64}(undef, N_prod)
    E_Vsprings = Array{Float64}(undef, N_prod)
    E_Vmodels = Array{Float64}(undef, N_prod)
    E_totals = Array{Float64}(undef, N_prod)
    temperatures = Array{Float64}(undef, N_prod)

    forces1 = Array{Float64}(undef, N_prod)
    forces2 = Array{Float64}(undef, N_prod)

    estim1 = estimator1(R, beta(T), int_prod, s)
    estim2 = estimator2(R, beta(T), int_prod, s)

    meter = Progress(N_prod, output=aux_output)

    for i in ProgressWrapper(1:N_prod, meter)
        step!(s, int_prod)

        constraints[i] = com_r(CN1, CN2, CJ, s)
        E = energy(int_prod, s)
        E_Ks[i] = pop!(E, "kinetic")
        E_Vsprings[i] = pop!(E, "spring")
        E_Vmodels[i] = pop!(E, "model")
        isempty(E) || error("Unexpected energy contribution.")
        E_totals[i] = E_Ks[i] + E_Vsprings[i] + E_Vmodels[i]
        temperatures[i] = get_temperature(int_prod)

        forces1[i] = estim1()
        forces2[i] = estim2()
    end

    data_table = []

    for (label, data) in [("xi (nm)", constraints),
                          ("K (kJ / mol)", E_Ks),
                          ("Vspring (kJ / mol)", E_Vsprings),
                          ("Vmodel (kJ / mol)", E_Vmodels),
                          ("E (kJ / mol)", E_totals),
                          ("T (K)", temperatures),
                          ("est1 (1 / nm)", forces1),
                          ("est2 (1 / nm)", forces2)]
        data_mean, data_stderr, data_actime = aggregate(mdp_prod.dt, data)
        push!(data_table, [label, data_mean, data_stderr, data_actime])
    end

    pretty_table(aux_output, permutedims(hcat(data_table...)),
                 ["", "Value", "Error", "Autocorrelation time (ps)"])

    s, aggregate(mdp_prod.dt, forces1), aggregate(mdp_prod.dt, forces2)
end


function run_restrained_equilibration(
        model_T::Type{<:WaterModel}, int_T::Type{<:Integrator}, N::Int,
        T::Float64, R::Float64, P::Int, mdp_equil::MolecularDynamicsParameters,
        bias_force_constant::Float64; qs_init::Maybe{Vector{Float64}}=nothing,
        output_skip::Int=100)
    is_constrained(int_T) && error("Integrator is constrained.")

    # Number of atoms in each water molecule.
    A = 3
    # Number of spatial dimensions.
    D = 3

    if !isnothing(qs_init)
        s = State(P, D, A, N; qs_init)
    else
        s = State(P, D, A, N; R)
    end

    bias = HarmonicRestraint(bias_force_constant, R, CN1, CN2, CJ)
    model = model_T(N)
    fc = ForceContainer(P, D, A, N)
    int_equil = int_T(model, mdp_equil.dt, beta(T), mdp_equil.gamma0, fc, s;
                      extra_forces=Dict("bias" => [bias]))

    for i in 1:num_steps(mdp_equil)
        step!(s, int_equil)

        i % output_skip == 0 || continue

        time = i * mdp_equil.dt
        xi = com_r(CN1, CN2, CJ, s)
        E = energy(int_equil, s)
        E_K = pop!(E, "kinetic")
        E_Vspring = pop!(E, "spring")
        E_Vmodel = pop!(E, "model")
        E_Vbias = pop!(E, "bias")
        isempty(E) || error("Unexpected energy contribution.")
        E_total = E_K + E_Vspring + E_Vmodel + E_Vbias
        T = get_temperature(int_equil)

        print("$(time) $(xi) $(E_K) $(E_Vspring) $(E_Vmodel) $(E_Vbias) ")
        println("$(E_total) $(T)")
    end

    s
end

function run_restrained_production(
        model_T::Type{<:WaterModel}, int_T::Type{<:Integrator}, N::Int,
        T::Float64, R::Float64, P::Int, mdp_equil::MolecularDynamicsParameters,
        mdp_prod::MolecularDynamicsParameters, bias_force_constant::Float64;
        qs_init::Maybe{Vector{Float64}}=nothing, output_skip::Int=100)
    is_constrained(int_T) && error("Integrator is constrained.")

    # Number of atoms in each water molecule.
    A = 3
    # Number of spatial dimensions.
    D = 3

    # Number of molecular dynamics steps for equilibration.
    N_equil = num_steps(mdp_equil)
    # Number of molecular dynamics steps for production.
    N_prod = num_steps(mdp_prod)

    if !isnothing(qs_init)
        s = State(P, D, A, N; qs_init)
    else
        s = State(P, D, A, N; R)
    end

    bias = HarmonicRestraint(bias_force_constant, R, CN1, CN2, CJ)
    model = model_T(N)
    fc = ForceContainer(P, D, A, N)
    int_equil = int_T(model, mdp_equil.dt, beta(T), mdp_equil.gamma0, fc, s;
                      extra_forces=Dict("bias" => [bias]))
    int_prod = int_T(model, mdp_prod.dt, beta(T), mdp_prod.gamma0, fc, s;
                     extra_forces=Dict("bias" => [bias]))

    for _ in 1:N_equil
        step!(s, int_equil)
    end

    for i in 1:N_prod
        step!(s, int_prod)

        i % output_skip == 0 || continue

        println(com_r(CN1, CN2, CJ, s))
    end

    s
end
