#!/usr/bin/env julia

using DelimitedFiles

using WaterMeanMeanForce

using ArgParse

s = ArgParseSettings()
s.autofix_names = true

models = join(list_models(), ", ")
integrators = join(list_integrators(), ", ")

@add_arg_table! s begin
    "--model"
        metavar = "M"
        range_tester = in(list_models())
        help = "water model ($(models))"
        required = true
    "--integrator"
        metavar = "I"
        range_tester = in(list_integrators())
        help = "molecular dynamics integrator ($(integrators))"
        required = true
    "-T"
        metavar = "T"
        help = "temperature (K)"
        arg_type = Float64
        required = true
    "-R"
        metavar = "R"
        help = "constraint/restraint distance (nm)"
        arg_type = Float64
        required = true
    "--num-links"
        metavar = "P"
        help = "number of Trotter links"
        arg_type = Int
        required = true
    "--time-step"
        metavar = "T"
        help = "time step (ps)"
        arg_type = Float64
        required = true
    "--centroid-friction"
        metavar = "G"
        help = "centroid friction (1 / ps)"
        arg_type = Float64
        required = true
    "--equil-duration"
        metavar = "D"
        help = "duration of equilibration (ps)"
        arg_type = Float64
        required = true
    "--umbrella-k"
        metavar = "K"
        help = "umbrella sampling force constant (kJ / nm^2 mol)"
        arg_type = Float64
    "--initial-configuration"
        metavar = "FILE"
        help = "path to initial configuration input file"
    "--final-configuration"
        metavar = "FILE"
        help = "path to final configuration output file"
end

c = parse_args(ARGS, s, as_symbols=true)

if !isnothing(c[:initial_configuration])
    qs_init = dropdims(readdlm(c[:initial_configuration]); dims=2)
else
    qs_init = nothing
end


mdp_equil = MolecularDynamicsParameters(c[:centroid_friction], c[:time_step],
                                        c[:equil_duration])

if isnothing(c[:umbrella_k])
    s = run_constrained_equilibration(get_model(c[:model]),
                                      get_integrator(c[:integrator]), 2, c[:T],
                                      c[:R], c[:num_links], mdp_equil; qs_init)
else
    s = run_restrained_equilibration(get_model(c[:model]),
                                     get_integrator(c[:integrator]), 2, c[:T],
                                     c[:R], c[:num_links], mdp_equil,
                                     c[:umbrella_k]; qs_init)
end

if !isnothing(c[:final_configuration])
    writedlm(c[:final_configuration], s.qs)
end
