#!/usr/bin/env julia

using WaterMeanMeanForce

using ArgParse

s = ArgParseSettings()
s.autofix_names = true

@add_arg_table! s begin
    "-T"
        metavar = "T"
        help = "comma-separated list of temperatures (K)"
        required = true
end

c = parse_args(ARGS, s, as_symbols=true)


for T_str in split(c[:T], ",")
    T = parse(Float64, T_str)
    beta, slope = harmonic_exact_slope(2, T)

    println("$(T) $(beta) $(slope)")
end
