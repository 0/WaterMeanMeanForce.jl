# Microscopic water models.

"""
    WaterModel

Potentials and forces describing water molecules and their interactions.
"""
abstract type WaterModel end


"""
    TermContainer

Container for potential and force terms.
"""
struct TermContainer
    "Intramolecular terms."
    terms_intra::Vector{Function}
    "Intermolecular terms."
    terms_inter::Vector{Function}
end


function build_model(model_T::Type{<:WaterModel},
                     potential_intra::Function, potential_inter::Function,
                     force_intra::Function, force_inter::Function,
                     N::Int)
    potential_intra_terms = []
    potential_inter_terms = []
    force_intra_terms = []
    force_inter_terms = []

    for n1 in 1:N
        push!(potential_intra_terms, potential_intra(n1))
        push!(force_intra_terms, force_intra(n1))

        for n2 in (n1+1):N
            push!(potential_inter_terms, potential_inter(n1, n2))
            push!(force_inter_terms, force_inter(n1, n2))
        end
    end

    potential = TermContainer(potential_intra_terms, potential_inter_terms)
    force = TermContainer(force_intra_terms, force_inter_terms)

    model_T(potential, force)
end


"""
    ForceContainer

Container for evaluated forces. Each value is the total force on that degree of
freedom from all terms.
"""
mutable struct ForceContainer
    "Forces (kJ / mol nm)."
    f::Array{Float64,4}
    "Whether the forces have been set."
    initialized::Bool
end

function ForceContainer(P::Int, D::Int, A::Int, N::Int)
    f = Array{Float64}(undef, P, D, A, N)
    initialized = false

    ForceContainer(f, initialized)
end


"""
    eval_pot(model::WaterModel, s::State)

Evaluate the `model` potential for `s`.
"""
function eval_pot(model::WaterModel, s::State{P,D,A,N}) where {P,D,A,N}
    result = 0.0 # kJ / mol

    for term_f in model.potential.terms_intra
        for j in 1:P
            result += term_f(j, s)
        end
    end

    for term_f in model.potential.terms_inter
        for j in 1:P
            result += term_f(j, s)
        end
    end

    result
end

"""
    eval_force!(result::ForceContainer, model::WaterModel, s::State)

Evaluate the `model` force for `s`, storing the result in `result`.
"""
function eval_force!(result::ForceContainer, model::WaterModel,
                     s::State{P,D,A,N}) where {P,D,A,N}
    for term_f! in model.force.terms_intra
        for j in 1:P
            term_f!(result, j, s)
        end
    end

    for term_f! in model.force.terms_inter
        for j in 1:P
            term_f!(result, j, s)
        end
    end

    nothing
end


const registered_models = Dict{String,Type{<:WaterModel}}()

list_models() = sort(collect(keys(registered_models)))

get_model(name::String) = registered_models[name]

include("model_harmonic.jl")
include("model_qspcfw.jl")
