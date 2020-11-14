# Center of mass constraint.

"""
    com(n::Int, j::Int, s::State)

Compute the center of mass position for molecule `n` at bead `j` in `s`.
"""
function com(n::Int, j::Int, s::State{P,D,A,N}) where {P,D,A,N}
    result = zeros(Float64, D) # nm

    for a in 1:A
        result .+= MS_FRAC[a] .* s.qs[j, :, a, n]
    end

    result
end

"""
    com_disp(n1::Int, n2::Int, j::Int, s::State)

Compute the center of mass displacement vector from molecule `n2` to `n1` at
bead `j` in `s`.
"""
com_disp(n1::Int, n2::Int, j::Int, s::State) = com(n1, j, s) .- com(n2, j, s)

"""
    com_r(n1::Int, n2::Int, j::Int, s::State)

Compute the center of mass distance between molecules `n1` and `n2` at bead `j`
in `s`.
"""
com_r(n1::Int, n2::Int, j::Int, s::State) = com_disp(n1, n2, j, s) |> norm

"""
    com_r_grad(n1::Int, n2::Int, j::Int, s::State)

Compute the gradient of the center of mass distance between molecules `n1` and
`n2` at bead `j` in `s`.
"""
function com_r_grad(n1::Int, n2::Int, j::Int,
                    s::State{P,D,A,N}) where {P,D,A,N}
    r = com_disp(n1, n2, j, s)
    r ./= norm(r)

    result = zeros(Float64, D, A, N)

    for a in 1:A
        result[:, a, n1] .= +MS_FRAC[a] .* r
        result[:, a, n2] .= -MS_FRAC[a] .* r
    end

    result
end
