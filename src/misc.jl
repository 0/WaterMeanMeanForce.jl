# Miscellaneous items.

Maybe{T} = Union{T,Nothing} where {T}

function acos_(x::Float64)
    if abs(x) > 1.0
        if 1.0 < x < 1.0 + 1e-12
            x = 1.0
        elseif -1.0 - 1e-12 < x < -1.0
            x = -1.0
        else
            error("Invalid cosine: $(x).")
        end
    end

    acos(x)
end

"Index of oxygen atom."
const AO = 1
"Index of first hydrogen atom."
const AH1 = 2
"Index of second hydrogen atom."
const AH2 = 3

"Water atomic masses (g / mol)."
const MS = [15.999, 1.008, 1.008]
"Fractional masses of atoms in a single water molecule."
const MS_FRAC = MS ./ sum(MS)

"First molecule of constraint."
const CN1 = 1
"Second molecule of constraint."
const CN2 = 2
"Bead index for constraint."
const CJ = 1

bead_prev(P::Int, j::Int) = mod(j - 1, 1:P)
bead_next(P::Int, j::Int) = mod(j + 1, 1:P)
