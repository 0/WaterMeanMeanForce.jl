# Additional force fields.

abstract type AdHocForce end


"""
    HarmonicRestraint

Harmonic restraint between the centers of mass of two molecules at a single
bead.
"""
struct HarmonicRestraint <: AdHocForce
    "Force constant (kJ / nm^2 mol)."
    k::Float64
    "Separation (nm)."
    R::Float64
    "Index of first molecule."
    n1::Int
    "Index of second molecule."
    n2::Int
    "Bead index."
    j::Int
end

function eval_pot(hr::HarmonicRestraint, s::State)
    0.5 * hr.k * (com_r(hr.n1, hr.n2, hr.j, s) - hr.R)^2
end

function eval_force!(result::ForceContainer, hr::HarmonicRestraint, s::State)
    factor = -hr.k * (com_r(hr.n1, hr.n2, hr.j, s) - hr.R) # kJ / nm mol
    result.f[hr.j, :, :, :] .+= factor .* com_r_grad(hr.n1, hr.n2, hr.j, s)

    nothing
end
