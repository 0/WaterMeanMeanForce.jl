# Physical units.

# mass:        g / mol
# length:      nm
# time:        ps
# momentum:    g nm / ps mol
# energy:      kJ / mol
# temperature: K
# charge:      e

const deg_to_rad = pi / 180.0
const to_kilo = 1e-3
const to_nano = 1e9
const to_pico = 1e12

# NIST Special Publication 811 (2008), Appendix B.8.
const A_to_nm = 0.1
const kcal_to_kJ = 4.184

# 2019 SI definitions from BIPM.
const e_to_C = 1.602_176_634e-19

# 2018 CODATA recommended values from NIST.
const N_A = 6.022_140_76e23 # 1 / mol
const planck = 6.626_070_15e-34 * to_kilo * to_pico * N_A # kJ ps / mol
const hbar = planck / 2pi
const kB = 1.380_649e-23 * to_kilo * N_A # kJ / mol K


"""
    beta(T::Float64)

Compute the reciprocal temperature for `T`.
"""
beta(T::Float64) = 1.0 / (kB * T) # mol / kJ
