from sympy import exp, sqrt, S
from symbols import *

# Symbolic number 1. Used to keep fractions symbolic.
one = S.One

# Uncorrected vibrational relaxation time of species s. The constant on pressure
# converts Pa to atm, since the original equation expects atmospheres.
tau_v_s_expr = 101325/p * exp(a_vib * (T**(-one/3) - b_vib) - 18.42)

# Total number density computed from mass density
n_expr = N_A * rho / M_bar

# Vibrational relaxation time of species s, with the Park 93 high temperature
# correction
tau_s_expr = tau_v_s + 1 / (n * sqrt(8 * k * T / (pi * m[s])) * sigma_v)

# Vibrational cross section of species s in Park 93 high temperature correction
sigma_v_expr = sigma_v_p * (50000 / T)**2
