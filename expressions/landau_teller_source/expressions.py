from sympy import exp, sqrt, S, Sum
from expressions.landau_teller_source.symbols import *

# Symbolic number 1. Used to keep fractions symbolic.
one = S.One

# Uncorrected vibrational relaxation time of species s. The constant on pressure
# converts Pa to atm, since the original equation expects atmospheres.
tau_v_s_expr = 101325/p * exp(a_vib * (T**(-one/3) - b_vib) - 18.42)

# Total number density computed from ideal gas law
n_expr = p / (k * T)

# Vibrational relaxation time of species s, with the Park 93 high temperature
# correction
tau_s_expr = tau_v_s + 1 / (n * sqrt(8 * k * T / (pi * m_bar)) * sigma_v)

# Mixture averaged species mass
m_bar_expr = Sum(Y[s] * m[s], (s, 0, ns-1))

# Vibrational cross section of species s in Park 93 high temperature correction
sigma_v_expr = sigma_v_p * (50000 / T)**2
