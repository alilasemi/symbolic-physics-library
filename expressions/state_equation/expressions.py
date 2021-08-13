from sympy import Sum
from physics.thermodynamics.symbols import *
from expressions.state_equation.symbols import *

# Ideal gas law, heavy pressure
p_h_expr = Sum(rho[s] / M[s], (s, 0, nh-1)) * R * T

# Ideal gas law, electron pressure
p_e_expr = Sum(rho[s] / M[s], (s, idx_e, idx_e)) * R * Tv

# Total pressure
p_expr = p_h + p_e
