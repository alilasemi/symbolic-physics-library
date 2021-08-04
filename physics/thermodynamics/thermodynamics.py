from sympy import Sum
from symbols import *


# Species energy per mass
e_s_expr = (H_RT[s] - 1) * R * T / M[s]

# Mixture energy
e_expr = Sum(e_s[s] * Y[s], (s, 0, ns-1))

# Mixture specific heat
cv_expr = Sum(cv_R[s] * R * Y[s] / M[s], (s, 0, ns-1))
