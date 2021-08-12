from sympy import Sum, Product, exp
from physics.chemistry.symbols import *

# Species mass production rate
wdot_expr = M[s] * Sum((beta[s, r] - alpha[s, r]) * (Rf[r] - Rb[r]), (r, 0, nr-1))

# Equilibrium constant
K_c_expr = (p_0 / (R * T)) ** nu \
        * exp( - Sum((beta[s, r] - alpha[s, r]) * G_RT[s], (s, 0, ns-1)))

# Forward rate
Rf_expr = kf * Product((rho[s] / M[s]) ** alpha[s, r], (s, 0, ns-1))
# Backward rate
Rb_expr = kb * Product((rho[s] / M[s]) ** beta[s, r], (s, 0, ns-1))

# Arrhenius forward rate coefficient
kf_arrhenius_expr = C * T**(eta) * exp(-theta / T)

# Three-body forward rate coefficient
kf_three_body_expr = kf * Sum(epsilon[s] * rho[s] / M[s], (s, 0, ns-1))

# Backward rate coefficient
kb_expr = kf / K_c
