import matplotlib.pyplot as plt
from matplotlib import rc
import sympy as sp
from sympy import exp, log, Sum, Product, diff, simplify

from symbols import *


# Species mass production rate
wdot_expr = M[s] * Sum((beta[s, r] - alpha[s, r]) * (Rf[r] - Rb[r]), (r, 0, nr-1))

# NASA9 Enthalpy fit
H_RT_expr = (-a[0] * T**(-2) + a[1] * log(T) / T + a[2] +
        a[3] * T / 2 + a[4] * T**2 / 3 + a[5] * T**3 / 4 +
        a[6] * T**4 / 5 + a[7] / T)

# NASA9 Entropy fit
S_R_expr = (-a[0] * T**(-2) / 2 - a[1] / T + a[2] * log(T) +
        a[3] * T + a[4] * T**2 / 2 + a[5] * T**3 / 3 +
        a[6] * T**4 / 4 + a[8])

# Resulting fit for Gibbs energy, G/RT
G_RT_expr = H_RT_expr - S_R_expr

# Resulting specific heat fits
cp_R_expr = simplify(diff(H_RT_expr * T, T))
cv_R_expr = cp_R_expr - 1

# Species energy per mass
e_s_expr = (H_RT[s] - 1) * R * T / M[s]

# Mixture energy
e_expr = Sum(e_s[s] * Y[s], (s, 0, ns-1))

# Mixture specific heat
cv_expr = Sum(cv_R[s] * R * Y[s] / M[s], (s, 0, ns-1))

# Equilibrium constant
K_c_expr = (p_0 / (R * T)) ** nu[r] \
        * exp( - Sum((beta[s, r] - alpha[s, r]) * G_RT[s], (s, 0, ns-1)))

# Forward rate
Rf_expr = kf * Product((rho[s] / M[s]) ** alpha[s, r], (s, 0, ns-1))
# Backward rate
Rb_expr = kb * Product((rho[s] / M[s]) ** beta[s, r], (s, 0, ns-1))

# Arrhenius forward rate coefficient
kf_arrhenius_expr = C[r] * T**(eta[r]) * exp(-theta[r] / T)

# Three-body forward rate coefficient
kf_three_body_expr = kf * Sum(epsilon[s] * rho[s] / M[s], (s, 0, ns-1))

# Backward rate coefficient
kb_expr = kf / K_c
