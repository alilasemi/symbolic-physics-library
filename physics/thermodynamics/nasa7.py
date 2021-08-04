from sympy import diff, simplify
from symbols import *


# NASA7 Enthalpy fit
#TODO: Add fits for NASA7
#H_RT_expr = (-a[0] * T**(-2) + a[1] * log(T) / T + a[2] +
#        a[3] * T / 2 + a[4] * T**2 / 3 + a[5] * T**3 / 4 +
#        a[6] * T**4 / 5 + a[7] / T)
#
## NASA7 Entropy fit
#S_R_expr = (-a[0] * T**(-2) / 2 - a[1] / T + a[2] * log(T) +
#        a[3] * T + a[4] * T**2 / 2 + a[5] * T**3 / 3 +
#        a[6] * T**4 / 4 + a[8])

# Resulting fit for Gibbs energy, G/RT
G_RT_expr = H_RT_expr - S_R_expr

# Resulting specific heat fits
cp_R_expr = simplify(diff(H_RT_expr * T, T))
cv_R_expr = cp_R_expr - 1
