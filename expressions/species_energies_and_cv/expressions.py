from sympy import exp, log, Sum, Product, diff, simplify
from symbols import *


# NASA9 Enthalpy fit
H_RT_expr = (-a[0] * T**(-2) + a[1] * log(T) / T + a[2] +
        a[3] * T / 2 + a[4] * T**2 / 3 + a[5] * T**3 / 4 +
        a[6] * T**4 / 5 + a[7] / T)
