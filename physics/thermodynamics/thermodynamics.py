from sympy import Sum, Piecewise
from physics.thermodynamics.symbols import *


# Species energy per mass
e_s_expr = (H_RT[s] - 1) * R * T / M[s]

# Species fully excited TR mass-specific heat
cv_tr_expr = (ndof/2) * R / M[s]
# Species fully excited TR energy per mass
e_s_tr_expr = cv_tr[s] * T
# Species fully excited VEE energy per mass
e_s_vee_expr = e_s[s] - e_s_tr[s]

# Mixture energy
e_expr = Sum(e_s[s] * Y[s], (s, 0, ns-1))

# Mixture specific heat
cv_expr = Sum(cv_R[s] * R * Y[s] / M[s], (s, 0, ns-1))

# Function to create piecewise expressions for the NASA fits. expression
# function of x. The coefficents are a, of size (n_ranges, n_coeffs).
def create_piecewise_expression_from_fit(expression, x, sym_a, a, temperature_ranges):
    n_ranges, n_coeffs = a.shape
    # TODO: Add out-of-bounds coefficients (currently just returns zero)
    return Piecewise(
            # The pieces are reversed so that high temperature comes first
            *reversed(
            # The default is set to 0
            [(0, True)] +
            # For the ith temp. range, the jth coefficient is plugged in
            [(expression.expression.subs([(sym_a[j], a[i, j])
                # A tuple is created by joining the expression with its
                # valid temperature range
                for j in range(n_coeffs)]), x < temperature_ranges[i + 1])
                for i in range(n_ranges)]))
