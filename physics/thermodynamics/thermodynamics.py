from sympy import Sum, Piecewise, diff
from physics.thermodynamics.symbols import *

def mass_average(expression):
    return Sum(expression[s] * Y[s], (s, 0, ns-1))

# == Species Quantities == #

# -- Translational-rotational component -- #
# Species fully excited TR mass-specific heat
cv_s_tr_expr = (ndof/2) * R / M[s]
# Species fully excited TR energy per mass
e_s_tr_expr = cv_s_tr[s] * T
# -- Vibrational-electronic-electron component -- #
# Species VEE energy per mass
e_s_vee_expr = e_s[s] - e_s_tr[s] - e_s_0[s]
# Species VEE mass-specific heat
cv_s_vee_expr = cv_s[s] - cv_s_tr[s]
# -- Combined energy from all components -- #
# Species energy per mass
e_s_expr = (H_RT[s] - 1) * R * T / M[s]
# Species mass-specific heat
cv_s_func = lambda e_s: diff(e_s, T)

# == Mixture Quantities == #
# Since these are mass-specific thermodynamic quantities, their mixture average
# is mass weighted.

# -- Translational-rotational component -- #
e_tr_expr  = mass_average(e_s_tr)
cv_tr_expr = mass_average(cv_s_tr)
# From other energies
e_tr_from_e_expr = e - e_vee - mass_average(e_s_0)
# -- Vibrational-electronic-electron component -- #
e_vee_expr  = mass_average(e_s_vee)
cv_vee_expr = mass_average(cv_s_vee)
# -- Combined energy from all components -- #
e_expr  = mass_average(e_s)
cv_expr = mass_average(cv_s)

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
                for j in range(n_coeffs)]), x > temperature_ranges[i])
                for i in range(n_ranges)]))
