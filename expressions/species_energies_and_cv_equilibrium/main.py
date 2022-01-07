import numpy as np
import pickle
import sympy as sp

from physics.thermodynamics.thermochemical_data import ThermochemicalData

from core.expression import Expression
import physics.constants as constants
import physics.thermodynamics.thermodynamics as exprs
import physics.thermodynamics.symbols as syms
import physics.thermodynamics.nasa9 as thermo_fit_exprs


def create(species):

    physics_file_name = 'physics.pkl'

    # -- Collect the necessary data -- #

    thermochem_data = ThermochemicalData('NASA9')
    ns = len(species)
    # TODO
    thermochem_data.data = {sp : thermochem_data[sp] for sp in species}
    ndofs = {'N2': 5, 'N2+': 5, 'N': 3, 'N+': 3, 'e-': 0}

    # -- Plug in -- #

    # For the combined energy of all modes
    e_s = np.empty(ns, dtype=object)
    e_s_0 = np.empty(ns)
    cv_s = np.empty(ns, dtype=object)
    for s, sp in enumerate(species):
        # e
        e_s[s] = Expression(exprs.e_s_expr)
        e_s[s].plug_in(syms.s, s)
        e_s[s].plug_in(syms.H_RT[s], thermo_fit_exprs.H_RT_expr)
        e_s[s].plug_in(syms.M[s], thermochem_data[sp].M)
        e_s[s].plug_in(syms.R, constants.R)
        # cv
        cv_s[s] = Expression(exprs.cv_s_func, e_s[s]).simplify()
        # Make piecewise expressions at the end
        # TODO: This is a mess, make this better
        e_s[s].expression = exprs.create_piecewise_expression_from_fit(e_s[s], syms.T, syms.a, thermochem_data[sp].a, thermochem_data[sp].temperatures)
        cv_s[s].expression = exprs.create_piecewise_expression_from_fit(cv_s[s], syms.T, syms.a, thermochem_data[sp].a, thermochem_data[sp].temperatures)
        # Formation energy
        e_s_0[s] = e_s[s].expression.subs(syms.T, thermochem_data.T_0)
    # Mixture averaged values
    e = Expression(exprs.e_expr)
    e.plug_in(syms.ns, ns).doit()
    e.plug_in(syms.e_s, e_s)
    cv = Expression(exprs.cv_expr)
    cv.plug_in(syms.ns, ns).doit()
    cv.plug_in(syms.cv_s, cv_s)

    # energy and cv will be bundled together, to allow common subexpression
    # elimination
    e_and_cv = [e, cv]

    # Save to file
    with open(physics_file_name, "wb") as physics_file:
        exprs_to_write = e, cv, e_and_cv, e_s, cv_s
        pickle.dump(exprs_to_write, physics_file,
                protocol=pickle.HIGHEST_PROTOCOL)
