import numpy as np
import pickle
import sympy as sp

from physics.thermodynamics.thermochemical_data import ThermochemicalData

from core.expression import Expression
import physics.constants as constants
import physics.thermodynamics.thermodynamics as exprs
import physics.thermodynamics.symbols as syms
import physics.thermodynamics.nasa9 as thermo_fit_exprs


def create():

    physics_file_name = 'physics.pkl'

    # -- Collect the necessary data -- #

    thermochem_data = ThermochemicalData('NASA9')
    # TODO: Un-hardcode
    species = ['N2', 'N2+', 'N', 'N+', 'e-']
    ns = len(species)
    thermochem_data.data = {sp : thermochem_data[sp] for sp in species}
    ndofs = {'N2': 5, 'N2+': 5, 'N': 3, 'N+': 3, 'e-': 0}

    # -- Plug in -- #

    # For the translational-rotational mode
    e_s_tr = np.empty(ns, dtype=object)
    cv_s_tr = np.empty(ns, dtype=object)
    for s, sp in enumerate(species):
        # cv
        cv_s_tr[s] = Expression(exprs.cv_s_tr_expr)
        cv_s_tr[s].plug_in(syms.s, s)
        cv_s_tr[s].plug_in(syms.ndof, ndofs[sp])
        cv_s_tr[s].plug_in(syms.M[s], thermochem_data[sp].M)
        cv_s_tr[s].plug_in(syms.R, constants.R)
        # e
        e_s_tr[s] = Expression(exprs.e_s_tr_expr)
        e_s_tr[s].plug_in(syms.cv_s_tr[syms.s], cv_s_tr[s])
    # Mixture averaged values
    e_tr = Expression(exprs.e_tr_expr)
    e_tr.plug_in(syms.ns, ns).doit()
    e_tr.plug_in(syms.e_s_tr, e_s_tr)
    cv_tr = Expression(exprs.cv_tr_expr)
    cv_tr.plug_in(syms.ns, ns).doit()
    cv_tr.plug_in(syms.cv_s_tr, cv_s_tr)

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
    # TR mode from other modes
    e_tr_from_e = Expression(exprs.e_tr_from_e_expr)
    e_tr_from_e.plug_in(syms.ns, ns).doit()
    e_tr_from_e.plug_in(syms.e_s_0, e_s_0)

    # For the vibrational-electronic-electron mode
    e_s_vee = np.empty(ns, dtype=object)
    cv_s_vee = np.empty(ns, dtype=object)
    for s, sp in enumerate(species):
        e_s_vee[s] = Expression(exprs.e_s_vee_expr)
        e_s_vee[s].plug_in(syms.e_s_tr[syms.s], e_s_tr[s])
        e_s_vee[s].plug_in(syms.e_s[syms.s], e_s[s])
        e_s_vee[s].plug_in(syms.e_s_0[syms.s], e_s_0[s])
        cv_s_vee[s] = Expression(exprs.cv_s_vee_expr)
        cv_s_vee[s].plug_in(syms.cv_s_tr[syms.s], cv_s_tr[s])
        cv_s_vee[s].plug_in(syms.cv_s[syms.s], cv_s[s])
    # Mixture averaged values
    e_vee = Expression(exprs.e_vee_expr)
    e_vee.plug_in(syms.ns, ns).doit()
    e_vee.plug_in(syms.e_s_vee, e_s_vee)
    cv_vee = Expression(exprs.cv_vee_expr)
    cv_vee.plug_in(syms.ns, ns).doit()
    cv_vee.plug_in(syms.cv_s_vee, cv_s_vee)

    # Save to file
    with open(physics_file_name, "wb") as physics_file:
        exprs_to_write = e, e_tr, e_tr_from_e, e_vee, cv, cv_tr, cv_vee, e_s, e_s_tr, e_s_vee, cv_s, cv_s_tr, cv_s_vee
        pickle.dump(exprs_to_write, physics_file,
                protocol=pickle.HIGHEST_PROTOCOL)
