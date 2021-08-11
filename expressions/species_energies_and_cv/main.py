import numpy as np
import pickle
import sympy as sp
import copy
from sympy.codegen.ast import Assignment

from physics.thermodynamics.thermochemical_data import ThermochemicalData

import expressions.species_energy.expressions as exprs
import physics.thermodynamics.thermodynamics as thermo_exprs
import expressions.species_energy.symbols as syms
import physics.thermodynamics.symbols as thermo_syms
import physics.thermodynamics.nasa9 as thermo_fit_exprs

from expression import Expression
import physics.constants as constants


def create():

    physics_file_name = 'physics.pkl'

    # Collect the necessary data
    thermochem_data = ThermochemicalData('NASA9')
    # TODO: Un-hardcode
    species = ['N2', 'N2+', 'N', 'N+', 'e-']
    ns = len(species)
    thermochem_data.data = {sp : thermochem_data[sp] for sp in species}
    ndofs = {'N2': 5, 'N2+': 5, 'N': 3, 'N+': 3, 'e-': 0}

    # -- Plug in -- #

    # Start with translational-rotational energy of species s
    #e_s_tr = Expression(thermo_exprs.e_s_tr_expr)
    #e_s_tr.plug_in(thermo_syms.cv_tr[syms.s], thermo_exprs.cv_tr_expr)

    # Plug in per-species data
    e_s_tr = np.empty(ns, dtype=object)
    cv_s_tr = np.empty(ns, dtype=object)
    for s, sp in enumerate(species):
        # cv
        cv_s_tr[s] = Expression(thermo_exprs.cv_s_tr_expr)
        cv_s_tr[s].plug_in(syms.s, s)
        cv_s_tr[s].plug_in(thermo_syms.ndof, ndofs[sp])
        cv_s_tr[s].plug_in(syms.M[s], thermochem_data[sp].M)
        cv_s_tr[s].plug_in(syms.R, constants.R)
        # e
        e_s_tr[s] = Expression(thermo_exprs.e_s_tr_expr)
        e_s_tr[s].plug_in(thermo_syms.cv_s_tr[syms.s], cv_s_tr[s])

    # Start with energy of species s
    # Plug in per-species data
    e_s = np.empty(ns, dtype=object)
    cv_s = np.empty(ns, dtype=object)
    # Plug in enthalpy fit
    for s, sp in enumerate(species):
        # e
        e_s[s] = Expression(thermo_exprs.e_s_expr)
        e_s[s].plug_in(syms.s, s)
        e_s[s].plug_in(syms.H_RT[s], thermo_fit_exprs.H_RT_expr)
        e_s[s].plug_in(syms.M[s], thermochem_data[sp].M)
        e_s[s].plug_in(syms.R, constants.R)
        # cv
        cv_s[s] = Expression(thermo_exprs.cv_s_func, e_s[s]).simplify()
        # Make piecewise expressions at the end
        # TODO: This is a mess, make this better
        e_s[s].expression = thermo_exprs.create_piecewise_expression_from_fit(e_s[s], syms.T, syms.a, thermochem_data[sp].a, thermochem_data[sp].temperatures)
        cv_s[s].expression = thermo_exprs.create_piecewise_expression_from_fit(cv_s[s], syms.T, syms.a, thermochem_data[sp].a, thermochem_data[sp].temperatures)

    # Plug in to the vibrational energy
    e_s_vee = np.empty(ns, dtype=object)
    cv_s_vee = np.empty(ns, dtype=object)
    for s, sp in enumerate(species):
        e_s_vee[s] = Expression(thermo_exprs.e_s_vee_expr)
        e_s_vee[s].plug_in(thermo_syms.e_s_tr[syms.s], e_s_tr[s])
        e_s_vee[s].plug_in(thermo_syms.e_s[syms.s], e_s[s])
        cv_s_vee[s] = Expression(thermo_exprs.cv_s_vee_expr)
        cv_s_vee[s].plug_in(thermo_syms.cv_s_tr[syms.s], cv_s_tr[s])
        cv_s_vee[s].plug_in(thermo_syms.cv_s[syms.s], cv_s[s])

    # Save to file
    with open(physics_file_name, "wb") as physics_file:
        exprs_to_write = e_s, e_s_tr, e_s_vee, cv_s, cv_s_tr, cv_s_vee
        pickle.dump(exprs_to_write, physics_file,
                protocol=pickle.HIGHEST_PROTOCOL)
