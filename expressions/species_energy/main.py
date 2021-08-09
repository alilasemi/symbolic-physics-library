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

    # -- Plug in expressions that are the same for all species -- #

    # Start with translational-rotational energy of species s
    e_s_tr = Expression(thermo_exprs.e_s_tr_expr)
    e_s_tr.plug_in(thermo_syms.cv_tr[syms.s], thermo_exprs.cv_tr_expr)

    # Plug in per-species data
    e_s_tr = np.array([copy.deepcopy(e_s_tr) for _ in species])
    for s, sp in enumerate(species):
        e_s_tr[s].plug_in(syms.s, s)
        e_s_tr[s].plug_in(thermo_syms.ndof, ndofs[sp])
        e_s_tr[s].plug_in(syms.M[s], thermochem_data[sp].M)
        e_s_tr[s].plug_in(syms.R, constants.R)

    # Start with energy of species s
    e_s = Expression(thermo_exprs.e_s_expr)
    # Plug in per-species data
    e_s = np.array([copy.deepcopy(e_s) for _ in species])
    # Plug in enthalpy fit
    for s, sp in enumerate(species):
        e_s[s].plug_in(syms.s, s)
        e_s[s].plug_in(syms.H_RT[s], thermo_fit_exprs.H_RT_expr)
        e_s[s].plug_in(syms.M[s], thermochem_data[sp].M)
        e_s[s].plug_in(syms.R, constants.R)
        # TODO: This is a mess, make this better
        e_s[s].expression = thermo_exprs.create_piecewise_expression_from_fit(e_s[s], syms.T, syms.a, thermochem_data[sp].a, thermochem_data[sp].temperatures)

    # Plug in to the vibrational energy
    e_s_vee = Expression(thermo_exprs.e_s_vee_expr)
    e_s_vee = np.array([copy.deepcopy(e_s_vee) for _ in species])
    for s, sp in enumerate(species):
        e_s_vee[s].plug_in(syms.s, s)
        e_s_vee[s].plug_in(thermo_syms.e_s_tr[s], e_s_tr[s])
        e_s_vee[s].plug_in(thermo_syms.e_s[s], e_s[s])

    # Save to file
    with open(physics_file_name, "wb") as physics_file:
        exprs_to_write = e_s_tr, e_s_vee
        pickle.dump(exprs_to_write, physics_file,
                protocol=pickle.HIGHEST_PROTOCOL)
