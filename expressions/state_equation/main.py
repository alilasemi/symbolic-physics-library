import numpy as np
import pickle
import sympy as sp

from physics.thermodynamics.thermochemical_data import ThermochemicalData

import expressions.state_equation.expressions as exprs
import expressions.state_equation.symbols as syms
import physics.thermodynamics.symbols as thermo_syms
from core.expression import Expression
import physics.constants as constants


def create():

    # Collect the necessary data
    thermochem_data = ThermochemicalData('NASA9')
    # TODO: Un-hardcode
    species = ['N2', 'N2+', 'N', 'N+', 'e-']
    heavies = ['N2', 'N2+', 'N', 'N+']
    idx_e = 4
    ns = len(species)
    nh = len(heavies)
    thermochem_data.data = {sp : thermochem_data[sp] for sp in species}
    M = np.array([thermochem_data[sp].M for sp in species])

    # -- Plug in -- #

    p_h = Expression(exprs.p_h_expr)
    p_h.plug_in(syms.nh, nh).doit()
    p_h.plug_in(thermo_syms.M, M)
    p_h.plug_in(thermo_syms.R, constants.R).simplify()

    p_e = Expression(exprs.p_e_expr)
    p_e.plug_in(syms.idx_e, idx_e).doit()
    p_e.plug_in(thermo_syms.M, M)
    p_e.plug_in(thermo_syms.R, constants.R)

    p = Expression(exprs.p_expr)
    p.plug_in(syms.p_h, p_h)
    p.plug_in(syms.p_e, p_e)

    # Save to file
    physics_file_name = 'physics.pkl'
    with open(physics_file_name, "wb") as physics_file:
        eqs_to_write = p
        pickle.dump(eqs_to_write, physics_file,
                protocol=pickle.HIGHEST_PROTOCOL)
