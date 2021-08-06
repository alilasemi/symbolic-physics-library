import math
import numpy as np
import pickle
import sympy as sp
from sympy.codegen.ast import Assignment

from physics.thermodynamics.thermochemical_data import ThermochemicalData

import expressions.landau_teller_source.expressions as exprs
import expressions.landau_teller_source.symbols as syms
from expression import Expression
import physics.constants as constants


def create():

    # Collect the necessary data
    thermochem_data = ThermochemicalData('NASA9')
    # TODO: Un-hardcode
    species = ['N2', 'N2+', 'N', 'N+', 'e-']
    vibrators = ['N2', 'N2+']
    heavies = ['N2', 'N2+', 'N', 'N+']
    ns = len(species)
    nv = len(vibrators)
    nh = len(heavies)
    thermochem_data.data = {sp : thermochem_data[sp] for sp in species}
    a_vib = np.array([[221., 221., 180., 180.], [221., 221., 180., 180.]])
    b_vib = np.array([[.029, .029, .0262, .0262], [.029, .029, .0262, .0262]])
    sigma_v_p = np.array([3e-21, 3e-21])

    # -- Plug in expressions that are the same for all species -- #

    # Start with relaxation time of species s
    tau_s = Expression(exprs.tau_s_expr)

    # Plug in total number density
    tau_s.plug_in(syms.n, exprs.n_expr)

    # Plug in mixture mass
    tau_s.plug_in(syms.m_bar, exprs.m_bar_expr)

    # Plug in number of species
    tau_s.plug_in(syms.ns, ns).doit()

    # Plug in masses of species
    tau_s.plug_in(syms.m, {s : thermochem_data[species[s]].m for s in range(ns)})

    # Plug in vibrational relaxation time of species s
    tau_s.plug_in(syms.tau_v_s, exprs.tau_v_s_expr)

    # Plug in vibrational cross section
    tau_s.plug_in(syms.sigma_v, exprs.sigma_v_expr)

    # Plug in constants
    tau_s.plug_in(syms.N_A, constants.N_A)
    tau_s.plug_in(syms.k, constants.k)
    tau_s.plug_in(syms.pi, math.pi)

    # -- Plug in species-specific quantities -- #

    # Start with the expression created so far
    tau_s = {s : {r : tau_s for r in range(nh)} for s in range(nv)}
    # Loop over vibrating species
    for s in range(nv):
        # Loop over heavy particles
        for r in range(nh):
            # Plug in a_vib
            tau_s[s][r].plug_in(syms.a_vib, a_vib[s, r])
            # Plug in b_vib
            tau_s[s][r].plug_in(syms.b_vib, b_vib[s, r])
            # Plug in sigma_v_p
            tau_s[s][r].plug_in(syms.sigma_v_p, sigma_v_p[s])
    breakpoint()


    # Create assignments to get equations
    #e_s_eq = [Assignment(syms.e_s[i], e_s[i].expression) for i in range(ns)]

    # Save to file
    with open(physics_file_name, "wb") as physics_file:
        eqs_to_write = e_s_eq
        pickle.dump(eqs_to_write, physics_file,
                protocol=pickle.HIGHEST_PROTOCOL)
