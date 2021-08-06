import math
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
    ns = len(species)
    thermochem_data.data = {sp : thermochem_data[sp] for sp in species}

    # Start with relaxation time of species s
    tau_s = Expression(exprs.tau_s_expr)

    # Plug in total number density
    tau_s.plug_in(syms.n, exprs.n_expr)

    # Plug in average molar mass
    tau_s.plug_in(syms.M_bar, exprs.M_bar_expr)

    # Plug in density
    tau_s.plug_in(syms.rho_t, exprs.rho_t_expr)

    # Number of species
    tau_s.plug_in(syms.ns, ns).doit()

    # Plug in molar masses
    tau_s.plug_in(syms.M, {s : thermochem_data[species[s]].M for s in range(ns)})

    # Plug in vibrational relaxation time of species s
    tau_s.plug_in(syms.tau_v_s, exprs.tau_v_s_expr)

    # Plug in vibrational cross section
    tau_s.plug_in(syms.sigma_v, exprs.sigma_v_expr)

    # Plug in constants
    tau_s.plug_in(syms.N_A, constants.N_A)
    tau_s.plug_in(syms.k, constants.k)
    tau_s.plug_in(syms.pi, math.pi)
    breakpoint()


    # Create assignments to get equations
    #e_s_eq = [Assignment(syms.e_s[i], e_s[i].expression) for i in range(ns)]

    # Save to file
    with open(physics_file_name, "wb") as physics_file:
        eqs_to_write = e_s_eq
        pickle.dump(eqs_to_write, physics_file,
                protocol=pickle.HIGHEST_PROTOCOL)
