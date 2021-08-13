import copy
import numpy as np
import pickle
import sympy as sp

from physics.thermodynamics.thermochemical_data import ThermochemicalData

import expressions.landau_teller_source.expressions as exprs
import expressions.landau_teller_source.symbols as syms
from core.expression import Expression
import physics.constants as constants


def create(e_s_vee, p):

    # Collect the necessary data
    thermochem_data = ThermochemicalData('NASA9')
    # TODO: Un-hardcode
    species = ['N2', 'N2+', 'N', 'N+', 'e-']
    vibrators = ['N2', 'N2+']
    heavies = ['N2', 'N2+', 'N', 'N+']
    vib_to_species = [0, 1]
    ns = len(species)
    nv = len(vibrators)
    nh = len(heavies)
    thermochem_data.data = {sp : thermochem_data[sp] for sp in species}
    M = np.array([thermochem_data[sp].M for sp in species])
    a_vib = np.array([[221., 221., 180., 180.], [221., 221., 180., 180.]])
    b_vib = np.array([[.029, .029, .0262, .0262], [.029, .029, .0262, .0262]])
    sigma_v_p = np.array([3e-21, 3e-21])

    # -- Plug in expressions that are the same for all species -- #

    # Start with relaxation time of species s
    tau = Expression(exprs.tau_expr)

    # Plug in total number density
    tau.plug_in(syms.n, exprs.n_expr)

    # Plug in mixture mass
    tau.plug_in(syms.m_bar, exprs.m_bar_expr)

    # Plug in number of species
    tau.plug_in(syms.ns, ns).doit()

    # Plug in masses of species
    tau.plug_in(syms.m, {s : thermochem_data[species[s]].m for s in range(ns)})

    # Plug in vibrational relaxation time of species s
    tau.plug_in(syms.tau_v_s, exprs.tau_v_s_expr)

    # Plug in vibrational cross section
    tau.plug_in(syms.sigma_v, exprs.sigma_v_expr)

    # Plug in constants
    tau.plug_in(syms.N_A, constants.N_A)
    tau.plug_in(syms.k, constants.k)
    tau.plug_in(syms.pi, np.pi)

    # -- Plug in species-specific quantities -- #

    # Start with the expression created so far
    tau = {s : {r : copy.deepcopy(tau) for r in range(nh)} for s in range(nv)}
    # Loop over vibrating species
    for s in range(nv):
        # Loop over heavy particles
        for r in range(nh):
            # Plug in a_vib
            tau[s][r].plug_in(syms.a_vib, a_vib[s, r])
            # Plug in b_vib
            tau[s][r].plug_in(syms.b_vib, b_vib[s, r])
            # Plug in sigma_v_p
            tau[s][r].plug_in(syms.sigma_v_p, sigma_v_p[s])

    # Loop over vibrating species
    Q_TV_s = np.empty(nv, dtype=object)
    for s in range(nv):
        # tau_s
        tau_s = Expression(exprs.tau_s_expr)
        tau_s.plug_in(syms.nh, nh).doit()
        tau_s.plug_in(syms.M, M)
        for r in range(nh):
            tau_s.plug_in(syms.tau[syms.s, r], tau[s][r])
        # Q_TV
        Q_TV_s[s] = Expression(exprs.Q_TV_s_expr)
        Q_TV_s[s].plug_in(syms.s, s)
        Q_TV_s[s].plug_in(syms.tau_s[s], tau_s)
        # TODO: unhack this
        Q_TV_s[s].plug_in(syms.e_s_vee[s], e_s_vee[vib_to_species[s]])
        # TODO: unhack this
        Q_TV_s[s].plug_in(syms.e_s_vee_star[s],
                e_s_vee[vib_to_species[s]].expression.subs(syms.T, syms.Tv))

    Q_TV = Expression(exprs.Q_TV_expr)
    Q_TV.plug_in(syms.nv, nv).doit()
    Q_TV.plug_in(syms.Q_TV_s, Q_TV_s)
    Q_TV.plug_in(syms.p, p)

    # Save to file
    physics_file_name = 'physics.pkl'
    with open(physics_file_name, "wb") as physics_file:
        eqs_to_write = Q_TV
        pickle.dump(eqs_to_write, physics_file,
                protocol=pickle.HIGHEST_PROTOCOL)
