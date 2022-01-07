import numpy as np
import pickle
import sympy as sp

from physics.thermodynamics.thermochemical_data import ThermochemicalData

from core.expression import Expression
import physics.constants as constants
import physics.chemistry.expressions as exprs
import physics.thermodynamics.nasa9 as thermo_fit_exprs
import physics.thermodynamics.thermodynamics as thermo_exprs
import physics.chemistry.symbols as syms
import physics.thermodynamics.symbols as thermo_syms


def create(kinetics_data):

    physics_file_name = 'physics.pkl'

    # Unpack
    species = kinetics_data.species
    ns = kinetics_data.ns
    nr = kinetics_data.nr
    alpha = kinetics_data.alpha
    beta = kinetics_data.beta
    epsilon = kinetics_data.epsilon
    nu = kinetics_data.nu
    three_body = kinetics_data.three_body
    C = kinetics_data.C
    eta = kinetics_data.eta
    theta = kinetics_data.theta

    # -- Collect the necessary data -- #

    thermochem_data = ThermochemicalData('NASA9')
    # TODO: Un-hardcode
    thermochem_data.data = {sp : thermochem_data[sp] for sp in species}
    M = np.array([thermochem_data[sp].M for sp in species])

    # -- Plug in -- #

    # Gibbs energy of each species
    G_RT = np.empty(ns, dtype=object)
    for s, sp in enumerate(species):
        G_RT[s] = Expression(thermo_fit_exprs.G_RT_expr)
        G_RT[s].expression = thermo_exprs.create_piecewise_expression_from_fit(G_RT[s], syms.T, thermo_syms.a, thermochem_data[sp].a, thermochem_data[sp].temperatures)

    Rf = np.empty(nr, dtype=object)
    Rb = np.empty(nr, dtype=object)
    # Reactions
    for r in range(nr):
        # kf
        if three_body[r]:
            kf = Expression(exprs.kf_three_body_expr)
            kf.plug_in(syms.kf, exprs.kf_arrhenius_expr)
            kf.plug_in(syms.ns, ns).doit()
            kf.plug_in(syms.M, M)
            kf.plug_in(syms.epsilon, epsilon[r])
        else:
            kf = Expression(exprs.kf_arrhenius_expr)
        kf.plug_in(syms.C, C[r])
        kf.plug_in(syms.eta, eta[r])
        kf.plug_in(syms.theta, theta[r])
        # K_c
        K_c = Expression(exprs.K_c_expr)
        K_c.plug_in(syms.R, constants.R)
        K_c.plug_in(syms.p_0, thermochem_data.p_0)
        K_c.plug_in(syms.nu, nu[r])
        K_c.plug_in(syms.ns, ns).doit()
        for s in range(ns):
            K_c.plug_in(syms.alpha[s, syms.r], alpha[s, r])
            K_c.plug_in(syms.beta[s, syms.r], beta[s, r])
        K_c.plug_in(syms.G_RT, G_RT)
        # TODO: This would be a good (optional) simplification spot
        kb = Expression(exprs.kb_expr)
        kb.plug_in(syms.kf, kf)
        kb.plug_in(syms.K_c, K_c)
        Rf[r] = Expression(exprs.Rf_expr)
        Rb[r] = Expression(exprs.Rb_expr)
        Rf[r].plug_in(syms.ns, ns).doit()
        Rb[r].plug_in(syms.ns, ns).doit()
        for s in range(ns):
            Rf[r].plug_in(syms.alpha[s, syms.r], alpha[s, r])
            Rb[r].plug_in(syms.beta[s, syms.r], beta[s, r])
        Rf[r].plug_in(syms.M, M)
        Rb[r].plug_in(syms.M, M)
        Rf[r].plug_in(syms.kf, kf)
        Rb[r].plug_in(syms.kb, kb)

    # Production
    wdot = [None] * ns
    for s in range(ns):
        wdot[s] = Expression(exprs.wdot_expr)
        wdot[s].plug_in(syms.M[syms.s], M[s])
        wdot[s].plug_in(syms.nr, nr).doit()
        for r in range(nr):
            wdot[s].plug_in(syms.alpha[syms.s, r], alpha[s, r])
            wdot[s].plug_in(syms.beta[syms.s, r], beta[s, r])
        wdot[s].plug_in(syms.Rf, Rf)
        wdot[s].plug_in(syms.Rb, Rb)

    # Save to file
    with open(physics_file_name, "wb") as physics_file:
        exprs_to_write = wdot
        pickle.dump(exprs_to_write, physics_file,
                protocol=pickle.HIGHEST_PROTOCOL)
