import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import pickle
import cantera as ct
import ctypes

from chemical_mechanism import ChemicalMechanism

def main():

    thermo_file_name = 'N2_ions.dat'
    chem = ChemicalMechanism(thermo_file_name)
    ns = chem.ns
    nr = chem.nr
    M = chem.M

    # Read data from DPLR
    dplr_state = np.loadtxt('reference_data/kinetics/reactor_n2_tceq.state.dat',
            skiprows=2)
    dplr_t = dplr_state[:, 0]
    nt = dplr_t.shape[0]
    dplr_T = dplr_state[:, -4]
    dplr_rho_t = dplr_state[:, -3]
    dplr_X = dplr_state[:, 1:ns+1] / dplr_state[:, [-2]] # n_s divided by n_tot
    dplr_Mbar = np.dot(dplr_X, M)
    # Convert mole fraction to mass fraction
    dplr_Y = np.empty_like(dplr_X)
    for i in range(dplr_state.shape[0]):
        dplr_Y[i] = dplr_X[i] * M / dplr_Mbar[i]
    # Convert mass fraction to partial density
    dplr_rho = dplr_Y * dplr_rho_t[:, np.newaxis]
    dplr_source = np.loadtxt('reference_data/kinetics/reactor_n2_tceq.source.dat',
            skiprows=2)[:, 1:ns+1]

    # Hack to use same composition and linear T
    #dplr_T = np.linspace(300, 20000, nt)
    #dplr_rho_t[:] = .00025
    #dplr_Y[:] = np.array([.29, .15, .21, .23, .12])
    #dplr_rho[:] = dplr_rho_t[0] * dplr_Y[0]

    # Load C library
    c_lib = ctypes.CDLL('./rxn_rates.so')
    # Set types
    c_lib.eval_spec_rates.argtypes = [ctypes.c_double, ctypes.c_void_p,
            ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
    c_lib.eval_spec_rates.restype = None
    # Call function for each T
    c_wdot = np.empty((ns, nt))
    Rf = np.empty((nr, nt))
    Rb = np.empty((nr, nt))
    for i in range(nt):
        c_wdot_i = np.empty(ns)
        Rf_i = np.empty(nr)
        Rb_i = np.empty(nr)
        rho = dplr_rho[i]
        c_lib.eval_spec_rates(
                dplr_T[i],
                rho.ctypes.data,
                Rf_i.ctypes.data,
                Rb_i.ctypes.data,
                c_wdot_i.ctypes.data)
        c_wdot[:, i] = c_wdot_i
        Rf[:, i] = Rf_i
        Rb[:, i] = Rb_i

    G_RT = np.empty((ns, nt))
    H_RT = np.empty((ns, nt))
    S_R = np.empty((ns, nt))
    kf = np.empty((nr, nt))
    kb = np.empty((nr, nt))
    wdot = np.empty((ns, nt))
    # Load gas and chemical mechanism with Cantera
    gas = ct.Solution('nitrogen.cti')
    # Loop over temperatures
    for i in range(nt):
        # Set temperature, density, and mass fractions
        gas.TDY = dplr_T[i], dplr_rho_t[i], dplr_Y[i]

        # Compute thermodynamics for each species
        G_RT[:, i] = gas.standard_gibbs_RT
        H_RT[:, i] = gas.standard_enthalpies_RT
        S_R[:, i] = gas.standard_entropies_R

        # Mass production rate of each species
        wdot[:, i] = gas.molecular_weights * gas.net_production_rates

    # Plot
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True)

    #fig = plt.figure(figsize=(7,7))
    #labels = ['N2', 'N2+', 'N', 'N+', 'e-']
    #for i in range(chem.ns):
    #    p = plt.plot(T, sp.lambdify('T', chem.species[i].G_RT)(T),
    #            lw=3, label=labels[i])
    #    color = p[0].get_color()
    #    plt.plot(T, G_RT[i], '--', color=color, lw=3, label=labels[i] + ', Cantera')
    #plt.xlabel('$T$', fontsize=20)
    #plt.ylabel('$\\frac{G^0_s}{RT}$', fontsize=20)
    #plt.tick_params(labelsize=20)
    #plt.grid(linestyle='--')
    #plt.legend(fontsize=14, ncol=2, loc='upper center')
    #plt.savefig('G_RT.png', bbox_inches='tight')

    ## Plot enthalpies
    #fig = plt.figure(figsize=(7,7))
    #labels = ['N2', 'N2+', 'N', 'N+', 'e-']
    #for i in range(chem.ns):
    #    p = plt.plot(T, sp.lambdify('T', chem.species[i].H_RT)(T),
    #            lw=3, label=labels[i])
    #    color = p[0].get_color()
    #    plt.plot(T, H_RT[i], '--', color=color, lw=3, label=labels[i] + ', Cantera')
    #plt.xlabel('$T$', fontsize=20)
    #plt.ylabel('$\\frac{H^0_s}{RT}$', fontsize=20)
    #plt.tick_params(labelsize=20)
    #plt.grid(linestyle='--')
    #plt.legend(fontsize=14, ncol=2, loc='upper center')
    #plt.savefig('H_RT.png', bbox_inches='tight')

    ## Plot entropies
    #fig = plt.figure(figsize=(7,7))
    #labels = ['N2', 'N2+', 'N', 'N+', 'e-']
    #for i in range(chem.ns):
    #    p = plt.plot(T, sp.lambdify('T', chem.species[i].S_R)(T),
    #            lw=3, label=labels[i])
    #    color = p[0].get_color()
    #    plt.plot(T, S_R[i], '--', color=color, lw=3, label=labels[i] + ', Cantera')
    #plt.xlabel('$T$', fontsize=20)
    #plt.ylabel('$\\frac{S^0_s}{R}$', fontsize=20)
    #plt.ylim([-8, 47])
    #plt.tick_params(labelsize=20)
    #plt.grid(linestyle='--')
    #plt.legend(fontsize=14, ncol=2, loc='lower center')
    #plt.savefig('S_R.png', bbox_inches='tight')

    # Plot kf
    #fig = plt.figure(figsize=(7,7))
    #labels = ['$N_2 + M$ --- $2N + M$',
    #        '$N + e^-$ --- $N^+ + 2e^-$',
    #        '$2N$ --- $N_2^+ + e^-$']
    #for i in range(chem.nr):
    #    #p = plt.semilogy(T, sp.lambdify('T', chem.reactions[i].kf)(T),
    #    #        lw=3, label=labels[i])
    #    #color = p[0].get_color()
    #    #plt.semilogy(T, kf[i], '--',
    #    #        color=color, lw=3, label=labels[i])
    #    plt.semilogy(T, (sp.lambdify('T', chem.reactions[i].kf)(T) - kf[i])/kf[i],
    #            lw=3, label=labels[i])
    #plt.xlabel('$T$', fontsize=20)
    #plt.ylabel('$\\frac{k_f - k_{f, \\textrm{ Cantera}}}{k_{f, \\textrm{ Cantera}}}$', fontsize=20)
    #plt.tick_params(labelsize=20)
    #plt.grid(linestyle='--')
    #plt.legend(fontsize=14, ncol=1, loc='center right')
    ##plt.ylim([1e-41, 1e15])
    #plt.ylim([1e-20, 1e4])
    #plt.savefig('kf.png', bbox_inches='tight')

    ## Plot kb
    #fig = plt.figure(figsize=(7,7))
    #for i in range(chem.nr):
    #    #p = plt.semilogy(T, sp.lambdify('T', chem.reactions[i].kb)(T),
    #    #        lw=3, label=labels[i])
    #    #color = p[0].get_color()
    #    #plt.semilogy(T, kb[i], '--',
    #    #        color=color, lw=3, label=labels[i])
    #    plt.semilogy(T, (sp.lambdify('T', chem.reactions[i].kb)(T) - kb[i])/kb[i],
    #            lw=3, label=labels[i])
    #plt.xlabel('$T$', fontsize=20)
    #plt.ylabel('$\\frac{k_b - k_{b, \\textrm{ Cantera}}}{k_{b, \\textrm{ Cantera}}}$', fontsize=20)
    #plt.tick_params(labelsize=20)
    #plt.grid(linestyle='--')
    #plt.legend(fontsize=14, ncol=1, loc='center right')
    ##plt.ylim([1e-89, 1e15])
    #plt.savefig('kb.png', bbox_inches='tight')

    # Plot wdot
    labels = ['N2', 'N2+', 'N', 'N+', 'e-']
    for i in range(chem.ns):
        fig = plt.figure(figsize=(7,7))
        plt.plot(dplr_T, c_wdot[i], lw=3, label=labels[i])
        plt.plot(dplr_T, wdot[i], lw=3, label=labels[i] + ', CT')
        plt.plot(dplr_T, dplr_source[:nt, i], lw=3, label=labels[i] + ', DPLR')
        plt.xlabel('$T$ (K)', fontsize=20)
        #plt.ylabel('$\\frac{\\dot{w} - \\dot{w}_{\\textrm{Cantera}}}{\\dot{w}_{\\textrm{Cantera}}}$', fontsize=20)
        plt.ylabel('$\\dot{w}$ (kg/m$^3$/s)', fontsize=20)
        #plt.ylabel('$\\dot{w}$', fontsize=20)
        plt.tick_params(labelsize=20)
        #plt.ylim([1e-3, 1e15])
        #plt.xlim([7500, 2e4])
        plt.grid(linestyle='--')
        plt.legend(fontsize=14, ncol=ns, loc='lower center')
        #plt.ylim([1e-89, 1e15])
        plt.savefig(f'wdot_{i}.png', bbox_inches='tight')

    plt.show()


if __name__ == "__main__":
    main()
