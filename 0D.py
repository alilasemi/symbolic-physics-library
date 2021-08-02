import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import pickle
import cantera as ct
import ctypes
import scipy
import h5py
import os

from chemical_mechanism import ChemicalMechanism

def main():

    dglegion_solution_dir = '/home/ali/projects/cfd_cases/0d_nitrogen/solutions/'
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
            skiprows=2)[:, 1:ns+1] * 1e6 # DPLR gives kg/cm^3/s, convert to m

    n_t = 2000
    t_final = nt * 1e-7
    rho_0 = dplr_rho[0]
    T_0 = dplr_T[0]
    e_in, _ = get_e_from_T(dplr_T[0], dplr_Y[0])
    t, rho = RK4(RHS, n_t, t_final, rho_0, T_0, args=(e_in,))

    # Get total density
    rho_t = np.sum(rho, axis=1, keepdims=True)
    # Get mass fraction
    Y = rho/rho_t
    # Get temperature
    T = np.empty_like(t)
    T_guess = T_0
    for i in range(t.size):
        T[i] = get_T_from_e(e_in, Y[i], T_guess)
        T_guess = T[i]

    # -- Read data from DG-Legion code -- #
    # Get list of all files in solution directory
    _, _, filenames = next(os.walk(dglegion_solution_dir))
    # Keep only the HDF5 files
    filenames = [name for name in filenames if name.endswith('.h5')]
    # Prepend the directory path
    filenames = [dglegion_solution_dir + name for name in filenames]
    # Iteration number of each file
    iters = [name.split('_000.h5')[0].split('_')[-1] for name in filenames]
    # Sort filenames by iteration number
    filenames = [name for _, name in sorted(zip(iters, filenames))]
    # Loop over each solution file
    dgl_t = np.empty(len(filenames))
    dgl_T = np.empty_like(dgl_t)
    dgl_Y = np.empty((dgl_t.size, ns))
    for i, filename in enumerate(filenames):
        # Open file
        h5file = h5py.File(filename, 'r')
        # Read time
        dgl_t[i] = h5file['Time'][()][0]
        # Read temperature
        dgl_T[i] = h5file['T'][()][0, 0]
        # Read mass fractions
        for s in range(ns):
            dgl_Y[i, s] = h5file[f'Y{s}'][()][0, 0]

    breakpoint()
    # Plot
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True)

    # Plot T
    fig = plt.figure(figsize=(7,7))
    plt.plot(dplr_t, dplr_T, lw=3, label='DPLR')
    plt.plot(t, T, lw=3, label='Python')
    plt.plot(dgl_t, dgl_T, lw=3, label='DG-Legion')
    plt.xlabel('$t$ (s)', fontsize=20)
    plt.ylabel('$T$ (K)', fontsize=20)
    plt.tick_params(labelsize=20)
    plt.grid(linestyle='--')
    plt.legend(fontsize=14, ncol=1, loc='upper right')
    plt.savefig(f'0D_T.png', bbox_inches='tight')

    # Plot Y
    labels = ['N2', 'N2+', 'N', 'N+', 'e-']
    for i in range(chem.ns):
        fig = plt.figure(figsize=(7,7))
        plt.plot(dplr_t, dplr_Y[:, i], lw=3, label=labels[i] + ', DPLR')
        plt.plot(t, Y[:, i], lw=3, label=labels[i] + ', Python')
        plt.plot(dgl_t, dgl_Y[:, i], lw=4, label=labels[i] + ', DG-Legion')
        plt.xlabel('$t$ (s)', fontsize=20)
        plt.ylabel('$\\rho$ (kg/m$^3$)', fontsize=20)
        plt.tick_params(labelsize=20)
        #plt.ylim([1e-3, 1e15])
        #plt.xlim([7500, 2e4])
        plt.grid(linestyle='--')
        plt.legend(fontsize=14, loc='lower center')
        #plt.ylim([1e-89, 1e15])
        plt.savefig(f'0D_sp{i}.png', bbox_inches='tight')

    plt.show()

def RHS(t, rho, e_in, T_guess):
    # Get total density
    rho_t = np.sum(rho)
    # Get mass fraction
    Y = rho/rho_t
    # Get temperature
    T = get_T_from_e(e_in, Y, T_guess)
    # Get mass production rate
    wdot = get_wdot(T, rho)
    return wdot, T

def get_T_from_e(e_in, Y, T_guess=300.):
    '''Use Newton iterations to invert the energy fit.'''

    # Solver settings
    max_iter = 300
    reltol = 1e-8
    T = T_guess # Initial guess

    # Newton iterations
    success = False
    for i in range(max_iter):
        e, cv = get_e_from_T(T, Y)
        delta = (e_in - e) / cv
        T += delta
        # Check for convergence
        if delta < T * reltol:
            success = True
            break
    # Error messages in case of failure
    if not success:
        print("Convergence failure in the Newton solve for temperature.")
        print("T      = ", T)
        print("delta  = ", delta)
        print("e_in   = ", e_in)
        print("e      = ", e)
        print("relerr = ", delta/T)
    return T

def get_e_from_T(T, Y):
    # Load C library
    c_lib = ctypes.CDLL('./energy.so')
    get_e_and_cv = c_lib.get_e_and_cv
    # Set types
    get_e_and_cv.argtypes = [ctypes.c_double, ctypes.c_void_p,
            ctypes.c_void_p, ctypes.c_void_p]
    get_e_and_cv.restype = None

    e = np.empty(1)
    cv = np.empty(1)
    get_e_and_cv(
            T,
            Y.ctypes.data,
            e.ctypes.data,
            cv.ctypes.data)
    return e[0], cv[0]

def get_wdot(T, rho):
    '''Compute the mass production rate.'''

    # Load C library
    c_lib = ctypes.CDLL('./rxn_rates.so')
    eval_spec_rates = c_lib.eval_spec_rates
    # Set types
    eval_spec_rates.argtypes = [ctypes.c_double, ctypes.c_void_p,
            ctypes.c_void_p]
    eval_spec_rates.restype = None

    wdot = np.empty(rho.size)
    c_lib.eval_spec_rates(
            T,
            rho.ctypes.data,
            wdot.ctypes.data)
    return wdot

def RK4(f, n_t, t_final, u_0, guess, args=()):

    t = np.linspace(0, t_final, n_t + 1)
    dt = t[1] - t[0]
    u = np.zeros((t.size, u_0.size))
    u[0] = u_0

    # For each time step
    for i in range(n_t):
        # Perform the four stages
        K1, guess = f(t[i], u[i],         *args, guess)
        K1 *= dt
        K2, guess = f(t[i], u[i] + .5*K1, *args, guess)
        K2 *= dt
        K3, guess = f(t[i], u[i] + .5*K2, *args, guess)
        K3 *= dt
        K4, guess = f(t[i], u[i] + K3,    *args, guess)
        K4 *= dt

        u[i+1] = u[i] + (K1 + 2*K2 + 2*K3 + K4) / 6

    return t, u


if __name__ == "__main__":
    main()
