import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import pickle
import ctypes
import scipy
import h5py
import os

from physics.thermodynamics.thermochemical_data import ThermochemicalData


def main():

    physics_file_name = 'physics.pkl'

    # Collect the necessary data
    thermochem_data = ThermochemicalData('NASA9')
    species = ['N2', 'N2+', 'N', 'N+', 'e-']
    ns = len(species)
    thermochem_data.data = {sp : thermochem_data[sp] for sp in species}
    M = np.array([thermochem_data[sp].M for sp in species])

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

    i_start = 25
    n_t = 2000
    t_final = nt * 1e-7 * 2e-8
    rho_0 = dplr_rho[i_start]
    T_0 = dplr_T[i_start]
    Tv_0 = T_0
    e_vee_in = get_e_vee_from_Tv(Tv_0, dplr_Y[i_start])
    e_tr_in = get_e_tr_from_T(T_0, dplr_Y[i_start])
    e_in = e_vee_in + e_tr_in
    t, rho = RK4(RHS, n_t, t_final, rho_0, Tv_0, args=(e_in, e_vee_in,))

    # Get total density
    rho_t = np.sum(rho, axis=1, keepdims=True)
    # Get mass fraction
    Y = rho/rho_t
    # Get temperature
    T = np.empty_like(t)
    Tv = np.empty_like(t)
    Tv_guess = Tv_0
    for i in range(t.size):
        Tv[i] = get_Tv_from_e_vee(e_vee_in, Y[i], Tv_guess)
        Tv_guess = Tv[i]
        T[i] = get_T_from_e_tr(e_tr_in, Y[i])

    # Plot
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True)

    # Plot T
    fig = plt.figure(figsize=(7,7))
    #plt.plot(dplr_t, dplr_T, lw=3, label='DPLR')
    plt.plot(t, T, 'r', lw=3, label='T')
    plt.plot(t, Tv, 'b', lw=3, label='Tv')
    plt.xlabel('$t$ (s)', fontsize=20)
    plt.ylabel('$T$, $T_v$ (K)', fontsize=20)
    plt.tick_params(labelsize=20)
    plt.grid(linestyle='--')
    plt.legend(fontsize=14, ncol=1, loc='upper right')
    plt.savefig(f'figs/0D_TTv.png', bbox_inches='tight')

    # Plot Y
    labels = ['N2', 'N2+', 'N', 'N+', 'e-']
    for i in range(ns):
        fig = plt.figure(figsize=(7,7))
        #plt.plot(dplr_t, dplr_Y[:, i], lw=3, label=labels[i] + ', DPLR')
        plt.plot(t, Y[:, i], lw=3, label=labels[i] + ', Python')
        #plt.plot(dgl_t, dgl_Y[:, i], lw=4, label=labels[i] + ', DG-Legion')
        plt.xlabel('$t$ (s)', fontsize=20)
        plt.ylabel('$Y$', fontsize=20)
        plt.tick_params(labelsize=20)
        #plt.ylim([1e-3, 1e15])
        #plt.xlim([7500, 2e4])
        plt.grid(linestyle='--')
        #plt.legend(fontsize=14, loc='lower center')
        plt.savefig(f'figs/0D_TTv_sp{i}.png', bbox_inches='tight')

    plt.show()

def RHS(t, rho, e_in, e_vee_in, Tv_guess):
    # Get total density
    rho_t = np.sum(rho)
    # Get mass fraction
    Y = rho/rho_t
    # Get TR temperature
    T = get_T_from_e_tr(e_in - e_vee_in, Y)
    # Get VEE temperature
    Tv = get_Tv_from_e_vee(e_vee_in, Y, Tv_guess)
    # Get mass production rate
    wdot = get_wdot(T, Tv, rho)
    return wdot, Tv

def get_Tv_from_e_vee(e_vee_in, Y, Tv_guess=300.):
    '''Use Newton iterations to invert the energy fit.'''

    # Solver settings
    max_iter = 300
    reltol = 1e-8
    Tv = Tv_guess # Initial guess

    # Newton iterations
    success = False
    for i in range(max_iter):
        e_vee = get_e_vee_from_Tv(Tv, Y)
        cv_vee = get_cv_vee_from_Tv(Tv, Y)
        delta = (e_vee_in - e_vee) / cv_vee
        Tv += delta
        # Check for convergence
        if delta < Tv * reltol:
            success = True
            break
    # Error messages in case of failure
    if not success:
        print("Convergence failure in the Newton solve for temperature.")
        print("Tv     = ", Tv)
        print("delta  = ", delta)
        print("e_in   = ", e_vee_in)
        print("e      = ", e_vee)
        print("relerr = ", delta/Tv)
    return Tv

def get_T_from_e_tr(e_tr, Y):
    '''No Newton iteration needed for T, since TR is assumed fully excited.'''
    # This is not a function of temperature...just plugged in 0 for now
    cv_tr = get_cv_tr_from_T(0, Y)
    return e_tr / cv_tr

def get_e_vee_from_Tv(Tv, Y):
    # Load C library
    c_lib = ctypes.CDLL('./generated/e_vee.so')
    func = c_lib.compute_e_vee
    # Set types
    func.argtypes = [ctypes.c_double, ctypes.c_void_p,
            ctypes.c_void_p]
    func.restype = None

    e_vee = np.empty(1)
    func(
            Tv,
            Y.ctypes.data,
            e_vee.ctypes.data)

    return e_vee[0]

def get_cv_vee_from_Tv(Tv, Y):
    # Load C library
    c_lib = ctypes.CDLL('./generated/cv_vee.so')
    func = c_lib.compute_cv_vee
    # Set types
    func.argtypes = [ctypes.c_double, ctypes.c_void_p,
            ctypes.c_void_p]
    func.restype = None

    cv_vee = np.empty(1)
    func(
            Tv,
            Y.ctypes.data,
            cv_vee.ctypes.data)

    return cv_vee[0]

def get_e_tr_from_T(T, Y):
    # Load C library
    c_lib = ctypes.CDLL('./generated/e_tr.so')
    func = c_lib.compute_e_tr
    # Set types
    func.argtypes = [ctypes.c_double, ctypes.c_void_p,
            ctypes.c_void_p]
    func.restype = None

    e_tr = np.empty(1)
    func(
            T,
            Y.ctypes.data,
            e_tr.ctypes.data)

    return e_tr[0]

def get_cv_tr_from_T(T, Y):
    # Load C library
    c_lib = ctypes.CDLL('./generated/cv_tr.so')
    func = c_lib.compute_cv_tr
    # Set types
    func.argtypes = [ctypes.c_double, ctypes.c_void_p,
            ctypes.c_void_p]
    func.restype = None

    cv_tr = np.empty(1)
    func(
            T,
            Y.ctypes.data,
            cv_tr.ctypes.data)

    return cv_tr[0]

def get_wdot(T, Tv, rho):
    '''Compute the mass production rate.'''

    # Load C library
    c_lib = ctypes.CDLL('./generated/wdot.so')
    func = c_lib.compute_wdot
    # Set types
    func.argtypes = [ctypes.c_double, ctypes.c_void_p,
            ctypes.c_void_p]
    func.restype = None

    wdot = np.empty(rho.size)
    func(
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
        breakpoint()
        K2, guess = f(t[i], u[i] + .5*K1, *args, guess)
        K2 *= dt
        K3, guess = f(t[i], u[i] + .5*K2, *args, guess)
        K3 *= dt
        K4, guess = f(t[i], u[i] + K3,    *args, guess)
        K4 *= dt

        u[i + 1] = u[i] + (K1 + 2*K2 + 2*K3 + K4) / 6

    return t, u


if __name__ == "__main__":
    main()
