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
    n_t = 200
    t_final = nt * 1e-7 * 1e-7
    rho_0 = dplr_rho[i_start]
    T_0 = dplr_T[i_start]
    Tv_0 = 1000
    e_vee_in = get_e_vee_from_Tv(Tv_0, dplr_Y[i_start])
    e_in = get_e_from_T(T_0, dplr_Y[i_start])
    u_0 = np.append(rho_0, e_vee_in)
    t, u = RK4(RHS, n_t, t_final, u_0, Tv_0, args=(e_in,))
    rho = u[:, :-1]
    e_vee = u[:, -1]

    # Get total density
    rho_t = np.sum(rho, axis=1, keepdims=True)
    # Get mass fraction
    Y = rho/rho_t
    # Get temperature
    T = np.empty_like(t)
    Tv = np.empty_like(t)
    Tv_guess = Tv_0
    for i in range(t.size):
        Tv[i] = get_Tv_from_e_vee(e_vee[i], Y[i], Tv_guess)
        Tv_guess = Tv[i]
        e_tr = get_e_tr_from_e(e_in, e_vee[i], Y[i])
        T[i] = get_T_from_e_tr(e_tr, Y[i])

    # Plot
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True)

    Q_TV = np.empty((5, 100))
    T_const = [1e4, 1.5e4, 2e4, 2.5e4, 3e4]
    for j in range(5):
        T = T_const[j] * np.ones(100)
        Tv = T - np.linspace(0, .9 * T_const[j], 100)
        rho = np.array([.0025, 0, 0, 0, 0])
        Y = np.array([1., 0, 0, 0, 0])
        for i in range(100):
            Q_TV[j, i] = get_Q_TV(T[i], Tv[i], rho, Y)
    # Plot LT source
    fig = plt.figure(figsize=(7,7))
    #plt.plot(dplr_t, dplr_T, lw=3, label='DPLR')
    for j in range(5):
        plt.plot(T - Tv, Q_TV[j], lw=3, label=f'T = {int(T_const[j])} K')
    plt.xlabel('$T - T_v$ (K)', fontsize=20)
    plt.ylabel('$Q_{T-V}$ (W/m$^3$)', fontsize=20)
    plt.tick_params(labelsize=20)
    plt.grid(linestyle='--')
    plt.legend(fontsize=14, ncol=1)#, loc='upper right')
    plt.savefig(f'figs/0D_Q_LT.pdf', bbox_inches='tight')
    plt.show()

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

def RHS(t, u, e_in, Tv_guess):
    # Unpack state
    rho = u[:-1]
    e_vee = u[-1]
    # Get total density
    rho_t = np.sum(rho)
    # Get mass fraction
    Y = rho/rho_t
    # Get TR temperature
    e_tr = get_e_tr_from_e(e_in, e_vee, Y)
    T = get_T_from_e_tr(e_tr, Y)
    # Get VEE temperature
    print(e_vee, Tv_guess)
    Tv = get_Tv_from_e_vee(e_vee, Y, Tv_guess)
    # Get species VEE energies
    e_s_vee = get_e_s_vee(Tv, len(rho))
    # Get mass production rate
    wdot = get_wdot(Tv, Tv, rho)
    # Get VEE source term
    Q_TV = get_Q_TV(T, Tv, rho, Y)
    # Combine
    udot = np.append(wdot, Q_TV)
    return udot, Tv

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

def get_e_tr_from_e(e, e_vee, Y):
    # Load C library
    c_lib = ctypes.CDLL('./generated/e_tr_from_e.so')
    func = c_lib.compute_e_tr_from_e
    # Set types
    func.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_void_p,
            ctypes.c_void_p]
    func.restype = None

    e_tr = np.empty(1)
    func(
            e,
            e_vee,
            Y.ctypes.data,
            e_tr.ctypes.data)

    return e_tr[0]

def get_e_from_T(T, Y):
    # Load C library
    c_lib = ctypes.CDLL('./generated/e.so')
    func = c_lib.compute_e
    # Set types
    func.argtypes = [ctypes.c_double, ctypes.c_void_p,
            ctypes.c_void_p]
    func.restype = None

    e = np.empty(1)
    func(
            T,
            Y.ctypes.data,
            e.ctypes.data)

    return e[0]

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

def get_e_s_vee(Tv, ns):
    # Load C library
    c_lib = ctypes.CDLL('./generated/e_s_vee.so')
    func = c_lib.compute_e_s_vee
    # Set types
    func.argtypes = [ctypes.c_double,
            ctypes.c_void_p]
    func.restype = None

    e_s_vee = np.empty(ns)
    func(
            Tv,
            e_s_vee.ctypes.data)

    return e_s_vee

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

def get_Q_TV(T, Tv, rho, Y):

    # Load C library
    c_lib = ctypes.CDLL('./generated/Q_TV.so')
    func = c_lib.compute_Q_TV
    # Set types
    func.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_void_p,
            ctypes.c_void_p, ctypes.c_void_p]
    func.restype = None

    Q_TV = np.empty(1)
    func(
            T,
            Tv,
            rho.ctypes.data,
            Y.ctypes.data,
            Q_TV.ctypes.data)
    return Q_TV[0]

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

        u[i + 1] = u[i] + (K1 + 2*K2 + 2*K3 + K4) / 6

    return t, u


if __name__ == "__main__":
    main()
