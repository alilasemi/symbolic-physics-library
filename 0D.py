import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import pickle
import cantera as ct
import ctypes
import scipy

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

    # Plot
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True)

    # Plot T
    fig = plt.figure(figsize=(7,7))
    plt.plot(t, T, lw=3, label='Generated')
    plt.plot(dplr_t, dplr_T, lw=3, label='DPLR')
    plt.xlabel('$t$ (s)', fontsize=20)
    plt.ylabel('$T$ (K)', fontsize=20)
    plt.tick_params(labelsize=20)
    plt.grid(linestyle='--')
    #plt.legend(fontsize=14, ncol=ns, loc='lower center')
    plt.savefig(f'0D_T.png', bbox_inches='tight')

    # labels = ['N2', 'N2+', 'N', 'N+', 'e-']
    # for i in range(chem.ns):
    #     fig = plt.figure(figsize=(7,7))
    #     plt.plot(dplr_T, c_wdot[i], lw=3, label=labels[i])
    #     plt.plot(dplr_T, wdot[i], lw=3, label=labels[i] + ', CT')
    #     plt.plot(dplr_T, dplr_source[:nt, i], lw=3, label=labels[i] + ', DPLR')
    #     plt.xlabel('$T$ (K)', fontsize=20)
    #     #plt.ylabel('$\\frac{\\dot{w} - \\dot{w}_{\\textrm{Cantera}}}{\\dot{w}_{\\textrm{Cantera}}}$', fontsize=20)
    #     plt.ylabel('$\\dot{w}$ (kg/m$^3$/s)', fontsize=20)
    #     #plt.ylabel('$\\dot{w}$', fontsize=20)
    #     plt.tick_params(labelsize=20)
    #     #plt.ylim([1e-3, 1e15])
    #     #plt.xlim([7500, 2e4])
    #     plt.grid(linestyle='--')
    #     plt.legend(fontsize=14, ncol=ns, loc='lower center')
    #     #plt.ylim([1e-89, 1e15])
    #     plt.savefig(f'wdot_{i}.png', bbox_inches='tight')

    plt.show()

def RHS(t, rho, e_in, T_guess):
    # Get total density
    rho_t = np.sum(rho)
    # Get mass fraction
    Y = rho/rho_t
    # Get temperature
    T = get_T_from_e(e_in, Y, T_guess)
    # Get mass production rate
    wdot = get_wdot(T, rho) * 1e-6
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
        print(T)
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
    get_e_cv = c_lib.get_e_cv
    # Set types
    get_e_cv.argtypes = [ctypes.c_double, ctypes.c_void_p,
            ctypes.c_void_p, ctypes.c_void_p]
    get_e_cv.restype = None

    e = np.empty(1)
    cv = np.empty(1)
    get_e_cv(
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
            ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p]
    eval_spec_rates.restype = None

    wdot = np.empty(rho.size)
    Rf = np.empty(3) # Not used?
    Rb = np.empty(3) # Not used?
    c_lib.eval_spec_rates(
            T,
            rho.ctypes.data,
            Rf.ctypes.data,
            Rb.ctypes.data,
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
