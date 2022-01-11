import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import pickle
import cantera as ct
import ctypes
import scipy
import h5py
import os
import sys
sys.path.append('../..')

def main():

    dglegion_solution_dir = '/home/ali/projects/cfd_cases/0d_nitrogen/solutions/'
    # TODO: Unhardcode
    ns = 5
    nr = 3
    M_N = .014007
    M = np.array([2*M_N, 2*M_N, M_N, M_N, 5.48579909e-7])

    # Run 0D simulation using Python ODE integrators
    n_t = 2000
    t_final = n_t * 1e-7
    rho_0 = np.array([2.5e-4, 0, 0, 0, 0])
    T_0 = 24000
    Y_0 = np.array([1., 0, 0, 0, 0])
    p_0 = 2.5e-4 * 296.8 * T_0
    e_in, _ = get_e_from_T(T_0, Y_0)
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

    # Create Cantera solution
    gas = ct.Solution('nitrogen.cti')
    # Set initial state
    gas.TPY = (T_0, p_0, Y_0)
    # Create 0D simulation
    reactor = ct.IdealGasReactor(gas)
    sim = ct.ReactorNet([reactor])
    sim.verbose = True
    # limit advance when temperature difference is exceeded
    delta_T_max = 20.
    reactor.set_advance_limit('temperature', delta_T_max)
    dt_max = 1.e-7
    t_end = 2000 * dt_max
    states = ct.SolutionArray(gas, extra=['t'])

    #print('{:10s} {:10s} {:10s} {:14s}'.format(
    #    't [s]', 'T [K]', 'P [Pa]', 'u [J/kg]'))
    while sim.time < t_end:
        sim.advance(sim.time + dt_max)
        states.append(reactor.thermo.state, t=sim.time)
        #print('{:10.3e} {:10.3f} {:10.3f} {:14.6f}'.format(
        #        sim.time, reactor.T, reactor.thermo.P, reactor.thermo.u))

    # Plot
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True)

    figsize = 7
    # Plot T
    fig = plt.figure(figsize=(figsize, figsize))
    plt.plot(states.t, states.T, lw=3, label='Cantera')
    plt.plot(t, T, lw=3, label='Sympy')
    plt.xlabel('$t$ (s)', fontsize=20)
    plt.ylabel('$T$ (K)', fontsize=20)
    plt.tick_params(labelsize=20)
    plt.grid(linestyle='--')
    plt.legend(fontsize=14, ncol=1, loc='upper right')
    plt.savefig(f'0D_T.png', bbox_inches='tight')

    # Plot Y
    labels = ['N2', 'N2+', 'N', 'N+', 'e-']
    for i in range(ns):
        fig = plt.figure(figsize=(figsize, figsize))
        plt.plot(states.t, states.Y[:, i], lw=3, label=labels[i] + ', Cantera')
        plt.plot(t, Y[:, i], lw=3, label=labels[i] + ', Sympy')
        plt.xlabel('$t$ (s)', fontsize=20)
        plt.ylabel('$Y$', fontsize=20)
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
    c_lib = ctypes.CDLL('./generated/e_and_cv.so')
    compute_e_and_cv = c_lib.compute_e_and_cv
    # Set types
    compute_e_and_cv.argtypes = [ctypes.c_double, ctypes.c_void_p,
            ctypes.c_void_p, ctypes.c_void_p]
    compute_e_and_cv.restype = None

    e = np.empty(1)
    cv = np.empty(1)
    compute_e_and_cv(
            T,
            Y.ctypes.data,
            e.ctypes.data,
            cv.ctypes.data)
    return e[0], cv[0]

def get_wdot(T, rho):
    '''Compute the mass production rate.'''

    # Load C library
    c_lib = ctypes.CDLL('./generated/wdot.so')
    compute_wdot = c_lib.compute_wdot
    # Set types
    compute_wdot.argtypes = [ctypes.c_double, ctypes.c_void_p,
            ctypes.c_void_p]
    compute_wdot.restype = None

    wdot = np.empty(rho.size)
    c_lib.compute_wdot(
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
