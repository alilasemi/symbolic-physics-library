import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import pickle
import sympy as sp

import physics.thermodynamics.symbols as syms

def main():

    physics_file_name = 'physics.pkl'
    species = ['N2', 'N2+', 'N', 'N+', 'e-']
    ns = len(species)

    with open(physics_file_name, "rb") as physics_file:
        e_s_expr = pickle.load(physics_file)

    nt = 100
    T = np.linspace(300, 20000-1, nt)
    e_s = np.empty((ns, nt))
    for s in range(ns):
        e_s[s] = sp.lambdify(syms.T, e_s_expr[s])(T)

    # Plot
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True)

    for s in range(ns):
        fig = plt.figure(figsize=(7,7))
        plt.plot(T, e_s[s], lw=3, label=species[s] + ', $e$')
        plt.xlabel('$T$', fontsize=20)
        plt.ylabel('$e_{s}$', fontsize=20)
        plt.tick_params(labelsize=20)
        plt.grid(linestyle='--')
        plt.legend(fontsize=14, ncol=2, loc='upper center')
        plt.savefig(f'e_s_{s}.png', bbox_inches='tight')

    plt.show()


if __name__ == "__main__":
    main()
