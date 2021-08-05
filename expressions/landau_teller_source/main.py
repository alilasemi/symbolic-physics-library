import numpy as np
import pickle
import sympy as sp
from sympy.codegen.ast import Assignment

import expressions as exprs
import symbols as syms


def main():

    # Start with relaxation time of species s
    tau_s = exprs.tau_s_expr

    # Plug in total number density
    tau_s = tau_s.subs(syms.n, exprs.n_expr)

    # Plug in vibrational relaxation time of species s
    tau_s = tau_s.subs(syms.tau_v_s, exprs.tau_v_s_expr)

    # Plug in vibrational cross section
    tau_s = tau_s.subs(syms.sigma_v, exprs.sigma_v_expr)
    breakpoint()


    # Create assignments to get equations
    #e_s_eq = [Assignment(syms.e_s[i], e_s[i].expression) for i in range(ns)]

    # Save to file
    with open(physics_file_name, "wb") as physics_file:
        eqs_to_write = e_s_eq
        pickle.dump(eqs_to_write, physics_file,
                protocol=pickle.HIGHEST_PROTOCOL)


if __name__ == "__main__":
    main()
