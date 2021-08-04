import numpy as np
import pickle
import sympy as sp
from sympy.codegen.ast import Assignment

import symbols as syms


def main():

    # Get thermo data somehow
    #thermo_data = ThermodynamicData(thermo_file_name)
    ns = chem.ns

    sizes = {syms.ns : ns}
    data = {syms.M : chem.M,
            syms.R : chem.R}

    # Thermo
    e_s = [None] * ns
    for s in range(ns):
        species_sizes = {**sizes, syms.s : s}
        H_RT = EnthalpyFit(chem.species[s].a,
                chem.species[s].temperature_ranges, species_sizes)
        # Species energy per mass
        e_s[s] = SpeciesEnergy(species_sizes, data, [H_RT[s]])

    # Create assignments to get equations
    e_s_eq = [Assignment(syms.e_s[i], e_s[i].expression) for i in range(ns)]

    # Save to file
    with open(physics_file_name, "wb") as physics_file:
        eqs_to_write = e_s_eq
        pickle.dump(eqs_to_write, physics_file,
                protocol=pickle.HIGHEST_PROTOCOL)


if __name__ == "__main__":
    main()
