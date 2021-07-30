import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import pickle
import sympy as sp
from sympy.codegen.ast import Assignment
import cantera as ct
import ctypes
import os

from chemical_mechanism import (ChemicalMechanism, MassProductionRate,
        ForwardRate, BackwardRate, ArrheniusRateConstant, BackwardRateConstant,
        EquilibriumConstant, GibbsEnergyFit, EnthalpyFit, EntropyFit, ThreeBodyRateConstant,
        SpeciesEnergy, MixtureEnergy, SpecificHeatFit, MixtureSpecificHeat)
from source_code import SourceCode
import symbols as syms

def main():

    # Thermo data file
    thermo_file_name = 'N2_ions.dat'
    # Physics file
    physics_file_name = 'physics.pkl'
    regenerate_physics = False
    # Code output file
    wdot_code_file_name = 'rxn_rates.c'
    e_and_cv_code_file_name = 'e_and_cv.c'
    e_s_code_file_name = 'e_s.c'

    chem = ChemicalMechanism(thermo_file_name)
    ns = chem.ns
    nr = chem.nr

    # If the physics file doens't exist, then it has to be regenerated
    if not os.path.isfile(physics_file_name):
        print(f'Physics file "{physics_file_name}" not found, regenerating '
                'physics...')
        regenerate_physics = True

    if regenerate_physics:
        sizes = {syms.ns : ns, syms.nr : nr}
        data = {syms.alpha : chem.alpha,
                syms.beta : chem.beta,
                syms.M : chem.M,
                syms.C : chem.C,
                syms.eta : chem.eta,
                syms.theta : chem.theta,
                syms.p_0 : chem.p_0,
                syms.R : chem.R,
                syms.nu : chem.nu,}

        # Thermo
        H_RT = [None] * ns
        S_R = [None] * ns
        G_RT = [None] * ns
        cv_R = [None] * ns
        e_s = [None] * ns
        for s in range(ns):
            species_sizes = {**sizes, syms.s : s}
            H_RT[s] = EnthalpyFit(chem.species[s].a,
                    chem.species[s].temperature_ranges, species_sizes)
            S_R[s] = EntropyFit(chem.species[s].a,
                    chem.species[s].temperature_ranges, species_sizes)
            G_RT[s] = GibbsEnergyFit(chem.species[s].a,
                    chem.species[s].temperature_ranges, species_sizes)
            cv_R[s] = SpecificHeatFit(chem.species[s].a, chem.species[s].temperature_ranges, species_sizes)
            # Species energy per mass
            e_s[s] = SpeciesEnergy(species_sizes, data, [H_RT[s]])

        # Mixture energy
        e = MixtureEnergy(sizes, data, e_s)
        # Mixture cv
        cv = MixtureSpecificHeat(sizes, data, cv_R)

        # Reactions
        Rf, Rb = [None] * nr, [None] * nr
        for r in range(nr):
            reaction_sizes = {**sizes, syms.r : r}
            reaction_data = {**data, syms.epsilon : chem.epsilon[r]}
            kf = chem.forward_rate[r](reaction_sizes, reaction_data)
            K_c = EquilibriumConstant(reaction_sizes, reaction_data, G_RT)
            kb = BackwardRateConstant(reaction_sizes, reaction_data, [kf, K_c])
            Rf[r] =  ForwardRate(reaction_sizes, reaction_data, [kf])
            Rb[r] = BackwardRate(reaction_sizes, reaction_data, [kb])
        # Production
        wdot = [None] * ns
        for s in range(ns):
            wdot[s] = MassProductionRate({**sizes, syms.s : s}, data, [Rf, Rb])

        # Create assignments to get equations
        wdot_eq = [Assignment(syms.wdot[i], wdot[i].expression) for i in range(ns)]
        e_eq = Assignment(syms.e, e.expression)
        cv_eq = Assignment(syms.cv, cv.expression)
        e_s_eq = [Assignment(syms.e_s[i], e_s[i].expression) for i in range(ns)]

        # Save to file
        with open(physics_file_name, "wb") as physics_file:
            eqs_to_write = (wdot_eq, e_eq, cv_eq, e_s_eq)
            pickle.dump(eqs_to_write, physics_file,
                    protocol=pickle.HIGHEST_PROTOCOL)
    else:
        with open(physics_file_name, "rb") as physics_file:
            wdot_eq, e_eq, cv_eq, e_s_eq = pickle.load(physics_file)

    # Generate C code
    src = SourceCode(wdot_code_file_name, 'wdot', wdot_eq)
    src = SourceCode(e_and_cv_code_file_name, 'e_and_cv', [e_eq, cv_eq],
            pointer=True)
    src = SourceCode(e_s_code_file_name, 'e_s', e_s_eq)

    print('Code generated - remember to compile.')
    print('Recommended commands:')
    print('gcc -fPIC -shared energy.h e_and_cv.c e_s.c -o energy.so')
    print('gcc -fPIC -shared rates.h rxn_rates.c -o rxn_rates.so')


if __name__ == "__main__":
    main()
