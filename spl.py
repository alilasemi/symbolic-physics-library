import pickle

import expressions.species_energies_and_cv.main as thermo_expression
import expressions.mass_production.main as chem_expression
from core.source_code import SourceCode

def main():

    chem_expression.create()
    breakpoint()

    thermo_expression.create()

    physics_file_name = 'physics.pkl'
    with open(physics_file_name, "rb") as physics_file:
        _, e_tr, e_vee, _, cv_tr, cv_vee, *_ = pickle.load(physics_file)

    SourceCode('generated/e_tr', 'e_tr', 'const double T, const double* __restrict__ Y, double* __restrict__ e_tr', e_tr, pointer=True)
    SourceCode('generated/e_vee', 'e_vee', 'const double T, const double* __restrict__ Y, double* __restrict__ e_vee', e_vee, pointer=True)
    SourceCode('generated/cv_tr', 'cv_tr', 'const double T, const double* __restrict__ Y, double* __restrict__ cv_tr', cv_tr, pointer=True)
    SourceCode('generated/cv_vee', 'cv_vee', 'const double T, const double* __restrict__ Y, double* __restrict__ cv_vee', cv_vee, pointer=True)



if __name__ == "__main__":
    main()
