import pickle

import expressions.state_equation.main as state_expression

# TODO
#import expressions.species_energies_and_cv.main as thermo_expression
import expressions.species_energies_and_cv_equilibrium.main as thermo_expression

import expressions.mass_production.main as chem_expression
import expressions.landau_teller_source.main as temp_relax_expression
from core.source_code import SourceCode
from physics.kinetics.kinetics_data import KineticsData

def main():

    # Inputs
    # TODO: Unhardcode
    # Paths
    kinetics_file_path = '/home/ali/projects/symbolic-physics-library/examples/nitrogen/N2-ion_park90.chem'
    # Thermal model
    thermal_model = 'equilibrium'

    # --- Data --- #
    # Kinetics data
    kinetics_data = KineticsData(kinetics_file_path)

    # --- Physics --- #

    # State equation
    state_expression.create()

    physics_file_name = 'physics.pkl'
    with open(physics_file_name, "rb") as physics_file:
        p = pickle.load(physics_file)

    # Thermodynamics
    thermo_expression.create(kinetics_data.species)

    physics_file_name = 'physics.pkl'
    with open(physics_file_name, "rb") as physics_file:
        # TODO
        #e, e_tr, e_tr_from_e, e_vee, _, cv_tr, cv_vee, _, _, e_s_vee, *_ = pickle.load(physics_file)
        _, _, e_and_cv, e_s, cv_s = pickle.load(physics_file)

    SourceCode('generated/e_and_cv', ['e', 'cv'], 'const double T, const double* __restrict__ Y, double* __restrict__ e, double* __restrict__ cv', e_and_cv, pointer=True)
    SourceCode('generated/e_s', 'e_s', 'const double T, double* __restrict__ e_s', e_s)
    SourceCode('generated/cv_s', 'cv_s', 'const double T, double* __restrict__ cv_s', cv_s)

    breakpoint()

    #SourceCode('generated/e', 'e', 'const double T, const double* __restrict__ Y, double* __restrict__ e', e, pointer=True)
    #SourceCode('generated/e_tr', 'e_tr', 'const double T, const double* __restrict__ Y, double* __restrict__ e_tr', e_tr, pointer=True)
    #SourceCode('generated/e_tr_from_e', 'e_tr_from_e', 'const double e, const double e_vee, const double* __restrict__ Y, double* __restrict__ e_tr_from_e', e_tr_from_e, pointer=True)
    #SourceCode('generated/e_vee', 'e_vee', 'const double T, const double* __restrict__ Y, double* __restrict__ e_vee', e_vee, pointer=True)
    #SourceCode('generated/e_s_vee', 'e_s_vee', 'const double T, double* __restrict__ e_s_vee', e_s_vee)
    #SourceCode('generated/cv_tr', 'cv_tr', 'const double T, const double* __restrict__ Y, double* __restrict__ cv_tr', cv_tr, pointer=True)
    #SourceCode('generated/cv_vee', 'cv_vee', 'const double T, const double* __restrict__ Y, double* __restrict__ cv_vee', cv_vee, pointer=True)

    # Temperature relaxation
    temp_relax_expression.create(e_s_vee, p)

    physics_file_name = 'physics.pkl'
    with open(physics_file_name, "rb") as physics_file:
        Q_TV = pickle.load(physics_file)

    SourceCode('generated/Q_TV', 'Q_TV', 'const double T, const double Tv, const double* __restrict__ rho, const double* __restrict__ Y, double* __restrict__ Q_TV', Q_TV, pointer=True)

    # Chemistry
    chem_expression.create(kinetics_data)

    physics_file_name = 'physics.pkl'
    with open(physics_file_name, "rb") as physics_file:
        wdot = pickle.load(physics_file)
    breakpoint()

    SourceCode('generated/wdot', 'wdot', 'const double T, const double* __restrict__ rho, double* __restrict__ wdot', wdot)


if __name__ == "__main__":
    main()
