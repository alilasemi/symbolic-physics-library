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
    # Paths
    kinetics_file_path = 'examples/nitrogen/N2-ion_park90.chem'

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
        _, _, e_and_cv, e_s, cv_s = pickle.load(physics_file)

    # Generate C code
    SourceCode('generated/e_and_cv', ['e', 'cv'], 'const double T, const double* __restrict__ Y, double* __restrict__ e, double* __restrict__ cv', e_and_cv, pointer=True)
    SourceCode('generated/e_s', 'e_s', 'const double T, double* __restrict__ e_s', e_s)
    SourceCode('generated/cv_s', 'cv_s', 'const double T, double* __restrict__ cv_s', cv_s)

    # Chemistry
    chem_expression.create(kinetics_data)

    physics_file_name = 'physics.pkl'
    with open(physics_file_name, "rb") as physics_file:
        wdot = pickle.load(physics_file)

    # Generate C code
    SourceCode('generated/wdot', 'wdot', 'const double T, const double* __restrict__ rho, double* __restrict__ wdot', wdot)

    # TODO: This section of the code (for the temperature relaxation term) is
    # not yet fully tested, so it is commented for now.

#    # Temperature relaxation
#    temp_relax_expression.create(e_s_vee, p)
#
#    physics_file_name = 'physics.pkl'
#    with open(physics_file_name, "rb") as physics_file:
#        Q_TV = pickle.load(physics_file)
#
#    SourceCode('generated/Q_TV', 'Q_TV', 'const double T, const double Tv, const double* __restrict__ rho, const double* __restrict__ Y, double* __restrict__ Q_TV', Q_TV, pointer=True)



if __name__ == "__main__":
    main()
