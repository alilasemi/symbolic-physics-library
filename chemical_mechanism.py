import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import sympy as sp
from sympy import exp, log

import symbols as syms
import physical_expressions

# This is a class to be inherited from for any mathematical expression that
# needs to be created. Each term in an expression must be defined as a
# subexpression, which is eventually plugged in.
class Expression:

    # Symbol representing this expression
    symbol = None

    # Sympy expression
    expression = None

    # All the sizes (for example, the bounds on sums or products)
    sizes = {}

    # All the data needed - basically, symbols to directly replace with numbers
    data = {}

    # All the expressions needed to evaluate this expression
    sub_expressions = []

    def __init__(self, sizes = {}, data = {}, sub_expressions = []):
        '''
        Set the given subexpressions and sizes, then plug them all in.
        '''
        self.sizes = sizes
        self.data = data
        self.sub_expressions = sub_expressions
        self.plug_in()

    def plug_in(self):
        '''
        Plug in all sizes, data, and sub expressions of this expression.
        '''
        # Plug in sizes first. This is so that sums and products are expanded.
        for symbol, size in self.sizes.items():
            # Plug in size and expand
            self.expression = self.expression.subs(symbol, size).doit()
            # Plug in to the symbol (useful for indexed symbols)
            self.symbol = self.symbol.subs(symbol, size)

        # Next, data is plugged in. This way, terms that become zero are
        # cancelled out as soon as possible.
        for symbol, value in self.data.items():
            # Plug in arrays elementwise
            if isinstance(value, np.ndarray):
                for i in range(value.size):
                    # Get the index of this element
                    index = np.unravel_index(i, value.shape)
                    # Plug in
                    self.expression = self.expression.subs(symbol[index],
                            value[index])
            # Plug in everything else directly
            else:
                self.expression = self.expression.subs(symbol, value)

        # Loop over subexpressions
        for sub_expression in self.sub_expressions:
            # If it's a list of expressions, plug in each one
            if isinstance(sub_expression, list):
                self.expression = self.expression.subs([(
                    sub_expression[i].symbol[i], sub_expression[i].expression)
                    for i in range(len(sub_expression))])
            # Otherwise, just plug it in
            else:
                self.expression = self.expression.subs(sub_expression.symbol,
                        sub_expression.expression)

class MassProductionRate(Expression):
    symbol = syms.wdot
    expression = physical_expressions.wdot_expr

class ForwardRate(Expression):
    symbol = syms.Rf
    expression = physical_expressions.Rf_expr

class BackwardRate(Expression):
    symbol = syms.Rb
    expression = physical_expressions.Rb_expr

class ArrheniusRateConstant(Expression):
    symbol = syms.kf
    expression = physical_expressions.kf_arrhenius_expr

class ThreeBodyRateConstant(Expression):
    symbol = syms.kf
    expression = physical_expressions.kf_three_body_expr
    sub_expressions = [ArrheniusRateConstant]

    def __init__(self, sizes = {}, data = {}, sub_expressions = []):
        # Plug in the Arrhenius forward rate first, do everything else after
        self.plug_in()
        super().__init__(sizes, data, sub_expressions)

class BackwardRateConstant(Expression):
    symbol = syms.kb
    expression = physical_expressions.kb_expr

class EquilibriumConstant(Expression):
    symbol = syms.K_c
    expression = physical_expressions.K_c_expr

class NASA9Fit(Expression):
    def __init__(self, a, temperature_ranges, sizes = {}, data = {}, sub_expressions = []):
        super().__init__(sizes, data, sub_expressions)
        self.create_piecewise_expression_from_fit(a, temperature_ranges)

    def create_piecewise_expression_from_fit(self, a, temperature_ranges):
        n_ranges, n_coeffs = a.shape
        # TODO: Add out-of-bounds coefficients (currently just returns zero)
        self.expression = sp.Piecewise(
                # The pieces are reversed so that high temperature comes first
                *reversed(
                # The default is set to 0
                [(0, True)] +
                # For the ith temp. range, the jth coefficient is plugged in
                [(self.expression.subs([(syms.a[j], a[i, j])
                    # A tuple is created by joining the expression with its
                    # valid temperature range
                    for j in range(n_coeffs)]), temperature_ranges[i])
                    for i in range(n_ranges)]))

class EnthalpyFit(NASA9Fit):
    symbol = syms.H_RT[syms.s]
    expression = physical_expressions.H_RT_expr

class EntropyFit(NASA9Fit):
    symbol = syms.S_R[syms.s]
    expression = physical_expressions.S_R_expr

class SpecificHeatFit(NASA9Fit):
    symbol = syms.cv_R[syms.s]
    expression = physical_expressions.cv_R_expr

class GibbsEnergyFit(NASA9Fit):
    symbol = syms.G_RT[syms.s]
    expression = physical_expressions.G_RT_expr

class SpeciesEnergy(Expression):
    symbol = syms.e_s[syms.s]
    expression = physical_expressions.e_s_expr

class MixtureEnergy(Expression):
    symbol = syms.e
    expression = physical_expressions.e_expr

class MixtureSpecificHeat(Expression):
    symbol = syms.cv
    expression = physical_expressions.cv_expr

class ChemicalMechanism:

    p_0 = 1e5 # 1 bar, converted to Pa
    R = 8.31446261815324 # Universal gas constant in SI units

    def __init__(self, thermo_file_name):
        self.thermo_file_name = thermo_file_name
        self.create_species()
        self.create_reactions()

    def create_species(self):
        # Species
        species = ['N2', 'N2+', 'N', 'N+', 'e-']
        self.ns = len(species)
        # Load thermo file into NASA9 objects
        self.species = load_thermo_data(self.thermo_file_name, species)
        self.M = np.array([s.M for s in self.species])

    def create_reactions(self):

        # TODO: Read these from DPLR chem files directly
        # Reactant stoichiometric coefficients
        self.alpha = np.array([
                [1, 0, 0],
                [0, 0, 0],
                [0, 1, 2],
                [0, 0, 0],
                [0, 1, 0],
                ])
        # Product stoichiometric coefficients
        self.beta = np.array([
                [0, 0, 0],
                [0, 0, 1],
                [2, 0, 0],
                [0, 1, 0],
                [0, 2, 1],
                ])
        # Net stoichiometric coefficient
        self.nu = np.sum(self.beta - self.alpha, axis=0)
        self.nr = self.alpha.shape[1]

        # Arrhenius parameters
        self.C = np.array([
            7.000e+18,
            2.500e+31,
            4.400e+04,
            ])
        self.eta = np.array([
            -1.60e+00,
            -3.82e+00,
             1.50e+00,
             ])
        self.theta = np.array([
            1.1320e+05,
            1.6860e+05,
            6.7500e+04,
            ])

        # Whether or not the third body collision partner, M, is present
        self.has_m = [True, False, False]

        self.forward_rate = [ThreeBodyRateConstant, ArrheniusRateConstant,
                ArrheniusRateConstant]

        # The units on C need to be in SI units:
        #     m^(3 a) * mol^(1 - a) * s,
        # where a is the sum of reactant stoichimetric coeffs (including M).
        # The DPLR input file assumes kilomoles, so this needs to be converted.
        kmol_to_mol = 1e3
        a = np.sum(self.alpha, axis=0) + self.has_m
        conversion = kmol_to_mol ** (1 - a)
        self.C /= conversion

        self.epsilon = [np.ones(self.ns, dtype=int) for r in range(self.nr)]
        # Make sure units equal when doing this normalization!
        self.epsilon[0] = np.array(
                [7.00e+18, 7.00e+18, 3.00e+19, 3.00e+19, 3.00e+21])/7.000e+18

class NASA9:
    '''
    H_RT
    '''
    species_name = ''
    M = 0;
    dHf = 0;
    n_ranges = 0;
    n_coeffs = 9;
    temperatures = np.array([], dtype=float)
    text = []

    def __init__(self, text):
        self.temperature_ranges = []
        # Store text from NASA9 file for this species
        self.text = text
        # Parse text
        self.parse_text()

    # Parse text from a NASA9 file into data structures.
    def parse_text(self):
        # The first line contains the species name
        self.species_name = self.text[0].split()[0]
        # The second line contains the number of temperature ranges, the
        # molecular weight, and the heat of formation
        self.n_ranges = int(self.text[1].split()[0])
        self.M        = float(self.text[1].split()[-2]) / 1000. # Convert to kg
        self.dHf      = float(self.text[1].split()[-1])
        self.temperatures = np.empty(self.n_ranges + 1)
        self.a = np.empty((self.n_ranges, self.n_coeffs))
        # Loop over temperature ranges
        for i in range(self.n_ranges):
            # Get text for this range
            range_text = self.text[2 + i*3 : 5 + i*3]
            # The first line in the range contains the temperatures
            self.temperatures[i : i + 2] = range_text[0].split('7')[0].split()
            self.temperature_ranges.append((syms.T > self.temperatures[i]))
            # The second line in the range contains the first five coefficients.
            # Must convert D to e for exponentials
            line = range_text[1].replace('D', 'e')
            for num in range(5):
                self.a[i, num] = float(line[num*16 : (num + 1)*16])
            # The third line in the range contains the last four coefficients.
            # Must convert D to e for exponentials
            line = range_text[2].replace('D', 'e')
            for num in range(2):
                self.a[i, 5 + num] = float(line[num*16 : (num + 1)*16])
            for num in range(2):
                self.a[i, 7 + num] = float(line[(num + 3)*16 : (num + 4)*16])

        # Plug the coefficients in to produce the piecewise fits
        self.H_RT = self.create_piecewise_expression_from_fit(
                physical_expressions.H_RT_expr)
        self.S_R = self.create_piecewise_expression_from_fit(
                physical_expressions.S_R_expr)
        self.G_RT = self.create_piecewise_expression_from_fit(
                physical_expressions.G_RT_expr)

    def create_piecewise_expression_from_fit(self, fit):
    # TODO: Add out-of-bounds coefficients (currently just returns zero)
        return sp.Piecewise(
                # The pieces are reversed so that high temperature comes first
                *reversed(
                # The default is set to 0
                [(0, True)] +
                # For the ith temp. range, the jth coefficient is plugged in
                [(fit.subs([(syms.a[j], self.a[i, j])
                    # A tuple is created by joining the expression with its
                    # valid temperature range
                    for j in range(self.n_coeffs)]), self.temperature_ranges[i])
                    for i in range(self.n_ranges)]))

# Load data from a NASA9 file into a list of NASA9 objects.
def load_thermo_data(thermo_file_name, species):
    thermo_data = []
    # Loop over species
    for s, species_name in enumerate(species):
        # Read the text for this species
        text = []
        with open(thermo_file_name, 'r') as thermo_file:
            relevant_text = False
            for line in thermo_file:
                if line.split()[0] == species_name:
                    text.append(line)
                    relevant_text = True
                elif relevant_text:
                    if line.startswith(' ') or line.startswith('-'):
                        text.append(line)
                    else:
                        break
        # Create NASA9 object and store
        thermo_data.append(NASA9(text))
    breakpoint()
    return thermo_data
