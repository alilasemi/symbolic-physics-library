from sympy import symbols, IndexedBase

# This file contains the definition of all symbols used for the expression
# created in this directory.

# Species index and total number of species
s = symbols('s')
ns = symbols('ns')

# Temperature
T = symbols('T', real=True, positive=True)
# Partial densities
rho = IndexedBase('rho', shape=(ns,))
# Mass fractions
Y = IndexedBase('Y', shape=(ns,))

# Universal gas constant
R = symbols('R')
# Molar mass
M = IndexedBase('M')

# Polynomial fit coefficients
a = IndexedBase('a', real=True)

# Nondimensional enthalpy entropy and Gibbs energy
H_RT = IndexedBase('H_RT')
S_R  = IndexedBase('S_R')
G_RT = IndexedBase('G_RT')
