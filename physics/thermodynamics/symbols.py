from sympy import symbols, IndexedBase

# This file contains the definition of all symbols used for the expressions
# created in this directory.

# Species index and total number of species
s = symbols('s')
ns = symbols('ns')

# Temperature
T = symbols('T', real=True, positive=True)
# Mass fractions
Y = IndexedBase('Y', shape=(ns,))

# Universal gas constant
R = symbols('R')
# Molar mass
M = IndexedBase('M')

# Polynomial fit coefficients
a = IndexedBase('a', real=True)

# Nondimensional enthalpy, entropy, and Gibbs energy
H_RT = IndexedBase('H_RT')
S_R  = IndexedBase('S_R')
G_RT = IndexedBase('G_RT')

# Mass-specific species energy
e_s = IndexedBase('e_s', shape=(ns,))
# Mass-specific total mixture energy
e = symbols('e')
# Nondimensional specific heat at constant volume for each species
cv_R = IndexedBase('cv_R')
# Mass-specific heat at constant volume for mixture
cv = symbols('cv')
