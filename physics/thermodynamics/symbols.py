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

# Species fully excited TR mass-specific heat
# Species fully excited TR energy per mass
e_s_tr  = IndexedBase('e_s_tr', shape=(ns,))
cv_s_tr = IndexedBase('cv_s_tr', shape=(ns,))
e_tr  = IndexedBase('e_tr', shape=(ns,))
cv_tr = IndexedBase('cv_tr', shape=(ns,))
# Species fully excited VEE energy per mass
e_s_vee  = IndexedBase('e_s_vee', shape=(ns,))
cv_s_vee = IndexedBase('cv_s_vee', shape=(ns,))
e_vee  = IndexedBase('e_vee', shape=(ns,))
cv_vee = IndexedBase('cv_vee', shape=(ns,))
# Combined energy from all components
e_s  = IndexedBase('e_s', shape=(ns,))
cv_s = IndexedBase('cv_s', shape=(ns,))
e  = IndexedBase('e', shape=(ns,))
cv = IndexedBase('cv', shape=(ns,))

# Number of translational-rotational degrees of freedom
ndof = symbols('ndof')
