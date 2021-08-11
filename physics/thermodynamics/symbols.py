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
e_s_tr  = IndexedBase('e_s_tr')
cv_s_tr = IndexedBase('cv_s_tr')
e_tr  = IndexedBase('e_tr')
cv_tr = IndexedBase('cv_tr')
# Species fully excited VEE energy per mass
e_s_vee  = IndexedBase('e_s_vee')
cv_s_vee = IndexedBase('cv_s_vee')
e_vee  = IndexedBase('e_vee')
cv_vee = IndexedBase('cv_vee')
# Combined energy from all components
e_s  = IndexedBase('e_s')
cv_s = IndexedBase('cv_s')
e  = IndexedBase('e')
cv = IndexedBase('cv')

# Number of translational-rotational degrees of freedom
ndof = symbols('ndof')
