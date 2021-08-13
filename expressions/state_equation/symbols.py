from sympy import symbols, IndexedBase

# This file contains the definition of all symbols used for the expression
# created in this directory.

# Total number of species
ns = symbols('ns')
# Heavies index and total number of heavies
s = symbols('s')
nh = symbols('nh')
# Index of electrons
idx_e = symbols('idx_e')

# Heavy pressure
p_h = symbols('p_h')
# Electron pressure
p_e = symbols('p_e')
# Pressure
p = symbols('p')
# Partial densities
rho = IndexedBase('rho', shape=(ns,))
