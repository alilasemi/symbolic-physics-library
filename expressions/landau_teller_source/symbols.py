from sympy import symbols, IndexedBase

# This file contains the definition of all symbols used for the expression
# created in this directory.

# Species index and total number of species
s = symbols('s')
ns = symbols('ns')
# Heavies index and total number of heavies
r = symbols('r')
nh = symbols('nh')
# Total number of vibrating species
nv = symbols('nv')

# TR Temperature
T = symbols('T', real=True, positive=True)
# VEE Temperature
Tv = symbols('Tv', real=True, positive=True)
# Total density
rho_t = symbols('rho_t')
# Partial densities
rho = IndexedBase('rho', shape=(ns,))
# Mole fractions
X = IndexedBase('X', shape=(ns,))
# Mass fractions
Y = IndexedBase('Y', shape=(ns,))

# Pressure
p = symbols('p')
# Molar mass of species s
M = IndexedBase('M')
# Mass of species s
m = IndexedBase('m')
# Mixture averaged molar mass
M_bar = symbols('M_bar')
# Mixture averaged species mass
m_bar = IndexedBase('m_bar')

# Boltzmann constant
k = symbols('k')
# Avagadro's number
N_A = symbols('N_A')
# pi
pi = symbols('pi')

# Vibrational relaxation parameters
a_vib, b_vib = symbols('a_vib b_vib')

# Mass-specific species vibrational energy
e_v_s = IndexedBase('e_s', shape=(ns,))

# Vibrational cross section of species s in Park 93 high temperature correction
sigma_v = IndexedBase('sigma_v')

# Vibrational cross section coefficient of species s in Park 93 high temperature
# correct
sigma_v_p = IndexedBase('sigma_v_p')

# Uncorrected vibrational relaxation time of species s
tau_v_s = IndexedBase('tau_v_s')

# Total number density computed from mass density
n = symbols('n')

# Vibrational relaxation time of species s, with the Park 93 high temperature
# correction
tau_s = IndexedBase('tau_s')
tau = IndexedBase('tau')

# VEE energy
e_s_vee = IndexedBase('e_s_vee')
# VEE energy evaluated at the TR temperature
e_s_vee_star = IndexedBase('e_s_vee_star')

# LT source term
Q_TV_s = IndexedBase('Q_TV_s')
