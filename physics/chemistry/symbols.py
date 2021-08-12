from sympy import symbols, IndexedBase

# -- Definition of symbols -- #

# Species index and total number of species
s = symbols('s')
ns = symbols('ns')
# Reactions index and total number of reactions
r = symbols('r')
nr = symbols('nr')

# Temperature
T = symbols('T', real=True, positive=True)
# Partial densities
rho = IndexedBase('rho', shape=(ns,))
# Mass fractions
Y = IndexedBase('Y', shape=(ns,))

# Standard pressure
p_0 = symbols('p_0')
# Universal gas constant
R = symbols('R')
# Molar mass
M = IndexedBase('M')

# Reactant and product stoichiometric coefficients
alpha = IndexedBase('alpha')
beta = IndexedBase('beta')
# Sum of stoichiometric coefficients
nu = symbols('nu')
# Rate coefficients
kf = symbols('kf')
kb = symbols('kb')
# Equilibrium constant
K_c = symbols('K_c')
# Forward and backward reaction rates
Rf = IndexedBase('Rf')
Rb = IndexedBase('Rb')
# Species mass production rates
wdot = IndexedBase('wdot', shape=(ns,))

# Nondimensional enthalpy entropy and Gibbs energy
H_RT = IndexedBase('H_RT')
S_R  = IndexedBase('S_R')
G_RT = IndexedBase('G_RT')

# Arrhenius fit parameters
C = symbols('C')
eta = symbols('eta')
theta = symbols('theta')
# Collision efficiencies
epsilon = IndexedBase('epsilon')
