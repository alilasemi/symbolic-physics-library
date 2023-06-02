# Symbolic Physics Library

This repository contains code which uses Sympy to generate C code to compute
thermodynamics and nonequilibrium chemistry terms. To get started, run:

    python spl.py

This will create `.c` files in the `generated` directory. For now, the following
files are generated:

- `cv_s.c`: computes the specific heat at constant volume for each species
- `e_and_cv.c`: computes the total internal energy and specific heat for the
  mixture
- `e_s.c`: computes the internal energy of each species
- `wdot.c`: computes the mass production rate of each species
