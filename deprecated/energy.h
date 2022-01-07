#ifndef ENERGY_H
#define ENERGY_H

void get_e_and_cv(const double T, const double* __restrict__ Y,
        double* __restrict__ e, double* __restrict__ cv);

void get_e_s(const double T, double* __restrict__ e_s);

#endif
