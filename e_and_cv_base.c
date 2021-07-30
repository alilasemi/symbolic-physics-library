#include "energy.h"
#include <math.h>

void get_e_and_cv(const double T, const double* __restrict__ Y,
        double* __restrict__ e, double* __restrict__ cv) {

    // --- Compute the mixture energy and specific heat at constant volume -- //
#pragma e_and_cv
}
