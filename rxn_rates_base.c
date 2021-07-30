#include "rates.h"
#include <math.h>

void eval_spec_rates(const double T, const double* __restrict__ rho,
        double* __restrict__ wdot) {

    // --- Compute the forward and backward rates for each reaction --- //
#pragma rates

    // --- Compute the species mass production rate for each species -- //
#pragma wdot
}
