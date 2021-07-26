#ifndef RATES_H
#define RATES_H

void eval_spec_rates(const double T, const double * __restrict__ Y,
        double * __restrict__ Rf,
        double * __restrict__ Rb,
        double * __restrict__ wdot);

#endif
