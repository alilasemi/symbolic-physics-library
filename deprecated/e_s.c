#include "energy.h"
#include <math.h>

void get_e_s(const double T, double* __restrict__ e_s) {

    // --- energy of each species -- //
    e_s[0] = 296.80220444419041*T*(((T > 6000.0) ? (
   2.8750777620000004e-16*pow(T, 4) - 2.4264885274999999e-11*pow(T, 3) + 8.2896777766666663e-7*pow(T, 2) - 0.015325460230000001*T + 202.02646350000001 - 642073.35400000005*log(T)/T + 4938707.04/T - 831013916.0/pow(T, 2)
)
: ((T > 1000.0) ? (
   2.1239087719999999e-16*pow(T, 4) - 4.8077637125000002e-12*pow(T, 3) + 4.9726889299999999e-8*pow(T, 2) - 0.00030698427499999999*T + 6.0669492199999997 - 2239.249073*log(T)/T + 12832.104149999999/T - 587712.40599999996/pow(T, 2)
)
: ((T > 200.0) ? (
   5.0394116180000006e-13*pow(T, 4) - 2.4064484049999999e-9*pow(T, 3) + 4.6154872966666664e-6*pow(T, 2) - 0.0042654572049999999*T + 6.0827383599999996 - 381.846182*log(T)/T + 710.84608600000001/T - 22103.714970000001/pow(T, 2)
)
: (
   0
)))) - 1);
    e_s[1] = 296.80801696295305*T*(((T > 6000.0) ? (
   -1.1243786180000001e-16*pow(T, 4) + 1.0361103075e-11*pow(T, 3) - 3.9168850833333332e-7*pow(T, 2) + 0.0078559664299999998*T - 96.035180499999996 + 313928.72340000002*log(T)/T - 2217361.8670000001/T + 371282977.0/pow(T, 2)
)
: ((T > 1000.0) ? (
   1.0823992940000001e-16*pow(T, 4) + 5.2562863625000003e-12*pow(T, 3) - 1.4538841033333331e-7*pow(T, 2) + 0.0015343385294999999*T - 2.8848863850000002 + 7058.8930300000002*log(T)/T + 134038.84830000001/T + 2845599.0019999999/pow(T, 2)
)
: ((T > 298.14999999999998) ? (
   3.2435120000000001e-13*pow(T, 4) - 1.4093262425000001e-9*pow(T, 3) + 2.2434921333333332e-6*pow(T, 2) - 0.0010661198905*T + 3.1649163699999998 + 269.62227030000003*log(T)/T + 179000.4424/T + 34740.474699999999/pow(T, 2)
)
: (
   0
)))) - 1);
    e_s[2] = 593.60440888838082*T*(((T > 6000.0) ? (
   2.5559720480000001e-17*pow(T, 4) - 2.7459192725000001e-12*pow(T, 3) + 1.2758574666666667e-7*pow(T, 2) - 0.0034239940649999998*T + 69.167827399999993 - 310757.49800000002*log(T)/T + 2550585.6179999998/T - 547518105.0/pow(T, 2)
)
: ((T > 1000.0) ? (
   -5.3544551420000006e-16*pow(T, 4) + 1.0031644699999999e-11*pow(T, 3) - 5.7650503333333327e-8*pow(T, 2) + 0.00014583600404999999*T + 2.3621882869999999 - 107.12315*log(T)/T + 56973.513299999999/T - 88765.013800000001/pow(T, 2)
)
: ((T > 200.0) ? (
   2.5 + 56104.637799999997/T
)
: (
   0
)))) - 1);
    e_s[3] = 593.62765941876523*T*(((T > 6000.0) ? (
   7.0798631860000008e-19*pow(T, 4) - 6.7285771575000001e-14*pow(T, 3) + 3.4082711866666664e-9*pow(T, 2) - 0.00010026967915*T + 4.9769866399999998 - 11131.652179999999*log(T)/T + 313628.46960000001/T - 16460921.48/pow(T, 2)
)
: ((T > 1000.0) ? (
   1.0092332558000001e-16*pow(T, 4) - 3.474585305e-12*pow(T, 3) + 4.5078343566666665e-8*pow(T, 2) - 0.00026441335949999999*T + 3.4773892900000001 - 855.79086099999995*log(T)/T + 231080.99840000001/T - 290497.03739999997/pow(T, 2)
)
: ((T > 298.14999999999998) ? (
   -8.8947019680000013e-16*pow(T, 4) + 4.6252783299999999e-12*pow(T, 3) - 1.0448158586666667e-8*pow(T, 2) + 1.368745378e-5*T + 2.4874888209999999 + 2.299958315*log(T)/T + 225628.47380000001/T - 5237.0792099999999/pow(T, 2)
)
: (
   0
)))) - 1);
    e_s[4] = 15156338.34320985*T*(((T > 298.14999999999998) ? (
   2.5 - 745.375/T
)
: (
   0
)) - 1);
}
