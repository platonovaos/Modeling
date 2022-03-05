#include "interpolation.h"

double interpolate(double x, vector<double> Xtable, vector<double> f)
{
    int x0 = 0, x1 = 1;

    for (int i = 0; i < Xtable.size(); i++) {
        x1 = i;

        if (x <= Xtable[i]) {
            break;
        }
    }

    x0 = x1 - 1;

    double res = f[x0] + (f[x1] - f[x0]) / (Xtable[x1] - Xtable[x0])
                * (x - Xtable[x0]);
    return res;
}


