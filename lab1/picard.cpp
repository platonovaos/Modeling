#include "picard.h"

yPicard initPicard(const int n)
{
    yPicard yp;
    yp.y1 = new double[n];
    yp.y2 = new double[n];
    yp.y3 = new double[n];
    yp.y4 = new double[n];

    return yp;
}

void deletePicard(yPicard yp)
{
    delete [] yp.y1;
    delete [] yp.y2;
    delete [] yp.y3;
    delete [] yp.y4;
}

double y1(const double x)
{
    return (pow(x, 3) / 3);
}

double y2(const double x)
{
    return (y1(x) + pow(x, 7) / 63);
}

double y3(const double x)
{
    return (y2(x) + 2 * pow(x, 11) / 2079 + pow(x, 15) / 59535);
}

double y4(const double x)
{
    return (y2(x) + 2 * pow(x, 11) / 2079 + 13 * pow(x, 15) / 218295 + 82 * pow(x, 19) / 37328445
            + 662 * pow(x, 23) / 10438212015 + 4 * pow(x, 27) / 3341878155 + pow(x, 31) / 109876902975);
}

void Picard(const int n, const double *x, yPicard yp)
{
    for (int i = 0; i < n; i++) {
        yp.y1[i] = y1(x[i]);
        yp.y2[i] = y2(x[i]);
        yp.y3[i] = y3(x[i]);
        yp.y4[i] = y4(x[i]);
    }
}
