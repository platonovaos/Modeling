#ifndef PICARD_H
#define PICARD_H

#include <cmath>

struct yPicard
{
    double *y1;
    double *y2;
    double *y3;
    double *y4;
};

yPicard initPicard(const int);
void deletePicard(yPicard);


void Picard(const int n, const double *x, yPicard yp);

#endif // PICARD_H
