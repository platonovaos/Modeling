#ifndef GETTABLES_H
#define GETTABLES_H

#include <vector>
using namespace std;

struct InParams
{
    double a[2], b[2], c[2];
    double m[2];
};

void getTable(InParams &ip);
void getParams(double &alpha0, double &alphaN, double &l, double &T0, double &R, double &F);

#endif // GETTABLES_H
