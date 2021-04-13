#include <iostream>
#include <cmath>
#include <iomanip>

#include "getTables.h"

using namespace std;

static vector<double> ITable, ToTable, mTable;
static vector<double> TTable, sigTable;

double grp = 0;
double gt0 = 0;

double interpolation(double Y, int nY, vector<double> tableY, vector<double> table)
{
    int imax = 0, imin = 0;
    for (int i = 0; i < nY; i++) {
        if (Y > tableY[i]) {
            imax = i;
        }
        else {
            imax = i;
            break;
        }
    }
    if (0 == imax) {
        imax = 1;
    }
    imin = imax - 1;

    double res = table[imin] + (table[imax] - table[imin]) / (tableY[imax] - tableY[imin]) * (Y - tableY[imin]);
    return res;
}

double Fint(double I, double Tw, double z)
{
    double t0 = interpolation(I, 9, ITable, ToTable);
    gt0 = t0;
    double m = interpolation(I, 9, ITable, mTable);
    double t = t0 + (Tw - t0) * pow(z, m);
    double sigma = interpolation(t, 11, TTable, sigTable);

    return sigma * z;
}

double iint(double I, double Tw)
{
    double a = 0, b = 1;
    double n = 100;
    double h = (b - a) / n;
    double s = (Fint(I, Tw, a) + Fint(I, Tw, b)) / 2;
    double x = 0;

    for (int i = 0; i < n - 1; i++) {
        x = x + h;
        s = s + Fint(I, Tw, x);
    }
    s = s * h;

    return s;
}

double Rp(double le, double R, double I, double Tw)
{
    return le / (2 * M_PI * R * R * iint(I, Tw));
}

double f(double I, double U, double le, double R, double Lk, double Rk, double Tw)
{
    grp = Rp(le, R, fabs(I), Tw);
    return (U - (Rk + grp) * I) / Lk;
}

double g(double I, double Ck)
{
    return -I / Ck;
}

double Inext(double I, double U, double le, double R, double Lk, double hn, double Rk, double Ck, double Tw)
{
    double k1 = f(I, U, le, R, Lk, Rk, Tw);
    double q1 = g(I, Ck);

    double k2 = f(I + hn * k1 / 2, U + hn * q1 / 2, le, R, Lk, Rk, Tw);
    double q2 = g(I + hn * k1 / 2, Ck);

    double k3 = f(I + hn * k2 / 2, U + hn * q2 / 2, le, R, Lk, Rk, Tw);
    double q3 = g(I + hn * k2 / 2, Ck);

    double k4 = f(I + hn * k3, U + hn * q3, le, R, Lk, Rk, Tw);
    double q4 = g(I + hn * k3, Ck);

    return I + hn * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
}

double  Unext(double I, double U, double le, double R, double Lk, double hn, double Rk, double Ck, double Tw)
{
    double k1 = f(I, U, le, R, Lk, Rk, Tw);
    double q1 = g(I, Ck);

    double k2 = f(I + hn * k1 / 2, U + hn * q1 / 2, le, R, Lk, Rk, Tw);
    double q2 = g(I + hn * k1 / 2, Ck);
    double k3 = f(I + hn * k2 / 2, U + hn * q2 / 2, le, R, Lk, Rk, Tw);
    double q3 = g(I + hn * k2 / 2, Ck);
    double k4 = f(I + hn * k3, U + hn * q3, le, R, Lk, Rk, Tw);
    double q4 = g(I + hn * k3, Ck);

    return U + hn * (q1 + 2 * q2 + 2 * q3 + q4) / 6;
}

int main()
{
    double R = 0, l = 0, Lk = 0, Ck = 0, Rk = 0,
                  Uco = 0, I0 = 0, Tw = 0;

    double Icur = 0, t = 0;
    double hn = 0.000001;

    getParams(R, l, Lk, Ck, Rk, Uco, I0, Tw);

    getTable1(ITable, ToTable, mTable);
    getTable2(TTable, sigTable);

    vector<double> IGraph;
    vector<double> TGraph;
    for (int i = 0; i < 12; i++) {
        Icur = Inext(Icur, Uco, l, R, Lk, hn, Rk, Ck, Tw);
        Uco = Unext(Icur, Uco, l, R, Lk, hn, Rk, Ck, Tw);
        t += hn;

        cout << "I = " << Icur << "\tU = " << Uco << endl;

        IGraph.push_back(Icur);
        TGraph.push_back(t);
    }
}
