#include <iostream>
#include <cmath>
#include <iomanip>

#include "getTables.h"
#include "interpolation.h"

using namespace std;

static vector<double> ITable, ToTable, mTable;
static vector<double> TTable, sigTable;

static vector<double> RpPlot, T0Plot;

double T(double z, double Tw, double I)
{
    double T0 = interpolate(I, ITable, ToTable);
    T0Plot.push_back(T0);

    double m = interpolate(I, ITable, mTable);

    double T = T0 + (Tw - T0) * pow(z, m);
    double sigma = interpolate(T, TTable, sigTable);

    return z * sigma;
}

double integrate(double I, double Tw)
{
    double a = 0, b = 1;
    double n = 100;
    double h = (b - a) / n;

    double res = (T(a, Tw, I) + T(b, Tw, I)) / 2;
    double x = a;

    for (int i = 0; i < n - 1; i++) {
        x += h;
        res += T(x, Tw, I);
    }
    res *= h;

    return res;
}

double Rp(double Le, double R, double I, double Tw)
{
    return Le / (2 * M_PI * pow(R, 2) * integrate(I, Tw));
}

double f(double y, double z, double Le, double R, double Lk, double Rk, double Tw)
{
    double resRp = Rp(Le, R, fabs(y), Tw);
    RpPlot.push_back(resRp);

    return (z - (Rk + resRp) * y) / Lk;
}

double g(double I, double Ck)
{
    return -I / Ck;
}

double RungeKutta4(double x, double y, double z, double R, double Lk, double hn, double Rk, double Ck, double Tw, int var)
{
    double k1 = hn * f(x, y, z, R, Lk, Rk, Tw);
    double q1 = hn * g(x, Ck);

    double k2 = hn * f(x + k1 / 2, y + q1 / 2, z, R, Lk, Rk, Tw);
    double q2 = hn * g(x + k1 / 2, Ck);

    double k3 = hn * f(x + k2 / 2, y + q2 / 2, z, R, Lk, Rk, Tw);
    double q3 = hn * g(x + k2 / 2, Ck);

    double k4 = hn * f(x + k3, y + q3, z, R, Lk, Rk, Tw);
    double q4 = hn * g(x + k3, Ck);

    double res = 0;
    if (var == 1) {
        res = x + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    }
    else {
        res = y + (q1 + 2 * q2 + 2 * q3 + q4) / 6;
    }

    return res;
}

int main()
{
    double R = 0, Le = 0, Lk = 0, Ck = 0, Rk = 0,
                  Uco = 0, I0 = 0, Tw = 0;

    double Icur = 0, t = 0;
    double h = 1e-6;

    getParams(R, Le, Lk, Ck, Rk, Uco, I0, Tw);

    getTable1(ITable, ToTable, mTable);
    getTable2(TTable, sigTable);

    vector<double> IPlot, UPlot;
    vector<double> tPlot;

    for (int i = 0; i < 12; i++) {
        Icur = RungeKutta4(Icur, Uco, Le, R, Lk, h, Rk, Ck, Tw, 1);
        Uco = RungeKutta4(Icur, Uco, Le, R, Lk, h, Rk, Ck, Tw, 2);
        t += h;

        tPlot.push_back(t);
        IPlot.push_back(Icur);
        UPlot.push_back(Uco);
    }

    return 0;
}
