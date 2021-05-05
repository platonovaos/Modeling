#include <iostream>
#include <cmath>
#include <iomanip>

#include "getTables.h"
#include "interpolation.h"

using namespace std;

struct Params
{
    double K; double M; double P;
};

struct Coefs
{
    vector<double> A; vector<double> B;
    vector<double> C; vector<double> D;
};

static vector<double> TLTable, LTable, TKTable, KTable;
static vector<double> TFPlot, TAlphaPlot, TF0Plot;
static double eps1 = 0.1, eps2 = 0.1;

double p(const double kf, const double T)
{
    return kf * interpolate(T, TKTable, KTable) * pow(T, 3);
}

double f(const double kf, const double T, const double T0)
{
    return kf * interpolate(T, TKTable, KTable) * pow(T0, 4);
}

void getParams0(Params &p0, const vector<double> T,
                 const double kf, const double T0, const double h, const double F0)
{
    double int0 = interpolate(T[0], TLTable, LTable);
    double int1 = interpolate(T[1], TLTable, LTable);
    double hi12 = (int0 + int1) / 2;

    double f0 = f(kf, T[0], T0);
    double f1 = f(kf, T[1], T0);
    double f12 = (f0 + f1) / 2;

    p0.K = hi12 / h;
    p0.M = -hi12 / h;
    p0.P = F0 + (h / 4) * (f0 + f12);
}

void getParamsN(Params &pN, const vector<double> T, const int N,
                  const double kf, const double T0, const double h, const double alpha)
{
    double intN1 = interpolate(T[N - 1], TLTable, LTable);
    double intN2 = interpolate(T[N - 2], TLTable, LTable);
    double hiN12 = (intN1 + intN2) / 2;

    double fN1 = f(kf, T[N - 1], T0);
    double fN2 = f(kf, T[N - 2], T0);
    double fN12 = (fN1 + fN2) / 2;

    pN.K = hiN12 / h;
    pN.M = -alpha - hiN12 / h;
    pN.P = -alpha * T0 - (h / 4) * (fN12 + fN1);
}

void getCoefs(vector<double> &T, const int N, Coefs &c,
              const double kf, const double T0, const double h)
{
    for (int i = 1; i < N; i++) {
        double lPrev = interpolate(T[i - 1], TLTable, LTable);
        double lCur = interpolate(T[i], TLTable, LTable);
        double lNext = interpolate(T[i + 1], TLTable, LTable);

        c.A.push_back((lPrev + lCur) / (2 * h));
        c.C.push_back((lCur + lNext) / (2 * h));
        c.B.push_back(c.A[i - 1] + c.C[i - 1] + p(kf, T[i]) * h);
        c.D.push_back(f(kf, T[i], T0) * h);
    }
}

double fS(const vector<double> T, const double idx, const double T0)
{
    return interpolate(T[idx], TKTable, KTable) * (pow(T[idx], 4) - pow(T0, 4));
}

double Simpson(const vector<double> T, const int N, const int h, const double T0)
{
    double res = 0;
    int halfN = (N % 2) ? (N + 1) / 2 : N / 2;

    for (int i = 0; i < halfN; i++) {
        res += h / 3 * (fS(T, 2 * i, T0) +
                        4 * fS(T, 2 * i + 1, T0) +
                        fS(T, 2 * i + 2, T0));
    }
    return res;
}

double f1(vector<double> T, const int N, const double F0, const double alpha, const double T0)
{
    return F0 - alpha * (T[N - 1] - T0);
}

double f2(const double kf, vector<double> T, const int N, const double h, const double T0)
{
    return kf * Simpson(T, N, h, T0);
}

bool isTemperature(const vector<double> T1, const vector<double> T2, const int N)
{
    bool res = true;
    for (int i = 0; i < N && res; i++) {
        double x = (T1[i] - T2[i]) / T1[i];
        if (abs(x) > eps1) {
            res = false;
        }
    }
    return res;
}

bool isBalanceEnergy(const vector<double> T, const int N,
                     const double kf, const double F0, const double alpha, const double T0, const double h)
{
    bool res = true;
    for (int i = 0; i < N && res; i++) {
        double ff = f1(T, N, F0, alpha, T0);
        double fs = f2(kf, T, N, h, T0);

        double x = (ff - fs) / ff;
        if (abs(x) > eps2) {
            res = false;
        }
    }
    return res;
}

vector<double> progonka(Coefs c, Params p0, Params pN)
{
    int n = c.A.size();
    vector<double> xi, eta;
    xi.push_back(-p0.M / p0.K);
    eta.push_back(p0.P / p0.K);

    for (int i = 0; i < n; i++) {
        double x = c.C[i] / (c.B[i] - c.A[i] * xi[i]);
        double e = (c.D[i] + c.A[i] * eta[i]) / (c.B[i] - c.A[i] * xi[i]);

        xi.push_back(x);
        eta.push_back(e);
    }

    vector<double> y;
    y.push_back((pN.P - pN.K * eta[n]) / (pN.M + pN.K * xi[n]));
    for (int i = n - 1; i >= 0; i--) {
        double yi = xi[i] * y[0] + eta[i];
        y.insert(y.begin(), yi);
    }

    return y;
}

int main()
{
    double Np = 0, l = 0, T0 = 0, sigma = 0, F0 = 0, alpha = 0;
    double h = 1e-3;

    getParams(Np, l, T0, sigma, F0, alpha);
    getTable1(TLTable, LTable);
    getTable2(TKTable, KTable);

    int N = static_cast<int>(l / h);
    double kf = 4 * Np * Np * sigma;

    vector<double> T;
    for (int i = 0; i < N; i++) {
        T.push_back(T0);
    }

    Params p0, pN;
    Coefs c;

    getParams0(p0, T, kf, T0, h, F0);
    getParamsN(pN, T, N, kf, T0, h, alpha);
    getCoefs(T, N, c, kf, T0, h);

    vector<double> Tres = progonka(c, p0, pN);

    int num = 0;
    for (int i = 0; i < 100 && isTemperature(T, Tres, N) && isBalanceEnergy(T, N, kf, F0, alpha, T0, h); i++) {
        T = Tres;

        getParams0(p0, T, kf, T0, h, F0);
        getParamsN(pN, T, N, kf, T0, h, alpha);
        getCoefs(T, N, c, kf, T0, h);

        Tres = add(T, alpha * (sub(progonka(c, p0, pN), T)));
        num = i;
    }

    plot(Tres, num);

    return 0;
}
