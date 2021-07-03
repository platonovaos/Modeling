#include <iostream>
#include <cmath>
#include <iomanip>

#include "getTables.h"
#include "interpolation.h"

using namespace std;

void plot(vector<double>, vector<double>);

struct Params
{
    double K, M, P;
};

struct Coefs
{
    vector<double> A, B, D, F;
};


static vector<double> XPlot, TPlot, TimePlot;
static double eps = 1e-4;


double k(const double T, const InParams p)
{
    return p.a[0] * (p.b[0] + p.c[0] * pow(T, p.m[0]));
}

double c(const double T, const InParams p)
{
    return p.a[1] + p.b[1] * pow(T, p.m[1]) - p.c[1] / pow(T, 2);
}

double alpha(const double x,
             const double alpha0, const double alphaN, const double l)
{
    double d = (alphaN * l) / (alphaN - alpha0);
    double c = - alpha0 * d;

    return c / (x - d);
}

double p(const double x,
         const double alpha0, const double alphaN, const double l, const double R)
{
    return 2 * alpha(x, alpha0, alphaN, l) / R;
}

double f(const double x,
         const double alpha0, const double alphaN, const double l, const double T0, const double R)
{
    return 2 * T0 * alpha(x, alpha0, alphaN, l) / R;
}


Params getParams0(const vector<double> TCap, const vector<double> T,
                  const InParams ip, const double h, const double tao,
                  const double alpha0, const double alphaN, const double T0,
                  const double F, const double l, const double R)
{
    double c0 = c(TCap[0], ip);
    double c1 = c(TCap[1], ip);
    double c12 = (c0 + c1) / 2;

    double hi0 = k(TCap[0], ip);
    double hi1 = k(TCap[1], ip);
    double hi12 = (hi0 + hi1) / 2;

    double p0 = p(0, alpha0, alphaN, l, R);
    double p12 = p(h / 2, alpha0, alphaN, l, R);

    double f0 = f(0, alpha0, alphaN, l, T0, R);
    double f12 = f(h/2, alpha0, alphaN, l, T0, R);

    Params p;
    p.K = h/8 * c12 + h/4 * c0 + tao/h * hi12 + h*tao/8 * p12 + h*tao/4 * p0;
    p.M = h/8 * c12 - tao/h * hi12 + h*tao/8 * p12;
    p.P = h/8 * c12 * (T[0] + T[1]) + h/4 * c0 * T[0] + F * tao + h*tao/4 * (f0 + f12);

    return p;
}

Params getParamsN(const vector<double> TCap, const vector<double> T, const double N,
                  const InParams ip, const double h, const double tao,
                  const double alpha0, const double alphaN, const double T0,
                  const double l, const double R)
{
    double cN = c(TCap[N - 1], ip);
    double cN1 = c(TCap[N - 2], ip);
    double cN12 = (cN1 + cN) / 2;

    double hiN = k(TCap[N - 1], ip);
    double hiN1 = k(TCap[N - 2], ip);
    double hiN12 = (hiN1 + hiN) / 2;

    double pN = p(1, alpha0, alphaN, l, R);
    double pN12 = p(1 - h / 2, alpha0, alphaN, l, R);

    double fN = f(1, alpha0, alphaN, l, T0, R);
    double fN12 = f(1 - h/2, alpha0, alphaN, l, T0, R);

    Params p;
    p.K = h/8 * cN12 + h/4 * cN + tao/h * hiN12 + tao*alphaN + h*tao/8 * pN12 + h*tao/4 * pN;
    p.M = h/8 * cN12 - tao/h * hiN12 + h*tao/8 * pN12;
    p.P = h/8 * cN12 * (T[N - 2] + T[N - 1]) + h/4 * cN * T[N - 1] + tao * alphaN * T0 + h*tao/4 * (fN + fN12);

    return p;
}

Coefs getCoefs(const vector<double> TCap, const vector<double> T, const double N,
               const InParams ip, const double h, const double tao,
               const double alpha0, const double alphaN, const double T0,
               const double l, const double R)
{
    Coefs cf;

    for (int i = 0; i < N; i++) {
        double hiPrev = k(TCap[i - 1], ip);
        double hiCur = k(TCap[i], ip);
        double hiNext = k(TCap[i + 1], ip);

        cf.A.push_back(tao/h * (hiPrev + hiCur) / 2);
        cf.D.push_back(tao/h * (hiCur + hiNext) / 2);
        cf.B.push_back(cf.A[i - 1] + cf.D[i - 1] + h * c(TCap[i], ip) + h*tao * p(i * h, alpha0, alphaN, l, R));
        cf.F.push_back(tao*h * f(i * h, alpha0, alphaN, l, T0, R) + h * c(T[i], ip) * T[i]);
    }

    return cf;

}

bool isTempBreak(const vector<double> T1, const vector<double> T2, const int N)
{
    bool exit = false;
    for (int i = 0; i < N && !exit; i++) {
        double x = (T1[i] - T2[i]) / T1[i];
        if (abs(x) < eps) {
            exit = true;
        }
    }
    return exit;
}

vector<double> getT(const Coefs c, const Params p0, const Params pN)
{
    vector<double> eps, eta;
    eps.push_back(-p0.M / p0.K);
    eta.push_back(p0.P / p0.K);

    double x = h;
    int n = 1;
    for (int i = 0; i < n; i++) {
        double ep = c.D[i] / (c.B[i] - c.A[i] * eps[i]);
        double et = (c.F[i] + c.A[i] * eta[i]) / (c.B[i] - c.A[i] * eps[i]);

        eps.push_back(ep);
        eta.push_back(et);
    }

    vector<double> T;
    T.push_back((pN.P - pN.K * eta[n]) / (pN.M + pN.K * eps[n]));
    for (int i = n - 1; i >= 0; i--) {
        double yi = eps[i] * T[i] + eta[i];
        T.insert(T.begin(), yi);
    }

    return T;
}

int main()
{
    InParams ip;
    double alpha0, alphaN, l, T0, R, F;
    double tao = 1;
    double h = 1e-2;

    getParams(alpha0, alphaN, l, T0, R, F);
    getTable(ip);

    int N = static_cast<int>(l / h);

    vector<double> TCap, T;
    for (int i = 0; i < N; i++) {
       TCap.push_back(0);
       T.push_back(T0);
    }

    Params p0, pN;
    Coefs cf;

    bool epsFlag = true;
    while (epsFlag) {
        p0 = getParams0(TCap, T, ip, h, tao, alpha0, alphaN, T0, F, l, R);
        pN = getParamsN(TCap, T, N, ip, h, tao, alpha0, alphaN, T0, l, R);
        T = getT(cf, p0, pN);
    }

    //plot(T, XPlot);
    //cout << "n = " << n "Here\n";

    return 0;
}
