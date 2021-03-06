#include "methods.h"

double f(const double x, const double y)
{
    return (x * x + y * y);
}

void Euler(const int n, const double h, double *x, double *y)
{
    y[0] = 0;
    for (int i = 0; i < n - 1; i++) {
        y[i + 1] = y[i] + h * f(x[i], y[i]);
    }
}

void RungeKutta(const int n, const double h, double *x, double *y)
{
    y[0] = 0;

    double alpha = 0.5;
    double k1 = 0, k2 = 0;

    for (int i = 0; i < n - 1; i++) {
        k1 = f(x[i], y[i]);
        k2 = f(x[i] + h / (2 * alpha), y[i] + h / (2 * alpha) * k1);

        y[i + 1] = y[i] + h * ((1 - alpha) * k1 + alpha * k2);
    }
}
