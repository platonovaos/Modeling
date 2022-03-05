#include <iostream>
#include <iomanip>

#include "picard.h"
#include "methods.h"

using namespace std;

const double h = 1e-6;

void tableOutput(const int numPoints, const double *x,
                const yPicard yp, const double *yEuler, const double *yRK);

int main()
{
    double x0 = 0, xn = 2.5;
    int n = static_cast<int>(ceil((xn - x0) / h) + 1);

    double *x = new double[n];
    for (int i = 0; i < n; i++) {
        x[i] = x0 + h * i;
    }

    yPicard yp = initPicard(n);
    Picard(n, x, yp);

    double *yEuler = new double[n];
    Euler(n, h, x, yEuler);

    double *yRK = new double[n];
    RungeKutta(n, h, x, yRK);

    tableOutput(n, x, yp, yEuler, yRK);

    delete [] yRK;
    delete [] yEuler;
    deletePicard(yp);
    delete [] x;

    return 0;
}

void tableOutput(const int numPoints, const double *x,
                const yPicard yp, const double *yEuler, const double *yRK)
{
    cout << "   X   |  Пикар 1 |  Пикар 2 |  Пикар 3 |  Пикар 4 |   Эйлер  | Рунге-Кутта " << endl;
    cout << "----------------------------------------------------------------------------" << endl;

    for (int i = 0; i < numPoints; i += 50000) {
        cout << fixed << setprecision(4) << x[i] << " |  ";

        cout << fixed << setprecision(4) << yp.y1[i] << "  |  ";
        cout << fixed << setprecision(4) << yp.y2[i] << "  |  ";
        cout << fixed << setprecision(4) << yp.y3[i] << "  |  ";
        cout << fixed << setprecision(4) << yp.y4[i] << "  |  ";

        cout << fixed << setprecision(4) << yEuler[i] << "  |  ";
        cout << fixed << setprecision(4) << yRK[i] << endl;
    }
    cout << endl;
}
