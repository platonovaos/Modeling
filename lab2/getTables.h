#ifndef GETTABLES_H
#define GETTABLES_H

#include <vector>
using namespace std;

void getTable1(vector<double> &ITable, vector<double> &ToTable, vector<double> &mTable);
void getTable2(vector<double> &TTable, vector<double> &sigTable);

void getParams(double &R, double &l, double &Lk, double &Ck, double &Rk,
               double &Uco, double &I0, double &Tw);

#endif // GETTABLES_H
