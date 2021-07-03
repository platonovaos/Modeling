#include "getTables.h"
#include <fstream>
#include <iostream>

using namespace std;

const string pathTable1 = "/home/platosha/Desktop/BMSTU/6sem/Modeling/lab2/lab2/table1.txt";
const string pathTable2 = "/home/platosha/Desktop/BMSTU/6sem/Modeling/lab2/lab2/table2.txt";
const string pathParams = "/home/platosha/Desktop/BMSTU/6sem/Modeling/lab2/lab2/params.txt";

void getTable1(vector<double> &ITable, vector<double> &ToTable, vector<double> &mTable)
{
    ifstream fin(pathTable1);

    if (fin.is_open()) {
        double I = 0, To = 0, m = 0;

        while (!fin.eof()) {
            fin >> I >> To >> m;

            ITable.push_back(I);
            ToTable.push_back(To);
            mTable.push_back(m);
        }
    }

    fin.close();

    ITable.pop_back();
    ToTable.pop_back();
    mTable.pop_back();
}

void getTable2(vector<double> &TTable, vector<double> &sigTable)
{
    ifstream fin(pathTable2);

    if (fin.is_open()) {
        double T = 0, sig = 0;

        while (!fin.eof()) {
            fin >> T >> sig;

            TTable.push_back(T);
            sigTable.push_back(sig);
        }
    }

    fin.close();

    TTable.pop_back();
    sigTable.pop_back();
}

void getParams(double &R, double &l, double &Lk, double &Ck, double &Rk,
               double &Uco, double &I0, double &Tw)
{
    ifstream fin(pathParams);

    if (fin.is_open()) {
        if (!fin.eof()) {
            fin >> R >> l >> Lk >> Ck >> Rk >> Uco >> I0 >> Tw;
        }
    }

    fin.close();
}
