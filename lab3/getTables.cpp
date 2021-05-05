#include "getTables.h"
#include <fstream>
#include <iostream>

using namespace std;

const string pathTable1 = "/home/platosha/Desktop/BMSTU/6sem/Modeling/lab3/lab3/table1.txt";
const string pathTable2 = "/home/platosha/Desktop/BMSTU/6sem/Modeling/lab3/lab3/table2.txt";
const string pathParams = "/home/platosha/Desktop/BMSTU/6sem/Modeling/lab3/lab3/params.txt";

void getTable1(vector<double> &TLTable, vector<double> &LTable)
{
    ifstream fin(pathTable1);

    if (fin.is_open()) {
        double T = 0, lambda = 0;

        while (!fin.eof()) {
            fin >> T >> lambda;

            TLTable.push_back(T);
            LTable.push_back(lambda);
        }
    }

    fin.close();

    TLTable.pop_back();
    LTable.pop_back();
}

void getTable2(vector<double> &TKTable, vector<double> &KTable)
{
    ifstream fin(pathTable2);

    if (fin.is_open()) {
        double T = 0, k = 0;

        while (!fin.eof()) {
            fin >> T >> k;

            TKTable.push_back(T);
            KTable.push_back(k);
        }
    }

    fin.close();

    TKTable.pop_back();
    KTable.pop_back();
}

void getParams(double &Np, double &l, double &T0, double &sigma, double &F0, double &alpha)
{
    ifstream fin(pathParams);

    if (fin.is_open()) {
        if (!fin.eof()) {
            fin >> Np >> l >> T0 >> sigma >> F0 >> alpha;
        }
    }

    fin.close();
}
