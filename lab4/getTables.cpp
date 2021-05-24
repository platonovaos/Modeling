#include "getTables.h"
#include <fstream>
#include <iostream>

using namespace std;

const string pathTable = "/home/platosha/Desktop/BMSTU/6sem/Modeling/lab4/lab4/table.txt";
const string pathParams = "/home/platosha/Desktop/BMSTU/6sem/Modeling/lab4/lab4/params.txt";

void getTable(InParams &ip)
{
    ifstream fin(pathTable);

    if (fin.is_open()) {
        double par1, par2;

        if (!fin.eof()) {
            fin >> par1 >> par2;
            ip.a[0] = par1;
            ip.a[1] = par2;

            fin >> par1 >> par2;
            ip.b[0] = par1;
            ip.b[1] = par2;

            fin >> par1 >> par2;
            ip.c[0] = par1;
            ip.c[1] = par2;

            fin >> par1 >> par2;
            ip.m[0] = par1;
            ip.m[1] = par2;
        }
    }

    fin.close();
}

void getParams(double &alpha0, double &alphaN, double &l, double &T0, double &R, double &F)
{
    ifstream fin(pathParams);

    if (fin.is_open()) {
        if (!fin.eof()) {
            fin >> alpha0 >> alphaN >> l >> T0 >> R >> F;
        }
    }

    fin.close();
}
