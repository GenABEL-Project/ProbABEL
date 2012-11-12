#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <stdio.h>
#include <getopt.h>
#if EIGEN
#include "mematrix.h"
#include "mematri1.h"
#else
#include "eigen_mematrix.h"
#include "eigen_mematrix.cpp"
#endif
#include "data.h"
#include "reg1.h"
#include "usage.h"

#define MAXITER 10
#define EPS 1.e-8

void print_dmatrix(double **matr, int nrow, int ncol)
{
    for (int j = 0; j < ncol; j++)
    {
        std::cout << "nc=" << j << ":";
        for (int i = 0; i < nrow; i++)
            std::cout << "\t" << matr[j][i];
        std::cout << "\n";
    }
}

int main(int argc, char * argv[])
{
}
