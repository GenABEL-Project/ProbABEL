/*
 *
 * Copyright (C) 2009--2014 Various members of the GenABEL team. See
 * the SVN commit logs for more details.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA  02110-1301, USA.
 *
 */


#include <stdio.h>
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include "eigen_mematrix.h"
#include "eigen_mematrix.cpp"
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
