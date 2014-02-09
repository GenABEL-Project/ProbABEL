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


#ifndef __EIGEN_MEMATRIX_H__
#define __EIGEN_MEMATRIX_H__
#include <Eigen/Dense>
#include <Eigen/LU>
#include <iostream>

using namespace Eigen;
using std::cout;
using std::cerr;


template<class DT> class mematrix {
 public:
    int nrow;
    int ncol;
    int nelements;
    Matrix<DT, Dynamic, Dynamic, RowMajor> data;

    mematrix()
    {
        nrow = ncol = nelements = 0;
        data.resize(1, 1);
    }
    mematrix(int nr, int nc);
    mematrix(const mematrix &M);
    ~mematrix()
    {
//        if (nelements > 0)
//            delete data;
    }

    mematrix & operator=(const mematrix &M);
    DT & operator[](int i);
//    mematrix operator+(DT toadd);
    mematrix operator+(const mematrix &M);
    mematrix operator-(DT toadd);
    mematrix operator-(const mematrix &M);
    mematrix operator*(DT toadd);
    mematrix operator*(const mematrix &M);
    mematrix operator*(const mematrix *M);

    void delete_column(const int delcol);
    void delete_row(const int delrow);

    void reinit(int nr, int nc);

    unsigned int getnrow(void)
    {
        return nrow;
    }
    unsigned int getncol(void)
    {
        return ncol;
    }
    DT get(int nr, int nc);
    void put(DT value, int nr, int nc);
    DT column_mean(int nc);
    DT column_sum(int nc);
    void print(void);
};

#endif
