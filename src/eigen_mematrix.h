/*
 *
 * Copyright (C) 2009--2015 Various members of the GenABEL team. See
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


    // Constructors and destructors
    mematrix()
    {
        nrow = ncol = nelements = 0;
        data.resize(1, 1);
    }
    mematrix(const int nr, const int nc);
    mematrix(const mematrix &M);
    // ~mematrix()


    // Operator overloading
    mematrix & operator=(const mematrix &M);
    const DT & operator[](const int i) const; /* Subscript for reading */
    DT & operator[](const int i);             /* Subscript for writing */
    mematrix operator+(const mematrix &M);
    mematrix operator-(const DT toadd);
    mematrix operator-(const mematrix &M);
    mematrix operator*(const DT toadd);
    mematrix operator*(const mematrix &M);
    mematrix operator*(const mematrix *M);


    // Other member functions
    void delete_column(const int delcol);
    void delete_row(const int delrow);

    void reinit(const int nr, const int nc);

    unsigned int getnrow(void) const
    {
        return nrow;
    }
    unsigned int getncol(void) const
    {
        return ncol;
    }
    DT get(const int nr, const int nc) const;
    void put(const DT value, const int nr, const int nc);
    DT column_mean(const int nc) const;
    DT column_sum(const int nc);
    void print(void) const;
};

#endif
