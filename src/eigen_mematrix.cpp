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


#ifndef EIGEN_MEMATRI1_H
#define EIGEN_MEMATRI1_H
#include "eigen_mematrix.h"
#include <iostream>
#include <string>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>

using namespace Eigen;
//
// constructors
//

template<class DT>
mematrix<DT>::mematrix(int nr, int nc)
{
    if (nr <= 0)
    {
        std::cerr << "mematrix(): nr <= 0\n";
        exit(1);
    }
    if (nc <= 0)
    {
        std::cerr << "mematrix(): nc <= 0\n";
        exit(1);
    }
    this->nrow = nr;
    this->ncol = nc;
    this->nelements = nr * nc;
    this->data.resize(nr, nc);
}

template<class DT>
mematrix<DT>::mematrix(const mematrix<DT> & M)
{
    ncol = M.ncol;
    nrow = M.nrow;
    nelements = M.nelements;
    data = M.data;
}

//
// operators
//
template<class DT>
mematrix<DT> &mematrix<DT>::operator=(const mematrix<DT> &M)
{
    if (this != &M)
    {
        ncol = M.ncol;
        nrow = M.nrow;
        nelements = M.nelements;
        data = M.data;
    }
    return *this;
}

template<class DT>
DT & mematrix<DT>::operator[](const int i)
{
    if (i < 0 || i >= (ncol * nrow))
    {
        std::cerr << "mematrix[]: " << i << " out of bounds (0,"
                  << nrow * ncol - 1 << ")\n";
        exit(1);
    }
    int column = i % ncol;
    int row = static_cast<int>( floor(static_cast<double>(i / ncol)) );

    return data(row, column);
}

// template<class DT>
// mematrix<DT> mematrix<DT>::operator+(DT toadd)
// {
//     mematrix<DT> temp(nrow, ncol);
//     for (int i = 0; i < nelements; i++)
//         temp.data[i] = data[i] + toadd;
//     return temp;
// }

template<class DT>
mematrix<DT> mematrix<DT>::operator+(const mematrix<DT> &M)
{
    if (ncol != M.ncol || nrow != M.nrow)
    {
        std::cerr << "mematrix+: matrices not equal in size ("
                  << nrow << "," << ncol << ") and ("
                  << M.nrow << "," << M.ncol << ")\n";
        exit(1);
    }
    mematrix<DT> temp;
    temp.data = data + M.data;
    temp.ncol = data.cols();
    temp.nrow = data.rows();
    temp.nelements = temp.nrow * temp.ncol;

    return temp;
}

template<class DT>
mematrix<DT> mematrix<DT>::operator-(const DT toadd)
{
    mematrix<DT> temp(nrow, ncol);
    temp.data = data.array() - toadd;
    return temp;
}

template<class DT>
mematrix<DT> mematrix<DT>::operator-(const mematrix<DT> &M)
{
    if (ncol != M.ncol || nrow != M.nrow)
    {
        std::cerr << "mematrix-: matrices not equal in size ("
                  << nrow << "," << ncol << ") and ("
                  << M.nrow << "," << M.ncol << ")\n";
        exit(1);
    }
    mematrix<DT> temp;
    temp.data = data - M.data;
    temp.ncol = temp.data.cols();
    temp.nrow = temp.data.rows();
    temp.nelements = temp.nrow * temp.ncol;
    return temp;
}

template<class DT>
mematrix<DT> mematrix<DT>::operator*(DT multiplyer)
{
//    MatrixXd add = MatrixXd::Constant(nrow, ncol, toadd);
    mematrix<DT> temp(nrow, ncol);
    temp.data = data * multiplyer;
    return temp;
}

template<class DT>
mematrix<DT> mematrix<DT>::operator*(const mematrix<DT> &M)
{
    if (ncol != M.nrow)
    {
        std::cerr << "mematrix*: ncol != nrow (" << nrow << ","
                  << ncol << ") and (" << M.nrow << "," << M.ncol << ")\n";
    }

    mematrix<DT> temp;
    temp.data = data * M.data;
    temp.ncol = temp.data.cols();
    temp.nrow = temp.data.rows();
    temp.nelements = temp.nrow * temp.ncol;
    // std::cerr << "mematrix*:  (" << nrow << "," << ncol << ") and ("
    //           << M.nrow << "," << M.ncol << "): result"
    //           << temp.nrow * temp.ncol << "\n";
    // std::cout.flush();

    return temp;
}

template<class DT>
mematrix<DT> mematrix<DT>::operator*(const mematrix<DT> *M)
{
    if (ncol != M->nrow)
    {
        std::cerr << "mematrix*: ncol != nrow (" << nrow << "," << ncol
                  << ") and (" << M->nrow << "," << M->ncol <<")\n";
    }
    mematrix<DT> temp;
    temp.data = data * M->data;
    temp.ncol = temp.data.cols();
    temp.nrow = temp.data.rows();
    temp.nelements = temp.nrow * temp.ncol;
//    std::cerr << "mematrix*:  (%d,%d) and (%d,%d):result%d\n", nrow,
//            ncol, M->nrow, M->ncol,temp.nrow * temp.ncol);

    return temp;
}

//
// operations
//
template<class DT>
void mematrix<DT>::reinit(int nr, int nc)
{
//    if (nelements > 0)
//        delete[] data;
    if (nr <= 0)
    {
        std::cerr << "mematrix(): number of rows less than 1\n";
        exit(1);
    }
    if (nc <= 0)
    {
        std::cerr << "mematrix(): number of columns less than 1\n";
        exit(1);
    }
    nrow = nr;
    ncol = nc;
    nelements = nr * nc;
//    std::cout << "[DEBUG] resizing" << std::endl;

    data.resize(nr, nc);
    data.setZero();
}

template<class DT>
DT mematrix<DT>::get(int nr, int nc)
{
#ifndef NDEBUG
    if (nc < 0 || nc > ncol -1)
    {
        std::cerr << "mematrix::get: column out of range: " << nc + 1
                  << " not between (1," << ncol << ")\n" << std::flush;
        exit(1);
    }
    if (nr < 0 || nr > nrow -1)
    {
        std::cerr << "mematrix::get: row out of range: " << nr + 1
                  << " not between (1," << nrow << ")\n" << std::flush;
        exit(1);
    }
#endif
    DT temp = data(nr, nc);
    return temp;
}

template<class DT>
void mematrix<DT>::put(DT value, int nr, int nc)
{
#ifndef NDEBUG
    if (nc < 0 || nc > ncol -1)
    {
        std::cerr << "mematrix::put: column out of range: " << nc + 1
                  << " not between (1," << ncol << ")\n" << std::flush;
        exit(1);
    }
    if (nr < 0 || nr > nrow -1)
    {
        std::cerr << "mematrix::put: row out of range: " << nr + 1
                  << " not between (1," << nrow << ")\n" << std::flush;
        exit(1);
    }
#endif
    data(nr, nc) = value;
}

template<class DT>
DT mematrix<DT>::column_mean(int nc)
{
    if (nc >= ncol || nc < 0)
    {
        std::cerr << "colmM bad column\n";
        exit(1);
    }

    return data.col(nc).mean();;
}

template<class DT>
mematrix<DT> column_sum(const mematrix<DT> &M)
{
    mematrix<DT> out;
    out.reinit(1, M.ncol);
    out.data = M.data.colwise().mean();
    return out;
}

template<class DT>
void mematrix<DT>::print(void)
{
    std::cout << "nrow=" << nrow << "; ncol=" << ncol
         << "; nelements=" << nelements << "\n";
    for (int i = 0; i < nrow; i++)
    {
        std:: cout << "nr=" << i << ":\t";
        for (int j = 0; j < ncol; j++)
        {
            printf("%1.9e\t", data.data()[i * ncol + j]);
        }
        std::cout << "\n";
    }
}

//
// other functions
//
template<class DT>
mematrix<DT> transpose(const mematrix<DT> &M)
{
    // cout << "[DEBUG TRANSPOSE PRE]nrow=" << M.nrow << "; ncol="
    //      << M.ncol << "; nelements=" << M.nelements;

    mematrix<DT> temp;
    temp.data = M.data.transpose();
    temp.ncol = M.nrow;
    temp.nrow = M.ncol;
    temp.nelements = M.nelements;
    // cout << "[DEBUG TRANSPOSE post]nrow=" << temp.nrow << "; ncol="
    //      << temp.ncol << "; nelements=" << temp.nelements;

    return temp;
}

template<class DT>
mematrix<DT> reorder(const mematrix<DT> &M, const mematrix<int> order)
{
    if (M.nrow != order.nrow)
    {
        std::cerr << "reorder: M & order have different # of rows\n";
        exit(1);
    }
    mematrix<DT> temp(M.nrow, M.ncol);

    for (int i = 0; i < temp.nrow; i++)
    {
        int source = order.data(i, 0);
        temp.data.row(source) = M.data.row(i);
    }
    return temp;
}

// template<class DT>
// mematrix<double> todouble(mematrix<DT> &M)
// {
//     mematrix<double> temp(M.nrow, M.ncol);
//     for (int i = 0; i < temp.nelements; i++)
//         temp.data[i] = double(M.data[i]);
//     return temp;
// }

template<class DT>
mematrix<DT> invert(const mematrix<DT> &M)
{
    if (M.ncol != M.nrow)
    {
        std::cerr << "invert: only square matrices possible\n";
        exit(1);
    }

    mematrix<DT> temp = M;

    temp.data = temp.data.inverse();
    return temp;
}

template<class DT>
mematrix<DT> productMatrDiag(const mematrix<DT> &M, const mematrix<DT> &D)
{
    // multiply all rows of M by value of first row of D
    if (M.ncol != D.nrow)
    {
        std::cerr << "productMatrDiag: wrong dimensions";
        exit(1);
    }
    mematrix<DT> temp = M;
    // make an array of the first row of D in the same way orientation as
    // M.data.row(i).array()
    Array<DT, Dynamic, Dynamic> row = D.data.block(0, 0, M.ncol, 1).transpose();

    for (int i = 0; i < temp.nrow; i++)
    {
        temp.data.row(i) = M.data.row(i).array() * row;
    }
    return temp;
}


template<class DT>
void mematrix<DT>::delete_column(const int delcol)
{
    if (delcol > ncol || delcol < 0)
    {
        fprintf(stderr, "mematrix::delete_column: column out of range\n");
        exit(1);
    }

    Matrix<DT, Dynamic, Dynamic, RowMajor> auxdata = data;
    data.resize(data.rows(), data.cols()-1);

    int rightColsSize = auxdata.cols() - delcol - 1;

    data.leftCols(delcol) = auxdata.leftCols(delcol);
    data.rightCols(rightColsSize) = auxdata.rightCols(rightColsSize);

    ncol--;
}


template<class DT>
void mematrix<DT>::delete_row(const int delrow)
{
    if (delrow > nrow || delrow < 0)
    {
        fprintf(stderr, "mematrix::delete_row: row out of range\n");
        exit(1);
    }

    Matrix<DT, Dynamic, Dynamic, RowMajor> auxdata = data;
    data.resize(data.rows()-1, data.cols());

    int bottomRowsSize = auxdata.rows() - delrow - 1;

    data.topRows(delrow) = auxdata.topRows(delrow);
    data.bottomRows(bottomRowsSize) = auxdata.bottomRows(bottomRowsSize);

    nrow--;
}



#endif
