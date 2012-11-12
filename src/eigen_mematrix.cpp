#ifndef EIGEN_MEMATRI1_H
#define EIGEN_MEMATRI1_H
#include "eigen_mematrix.h"
#include <Eigen/Dense>
#include <Eigen/LU>
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
        fprintf(stderr, "mematrix(): nr <= 0\n");
        exit(1);
    }
    if (nc <= 0)
    {
        fprintf(stderr, "mematrix(): nc <= 0\n");
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
        fprintf(stderr, "mematrix[]: %d out of bounds (0,%d)\n", i,
                nrow * ncol - 1);
        exit(1);
    }
    int column = i % ncol;
    int row = (int) floor((double) i / ncol);

    return data(row, column);
}

//template<class DT>
//mematrix<DT> mematrix<DT>::operator+(DT toadd)
//{
//    mematrix<DT> temp(nrow, ncol);
//    for (int i = 0; i < nelements; i++)
//        temp.data[i] = data[i] + toadd;
//    return temp;
//}
template<class DT>
mematrix<DT> mematrix<DT>::operator+(const mematrix<DT> &M)
{
    if (ncol != M.ncol || nrow != M.nrow)
    {
        fprintf(stderr,
                "mematrix+: matrices not equal in size (%d,%d) and (%d,%d)",
                nrow, ncol, M.nrow, M.ncol);
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
        fprintf(stderr,
                "mematrix-: matrices not equal in size (%d,%d) and (%d,%d)",
                nrow, ncol, M.nrow, M.ncol);
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
        fprintf(stderr, "mematrix*: ncol != nrow (%d,%d) and (%d,%d)", nrow,
                ncol, M.nrow, M.ncol);
    }

    mematrix<DT> temp;
    temp.data = data * M.data;
    temp.ncol = temp.data.cols();
    temp.nrow = temp.data.rows();
    temp.nelements = temp.nrow * temp.ncol;
//    fprintf(stderr, "mematrix*:  (%d,%d) and (%d,%d):result%d\n", nrow,
//            ncol, M.nrow, M.ncol,temp.nrow * temp.ncol);
//	std::cout.flush();

    return temp;
}

template<class DT>
mematrix<DT> mematrix<DT>::operator*(const mematrix<DT> *M)
{
    if (ncol != M->nrow)
    {
        fprintf(stderr, "mematrix*: ncol != nrow (%d,%d) and (%d,%d)", nrow,
                ncol, M->nrow, M->ncol);
    }
    mematrix<DT> temp;
    temp.data = data * M->data;
    temp.ncol = temp.data.cols();
    temp.nrow = temp.data.rows();
    temp.nelements = temp.nrow * temp.ncol;
//    fprintf(stderr, "mematrix*:  (%d,%d) and (%d,%d):result%d\n", nrow,
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
        fprintf(stderr, "mematrix(): number of rows smaller then 1\n");
        exit(1);
    }
    if (nc <= 0)
    {
        fprintf(stderr, "mematrix(): number of columns smaller then 1\n");
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
#if !NDEBUG
    if (nc < 0 || nc > ncol)
    {
        fprintf(stderr,
                "mematrix::get: column out of range: %d not in (0,%d)\n", nc,
                ncol);
        exit(1);
    }
    if (nr < 0 || nr > nrow)
    {
        printf("mematrix::get: row out of range: %d not in (0,%d)\n", nr, nrow);
        exit(1);
    }
#endif
    DT temp = data(nr, nc);
    return temp;
}
template<class DT>
void mematrix<DT>::put(DT value, int nr, int nc)
{
#if !NDEBUG
    if (nc < 0 || nc > ncol)
    {
        fprintf(stderr,
                "mematrix::put: column out of range: %d not in (0,%d)\n", nc,
                ncol);
        exit(1);
    }
    if (nr < 0 || nr > nrow)
    {
        printf("mematrix::put: row out of range: %d not in (0,%d)\n", nr, nrow);
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
        fprintf(stderr, "colmM bad column\n");
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
    cout << "nrow=" << nrow << "; ncol=" << ncol << "; nelements=" << nelements
            << "\n";
    for (int i = 0; i < nrow; i++)
    {
        cout << "nr=" << i << ":\t";
        for (int j = 0; j < ncol; j++)
            cout << data.data()[i * ncol + j] << "\t";
        cout << "\n";
    }
}
// other functions
//

template<class DT>
mematrix<DT> transpose(const mematrix<DT> &M)
{
//cout << "[DEBUG TRANSPOSE PRE]nrow=" << M.nrow << "; ncol=" << M.ncol << "; nelements=" << M.nelements;

    mematrix<DT> temp;
    temp.data = M.data.transpose();
    temp.ncol = M.nrow;
    temp.nrow = M.ncol;
    temp.nelements = M.nelements;
//cout << "[DEBUG TRANSPOSE post]nrow=" << temp.nrow << "; ncol=" << temp.ncol << "; nelements=" << temp.nelements;

    return temp;
}

template<class DT>
mematrix<DT> reorder(const mematrix<DT> &M, const mematrix<int> order)
{
    if (M.nrow != order.nrow)
    {
        fprintf(stderr, "reorder: M & order have differet # of rows\n");
        exit(1);
    }
    mematrix<DT> temp(M.nrow, M.ncol);
//FIXME(maarten): commented out to get compilation running
//    for (int i = 0; i < temp.nrow; i++)
//        for (int j = 0; j < temp.ncol; j++)
//            temp.data[order[i] * temp.ncol + j] = M.data[i * M.ncol + j];
    return temp;
}
//
//
//template<class DT>
//mematrix<double> todouble(mematrix<DT> &M)
//{
//    mematrix<double> temp(M.nrow, M.ncol);
//    for (int i = 0; i < temp.nelements; i++)
//        temp.data[i] = double(M.data[i]);
//    return temp;
//}
//

//

template<class DT>
mematrix<DT> invert(const mematrix<DT> &M)
{
    if (M.ncol != M.nrow)
    {
        fprintf(stderr, "invert: only square matrices possible\n");
        exit(1);
    }

    mematrix<DT> temp = M;

    temp.data = temp.data.inverse();
    return temp;
}

template<class DT>
mematrix<DT> productMatrDiag(const mematrix<DT> &M, const mematrix<DT> &D)
{
    //multiply all rows of M by value of first row of D
    if (M.ncol != D.nrow)
    {
        fprintf(stderr, "productMatrDiag: wrong dimenstions");
        exit(1);
    }
    mematrix<DT> temp = M;
    //make a array of the first row of D in the same way orientation as M.data.row(i).array()
    Array<DT, Dynamic, Dynamic> row = D.data.block(0, 0, M.ncol, 1).transpose();

    for (int i = 0; i < temp.nrow; i++)
    {
        temp.data.row(i) = M.data.row(i).array() * row;
    }
    return temp;
}

#endif
