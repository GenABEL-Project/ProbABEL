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


#ifndef MEMATRI1_H
#define MEMATRI1_H

#include <cstdlib>
#include <string>
#include <cstdarg>
#include <cstdio>

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
    nrow = nr;
    ncol = nc;
    nelements = nr * nc;
    data = new (nothrow) DT[ncol * nrow];
    if (!data)
    {
        fprintf(stderr, "mematrix(nr,nc): cannot allocate memory (%d,%d)\n",
                nrow, ncol);
        exit(1);
    }
}

template<class DT>
mematrix<DT>::mematrix(const mematrix<DT> & M)
{
    ncol = M.ncol;
    nrow = M.nrow;
    nelements = M.nelements;
    data = new (nothrow) DT[M.ncol * M.nrow];
    if (!data)
    {
        fprintf(stderr,
                "mematrix const(mematrix): cannot allocate memory (%d,%d)\n",
                M.nrow, M.ncol);
        exit(1);
    }
    //	std::cerr << "mematrix const(mematrix): can allocate memory ("
    //            << M.nrow << "," << M.ncol << ")\n";
    for (int i = 0; i < M.ncol * M.nrow; i++)
        data[i] = M.data[i];
}

//
// operators
//
template<class DT>
mematrix<DT> &mematrix<DT>::operator=(const mematrix<DT> &M)
{
    if (this != &M)
    {
        if (data != NULL)
            delete[] data;
        data = new (nothrow) DT[M.ncol * M.nrow];
        if (!data)
        {
            fprintf(stderr, "mematrix=: cannot allocate memory (%d,%d)\n",
                    M.nrow, M.ncol);
            delete[] data;
            exit(1);
        }
        ncol = M.ncol;
        nrow = M.nrow;
        nelements = M.nelements;
        for (int i = 0; i < M.ncol * M.nrow; i++)
        {
            data[i] = M.data[i];
        }
    }
    return *this;
}

template<class DT>
DT &mematrix<DT>::operator[](int i)
{
    if (i < 0 || i >= (ncol * nrow))
    {
        fprintf(stderr, "mematrix[]: %d out of bounds (0,%d)\n", i,
                nrow * ncol - 1);
        exit(1);
    }
    return data[i];
}

template<class DT>
mematrix<DT> mematrix<DT>::operator+(DT toadd)
{
    mematrix<DT> temp(nrow, ncol);
    for (int i = 0; i < nelements; i++)
        temp.data[i] = data[i] + toadd;
    return temp;
}

template<class DT>
mematrix<DT> mematrix<DT>::operator+(mematrix<DT> &M)
{
    if (ncol != M.ncol || nrow != M.nrow)
    {
        fprintf(stderr,
                "mematrix+: matrices not equal in size (%d,%d) and (%d,%d)",
                nrow, ncol, M.nrow, M.ncol);
        exit(1);
    }
    mematrix<DT> temp(nrow, ncol);
    for (int i = 0; i < nelements; i++)
        temp.data[i] = data[i] + M.data[i];
    return temp;
}

template<class DT>
mematrix<DT> mematrix<DT>::operator-(DT toadd)
{
    mematrix<DT> temp(nrow, ncol);
    for (int i = 0; i < nelements; i++)
        temp.data[i] = data[i] - toadd;
    return temp;
}

template<class DT>
mematrix<DT> mematrix<DT>::operator-(mematrix<DT> &M)
{
    if (ncol != M.ncol || nrow != M.nrow)
    {
        fprintf(stderr,
                "mematrix-: matrices not equal in size (%d,%d) and (%d,%d)",
                nrow, ncol, M.nrow, M.ncol);
        exit(1);
    }
    mematrix<DT> temp(nrow, ncol);
    for (int i = 0; i < nelements; i++)
        temp.data[i] = data[i] - M.data[i];
    return temp;
}

template<class DT>
mematrix<DT> mematrix<DT>::operator*(DT toadd)
{
    // A che naschet std::string vmesto DT? Maksim.
    mematrix<DT> temp(nrow, ncol);
    for (int i = 0; i < nelements; i++)
        temp.data[i] = data[i] * toadd;
    return temp;
}

template<class DT>
mematrix<DT> mematrix<DT>::operator*(mematrix<DT> &M)
{
    if (ncol != M.nrow)
    {
        fprintf(stderr, "mematrix*: ncol != nrow (%d,%d) and (%d,%d)", nrow,
                ncol, M.nrow, M.ncol);
        exit(1);
    }
    mematrix<DT> temp(nrow, M.ncol);
    for (int j = 0; j < temp.nrow; j++)
    {
        for (int i = 0; i < temp.ncol; i++)
        {
            DT sum = 0;
            for (int j1 = 0; j1 < ncol; j1++)
                sum += data[j * ncol + j1] * M.data[j1 * M.ncol + i];
            temp[j * temp.ncol + i] = sum;
        }
    }
    return temp;
}

template<class DT>
mematrix<DT> mematrix<DT>::operator*(mematrix<DT> *M)
{
    if (ncol != M->nrow)
    {
        fprintf(stderr, "mematrix*: ncol != nrow (%d,%d) and (%d,%d)", nrow,
                ncol, M->nrow, M->ncol);
        exit(1);
    }
    mematrix<DT> temp(nrow, M->ncol);
    for (int j = 0; j < temp.nrow; j++)
    {
        for (int i = 0; i < temp.ncol; i++)
        {
            DT sum = 0;
            for (int j1 = 0; j1 < ncol; j1++)
                sum += data[j * ncol + j1] * M->data[j1 * M->ncol + i];
            temp[j * temp.ncol + i] = sum;
        }
    }
    return temp;
}

//
// operations
//
template<class DT>
void mematrix<DT>::reinit(int nr, int nc)
{
    if (nelements > 0)
        delete[] data;
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
    data = new (nothrow) DT[ncol * nrow];
    if (!data)
    {
        fprintf(stderr, "mematrix(nr,nc): cannot allocate memory (%d,%d)\n",
                nrow, ncol);
        exit(1);
    }
}

template<class DT>
DT mematrix<DT>::get(int nr, int nc)
{
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
    DT temp = data[nr * ncol + nc];
    return temp;
}

template<class DT>
void mematrix<DT>::put(DT value, int nr, int nc)
{
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
    data[nr * ncol + nc] = value;
}

template<class DT>
DT mematrix<DT>::column_mean(int nc)
{
    if (nc >= ncol || nc < 0)
    {
        fprintf(stderr, "colmM bad column\n");
        exit(1);
    }
    DT out = 0.0;
    for (int i = 0; i < nrow; i++)
        out += DT(data[i * ncol + nc]);
    out /= DT(nrow);
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
        {
            printf("%e\t", data[i * ncol + j]);
        }
        cout << "\n";
    }
}

template<class DT>
void mematrix<DT>::delete_column(const int delcol)
{
    if (delcol > ncol || delcol < 0)
    {
        fprintf(stderr, "mematrix::delete_column: column out of range\n");
        exit(1);
    }
    mematrix<DT> temp = *this;
    if (nelements > 0)
        delete[] data;
    ncol--;
    nelements = ncol * nrow;
    data = new (nothrow) DT[ncol * nrow];
    if (!data)
    {
        fprintf(stderr,
                "mematrix::delete_column: cannot allocate memory (%d,%d)\n",
                nrow, ncol);
        delete[] data;
        exit(1);
    }
    int newcol = 0;
    for (int nr = 0; nr < temp.nrow; nr++)
    {
        newcol = 0;
        for (int nc = 0; nc < temp.ncol; nc++)
            if (nc != delcol)
                data[nr * ncol + (newcol++)] = temp[nr * temp.ncol + nc];
    }
}

template<class DT>
void mematrix<DT>::delete_row(const int delrow)
{
    if (delrow > nrow || delrow < 0)
    {
        fprintf(stderr, "mematrix::delete_row: row out of range\n");
        exit(1);
    }
    mematrix<DT> temp = *this;
    if (nelements > 0)
        delete[] data;
    nrow--;
    nelements = ncol * nrow;
    data = new (nothrow) DT[ncol * nrow];
    if (!data)
    {
        fprintf(stderr,
                "mematrix::delete_row: cannot allocate memory (%d,%d)\n", nrow,
                ncol);
        delete[] data;
        exit(1);
    }
    int newrow = 0;
    for (int nc = 0; nc < temp.ncol; nc++)
    {
        newrow = 0;
        for (int nr = 0; nr < temp.nrow; nr++)
            if (nr != delrow)
                data[nr * ncol + (newrow++)] = temp[nr * temp.ncol + nc];
    }
}

//
// other functions
//
template<class DT>
mematrix<DT> column_sum(mematrix<DT> &M)
{
    mematrix<DT> out;
    out.reinit(1, M.ncol);
    for (int j = 0; j < M.ncol; j++)
    {
        DT sum = 0;
        for (int i = 0; i < M.nrow; i++)
            sum = sum + DT(M.data[i * M.ncol + j]);
        out.put(sum, 0, j);
    }
    return out;
}

template<class DT>
mematrix<DT> column_mean(mematrix<DT> &M)
{
    mematrix<DT> out;
    out.reinit(1, M.ncol);
    for (int j = 0; j < M.ncol; j++)
    {
        DT sum = 0;
        for (int i = 0; i < M.nrow; i++)
            sum = sum + DT(M.data[i * M.ncol + j]);
        sum /= DT(M.nrow);
        out.put(sum, 0, j);
    }
    return out;
}

template<class DT>
mematrix<DT> transpose(mematrix<DT> &M)
{
    mematrix<DT> temp(M.ncol, M.nrow);
    for (int i = 0; i < temp.nrow; i++)
        for (int j = 0; j < temp.ncol; j++)
            temp.data[i * temp.ncol + j] = M.data[j * M.ncol + i];
    return temp;
}

template<class DT>
mematrix<DT> reorder(mematrix<DT> &M, mematrix<int> order)
{
    if (M.nrow != order.nrow)
    {
        std::cerr << "reorder: M & order have different # of rows\n";
        exit(1);
    }
    mematrix<DT> temp(M.nrow, M.ncol);
    for (int i = 0; i < temp.nrow; i++)
        for (int j = 0; j < temp.ncol; j++)
            temp.data[order[i] * temp.ncol + j] = M.data[i * M.ncol + j];
    return temp;
}

template<class DT>
mematrix<DT> productMatrDiag(mematrix<DT> &M, mematrix<DT> &D)
{
    //multiply all rows of M by value of first row of D
    if (M.ncol != D.nrow)
    {
        fprintf(stderr, "productMatrDiag: wrong dimenstions");
        exit(1);
    }
    mematrix<DT> temp(M.nrow, M.ncol);
    for (int i = 0; i < temp.nrow; i++){
        for (int j = 0; j < temp.ncol; j++){
            temp.data[i * temp.ncol + j] = M.data[i * M.ncol + j] * D.data[j];
    //			temp.put(M.get(i,j)*D.get(j,0),i,j);
        }
    }
    return temp;
}

template<class DT>
mematrix<double> todouble(mematrix<DT> &M)
{
    mematrix<double> temp(M.nrow, M.ncol);
    for (int i = 0; i < temp.nelements; i++)
        temp.data[i] = double(M.data[i]);
    return temp;
}

template<class DT>
mematrix<DT> productXbySymM(mematrix<DT> &X, mematrix<DT> &M)
{
    if (M.ncol < 1 || M.nrow < 1 || X.ncol < 1 || X.nrow < 1)
    {
        fprintf(stderr,
                "productXbySymM: M.ncol<1 || M.nrow<1 || X.ncol<1 || X.nrow < 1\n");
        exit(1);
    }
    if (M.ncol != M.nrow)
    {
        fprintf(stderr, "productXbySymM: M.ncol != M.nrow\n");
        exit(1);
    }
    if (M.ncol != X.ncol)
    {
        fprintf(stderr, "productXbySymM: M.ncol != X.ncol\n");
        exit(1);
    }
    if (M.ncol != X.ncol)
    {
        fprintf(stderr, "productXbySymM: M.ncol != X.ncol\n");
        exit(1);
    }

    mematrix<DT> out(X.nrow, X.ncol);
    int i, j, k;

    double temp1, temp2, value1, value2; // not good should be of <DT>!
    for (k = 0; k < X.nrow; k++)
    {
        temp1 = 0.;
        for (i = 0; i < X.ncol; i++)
        {
            temp1 = X.get(k, i);
            temp2 = 0.;
            for (j = (i + 1); j < X.ncol; j++)
            {
                value1 = out.get(k, j) + temp1 * M.get(i, j);
                out.put(value1, k, j);
                temp2 += M.get(i, j) * X.get(k, j);
            }
            value2 = out.get(k, i) + temp2 + M.get(i, i) * X.get(k, i);
            out.put(value2, k, i);
        }
    }

    return out;
}

// written by Mike Dinolfo 12/98
// modified Yurii Aulchenko 2008-04-22
template<class DT>
mematrix<DT> invert(mematrix<DT> &M)
{
    if (M.ncol != M.nrow)
    {
        fprintf(stderr, "invert: only square matrices possible\n");
        exit(1);
    }
    if (M.ncol == 1)
    {
        mematrix<DT> temp(1, 1);
        temp[0] = 1. / M[0];
    }
    /*
     for (int i=0;i<M.ncol;i++)
     if (M.data[i*M.ncol+i]==0)
     {
     fprintf(stderr,"invert: zero elements in diagonal\n");
     mematrix<DT> temp = M;
     for (int i = 0; i < M.ncol; i++)
     for (int j = 0; j < M.ncol; j++)
     temp.put(NAN,i,j);
     return temp;
     //exit(1);
     }
     */
    int actualsize = M.ncol;
    int maxsize = M.ncol;
    mematrix<DT> temp = M;
    for (int i = 1; i < actualsize; i++)
        temp.data[i] /= temp.data[0]; // normalize row 0
    for (int i = 1; i < actualsize; i++)
    {
        for (int j = i; j < actualsize; j++)
        { // do a column of L
            DT sum = 0.0;
            for (int k = 0; k < i; k++)
                sum += temp.data[j * maxsize + k] * temp.data[k * maxsize + i];
            temp.data[j * maxsize + i] -= sum;
        }
        if (i == actualsize - 1)
            continue;
        for (int j = i + 1; j < actualsize; j++)
        { // do a row of U
            DT sum = 0.0;
            for (int k = 0; k < i; k++)
                sum += temp.data[i * maxsize + k] * temp.data[k * maxsize + j];
            temp.data[i * maxsize + j] = (temp.data[i * maxsize + j] - sum)
                    / temp.data[i * maxsize + i];
        }
    }
    for (int i = 0; i < actualsize; i++) // invert L
        for (int j = i; j < actualsize; j++)
        {
            DT x = 1.0;
            if (i != j)
            {
                x = 0.0;
                for (int k = i; k < j; k++)
                    x -= temp.data[j * maxsize + k]
                            * temp.data[k * maxsize + i];
            }
            temp.data[j * maxsize + i] = x / temp.data[j * maxsize + j];
        }
    for (int i = 0; i < actualsize; i++) // invert U
        for (int j = i; j < actualsize; j++)
        {
            if (i == j)
                continue;
            DT sum = 0.0;
            for (int k = i; k < j; k++)
                sum += temp.data[k * maxsize + j]
                        * ((i == k) ? 1.0 : temp.data[i * maxsize + k]);
            temp.data[i * maxsize + j] = -sum;
        }
    for (int i = 0; i < actualsize; i++) // final inversion
        for (int j = 0; j < actualsize; j++)
        {
            DT sum = 0.0;
            for (int k = ((i > j) ? i : j); k < actualsize; k++)
                sum += ((j == k) ? 1.0 : temp.data[j * maxsize + k])
                        * temp.data[k * maxsize + i];
            temp.data[j * maxsize + i] = sum;
        }
    return temp;
}

//_________Maksim____________
template<class DT>
DT var(mematrix<DT> &M)
{
    DT sum = 0;
    for (int i = 0; i < M.nelements; i++)
    {
        sum += M.data[i];
    }
    DT mean = sum / M.nelements;

    DT sum2 = 0;
    for (int i = 0; i < M.nelements; i++)
    {
        sum2 += pow(M.data[i] - mean, 2);
    }

    return sum2 / (M.nelements - 1);
}
//_________Maksim____________
#endif /* MEMATRI1_H */
