/*
 * cholesky.cpp
 *
 *  Created on: Mar 15, 2012
 *      Author: mkooyman
 */
#include <string>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>

#if EIGEN
#include "eigen_mematrix.h"
#include "eigen_mematrix.cpp"
#else
#include "mematrix.h"
#include "mematri1.h"
#endif


/*  SCCS @(#)cholesky2.c    5.2 10/27/98
 ** subroutine to do Cholesky decompostion on a matrix: C = FDF'
 **   where F is lower triangular with 1's on the diagonal, and D is diagonal
 **
 ** arguments are:
 **     n         the size of the matrix to be factored
 **     **matrix  a ragged array containing an n by n submatrix to be factored
 **     toler     the threshold value for detecting "singularity"
 **
 **  The factorization is returned in the lower triangle, D occupies the
 **    diagonal and the upper triangle is left undisturbed.
 **    The lower triangle need not be filled in at the start.
 **
 **  Return value:  the rank of the matrix (non-negative definite), or -rank
 **     it not SPD or NND
 **
 **  If a column is deemed to be redundant, then that diagonal is set to zero.
 **
 **   Terry Therneau
 */

// modified Yurii Aulchenko 2008-05-20
int cholesky2_mm(mematrix<double> &matrix, double toler)
{
    if (matrix.ncol != matrix.nrow)
    {
        fprintf(stderr, "cholesky2_mm: matrix should be square\n");
        exit(1);
    }
    int n = matrix.ncol;
    double temp;
    int i, j, k;
    double eps, pivot;
    int rank;
    int nonneg;

    nonneg = 1;
    eps = 0;
    for (i = 0; i < n; i++)
    {
        if (matrix[i * n + i] > eps)
            eps = matrix[i * n + i];
        for (j = (i + 1); j < n; j++)
            matrix[j * n + i] = matrix[i * n + j];
    }
    eps *= toler;

    rank = 0;
    for (i = 0; i < n; i++)
    {
        pivot = matrix[i * n + i];
        if (pivot < eps)
        {
            matrix[i * n + i] = 0;
            if (pivot < -8 * eps)
                nonneg = -1;
        }
        else
        {
            rank++;
            for (j = (i + 1); j < n; j++)
            {
                temp = matrix[j * n + i] / pivot;
                matrix[j * n + i] = temp;
                matrix[j * n + j] -= temp * temp * pivot;
                for (k = (j + 1); k < n; k++)
                    matrix[k * n + j] -= temp * matrix[k * n + i];
            }
        }
    }
    return (rank * nonneg);
}

void chinv2_mm(mematrix<double> &matrix)
{
    if (matrix.ncol != matrix.nrow)
    {
        fprintf(stderr, "cholesky2_mm: matrix should be square\n");
        exit(1);
    }

    int n = matrix.ncol;
    register double temp;
    register int i, j, k;

    /*
     ** invert the cholesky in the lower triangle
     **   take full advantage of the cholesky's diagonal of 1's
     */
    for (i = 0; i < n; i++)
    {
        if (matrix[i * n + i] > 0)
        {
            matrix[i * n + i] = 1 / matrix[i * n + i]; /*this line inverts D */
            for (j = (i + 1); j < n; j++)
            {
                matrix[j * n + i] = -matrix[j * n + i];
                for (k = 0; k < i; k++) /*sweep operator */
                    matrix[j * n + k] += matrix[j * n + i] * matrix[i * n + k];
            }
        }
    }

    /*
     ** lower triangle now contains inverse of cholesky
     ** calculate F'DF (inverse of cholesky decomp process) to get inverse
     **   of original matrix
     */
    for (i = 0; i < n; i++)
    {
        if (matrix[i * n + i] == 0)
        { /* singular row */
            for (j = 0; j < i; j++)
                matrix[j * n + i] = 0;
            for (j = i; j < n; j++)
                matrix[i * n + j] = 0;
        }
        else
        {
            for (j = (i + 1); j < n; j++)
            {
                temp = matrix[j * n + i] * matrix[j * n + j];
                if (j != i)
                    matrix[i * n + j] = temp;
                for (k = i; k < j; k++)
                    matrix[i * n + k] += temp * matrix[j * n + k];
            }
        }
    }
    // ugly fix to return only inverse
    for (int col = 1; col < n; col++)
        for (int row = 0; row < col; row++)
            matrix[col * n + row] = matrix[row * n + col];
}
