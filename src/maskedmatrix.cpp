/*
 * maskedmatrix.cpp
 *
 *  Created on: May 22, 2012
 *      Author: mkooyman
 *
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


#include "maskedmatrix.h"
#if EIGEN
#include "eigen_mematrix.h"
#include "eigen_mematrix.cpp"
#else
#include "mematrix.h"
#include "mematri1.h"
#endif

masked_matrix::masked_matrix()
{
    length_of_mask = 0;
    masked_data = NULL;
    mask_of_old = NULL;
}

masked_matrix::masked_matrix(mematrix<double> M) : matrix_original(M)
{
//    matrix_original = M;
    masked_data = &matrix_original;
    mask_of_old = new unsigned short int[M.nrow];
    for (int i = 0; i < M.nrow; i++)
    {
        mask_of_old[i] = 0;
    }
    //TODO:set length of mask for all types
    length_of_mask = M.nrow;
    //TODO:set type (row,column,symmetric)
    //type="symmetric";
}

void masked_matrix::set_matrix(const mematrix<double> &M)
{
    matrix_original = M;
    masked_data = &matrix_original;
    mask_of_old = new unsigned short int[M.nrow];
    for (int i = 0; i < M.nrow; i++)
    {
        mask_of_old[i] = 0;
    }
    //TODO:set length of mask for all types
    length_of_mask = M.nrow;
    //TODO:set type (row,column,symmetric)
    //type="symmetric";
}

masked_matrix::~masked_matrix()
{
    delete[] mask_of_old;
}

void masked_matrix::update_mask(short unsigned int *newmask)
{
    //find length of masked matrix
    int nmeasured = 0;
    for (int i = 0; i < length_of_mask; i++)
    {
        if (newmask[i] == 0)
        {
            nmeasured++;
        }
    }

    //Check update mask is the same as original matrix
    if (nmeasured == length_of_mask)
    {
        //masked matrix is the same as original matrix
        masked_data = &matrix_original;
    }
    else
    {
        //Check update mask is the same as old matrix
        if (is_equal_array(newmask, mask_of_old, length_of_mask))
        {
            //new mask is the same as old matrix
            masked_data = &matrix_masked_data;
        }
        else
        {
            // new mask differs from old matrix and create new.
            // mask_of_old = newmask;
            //TODO(maarten): there must be a smarter way to copy these values
            for (int i = 0; i < length_of_mask; i++)
            {
                mask_of_old[i] = newmask[i];
            }
            mask_symmetric(nmeasured);
            masked_data = &matrix_masked_data;
        }
    }
}

void masked_matrix::mask_symmetric(int nmeasured)
{
    // Mask a symmetric matrix: this matrix is always a square matrix and will
    // mask rows and columns. The result is always a square matrix
    matrix_masked_data.reinit(nmeasured, nmeasured);
    int i1 = 0, j1 = 0;

    for (int i = 0; i < length_of_mask; i++)
        if (mask_of_old[i] == 0)
        {
            j1 = 0;
            for (int j = 0; j < length_of_mask; j++)
                if (mask_of_old[j] == 0)
                {
                    //std::cout << "val" << i1 << " " << j1 << "\n";
                    matrix_masked_data.put(matrix_original.get(i, j), i1, j1);
                    j1++;
                }
            i1++;
        }
}

bool masked_matrix::is_equal_array(unsigned short int *a, unsigned short int *b,
                                   int size)
{
    for (int i = 0; i < size; i++)
    {
        if (a[i] != b[i])
        {
            return false;
        }
    }
    return true;
}
