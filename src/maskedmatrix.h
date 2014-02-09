/*
 * maskedmatrix.h
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


#ifndef MASKEDMATRIX_H_
#define MASKEDMATRIX_H_

#if EIGEN
#include "eigen_mematrix.h"
#include "eigen_mematrix.cpp"
#else
#include "mematrix.h"
#include "mematri1.h"
#endif

class masked_matrix {
 public:
    masked_matrix();
    masked_matrix(mematrix<double> M);
    void set_matrix(const mematrix<double> &M);
    ~masked_matrix();
    void update_mask(short unsigned int *newmask);
//    mematrix<double>* get_matrix();
    mematrix<double> matrix_original;
    mematrix<double> *masked_data;
    int length_of_mask;

 private:
    mematrix<double> matrix_masked_data;
    unsigned short int *mask_of_old;
    void mask_symmetric(int nmeasured);
    bool is_equal_array(unsigned short int *a, unsigned short int *b, int size);
};

#endif /* MASKEDMATRIX_H_ */
