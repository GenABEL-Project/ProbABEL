/*
 * maskedmatrix.h
 *
 *  Created on: May 22, 2012
 *      Author: mkooyman
 *
 *
 * Copyright (C) 2009--2016 Various members of the GenABEL team. See
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

#include "eigen_mematrix.h"
#include "eigen_mematrix.cpp"


class masked_matrix {
 public:
    masked_matrix();
    masked_matrix(mematrix<double> M);
    // ~masked_matrix();


    mematrix<double> matrix_original;
    mematrix<double> *masked_data;
    int length_of_mask;


    void set_matrix(const mematrix<double> &M);
    void update_mask(std::vector<bool> newmask);


 private:
    mematrix<double> matrix_masked_data;
    std::vector<bool> mask_of_old;

    void mask_symmetric(int nmeasured);
};

#endif // MASKEDMATRIX_H_
