/*
 * maskedmatrix.h
 *
 *  Created on: May 22, 2012
 *      Author: mkooyman
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
