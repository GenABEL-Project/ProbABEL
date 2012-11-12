/*
 * cholesky.h
 *
 *  Created on: Mar 15, 2012
 *      Author: mkooyman
 */

#ifndef CHOLESKY_H_
#define CHOLESKY_H_

#if EIGEN
#include "eigen_mematrix.h"
#include "eigen_mematrix.cpp"
#else
#include "mematrix.h"
#endif

int cholesky2_mm(mematrix<double> &matrix, double toler);
void chinv2_mm(mematrix<double> &matrix);

#endif /* CHOLESKY_H_ */
