/*
 * coxph_reg.h
 *
 *  Created on: Oct 6, 2012
 *      Author: mkooyman
 */

#ifndef COXPH_REG_H_
#define COXPH_REG_H_
#include "coxph_data.h"

#if EIGEN
#include "eigen_mematrix.h"
#else
#include "mematrix.h"
#endif

class coxph_reg {
public:
	mematrix<double> beta;
	mematrix<double> sebeta;
	mematrix<double> residuals;
	double sigma2;
	double loglik;
	double chi2_score;
	int niter;

	coxph_reg(const coxph_data &cdatain);
	virtual ~coxph_reg();
	void coxph_reg::estimate( coxph_data &cdatain, int verbose, int maxiter,
			double eps, double tol_chol, int model, int interaction, int ngpreds,
			bool iscox, int nullmodel = 0);
};

#endif /* COXPH_REG_H_ */
