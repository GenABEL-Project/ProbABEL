/*
 * coxph_data.h
 *
 *  Created on: Mar 31, 2012
 *      Author: mkooyman
 */

#ifndef COXPH_DATA_H_
#define COXPH_DATA_H_

#if EIGEN
#include "eigen_mematrix.h"
#include "eigen_mematrix.cpp"
#else
#include "mematrix.h"
#include "mematri1.h"
#endif

#include "reg1.h"
#include "gendata.h"
#include "phedata.h"

class coxph_data {
public:
    coxph_data get_unmasked_data();

    coxph_data()
    {
        nids = 0;
        ncov = 0;
        ngpreds = 0;
        masked_data = NULL;
    }

    coxph_data(const coxph_data &obj);
    coxph_data(phedata &phed, gendata &gend, int snpnum);
    void update_snp(gendata &gend, int snpnum);
    ~coxph_data();

    int nids;
    int ncov;
    int ngpreds;
    mematrix<double> weights;
    mematrix<double> stime;
    mematrix<int>    sstat;
    mematrix<double> offset;
    mematrix<int>    strata;
    mematrix<double> X;
    mematrix<int>    order;
    unsigned short int * masked_data;
};

class coxph_reg {
public:
    mematrix<double> beta;
    mematrix<double> sebeta;
    mematrix<double> residuals;
    double sigma2;
    double loglik;
    double chi2_score;
    int niter;

    coxph_reg(coxph_data &cdatain);
    void estimate(coxph_data &cdatain, int verbose, int maxiter, double eps,
                  double tol_chol, int model, int interaction, int ngpreds,
                  bool iscox, int nullmodel);
};

#endif /* COXPH_DATA_H_ */
