//==============================================================================
//
//           Filename:  src/reg1.h
//
//        Description:  ProbABEL
//
//            Version:  0.1-3
//            Created:  ---
//           Revision:  none
//  last modification:  11-Jan-2009
//
//             Author:  Yurii S. Aulchenko
//                        modified by:  Maksim V. Struchalin, 11-Jan-2009
//
// Modified by Han Chen (hanchen@bu.edu) on Nov 9, 2009
// based on src/reg1.h version 0.2 as of Oct 19, 2009
//
//            Company:  ErasmusMC,
//                      Epidemiology & Biostatistics Department,
//                      Rotterdam,
//                      The Netherlands.
//              Email:  i.aoultchenko@erasmusmc.nl, m.struchalin@erasmusmc.nl
//
//==============================================================================
//#include "coxph_data.h"
#ifndef REG1_H_
#define REG1_H_
#include <cmath>
#include "cholesky.h"
#include "regdata.h"
#include "maskedmatrix.h"


mematrix<double> apply_model(mematrix<double>& X, int model, int interaction,
        int ngpreds, bool is_interaction_excluded, bool iscox = false,
        int nullmodel = 0);

mematrix<double> t_apply_model(mematrix<double>& X, int model, int interaction,
        int ngpreds, bool iscox, int nullmodel = 0);

class base_reg {
 public:
    mematrix<double> beta;
    mematrix<double> sebeta;
    //Han Chen
    mematrix<double> covariance;
    //Oct 26, 2009
    mematrix<double> residuals;
    double sigma2;
    double loglik;
    double chi2_score;
    regdata reg_data;

    void base_score(mematrix<double>& resid,  int verbose,
            double tol_chol, int model, int interaction, int ngpreds,
            const masked_matrix& invvarmatrix, int nullmodel);
};

class linear_reg: public base_reg {
 public:
    linear_reg(regdata& rdatain);
    ~linear_reg()
    {
        delete [] reg_data.masked_data ;
        //		delete beta;
        //		delete sebeta;
        //		delete residuals;
    }

    void estimate( int verbose, double tol_chol, int model,
                  int interaction, int ngpreds,
                  masked_matrix& invvarmatrixin,
                  int robust, int nullmodel = 0);

    void score(mematrix<double>& resid,  int verbose,
               double tol_chol, int model, int interaction, int ngpreds,
               const masked_matrix& invvarmatrix, int nullmodel = 0);
};

class logistic_reg: public base_reg {
 public:
    int niter;

    logistic_reg(regdata& rdatain);
    ~logistic_reg()
    {
        delete [] reg_data.masked_data ;
        //		delete beta;
        //		delete sebeta;
    }

    void estimate( int verbose, int maxiter, double eps,
                  double tol_chol, int model, int interaction, int ngpreds,
                  masked_matrix& invvarmatrixin, int robust,
                  int nullmodel = 0);
    // just a stupid copy from linear_reg
    void score(mematrix<double>& resid,  int verbose,
               double tol_chol, int model, int interaction, int ngpreds,
               masked_matrix& invvarmatrix, int nullmodel = 0);
};

#endif
