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


//#include "coxph_data.h"
#ifndef REG1_H_
#define REG1_H_
#include <cmath>
#include "regdata.h"
#include "maskedmatrix.h"


mematrix<double> apply_model(const mematrix<double>& X,
                             const int model,
                             const int interaction,
                             const int ngpreds,
                             const bool is_interaction_excluded,
                             const bool iscox = false,
                             const int nullmodel = 0);

mematrix<double> t_apply_model(const mematrix<double>& X,
                               const int model,
                               const int interaction,
                               const int ngpreds,
                               const bool iscox,
                               const int nullmodel = 0);


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

    void base_score(const mematrix<double>& resid,
                    const int model,
                    const int interaction, const int ngpreds,
                    const masked_matrix& invvarmatrix,
                    int nullmodel);
};


class linear_reg: public base_reg {
 public:
    linear_reg(const regdata& rdatain);

    void estimate(const int verbose, const int model,
                  const int interaction, const int ngpreds,
                  masked_matrix& invvarmatrixin,
                  const int robust, const int nullmodel = 0);

    void score(const mematrix<double>& resid,
               const int model, const int interaction,
               const int ngpreds, const masked_matrix& invvarmatrix,
               int nullmodel = 0);

 private:
    void mmscore_regression(const mematrix<double>& X,
                            const masked_matrix& W_masked,
                            LDLT<MatrixXd>& Ch);
    void logLikelihood(const mematrix<double>& X);
    void LeastSquaredRegression(const mematrix<double> & X,
                                LDLT<MatrixXd>& Ch);
    void RobustSEandCovariance(const mematrix<double> & X,
                               mematrix <double> robust_sigma2,
                               const MatrixXd tXX_inv,
                               const int offset);
    void PlainSEandCovariance(const double sigma2_internal,
                              const MatrixXd & tXX_inv,
                              const int offset);
};


class logistic_reg: public base_reg {
 public:
    int niter;

    logistic_reg(const regdata& rdatain);

    void estimate(const int verbose,
                  const int model, const int interaction, const int ngpreds,
                  masked_matrix& invvarmatrixin, const int robust,
                  const int nullmodel = 0);

    // just a stupid copy from linear_reg
    void score(const mematrix<double>& resid,
               const int model, const int interaction,
               const int ngpreds, const masked_matrix& invvarmatrix,
               int nullmodel = 0);

 private:
    /**
     * \brief Constant that contains the maximum number of iterations
     * done during the regression routine.
     */
    static const int MAXITER = 20;

    /**
     * \brief Constant containing the tolerance for convergence.
     *
     * Iteration continues until the percent change in loglikelihood
     * is <= EPS.
     */
    static const double EPS = 1e-8;
};

#endif // REG1_H_
