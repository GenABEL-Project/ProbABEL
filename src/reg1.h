/**
 * \file reg1.h
 * \author Yurii S. Aulchenko
 * \author M. Kooyman
 * \author L.C. Karssen
 * \author Maksim V. Struchalin
 * \author Han Chen (hanchen@bu.edu)
 *
 * \brief Describes various classes containing regression objects.
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

/**
 * \brief The base_reg class defines a basic regression object.
 *
 * It contains all elements like a design matrix \f$X\f$, the vector
 * of regression coefficients \f$\beta\f$, that are shared by more
 * specialised classes like linear_reg and logistic_reg.
 */
class base_reg {
 public:
    /**
     * \brief The vector containing the regression coefficients.
     *
     * Note that the SNP coefficients are the last elements of the
     * vector. The mean and coefficients for other covariates come
     * first.
     *
     */
    mematrix<double> beta;

    /**
     * \brief The vector containing the standard errors of the
     * regression coefficients.
     *
     */
    mematrix<double> sebeta;

    //Han Chen
    mematrix<double> covariance;
    //Oct 26, 2009
    mematrix<double> residuals;

    /**
     * \brief The MLE of the residual variance.
     *
     * See the equation following Eq.(1) in the ProbABEL paper
     * (Aulchenko et al. 2010).
     *
     */
    double sigma2;

    /**
     * \brief The loglikelihood of the model
     *
     */
    double loglik;

    double chi2_score;
    regdata reg_data;

    void base_score(const mematrix<double>& resid,
                    const int model,
                    const int interaction, const int ngpreds,
                    const masked_matrix& invvarmatrix,
                    int nullmodel);
};


/**
 * \brief An extension of the  base_reg class specialised in linear
 * regression.
 *
 * This class contains functions for estimation of the model
 * parameters, loglikelihood, etc. using linear regression.
 */
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


/**
 * \brief An extension of the  base_reg class specialised in logistic
 * regression.
 *
 * This class contains functions for estimation of the model
 * parameters, loglikelihood, etc. using logistic regression.
 */
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
