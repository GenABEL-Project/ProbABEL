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

    void base_score(mematrix<double>& resid, regdata& rdata, int verbose,
            double tol_chol, int model, int interaction, int ngpreds,
            const masked_matrix& invvarmatrix, int nullmodel);
};

class linear_reg: public base_reg {
 public:
    linear_reg(regdata& rdatain);
    ~linear_reg()
    {
        //		delete beta;
        //		delete sebeta;
        //		delete residuals;
    }

    void estimate(regdata& rdatain, int verbose, double tol_chol, int model,
                  int interaction, int ngpreds,
                  masked_matrix& invvarmatrixin,
                  int robust, int nullmodel = 0);

    void score(mematrix<double>& resid, regdata& rdatain, int verbose,
               double tol_chol, int model, int interaction, int ngpreds,
               const masked_matrix& invvarmatrix, int nullmodel = 0);
};

class logistic_reg: public base_reg {
 public:
    int niter;

    logistic_reg(regdata& rdatain);
    ~logistic_reg()
    {
        //		delete beta;
        //		delete sebeta;
    }

    void estimate(regdata& rdatain, int verbose, int maxiter, double eps,
                  double tol_chol, int model, int interaction, int ngpreds,
                  masked_matrix& invvarmatrixin, int robust,
                  int nullmodel = 0);
    // just a stupid copy from linear_reg
    void score(mematrix<double>& resid, regdata& rdata, int verbose,
               double tol_chol, int model, int interaction, int ngpreds,
               masked_matrix& invvarmatrix, int nullmodel = 0);
};

#endif
