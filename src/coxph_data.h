/**
 * \file coxph_data.h
 * \author Yurii S. Aulchenko
 * \author M. Kooyman
 * \author L.C. Karssen
 * \author Maksim V. Struchalin
 *
 * \brief Describes the coxph_data and coxph_reg classes used in for
 * Cox Proportional Hazards regression. The coxph_data class is
 * similar to the regdata class.
 *
 *  Created on: Mar 31, 2012
 *
 *
 * Copyright (C) 2009--2015 Various members of the GenABEL team. See
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


#ifndef COXPH_DATA_H_
#define COXPH_DATA_H_

#include <vector>

#include "eigen_mematrix.h"
#include "gendata.h"
#include "mlinfo.h"
#include "phedata.h"
#include "reg1.h"

/**
 * \brief A coxph_data object contains the data used for Cox PH
 * regression.
 *
 * The Cox regression code in coxfit2.c requires all data to be sorted
 * based on follow-up time, which is done in the constructor.
 *
 * A similar class for linear and logistic regression can be found in
 * the regdata class.
 */
class coxph_data {
 public:
    int nids;                   /**< \brief Number of IDs/samples.  */

    /**
     * \brief Number of covariates + possible nr of genomic
     * predictors.
     *
     * If snpnum >=0 then this equals the number of covariates + the
     * number of regdata::ngpreds.
     */
    int ncov;

   /**
     * \brief Number of genomic predictors, 1 for dosage data, 2 for
     * probability data.
     */
    int ngpreds;

    /**
     * \brief Number of non-masked genotypes.
     */
    unsigned int gcount;

    /**
     * \brief Allele frequency.
     *
     * Calculation is only based on non-masked SNPs.
     */
    double freq;

    /**
     * \brief The design matrix.
     *
     * The matrix dimensions are coxph_data::nids x coxph_data::ncov.
     * After filling the matrix is it is time-sorted based on
     * coxph_data::order. Note that coxfit2() expects the data to be
     * in column-major order, so at the end of the constructor the
     * matrix is transposed!
     */
    mematrix<double> X;

    /**
     * \brief A vector that contains the follow-up times for each
     * individual, as read from the phenotype file.
     *
     * The vector is coxph_data::nids long. After reading the data
     * from the phenotype file this vector is time-sorted based on
     * coxph_data::order.
     */
    mematrix<double> stime;

    /**
     * \brief A vector containing the survival status of each
     * individual, as read from the phenotype file.
     *
     * Values should be 0 or 1. The vector is coxph_data::nids
     * long. After reading the data from the phenotype file this
     * vector is time-sorted based on coxph_data::order.
     */
    mematrix<int> sstat;

    /**
     * \brief A vector containing case weights. These are set to 1.0
     * in the constructor.
     *
     * The vector is coxph_data::nids long. After reading the data
     * from the phenotype file this vector is time-sorted based on
     * coxph_data::order.
     */
    mematrix<double> weights;

    /**
     * \brief A vector containing an offset for the linear
     * predictor. Set to 0.0 in the constructor. This is a mandatory
     * input for coxfit2().
     *
     * The vector is coxph_data::nids long. After reading the data
     * from the phenotype file this vector is time-sorted based on
     * coxph_data::order.
     */
    mematrix<double> offset;

    /**
     * \brief A vector containing strata. See coxfit2() for more
     * details. All strata are set to 0 in the constructor.
     *
     * The vector is coxph_data::nids long. After reading the data
     * from the phenotype file this vector is time-sorted based on
     * coxph_data::order.
     */
    mematrix<int> strata;

    /**
     * \brief A vector used to store the follow-up times in ascending
     * order. It is used to reorder the other vectors and
     * coxph_data::X.
     *
     * The vector is coxph_data::nids long.
     */
    mematrix<int> order;

    /**
     * \brief A vector that contains ones/zeros to indicate which data points
     * should be omitted because no genetic data is present.
     *
     * The vector is regdata::nids long. A value of 1 or 'true' means
     * that that ID/sample will be masked because the SNP data is NA
     * for that ID.
     */
    std::vector<bool> masked_data;


    // Constructors and destructors
    coxph_data();
    coxph_data(const coxph_data &obj);
    coxph_data(const phedata &phed, const gendata &gend, const int snpnum);
    // ~coxph_data();


    // Member functions
    coxph_data get_unmasked_data() const;
    void update_snp(const gendata *gend,
                    const int snpnum,
                    mlinfo &snpinfo,
                    const bool flipMAF);
    void remove_snp_from_X();
};


class coxph_reg {
 public:
    // Variables
    mematrix<double> beta;
    mematrix<double> sebeta;
    mematrix<double> residuals;
    double sigma2;
    double loglik;
    double chi2_score;
    int niter;


    // Constructors and destructors
    coxph_reg(const coxph_data &cdatain);


    // Member functions
    void estimate(const coxph_data &cdatain, const int model,
                  const std::vector<std::string> &modelNames,
                  const int interaction, const int ngpreds, const bool iscox,
                  const int nullmodel, const mlinfo &snpinfo, const int cursnp);

 private:
    /**
     * \brief Constant that contains the maximum number of iterations
     * done during the regression routine.
     */
    static const int MAXITER = 20;

    /**
     * \brief Constant containing the tolerance for
     * convergence.
     *
     * Iteration continues until the percent change in loglikelihood
     * is <= EPS.
     */
    static const double EPS = 1e-8;

    /**
     * \brief Constant containing the precision for the Cholesky
     * decomposition.
     */
    static const double CHOLTOL = 1.5e-12;
};

#endif /* COXPH_DATA_H_ */
