/**
 * \file regdata.h
 * \author Yurii S. Aulchenko
 * \author M. Kooyman
 * \author L.C. Karssen
 * \author Maksim V. Struchalin
 *
 * \brief Describes the regdata class containing a linear or logistic
 * regression object.
 *
 *  Created on: Mar 29, 2012
 *      Author: mkooyman
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


#ifndef REGDATA_H_
#define REGDATA_H_

#include "eigen_mematrix.h"
#include "eigen_mematrix.cpp"
#include "gendata.h"
#include "phedata.h"
#include "mlinfo.h"


/**
 * \brief A regdata object contains data used for linear or
 * logistic regression.
 *
 * A similar object for CoxPH regression can be found in the
 * coxph_data class.
 *
 */
class regdata {
 public:
    int nids;                   /**< Number of IDs/samples. */

    /**
     * \brief Number of covariates + possible nr of genomic predictors.
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
      * \brief Number of outcomes, taken from phedata::noutcomes.
      */
    int noutcomes;

    /**
     * \brief A vector that contains ones/zeros to indicate which data
     * points should be omitted because no genetic data is present.
     *
     * The vector is regdata::nids long. A value of 1 or 'true' means
     * that that ID/sample will be masked because the SNP data is NA
     * for that ID.
     */
    std::vector<bool> masked_data;

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
     * The matrix dimensions are regdata::nids x (regdata::ncov + 1)
     */
    mematrix<double> X;

    /**
     * \brief Matrix containing the phenotype data.
     *
     * The matrix dimensions are regdata::nids x
     * regdata::noutcomes.
     */
    mematrix<double> Y;


    // Constructors and destructors
    regdata();
    regdata(const regdata &obj);
    regdata(const phedata &phed, const gendata &gend, const int snpnum);
    //~regdata();


    // Member functions.
    void update_snp(const gendata *gend,
                    const int snpnum,
                    const mlinfo &snpinfo,
                    const bool flipMAF);
    void remove_snp_from_X();
    regdata get_unmasked_data() const;
    mematrix<double> extract_genotypes();

 private:
};

#endif /* REGDATA_H_ */
