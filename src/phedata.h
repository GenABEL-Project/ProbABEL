/**
 * \file phedata.h
 * \author Y.S. Aulchenko
 * \author L.C. Karssen
 * \author M. Kooyman
 * \author Maksim V. Struchalin
 *
 * \brief Describes the phdata class containing the phenotype data.
 *
 *  Created on: Mar 6, 2012
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


#ifndef PHEDATA_H_
#define PHEDATA_H_

#include "eigen_mematrix.h"
#include "eigen_mematrix.cpp"


/**
 * \brief A phedata object contains the phenotype data used in the
 * regression analysis.
 *
 * Typical data contained in this object are the values for outcome(s)
 * and covariates.
 */
class phedata {
 public:
    phedata()
    {
        is_interaction_excluded = 0;
        nids_all                = 0;
        nids                    = 0;
        noutcomes               = 0;
        ncov                    = 0;
        n_model_terms           = 0;
        allmeasured             = NULL;
        idnames                 = NULL;
        model_terms             = NULL;
    }

    phedata(const char * fname, const int noutc, const int npeople,
            const int interaction, const bool iscox);

    ~phedata();


    void setphedata(const char * fname, const int noutc, const int npeople,
                    const int interaction, const bool iscox);
    void set_is_interaction_excluded(const bool inter_excluded);


    bool is_interaction_excluded;

    /**
     * \brief The total number of individuals in the phenotype file,
     * or, if the <tt>-n/\--ids</tt> command line option is given, the
     * number specified with that option.
     *
     * After reading the phenotype data, this variable contains the
     * total number of individuals found in the phenotype file (or
     * maximised by the number specified by the <tt>-n/\--nids</tt>
     * option). This will be the maximum number of individuals used in
     * the analysis (some may still be dropped if their phenotype was
     * not set).
     *
     * \sa cmdvars::getNpeople()
     */
    int nids_all;

    /**
     * \brief The number of individuals for which phenotype data is
     * complete.
     *
     * This is the number of people (<tt>#nids <= #nids_all</tt>) used
     * in the analysis. Any people with incomplete outcome or
     * covariate data are excluded.
     */
    int nids;

    /**
     * \brief The number of outcomes in the phenotype file.
     *
     * Usually this is 1, e.g. height (palinear) or affection status
     * (palogist). For Cox PH regression (pacoxph) this is 2: the
     * followup time and the survival status.
     */
    int noutcomes;

    /**
     * \brief The number of covariates in the phenotype file.
     *
     * This variable contains the total number of covariates as found
     * in the phenotype file (so excluding any genotypes or
     * interactions).
     */
    int ncov;

    /**
     * \brief The number of terms in the statistical model used for
     * the regression.
     *
     * For example, for a regression formula <tt>Y ~ mu + sex +
     * age</tt> #n_model_terms is equal to 3. If interaction terms
     * are requested, these are added as well.
     */
    int n_model_terms;

    /**
     * \brief An array of length #nids_all indicating if the
     * phenotype and covariate data for the given individual is
     * present or not.
     *
     * An element in the array is set to 1 if none of the columns in
     * the phenotype file for that individual has starts with an 'N'
     * or 'n'. This detects NA's in the phenotype/covariate data. If
     * one of the values is NA, the #allmeasured row corresponding
     * to this individual is set to 0.
     */
    unsigned short int * allmeasured;

    /**
     * \brief An array of strings containing the individuals IDs for
     * those individuals that are included in the analysis.
     *
     * Length is #nids, so people with incomplete phenotype
     * information (see <tt>allmeasured</tt>) are not in this array.
     */
    std::string * idnames;

    /**
     * \brief The statistical model to be used for the regression.
     *
     * Depending on the covariates in the phenotype file and the
     * interactions specified, this string contains a textual
     * representation of the (full) statistical model that will be
     * analysed. For example <tt>Y ~ mu + sex + age + SNP_A1</tt>.
     */
    std::string model;

    /**
     * \brief A list of strings containing each of the model terms
     * used in the textual representation of the statistical model.
     *
     */
    std::string * model_terms;

    /**
     * \brief The matrix with covariate data.
     *
     * After reading the phenotype data, this matrix contains the
     * values of the covariates for each individual with complete
     * phenotype data. The number of rows is equal to #nids, the
     * number of columns is equal to #ncov.
     */
    mematrix<double> X;

    /**
     * \brief The vector/matrix contain the phenotype values.
     *
     * The number of columns is equal to the number of outcomes
     * (#noutcomes). The number of rows is equal to the number of
     * individuals with complete phenotype data.
     */
    mematrix<double> Y;
};

#endif /* PHEDATA_H_ */
