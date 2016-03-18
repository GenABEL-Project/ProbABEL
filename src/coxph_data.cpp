/**
 * \file coxph_data.cpp
 * \author Yurii S. Aulchenko
 * \author L.C. Karssen
 * \author M. Kooyman
 * \author Maksim V. Struchalin
 *
 * \brief Describes functions of the coxphdata class for Cox PH
 * regression objects
 *
 *  Created on: Mar 31, 2012
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


#include "coxph_data.h"
#include <iostream>
#include <algorithm>
#include <cmath>
extern "C" {
#include "survproto.h"
}

#include "fvlib/AbstractMatrix.h"
#include "fvlib/CastUtils.h"
#include "fvlib/const.h"
#include "fvlib/convert_util.h"
#include "fvlib/FileVector.h"
#include "fvlib/frutil.h"
#include "fvlib/frversion.h"
#include "fvlib/Logger.h"
#include "fvlib/Transposer.h"

/**
 * Definition of the compare function used for the sort of times. Used
 * by the qsort() function.
 */
int cmpfun(const void *a, const void *b)
{
    double el1 = *(double*) a;
    double el2 = *(double*) b;
    if (el1 > el2)
    {
        return 1;
    }
    if (el1 < el2)
    {
        return -1;
    }
    if (el1 == el2)
    {
        return 0;
    }

    // You should never come here...
    return -9;
}


/**
 * Constructor. Initialises all values to zero.
 */
coxph_data::coxph_data()
{
    nids    = 0;
    ncov    = 0;
    ngpreds = 0;
    gcount  = 0;
    freq    = 0;
}


/**
 * Copy constructor. Creates a coxph_data object by copying the values
 * of another one.
 *
 * \param obj Reference to the coxph_data object to be copied to the new
 * object.
 */
coxph_data::coxph_data(const coxph_data &obj) : X(obj.X),
                                                stime(obj.stime),
                                                sstat(obj.sstat),
                                                weights(obj.weights),
                                                offset(obj.offset),
                                                strata(obj.strata),
                                                order(obj.order),
                                                masked_data(obj.masked_data)
{
    nids        = obj.nids;
    ncov        = obj.ncov;
    ngpreds     = obj.ngpreds;
    gcount      = 0;
    freq        = 0;
}


/**
 * Constructor that fills a coxph_data object with phenotype and genotype
 * data.
 *
 * @param phed Reference to a phedata object with phenotype data
 * @param gend Reference to a gendata object with genotype data
 * @param snpnum The number of the SNP in the genotype data object to
 * be added to the design matrix regdata::X. When set to a number < 0
 * no SNP data is added to the design matrix (e.g. when calculating
 * the null model).
 */
coxph_data::coxph_data(const phedata &phed, const gendata &gend,
                       const int snpnum)
{
    freq        = 0;
    gcount      = 0;
    nids        = gend.nids;
    masked_data = std::vector<bool>(nids, false);

    ngpreds = gend.ngpreds;
    if (snpnum >= 0)
    {
        ncov = phed.ncov + ngpreds;
    }
    else
    {
        ncov = phed.ncov;
    }

    if (phed.noutcomes != 2)
    {
        std::cerr << "coxph_data: number of outcomes should be 2 (now: "
                  << phed.noutcomes << ")\n";
        exit(1);
    }

    X.reinit(nids, (ncov + 1)); // Note: ncov takes ngpreds into
                                // account, see above!
    stime.reinit(nids, 1);
    sstat.reinit(nids, 1);
    weights.reinit(nids, 1);
    offset.reinit(nids, 1);
    strata.reinit(nids, 1);
    order.reinit(nids, 1);

    for (int i = 0; i < nids; i++)
    {
        stime[i] = (phed.Y).get(i, 0);
        sstat[i] = static_cast<int>((phed.Y).get(i, 1));
        if (sstat[i] != 1 && sstat[i] != 0)
        {
            std::cerr << "coxph_data: status not 0/1 "
                      <<"(correct order: id, fuptime, status ...)"
                      << endl;
            exit(1);
        }
    }

    // Add a column with a constant (=1) to the X matrix (the mean)
    for (int i = 0; i < nids; i++)
    {
        X.put(1., i, 0);
    }

    // Insert the covariate data into X (note we use phed.ncov and not
    // ncov, which includes ngpreds is not computing the null model!)
    for (int j = 1; j <= phed.ncov; j++)
    {
        for (int i = 0; i < nids; i++)
        {
            X.put((phed.X).get(i, j - 1), i, j);
        }
    }

    // Insert the genotype data into X
    if (snpnum > 0)
    {
        for (int j = 0; j < ngpreds; j++)
        {
            double *snpdata = new double[nids];
            gend.get_var(snpnum * ngpreds + j, snpdata);
            for (int i = 0; i < nids; i++)
            {
                X.put(snpdata[i], i, (ncov - ngpreds + j));
            }
            delete[] snpdata;
        }
    }

    for (int i = 0; i < nids; i++)
    {
        weights[i] = 1.0;
        offset[i] = 0.0;
        strata[i] = 0;
    }

    // sort by time
    double *tmptime = new double[nids];
    int *passed_sorted = new int[nids];
    std::fill(passed_sorted, passed_sorted + nids, 0);


    for (int i = 0; i < nids; i++)
    {
        tmptime[i] = stime[i];
    }

    qsort(tmptime, nids, sizeof(double), cmpfun);

    for (int i = 0; i < nids; i++)
    {
        int passed = 0;
        for (int j = 0; j < nids; j++)
        {
            if (tmptime[j] == stime[i])
            {
                if (!passed_sorted[j])
                {
                    order[i] = j;
                    passed_sorted[j] = 1;
                    passed = 1;
                    break;
                }
            }
        }
        if (passed != 1)
        {
            std::cerr << "cannot recover element " << i << "\n";
            exit(1);
        }
    }

    stime   = reorder(stime, order);
    sstat   = reorder(sstat, order);
    weights = reorder(weights, order);
    strata  = reorder(strata, order);
    offset  = reorder(offset, order);
    X       = reorder(X, order);

    // The coxfit2() function expects data in column major order.
    X = transpose(X);

    // X.print();
    // offset.print();
    // weights.print();
    // stime.print();
    // sstat.print();

    delete[] tmptime;
    delete[] passed_sorted;
}


/**
 * \brief Update the SNP dosages/probabilities in the design matrix
 * coxph_data::X.
 *
 * Adds the genetic information for a new SNP to the design
 * matrix.
 *
 * \param gend Object that contains the genetic data from which the
 * dosages/probabilities will be added to the design matrix.
 * \param snpnum Number of the SNP for which the dosage/probability
 * data will be extracted from the gend object.
 * \param snpinfo An object of \ref mlinfo class that contains
 * information as read from the info file.
 * \param flipMAF A boolean indicating whether the reference and
 * coding alleles should be flipped based on Minor Allele Frequency
 * (as read from the \a snpinfo object). See also cmdvars::flipMAF.
*/
void coxph_data::update_snp(const gendata *gend,
                            const int snpnum,
                            mlinfo &snpinfo,
                            const bool flipMAF) {
    /*
     * This is the main part of the fix of bug #1846
     * (C) of the fix:
     *   UMC St Radboud Nijmegen,
     *   Dept of Epidemiology & Biostatistics,
     *   led by Prof. B. Kiemeney
     */
    /**
     * Note this sorts by "order"!!!
     * Here we deal with transposed X, hence last two arguments are swapped
     * compared to regdata::update_snp().
     * Also, the starting column-1 is not necessary for cox X therefore
     * 'ncov-j' changes to 'ncov-j-1'
     */

    // Reset counter for frequency since it is a new SNP
    gcount = 0;
    freq   = 0.0;

    for (int j = 0; j < ngpreds; j++)
    {
        // Row in X that  (note X is transposed when doing Cox
        // regression) contains the SNP data we are updating
        int snprow = ncov - j;
        // Double check: The number of rows in the X matrix should be
        // equal to ncov + 1 for mu.
        assert(snprow == X.nrow - 1 - j);

        double *snpdata = new double[nids];
        masked_data = std::vector<bool>(nids, false);

        gend->get_var(snpnum * ngpreds + j, snpdata);

        double *PA1A2 = new double[nids];
        if (flipMAF && ngpreds == 2 && j == 0) {
            // Read next genotype column because we need
            // P_A1A2 as well to calculate P_A2A2
            gend->get_var(snpnum * ngpreds + j + 1, PA1A2);
        }

        for (int i = 0; i < nids; i++) {
            if (flipMAF && (snpinfo.Freq1[snpnum] > 0.5)){
                // Flip the allele coding.
                snpinfo.allelesFlipped[snpnum] = true;
                // For dosage data: dosage is dosage_A1 (MaCH tutorial)
                //  so: dosage_A2 = 2 - dosage_A1
                //
                // For probability data:
                //   P_1 == P_A1A1  (MaCH tutorial)
                //   P_2 == P_A1A2  (MaCH tutorial)
                // and d_A1 = P_A1A2 + 2 P_A1A1, so when flipping:
                //   P_1 = P_A2A2 = 1 - P_A1A1 - P_A1A2== 1 - P_1 - P_2
                //   P_2 = P_A2A1 = P_A1A2 = P_2
                if (ngpreds == 1)
                {
                    // Dosage data
                    X.put(2 - snpdata[i], snprow, order[i]);
                }
                else if (ngpreds == 2 && j == 0) {
                    // Probability data, first probability (= P_A1A1)
                    X.put(1 - snpdata[i] - PA1A2[i], snprow, order[i]);
                }
                else if (ngpreds == 2 && j == 1)
                {
                    // Probability data, second probability
                    X.put(snpdata[i], snprow, order[i]);
                }
                else {
                    // You should never come here...
                    std::cerr << "Error: "
                              << "ngpreds != 1 or 2 while reading genetic data"
                              << std::endl;
                    exit(1);
                }
            } // end if (flipMAF && Freq1 > 0.5)
            else
            {
                // No flipping needed, simply copy the genetic data to the
                // snpdata array.
                X.put(snpdata[i], snprow, order[i]);
            }

            if (std::isnan(snpdata[i])) {
                masked_data[order[i]] = true;
                // SNP not masked
            } else {
                // check for first predictor
                if (j == 0) {
                    gcount++;
                    if (ngpreds == 1) {
                        freq += snpdata[i] * 0.5;
                    } else if (ngpreds == 2) {
                        freq += snpdata[i];
                    }
                } else if (j == 1) {
                    // Add second genotype in two predictor data form
                    freq += snpdata[i] * 0.5;
                }
            }    // end std::isnan(snpdata[i]) snp
        }    // end for loop: i=0 to nids

        delete[] snpdata;
    }    // End for loop: j=0 to ngpreds

    freq /= static_cast<double>(gcount);  // Allele frequency
}


/**
 * \brief Remove SNP information from the design matrix coxph_data::X.
 *
 * coxph_data::update_snp() adds SNP information to the design
 * matrix. This function allows you to strip that information from X
 * again. This is used for example when calculating the null model.
 *
 */
void coxph_data::remove_snp_from_X()
{
    if (ngpreds == 1)
    {
        X.delete_row(X.nrow -1);
    }
    else if (ngpreds == 2)
    {
        X.delete_row(X.nrow -1);
        X.delete_row(X.nrow -1);
    }
    else
    {
        cerr << "Error: ngpreds should be 1 or 2. "
             << "You should never come here!\n";
    }
}


/**
 * \brief Create a new coxph_data object that contains only the
 * non-masked data.
 *
 * The non-masked data is extracted according to the data in the
 * regdata::masked_data array. The resulting regdata::nids corresponds
 * to the number of IDs for which genotype data is present.
 *
 * Note that the regdata::masked_data array of the new object should
 * contain only zeros (i.e. not masked).
 *
 * @return A new regdata object containing only the rows from
 * regdata::X and regdata::Y for which genotype data is present.
 */
coxph_data coxph_data::get_unmasked_data() const
{
    coxph_data to;

    int nmeasured = std::count(masked_data.begin(), masked_data.end(), 0);
    to.nids = nmeasured;
    to.ncov = ncov;
    to.ngpreds = ngpreds;
    int dim1X = X.nrow;
    (to.weights).reinit(to.nids, 1);
    (to.stime).reinit(to.nids, 1);
    (to.sstat).reinit(to.nids, 1);
    (to.offset).reinit(to.nids, 1);
    (to.strata).reinit(to.nids, 1);
    (to.order).reinit(to.nids, 1);
    (to.X).reinit(dim1X, to.nids);

    int j = 0;
    for (int i = 0; i < nids; i++)
    {
        if (masked_data[i] == 0)
        {
            (to.weights).put(weights.get(i, 0), j, 0);
            (to.stime).put(stime.get(i, 0), j, 0);
            (to.sstat).put(sstat.get(i, 0), j, 0);
            (to.offset).put(offset.get(i, 0), j, 0);
            (to.strata).put(strata.get(i, 0), j, 0);
            (to.order).put(order.get(i, 0), j, 0);
            for (int nc = 0; nc < dim1X; nc++)
            {
                (to.X).put(X.get(nc, i), nc, j);
            }
            j++;
        }
    }

    to.masked_data = masked_data;
    return (to);
}


coxph_reg::coxph_reg(const coxph_data &cdatain)
{
    coxph_data cdata = cdatain.get_unmasked_data();
    beta.reinit(cdata.X.nrow, 1);
    sebeta.reinit(cdata.X.nrow, 1);
    loglik = - INFINITY;
    sigma2 = -1.;
    chi2_score = -1.;
    niter = 0;
}


void coxph_reg::estimate(const coxph_data &cdatain,
                         const int model,
                         const int interaction, const int ngpreds,
                         const bool iscox, const int nullmodel,
                         const mlinfo &snpinfo, const int cursnp)
{
    coxph_data cdata = cdatain.get_unmasked_data();

    mematrix<double> X = t_apply_model(cdata.X, model, interaction, ngpreds,
                                       iscox, nullmodel);

    int length_beta = X.nrow;
    beta.reinit(length_beta, 1);
    sebeta.reinit(length_beta, 1);
    mematrix<double> newoffset = cdata.offset -
        (cdata.offset).column_mean(0);
    mematrix<double> means(X.nrow, 1);

    for (int i = 0; i < X.nrow; i++)
    {
        beta[i] = 0.;
    }

    mematrix<double> u(X.nrow, 1);
    mematrix<double> imat(X.nrow, X.nrow);

    double *work = new double[X.ncol * 2 +
                              2 * (X.nrow) * (X.nrow) +
                              3 * (X.nrow)];
    double loglik_int[2];
    int flag;

    // Use Efron method of handling ties (for Breslow: 0.0), like in
    // R's coxph()
    double sctest = 1.0;

    // Set the maximum number of iterations that coxfit2() will run to
    // the default value from the class definition.
    int maxiterinput = MAXITER;
    // Make separate variables epsinput and tolcholinput that are not
    // const to send to coxfit2(), this way we won't have to alter
    // that function (which is a good thing: we want to keep it as
    // pristine as possible because it is copied from the R survival
    // package).
    double epsinput = EPS;
    double tolcholinput = CHOLTOL;

    coxfit2(&maxiterinput, &cdata.nids, &X.nrow, cdata.stime.data.data(),
            cdata.sstat.data.data(), X.data.data(), newoffset.data.data(),
            cdata.weights.data.data(), cdata.strata.data.data(),
            means.data.data(), beta.data.data(), u.data.data(),
            imat.data.data(), loglik_int, &flag, work, &epsinput,
            &tolcholinput, &sctest);

    // After coxfit2() maxiterinput contains the actual number of
    // iterations that were used. Store it in niter.
    niter = maxiterinput;


    // Check the results of the Cox fit; mirrored from the same checks
    // in coxph.fit.S and coxph.R from the R survival package.

    bool setToNAN = false;

    // Based on coxph.fit.S lines with 'which.sing' and coxph.R line
    // with if(any(is.NA(coefficients))). These lines set coefficients
    // to NA if flag < nvar (with nvar = ncol(x)) and MAXITER >
    // 0. coxph.R then checks for any NAs in the coefficients and
    // outputs the warning message if NAs were found.
    if (flag < X.nrow)
    {
        int which_sing = 0;
        MatrixXd imateigen = imat.data;
        VectorXd imatdiag = imateigen.diagonal();

        // Start at i=1 to exclude the beta coefficient for the
        // (constant) mean from the check.
        for (int i=1; i < imatdiag.size(); i++)
        {
            if (imatdiag[i] == 0)
                {
                    which_sing = i;
                    setToNAN = true;
                    std::cerr << "Warning for " << snpinfo.name[cursnp]
                              << ": X matrix deemed to be singular (variable "
                              << which_sing + 1 << ")" << std::endl;
                }
        }
    }

    if (niter >= MAXITER)
    {
        cerr << "Warning for " << snpinfo.name[cursnp]
             << ": nr of iterations > the maximum (" << MAXITER << "): "
             << niter << endl;
    }

    if (flag == 1000)
    {
        cerr << "Warning for " << snpinfo.name[cursnp]
             << ": Cox regression ran out of iterations and did not converge,"
             << " setting beta and se to 'NaN'\n";
        setToNAN = true;
    } else {
        VectorXd ueigen = u.data;
        MatrixXd imateigen = imat.data;
        VectorXd infs = ueigen.transpose() * imateigen;
        infs = infs.cwiseAbs();
        VectorXd betaeigen = beta.data;
        bool problems = false;

        assert(betaeigen.size() == infs.size());

        // We check the beta's for all coefficients
        // (incl. covariates), maybe stick to only checking the SNP
        // coefficient?
        for (int i = 0; i < infs.size(); i++) {
            if (infs[i] > EPS &&
                infs[i] > sqrt(EPS) * abs(betaeigen[i])) {
                problems = true;
            }
        }

        if (problems) {
            cerr << "Warning for " << snpinfo.name[cursnp]
                 << ": beta may be infinite,"
                 << " setting beta and se to 'NaN'\n";

            setToNAN = true;
        }
    }

    for (int i = 0; i < X.nrow; i++)
    {
        if (setToNAN)
        {
            // Cox regression failed
            sebeta[i] = NAN;
            beta[i]   = NAN;
            loglik    = NAN;
        } else {
            sebeta[i] = sqrt(imat.get(i, i));
            loglik = loglik_int[1];
        }
    }

    delete[] work;
}
