/**
 * \file regdata.cpp
 * \author Yurii S. Aulchenko
 * \author M. Kooyman
 * \author L.C. Karssen
 * \author Maksim V. Struchalin
 *
 * \brief Describes functions of the regdata class containing a
 * linear or logistic regression object.
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


#include "fvlib/AbstractMatrix.h"
#include "fvlib/CastUtils.h"
#include "fvlib/const.h"
#include "fvlib/convert_util.h"
#include "fvlib/FileVector.h"
#include "fvlib/frutil.h"
#include "fvlib/frversion.h"
#include "fvlib/Logger.h"
#include "fvlib/Transposer.h"

#include <algorithm> // STL algoritms
#include "regdata.h"


/**
 * Constructor. Initialises all values to zero.
 *
 */
regdata::regdata()
{
    nids                    = 0;
    ncov                    = 0;
    ngpreds                 = 0;
    noutcomes               = 0;
    gcount                  = 0;
    freq                    = 0;
}


/**
 * Copy constructor. Creates a regdata object by copying the values of
 * another one.
 *
 * \param obj Reference to the regdata object to be copied to the new
 * object
 */
regdata::regdata(const regdata &obj) : masked_data(obj.masked_data),
                                       X(obj.X), Y(obj.Y)
{
    nids = obj.nids;
    ncov = obj.ncov;
    ngpreds = obj.ngpreds;
    noutcomes = obj.noutcomes;
    gcount = obj.gcount;
    freq = obj.freq;
}


/**
 * Constructor that fills a regdata object with phenotype and genotype
 * data.
 *
 * \param phed Reference to a phedata object with phenotype data
 * \param gend Reference to a gendata object with genotype data
 * \param snpnum The number of the SNP in the genotype data object to
 * be added to the design matrix regdata::X. When set to a number < 0
 * no SNP data is added to the design matrix (e.g. when calculating
 * the null model).
 */
regdata::regdata(const phedata &phed,
                 const gendata &gend,
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

    noutcomes = phed.noutcomes;
    X.reinit(nids, (ncov + 1));
    Y.reinit(nids, noutcomes);

    for (int i = 0; i < nids; i++)
    {
        X.put(1., i, 0);
        Y.put((phed.Y).get(i, 0), i, 0);
    }

    for (int j = 1; j <= phed.ncov; j++)
    {
        for (int i = 0; i < nids; i++)
        {
            X.put((phed.X).get(i, j - 1), i, j);
        }
    }

    if (snpnum > 0)
        for (int j = 0; j < ngpreds; j++)
        {
            double *snpdata = new double[nids];
            gend.get_var(snpnum * ngpreds + j, snpdata);
            for (int i = 0; i < nids; i++)
            {
                X.put(snpdata[i], i, (ncov - ngpreds + 1 + j));
            }
            delete[] snpdata;
        }
        // for (int i=0;i<nids;i++)
        //     for (int j=0;j<ngpreds;j++)
        //       X.put((gend.G).get(i,(snpnum*ngpreds+j)),i,(ncov-ngpreds+1+j));
}


/**
 * \brief Update the SNP dosages/probabilities in the design matrix
 * regdata::X.
 *
 * Adds the genetic information for a new SNP to the design
 * matrix. NOTE: For probability data, the order of the two
 * probabilities is reversed compared to the way they are stored in
 * the input file. Mach stores the probabilities as \f$P_{A_1A_1}\f$
 * \f$P_{A_1A_2}\f$.
 *
 * @param gend Object that contains the genetic data from which the
 * dosages/probabilities will be added to the design matrix.
 * @param snpnum Number of the SNP for which the dosage/probability
 * data will be extracted from the gend object.
 */
void regdata::update_snp(const gendata *gend, const int snpnum)
{
    // Reset counter for frequency since it is a new SNP
    gcount = 0;
    freq   = 0.0;

    // Add genotype data (dosage or probabilities) to the design
    // matrix X. Start filling from the last column, so for
    // probability data (ngpreds==2) the order of the two
    // probabilities is reversed.
    for (int j = 0; j < ngpreds; j++)
    {
        double *snpdata = new double[nids];
        masked_data = std::vector<bool>(nids, false);

        gend->get_var(snpnum * ngpreds + j, snpdata);

        for (int i = 0; i < nids; i++) {
            X.put(snpdata[i], i, (ncov - j));
            if (std::isnan(snpdata[i])) {
                masked_data[i] = true;
            } else {
                // SNP not masked
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
            }  // End if std::isnan(snpdata[i]) snp
        }  // End for loop: i = 0 to nids

        delete[] snpdata;
    }  // End for loop: j = 0 to ngpreds

    freq /= static_cast<double>(gcount); // Allele frequency
}


/**
 * \brief Remove SNP information from the design matrix regdata::X.
 *
 * \ref update_snp adds SNP information to the design matrix. This
 * function allows you to strip that information from X again.
 * This is used for example when calculating the null model.
 */
void regdata::remove_snp_from_X()
{
    if (ngpreds == 1)
    {
        X.delete_column(X.ncol -1);
    }
    else if (ngpreds == 2)
    {
        X.delete_column(X.ncol -1);
        X.delete_column(X.ncol -1);
    }
    else
    {
        cerr << "Error: ngpreds should be 1 or 2. "
             << "You should never come here!\n";
    }
}


/**
 * \brief Create a new regdata object that contains only the
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
regdata regdata::get_unmasked_data() const
{
    regdata to;
    int nmeasured = std::count(masked_data.begin(), masked_data.end(), 0);

    to.nids                    = nmeasured;
    to.ncov                    = ncov;
    to.ngpreds                 = ngpreds;
    to.noutcomes               = noutcomes;
    int dim2Y = Y.ncol;
    int dim2X = X.ncol;
    (to.X).reinit(to.nids, dim2X);
    (to.Y).reinit(to.nids, dim2Y);

    int j = 0;
    for (int i = 0; i < nids; i++)
    {
        if (masked_data[i] == 0)
        {
            for (int nc = 0; nc < dim2X; nc++)
            {
                (to.X).put(X.get(i, nc), j, nc);
            }

            for (int nc = 0; nc < dim2Y; nc++)
            {
                (to.Y).put(Y.get(i, nc), j, nc);
            }
            j++;
        }
    }

    to.masked_data = masked_data;
    return (to);
}


/**
 * Extracts the genotype data from the design matrix regdata::X of a
 * regdata object.
 *
 * @return A new mematrix object of dimensions regdata::X.nrow x
 * ngpreds that contains the genotype data.
 */
mematrix<double> regdata::extract_genotypes(void)
{
    mematrix<double> out;
    out.reinit(X.nrow, ngpreds);
    for (int i = 0; i < X.nrow; i++)
    {
        for (int j = 0; j < ngpreds; j++)
        {
            out[i * ngpreds + j] = X.get(i, (ncov - ngpreds + 1 + j));
        }
    }
    return out;
}
