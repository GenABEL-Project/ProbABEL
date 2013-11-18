/*
 * regdata.cpp
 *
 *  Created on: Mar 29, 2012
 *      Author: mkooyman
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

#include "regdata.h"

regdata::regdata()
{
    nids                    = 0;
    ncov                    = 0;
    ngpreds                 = 0;
    noutcomes               = 0;
    is_interaction_excluded = false;
    masked_data             = NULL;
}


regdata::regdata(const regdata &obj) : X(obj.X), Y(obj.Y)
{
    nids = obj.nids;
    ncov = obj.ncov;
    ngpreds = obj.ngpreds;
    noutcomes = obj.noutcomes;
    is_interaction_excluded = obj.is_interaction_excluded;
    masked_data = new unsigned short int[nids];

    for (int i = 0; i < nids; i++)
    {
        masked_data[i] = obj.masked_data[i];
    }
}

regdata::regdata(phedata &phed, gendata &gend, int snpnum,
                 bool ext_is_interaction_excluded)
{
    nids = gend.nids;
    masked_data = new unsigned short int[nids];
    for (int i = 0; i < nids; i++)
    {
        masked_data[i] = 0;
    }

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
    is_interaction_excluded = ext_is_interaction_excluded;
}

void regdata::update_snp(gendata &gend, int snpnum)
{
    // Add genotypic data (dosage or probabilities) to the design
    // matrix X.

    for (int j = 0; j < ngpreds; j++)
    {
        double *snpdata = new double[nids];
        for (int i = 0; i < nids; i++)
        {
            masked_data[i] = 0;
        }

        gend.get_var(snpnum * ngpreds + j, snpdata);

        for (int i = 0; i < nids; i++)
        {
            X.put(snpdata[i], i, (ncov - j));
            if (std::isnan(snpdata[i]))
            {
                masked_data[i] = 1;
            }
        }
        delete[] snpdata;
    }
}

void regdata::remove_snp_from_X()
{
    // update_snp() adds SNP information to the design matrix. This
    // function allows you to strip that information from X again.
    // This is used for example when calculating the null model.

    if(ngpreds == 1)
    {
        X.delete_column(X.ncol -1);
    }
    else if(ngpreds == 2)
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

regdata::~regdata()
{
    delete[] regdata::masked_data;
    //      delete X;
    //      delete Y;
}

regdata regdata::get_unmasked_data()
{
    regdata to; // = regdata(*this);
    int nmeasured = 0;
    for (int i = 0; i < nids; i++)
    {
        if (masked_data[i] == 0)
            nmeasured++;
    }

    to.nids                    = nmeasured;
    to.ncov                    = ncov;
    to.ngpreds                 = ngpreds;
    to.noutcomes               = noutcomes;
    to.is_interaction_excluded = is_interaction_excluded;
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

    //delete [] to.masked_data;
    to.masked_data = new unsigned short int[to.nids];
    for (int i = 0; i < to.nids; i++)
    {
        to.masked_data[i] = 0;
    }
    // std::cout << "get_unmasked: " << to.nids << " "
    //           << dim2X << " " << dim2Y << "\n";
    return (to);
}

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
