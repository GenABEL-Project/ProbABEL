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

regdata::regdata() {
}
;

regdata::regdata(const regdata &obj) {
    nids = obj.nids;
    ncov = obj.ncov;
    ngpreds = obj.ngpreds;
    noutcomes = obj.noutcomes;
    X = obj.X;
    Y = obj.Y;
    masked_data = new unsigned short int[nids];
    for (int i = 0; i < nids; i++) {
        masked_data[i] = 0;
    }
}
regdata::regdata(phedata &phed, gendata &gend, int snpnum) {
    nids = gend.nids;
    masked_data = new unsigned short int[nids];
    for (int i = 0; i < nids; i++) {
        masked_data[i] = 0;
    }

    ngpreds = gend.ngpreds;
    if (snpnum >= 0) {
        ncov = phed.ncov + ngpreds;
    } else {
        ncov = phed.ncov;
    }
    noutcomes = phed.noutcomes;
    X.reinit(nids, (ncov + 1));
    Y.reinit(nids, noutcomes);
    for (int i = 0; i < nids; i++) {
        X.put(1., i, 0);
        Y.put((phed.Y).get(i, 0), i, 0);
    }
    for (int j = 1; j <= phed.ncov; j++)
        for (int i = 0; i < nids; i++)
            X.put((phed.X).get(i, j - 1), i, j);
    if (snpnum > 0)
        for (int j = 0; j < ngpreds; j++) {
            float snpdata[nids];
            gend.get_var(snpnum * ngpreds + j, snpdata);
            for (int i = 0; i < nids; i++)
                X.put(snpdata[i], i, (ncov - ngpreds + 1 + j));
        }
    //          for (int i=0;i<nids;i++)
    //              for (int j=0;j<ngpreds;j++)
    //                  X.put((gend.G).get(i,(snpnum*ngpreds+j)),i,(ncov-ngpreds+1+j));
}
void regdata::update_snp(gendata &gend, int snpnum) {
    for (int j = 0; j < ngpreds; j++) {
        float snpdata[nids];
        for (int i = 0; i < nids; i++)
            masked_data[i] = 0;
        gend.get_var(snpnum * ngpreds + j, snpdata);
        for (int i = 0; i < nids; i++) {
            X.put(snpdata[i], i, (ncov + 1 - j - 1));
            if (isnan(snpdata[i]))
                masked_data[i] = 1;
        }
    }
}
regdata::~regdata() {
    delete[] regdata::masked_data;
    //      delete X;
    //      delete Y;
}

regdata regdata::get_unmasked_data() {
    regdata to; // = regdata(*this);
    int nmeasured = 0;
    for (int i = 0; i < nids; i++)
        if (masked_data[i] == 0)
            nmeasured++;
    to.nids = nmeasured;
    //cout << to.nids << " in get_unmasked_data\n";
    to.ncov = ncov;
    to.ngpreds = ngpreds;
    to.noutcomes = noutcomes;
    int dim2Y = Y.ncol;
    int dim2X = X.ncol;
    (to.X).reinit(to.nids, dim2X);
    (to.Y).reinit(to.nids, dim2Y);

    int j = 0;
    for (int i = 0; i < nids; i++) {
        if (masked_data[i] == 0) {
            for (int nc = 0; nc < dim2X; nc++) {
                (to.X).put(X.get(i, nc), j, nc);
            }

            for (int nc = 0; nc < dim2Y; nc++) {
                (to.Y).put(Y.get(i, nc), j, nc);
            }
            j++;
        }
    }

    //delete [] to.masked_data;
    to.masked_data = new unsigned short int[to.nids];
    for (int i = 0; i < to.nids; i++)
        to.masked_data[i] = 0;
    //fprintf(stdout,"get_unmasked: %i %i %i\n",to.nids,dim2X,dim2Y);
    return (to);
}

mematrix<double> regdata::extract_genotypes(void) {
    mematrix<double> out;
    out.reinit(X.nrow, ngpreds);
    for (int i = 0; i < X.nrow; i++)
        for (int j = 0; j < ngpreds; j++)
            out[i * ngpreds + j] = X.get(i, (ncov - ngpreds + 1 + j));
    return out;
}

