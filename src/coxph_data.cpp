/*
 * coxph_data.cpp
 *
 *  Created on: Mar 31, 2012
 *      Author: mkooyman
 */
#include "coxph_data.h"
#include <iostream>
#include <cmath>
extern "C" {
#include "survproto.h"
}

//  #include "reg1.h"
#include "fvlib/AbstractMatrix.h"
#include "fvlib/CastUtils.h"
#include "fvlib/const.h"
#include "fvlib/convert_util.h"
#include "fvlib/FileVector.h"
#include "fvlib/frutil.h"
#include "fvlib/frversion.h"
#include "fvlib/Logger.h"
#include "fvlib/Transposer.h"

// compare for sort of times
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

coxph_data::coxph_data(const coxph_data &obj) : sstat(obj.sstat)
{
    nids = obj.nids;
    ncov = obj.ncov;
    ngpreds = obj.ngpreds;
    weights = obj.weights;
    stime = obj.stime;
    offset = obj.offset;
    strata = obj.strata;
    X = obj.X;
    order = obj.order;
    masked_data = new unsigned short int[nids];
    for (int i = 0; i < nids; i++)
        masked_data[i] = 0;
}

coxph_data::coxph_data(phedata &phed, gendata &gend, int snpnum)
{
    nids = gend.nids;
    masked_data = new unsigned short int[nids];
    for (int i = 0; i < nids; i++)
        masked_data[i] = 0;

    ngpreds = gend.ngpreds;
    if (snpnum >= 0)
        ncov = phed.ncov + ngpreds;
    else
        ncov = phed.ncov;

    if (phed.noutcomes != 2)
    {
        std::cerr << "coxph_data: number of outcomes should be 2 (now: "
                  << phed.noutcomes << ")\n";
        exit(1);
    }

    //      X.reinit(nids,(ncov+1));
    X.reinit(nids, ncov);
    stime.reinit(nids, 1);
    sstat.reinit(nids, 1);
    weights.reinit(nids, 1);
    offset.reinit(nids, 1);
    strata.reinit(nids, 1);
    order.reinit(nids, 1);
    for (int i = 0; i < nids; i++)
    {
        //          X.put(1.,i,0);
        stime[i] = (phed.Y).get(i, 0);
        sstat[i] = int((phed.Y).get(i, 1));
        if (sstat[i] != 1 && sstat[i] != 0)
        {
          std::cerr << "coxph_data: status not 0/1 (correct order: id, fuptime, status ...)"
                    << endl;
            exit(1);
        }
    }

    for (int j = 0; j < phed.ncov; j++)
        for (int i = 0; i < nids; i++)
            X.put((phed.X).get(i, j), i, j);

    if (snpnum > 0)
        for (int j = 0; j < ngpreds; j++)
        {
            float snpdata[nids];
            gend.get_var(snpnum * ngpreds + j, snpdata);
            for (int i = 0; i < nids; i++)
                X.put(snpdata[i], i, (ncov - ngpreds + j));
        }

    //          for (int i=0;i<nids;i++)
    //              for (int j=0;j<ngpreds;j++)
    //                  X.put((gend.G).get(i,(snpnum*ngpreds+j)),i,(ncov-ngpreds+j));

    for (int i = 0; i < nids; i++)
    {
        weights[i] = 1.0;
        offset[i] = 0.0;
        strata[i] = 0;
    }
    // sort by time
    double tmptime[nids];
    int passed_sorted[nids];
    for (int i = 0; i < nids; i++)
    {
        tmptime[i] = stime[i];
        passed_sorted[i] = 0;
    }
    qsort(tmptime, nids, sizeof(double), cmpfun);
    for (int i = 0; i < nids; i++)
    {
        int passed = 0;
        for (int j = 0; j < nids; j++)
            if (tmptime[j] == stime[i])
                if (!passed_sorted[j])
                {
                    order[i] = j;
                    passed_sorted[j] = 1;
                    passed = 1;
                    break;
                }
        if (passed != 1)
        {
            std::cerr << "cannot recover element " << i << "\n";
            exit(1);
        }
    }
    stime = reorder(stime, order);
    sstat = reorder(sstat, order);
    weights = reorder(weights, order);
    strata = reorder(strata, order);
    offset = reorder(offset, order);
    X = reorder(X, order);
    X = transpose(X);
    //      X.print();
    //      offset.print();
    //      weights.print();
    //      stime.print();
    //      sstat.print();
}
void coxph_data::update_snp(gendata &gend, int snpnum)
{
  /**
   * This is the main part of the fix of bug #1846
   * (C) of the fix:
   *   UMC St Radboud Nijmegen,
   *   Dept of Epidemiology & Biostatistics,
   *   led by Prof. B. Kiemeney
   *
   * Note this sorts by "order"!!!
   * Here we deal with transposed X, hence last two arguments are swapped
   * compared to the other 'update_snp'
   * Also, the starting column-1 is not necessary for cox X therefore
   * 'ncov-j' changes to 'ncov-j-1'
   **/

    for (int j = 0; j < ngpreds; j++)
    {
        float snpdata[nids];
        for (int i = 0; i < nids; i++)
            masked_data[i] = 0;

        gend.get_var(snpnum * ngpreds + j, snpdata);

        for (int i = 0; i < nids; i++)
        {
            X.put(snpdata[i], (ncov - j - 1), order[i]);
            if (isnan(snpdata[i]))
                masked_data[order[i]] = 1;
        }
    }
}

coxph_data::~coxph_data()
{
    delete[] coxph_data::masked_data;
    //      delete X;
    //      delete sstat;
    //      delete stime;
    //      delete weights;
    //      delete offset;
    //      delete strata;
    //      delete order;
}

coxph_data coxph_data::get_unmasked_data()
{
//      std::cout << " !!! in get_unmasked_data !!! ";
    coxph_data to; // = coxph_data(*this);
    // filter missing data

    int nmeasured = 0;
    for (int i = 0; i < nids; i++)
        if (masked_data[i] == 0)
            nmeasured++;
    to.nids = nmeasured;
//      std::cout << "nmeasured=" << nmeasured << "\n";
    to.ncov = ncov;
//      std::cout << "ncov=" << ncov << "\n";
    to.ngpreds = ngpreds;
    int dim1X = X.nrow;
//      std::cout << "X.ncol=" << X.ncol << "\n";
//      std::cout << "X.nrow=" << X.nrow << "\n";
    (to.weights).reinit(to.nids, 1);
    (to.stime).reinit(to.nids, 1);
    (to.sstat).reinit(to.nids, 1);
    (to.offset).reinit(to.nids, 1);
    (to.strata).reinit(to.nids, 1);
    (to.order).reinit(to.nids, 1);
    (to.X).reinit(dim1X, to.nids);

//      std::cout << "(to.X).ncol=" << (to.X).ncol << "\n";
//      std::cout << "(to.X).nrow=" << (to.X).nrow << "\n";
//      std::cout << " !!! just before cycle !!! ";
    int j = 0;
    for (int i = 0; i < nids; i++)
    {
//          std::cout << nids << " " << i << " " << masked_data[i] << "\n";
        if (masked_data[i] == 0)
        {
            (to.weights).put(weights.get(i, 1), j, 1);
            (to.stime).put(stime.get(i, 1), j, 1);
            (to.sstat).put(sstat.get(i, 1), j, 1);
            (to.offset).put(offset.get(i, 1), j, 1);
            (to.strata).put(strata.get(i, 1), j, 1);
            (to.order).put(order.get(i, 1), j, 1);
            for (int nc = 0; nc < dim1X; nc++)
                (to.X).put(X.get(nc, i), nc, j);
            j++;
        }
    }
//      std::cout << " !!! just after cycle !!! ";

    //delete [] to.masked_data;
    to.masked_data = new unsigned short int[to.nids];
    for (int i = 0; i < to.nids; i++)
        to.masked_data[i] = 0;
    // std::cout << "get_unmasked: " << to.nids << " "
    //           << dim2X << " " << dim2Y << "\n";
    return (to);
}

coxph_reg::coxph_reg(coxph_data &cdatain)
{
    coxph_data cdata = cdatain.get_unmasked_data();
    beta.reinit(cdata.X.nrow, 1);
    sebeta.reinit(cdata.X.nrow, 1);
    loglik = -9.999e+32;
    sigma2 = -1.;
    chi2_score = -1.;
    niter = 0;
}

void coxph_reg::estimate(coxph_data &cdatain, int verbose, int maxiter,
                         double eps, double tol_chol, int model,
                         int interaction, int ngpreds, bool iscox,
                         int nullmodel)
{
    coxph_data cdata = cdatain.get_unmasked_data();
    mematrix<double> X = t_apply_model(cdata.X, model, interaction, ngpreds,
                                       iscox, nullmodel);
    //      X.print();
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
    double work[X.ncol * 2 + 2 * (X.nrow) * (X.nrow) + 3 * (X.nrow)];
    double loglik_int[2];
    int flag;
    double sctest = 1.0;

    //TODO(maarten): this function works only  in combination with eigen remove .data() from arguments to make is available under old matrix
    coxfit2(&maxiter, &cdata.nids, &X.nrow, cdata.stime.data.data(),
            cdata.sstat.data.data(), X.data.data(), newoffset.data.data(),
            cdata.weights.data.data(), cdata.strata.data.data(),
            means.data.data(), beta.data.data(), u.data.data(),
            imat.data.data(), loglik_int, &flag, work, &eps, &tol_chol,
            &sctest);

    for (int i = 0; i < X.nrow; i++)
    {
        sebeta[i] = sqrt(imat.get(i, i));
    }
    loglik = loglik_int[1];
    niter = maxiter;
}
