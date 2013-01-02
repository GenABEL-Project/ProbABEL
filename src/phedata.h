/*
 * phedata.h
 *
 *  Created on: Mar 6, 2012
 *      Author: mkooyman
 */

#ifndef PHEDATA_H_
#define PHEDATA_H_

#if EIGEN
#include "eigen_mematrix.h"
#include "eigen_mematrix.cpp"
#else
#include "mematrix.h"
#include "mematri1.h"
#endif

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
        allmeasured = NULL;
        idnames     = NULL;
        model_terms = NULL;
    }

    phedata(char * fname, int noutc, int npeople, int interaction, bool iscox);

    void setphedata(char * fname, int noutc, int npeople, int interaction,
                    bool iscox);

    void set_is_interaction_excluded(bool inter_excluded);

    bool is_interaction_excluded;
    int nids_all;
    int nids;
    int noutcomes;
    int ncov;
    int n_model_terms;
    unsigned short int * allmeasured;
    std::string * idnames;
    std::string model;
    std::string * model_terms;
    mematrix<double> X;       /* Will contain the values of the covariate(s) */
    mematrix<double> Y;       /* Will contain the values of the outcome(s) */

    ~phedata();
};

#endif /* PHEDATA_H_ */
