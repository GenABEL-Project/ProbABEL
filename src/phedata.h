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
        nids_all = 0;
        nids = 0;
        noutcomes = 0;
        ncov = 0;
        allmeasured = NULL;
        idnames = NULL;
        is_interaction_excluded = 0;
        n_model_terms = 0;
        model_terms = NULL;

    }
    phedata(char * fname, int noutc, int npeople, int interaction, bool iscox);

    void setphedata(char * fname, int noutc, int npeople, int interaction,
            bool iscox);
    bool is_interaction_excluded;

    int nids_all;
    int nids;
    int noutcomes;
    int ncov;
    unsigned short int * allmeasured;
    mematrix<double> X;       /* Will contain the values of the covariate(s) */
    mematrix<double> Y;       /* Will contain the values of the outcome(s) */
    std::string * idnames;
    std::string model;
    std::string * model_terms;
    int n_model_terms;
    void set_is_interaction_excluded(bool inter_excluded);
    ~phedata();

};
#endif /* PHEDATA_H_ */
