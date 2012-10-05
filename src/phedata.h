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
#else
#include "mematrix.h"
#endif

class phedata
{

public:
    phedata()
    {

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
    mematrix<double> X;
    mematrix<double> Y;
    std::string * idnames;
    std::string model;
    std::string * model_terms;
    int n_model_terms;
    void set_is_interaction_excluded(bool inter_excluded);
    ~phedata();

};
#endif /* PHEDATA_H_ */