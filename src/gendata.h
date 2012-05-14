/*
 * gendata.h
 *
 *  Created on: Mar 8, 2012
 *      Author: mkooyman
 */

#ifndef GENDATA_H_
#define GENDATA_H_
#include <string>
#include "fvlib/FileVector.h"

#if EIGEN
#include "eigen_mematrix.h"
#else
#include "mematrix.h"
#endif

class gendata
{
public:
    int nsnps;
    int nids;
    int ngpreds;
    gendata();
    void re_gendata(char * fname, int insnps, int ingpreds, int npeople,
            int nmeasured, unsigned short int * allmeasured, int skipd,
            std::string * idnames);
    void re_gendata(string filename, int insnps, int ingpreds, int npeople,
            int nmeasured, unsigned short int * allmeasured,
            std::string * idnames);
    void get_var(int var, float * data);
    ~gendata();

    // MAKE THAT PRIVATE, ACCESS THROUGH GET_SNP
    // ANOTHER PRIVATE OBJECT IS A POINTER TO DATABELBASECPP
    // UPDATE SNP, ALL REGRESSION METHODS: ACCOUNT FOR MISSING
private:
    mematrix<float> G;
    AbstractMatrix * DAG;
    unsigned short int * DAGmask;
    //	mematrix<double> G;
};

#endif /* GENDATA_H_ */
