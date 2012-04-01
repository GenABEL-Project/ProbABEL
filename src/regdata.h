/*
 * regdata.h
 *
 *  Created on: Mar 29, 2012
 *      Author: mkooyman
 */

#ifndef REGDATA_H_
#define REGDATA_H_

#include "mematrix.h"
#include "mematri1.h"
#include "gendata.h"
#include "phedata.h"

class regdata
{
public:
    int nids;
    int ncov;
    int ngpreds;
    int noutcomes;
    unsigned short int * masked_data;
    mematrix<double> X;
    mematrix<double> Y;
    regdata();
    regdata(const regdata &obj);
    regdata(phedata &phed, gendata &gend, int snpnum);
    mematrix<double> extract_genotypes();
    void update_snp(gendata &gend, int snpnum);
    regdata get_unmasked_data();
    ~regdata();
private:

};

#endif /* REGDATA_H_ */
