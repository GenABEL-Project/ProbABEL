/*
 * coxph_data.h
 *
 *  Created on: Mar 31, 2012
 *      Author: mkooyman
 */

#ifndef COXPH_DATA_H_
#define COXPH_DATA_H_

class coxph_data
{
public:
    int nids;
    int ncov;
    int ngpreds;
    mematrix<double> weights;
    mematrix<double> stime;
    mematrix<int> sstat;
    mematrix<double> offset;
    mematrix<int> strata;
    mematrix<double> X;
    mematrix<int> order;
    unsigned short int * masked_data;
    coxph_data get_unmasked_data();
    coxph_data()
    {
    }
    coxph_data(const coxph_data &obj);
    coxph_data(phedata &phed, gendata &gend, int snpnum);
    void update_snp(gendata &gend, int snpnum);
    ~coxph_data();

};

#endif /* COXPH_DATA_H_ */
