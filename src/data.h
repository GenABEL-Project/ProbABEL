/*
 * data.h
 *
 *  Created on: Mar 8, 2012
 *      Author: mkooyman
 */

#ifndef DATA_H_
#define DATA_H_
#include <string>

extern bool is_interaction_excluded;

unsigned int Nmeasured(char * fname, int nphenocols, int npeople);
#include "phedata.h"
#include "gendata.h"

class mlinfo {
public:
    int nsnps;
    std::string * name;
    std::string * A1;
    std::string * A2;
    double * Freq1;
    double * MAF;
    double * Quality;
    double * Rsq;
    std::string * map;
    mlinfo()
    {
        Freq1 = NULL;
        MAF = NULL;
        Quality = NULL;
        Rsq = NULL;
        nsnps = 0;
        A1 = NULL;
        A2 = NULL;
        name = NULL;
        map = NULL;
    }
    mlinfo(char * filename, char * mapname);
    ~mlinfo();
};

class InvSigma {
private:
    static const unsigned MAXIMUM_PEOPLE_AMOUNT = 1000000;
    unsigned npeople; //amount of people
    std::string filename;
    mematrix<double> matrix; //file is stored here

public:
    InvSigma(const char * filename_, phedata * phe);
    mematrix<double> & get_matrix();
    ~InvSigma();
};

#endif /* DATA_H_ */
